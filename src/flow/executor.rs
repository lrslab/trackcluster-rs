//! Bounded per-gene worker execution with centralized panic and progress handling.

use std::panic::{self, AssertUnwindSafe};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};
use std::time::{Duration, Instant};

use crate::flow::full::{
    BatchRunOptions, DownsampleRecord, GeneFailureKind, GeneOutcome, GeneSkipReason,
    ProcessGeneResult, ResumeDecision,
};

/// Complete, deterministic accounting from one bounded worker run.
pub(super) struct ExecutionReport {
    pub(super) processed: usize,
    pub(super) skipped: usize,
    pub(super) skipped_completed_outputs: usize,
    pub(super) skipped_empty_reads: usize,
    pub(super) skipped_no_usable_reads: usize,
    pub(super) rejected_read_tracks: usize,
    pub(super) genes_with_rejected_reads: usize,
    pub(super) prepare_rejected_read_tracks: usize,
    pub(super) errors: usize,
    pub(super) failed_missing_inputs: usize,
    pub(super) failed_processing: usize,
    pub(super) failed_panics: usize,
    pub(super) elapsed: Duration,
    pub(super) worker_count: usize,
    pub(super) error_lines: Vec<String>,
    pub(super) downsample_records: Vec<DownsampleRecord>,
    pub(super) resume_decisions: Vec<ResumeDecision>,
    /// Genes with a complete, mergeable output set from this execution.
    pub(super) mergeable_genes: Vec<String>,
}

impl ExecutionReport {
    pub(super) fn failed_gene_count(&self) -> usize {
        self.resume_decisions
            .iter()
            .filter(|decision| decision.action == "fail")
            .count()
    }

    pub(super) fn infrastructure_error_count(&self) -> usize {
        self.errors.saturating_sub(self.failed_gene_count())
    }
}

#[derive(Debug, Default)]
struct WorkerState {
    gene: Option<String>,
    started_at: Option<Instant>,
}

fn panic_message(payload: Box<dyn std::any::Any + Send>) -> String {
    payload
        .downcast_ref::<&str>()
        .copied()
        .map(str::to_owned)
        .or_else(|| payload.downcast_ref::<String>().cloned())
        .unwrap_or_else(|| "<non-string panic>".to_owned())
}

fn cloned_records<T: Clone>(records: &Arc<Mutex<Vec<T>>>) -> Vec<T> {
    records
        .lock()
        .unwrap_or_else(|poisoned| poisoned.into_inner())
        .clone()
}

/// Execute a finite gene list using no more than the configured number of workers.
///
/// Every task panic becomes a gene-level failure and every unexpected worker panic becomes a
/// synthetic error record, so the caller has one error-propagation path.
pub(super) fn execute_genes(
    genes: Vec<String>,
    args: &BatchRunOptions,
    process_gene: fn(&str, &BatchRunOptions) -> ProcessGeneResult,
) -> ExecutionReport {
    let total = genes.len();
    let worker_count = args.runtime.worker_count(total);
    let started = Instant::now();

    let processed = Arc::new(AtomicUsize::new(0));
    let skipped = Arc::new(AtomicUsize::new(0));
    let skipped_completed_outputs = Arc::new(AtomicUsize::new(0));
    let skipped_empty_reads = Arc::new(AtomicUsize::new(0));
    let skipped_no_usable_reads = Arc::new(AtomicUsize::new(0));
    let rejected_read_tracks = Arc::new(AtomicUsize::new(0));
    let genes_with_rejected_reads = Arc::new(AtomicUsize::new(0));
    let errors = Arc::new(AtomicUsize::new(0));
    let failed_missing_inputs = Arc::new(AtomicUsize::new(0));
    let failed_processing = Arc::new(AtomicUsize::new(0));
    let failed_panics = Arc::new(AtomicUsize::new(0));
    let done = Arc::new(AtomicUsize::new(0));
    let error_lines = Arc::new(Mutex::new(Vec::new()));
    let downsample_records = Arc::new(Mutex::new(Vec::new()));
    let resume_decisions = Arc::new(Mutex::new(Vec::new()));
    let mergeable_genes = Arc::new(Mutex::new(Vec::new()));
    let queue = Arc::new(Mutex::new(genes));
    let worker_states: Arc<Vec<Mutex<WorkerState>>> = Arc::new(
        (0..worker_count)
            .map(|_| Mutex::new(WorkerState::default()))
            .collect(),
    );

    let (heartbeat_stop_tx, heartbeat_handle) = if args.runtime.heartbeat_seconds > 0 {
        use std::sync::mpsc::{self, RecvTimeoutError};
        let (stop_tx, stop_rx) = mpsc::channel::<()>();
        let processed = Arc::clone(&processed);
        let skipped = Arc::clone(&skipped);
        let errors = Arc::clone(&errors);
        let done = Arc::clone(&done);
        let queue = Arc::clone(&queue);
        let worker_states = Arc::clone(&worker_states);
        let heartbeat_seconds = args.runtime.heartbeat_seconds;
        let heartbeat_top = args.runtime.heartbeat_top.max(1);

        let handle = std::thread::spawn(move || {
            let mut last_done = done.load(Ordering::Relaxed);
            loop {
                match stop_rx.recv_timeout(Duration::from_secs(heartbeat_seconds)) {
                    Ok(()) | Err(RecvTimeoutError::Disconnected) => break,
                    Err(RecvTimeoutError::Timeout) => {}
                }

                let done_now = done.load(Ordering::Relaxed);
                let processed_now = processed.load(Ordering::Relaxed);
                let skipped_now = skipped.load(Ordering::Relaxed);
                let errors_now = errors.load(Ordering::Relaxed);
                let queue_remaining = queue.lock().map(|guard| guard.len()).unwrap_or_default();
                eprintln!(
                    "heartbeat: {done_now}/{total} (processed={processed_now}, skipped={skipped_now}, errors={errors_now}) queue_remaining={queue_remaining} elapsed={:?}",
                    started.elapsed()
                );

                if done_now == last_done && done_now < total {
                    let mut inflight = Vec::new();
                    for state_lock in worker_states.iter() {
                        let Ok(state) = state_lock.lock() else {
                            continue;
                        };
                        let (Some(gene), Some(started_at)) =
                            (state.gene.as_ref(), state.started_at)
                        else {
                            continue;
                        };
                        inflight.push((gene.clone(), started_at.elapsed()));
                    }
                    if inflight.is_empty() {
                        eprintln!("heartbeat: no in-flight genes (all workers idle?)");
                    } else {
                        inflight.sort_by(|left, right| right.1.cmp(&left.1));
                        let mut line = String::from("in_flight(top):");
                        for (gene, duration) in inflight.into_iter().take(heartbeat_top) {
                            line.push(' ');
                            line.push_str(&format!("{gene}={:.1}s", duration.as_secs_f64()));
                        }
                        eprintln!("{line}");
                    }
                }
                last_done = done_now;
                if done_now >= total {
                    break;
                }
            }
        });
        (Some(stop_tx), Some(handle))
    } else {
        (None, None)
    };

    let mut handles = Vec::with_capacity(worker_count);
    for worker_idx in 0..worker_count {
        let queue = Arc::clone(&queue);
        let processed = Arc::clone(&processed);
        let skipped = Arc::clone(&skipped);
        let skipped_completed_outputs = Arc::clone(&skipped_completed_outputs);
        let skipped_empty_reads = Arc::clone(&skipped_empty_reads);
        let skipped_no_usable_reads = Arc::clone(&skipped_no_usable_reads);
        let rejected_read_tracks = Arc::clone(&rejected_read_tracks);
        let genes_with_rejected_reads = Arc::clone(&genes_with_rejected_reads);
        let errors = Arc::clone(&errors);
        let failed_missing_inputs = Arc::clone(&failed_missing_inputs);
        let failed_processing = Arc::clone(&failed_processing);
        let failed_panics = Arc::clone(&failed_panics);
        let done = Arc::clone(&done);
        let error_lines = Arc::clone(&error_lines);
        let downsample_records = Arc::clone(&downsample_records);
        let resume_decisions = Arc::clone(&resume_decisions);
        let mergeable_genes = Arc::clone(&mergeable_genes);
        let worker_states = Arc::clone(&worker_states);
        let args = args.clone();

        handles.push(std::thread::spawn(move || loop {
            let gene = queue
                .lock()
                .unwrap_or_else(|poisoned| poisoned.into_inner())
                .pop();
            let Some(gene) = gene else {
                break;
            };
            if let Ok(mut state) = worker_states[worker_idx].lock() {
                state.gene = Some(gene.clone());
                state.started_at = Some(Instant::now());
            }

            let result = match panic::catch_unwind(AssertUnwindSafe(|| process_gene(&gene, &args)))
            {
                Ok(result) => result,
                Err(payload) => {
                    ProcessGeneResult::failed(GeneFailureKind::Panic, panic_message(payload))
                }
            };
            if let Ok(mut state) = worker_states[worker_idx].lock() {
                state.gene = None;
                state.started_at = None;
            }

            if result.rejected_read_tracks > 0 {
                rejected_read_tracks.fetch_add(result.rejected_read_tracks, Ordering::Relaxed);
                genes_with_rejected_reads.fetch_add(1, Ordering::Relaxed);
            }
            if result.all_reads_rejected {
                skipped_no_usable_reads.fetch_add(1, Ordering::Relaxed);
            }

            let (resume_action, mergeable) = match &result.outcome {
                GeneOutcome::Processed => {
                    processed.fetch_add(1, Ordering::Relaxed);
                    ("rebuild", true)
                }
                GeneOutcome::Skipped(reason) => {
                    skipped.fetch_add(1, Ordering::Relaxed);
                    match reason {
                        GeneSkipReason::CompletedOutputs => {
                            skipped_completed_outputs.fetch_add(1, Ordering::Relaxed);
                            ("reuse", true)
                        }
                        GeneSkipReason::EmptyReads => {
                            skipped_empty_reads.fetch_add(1, Ordering::Relaxed);
                            ("skip", true)
                        }
                        GeneSkipReason::NoUsableReads => ("skip", true),
                    }
                }
                GeneOutcome::Failed(failure) => {
                    errors.fetch_add(1, Ordering::Relaxed);
                    match failure.kind {
                        GeneFailureKind::MissingInputs => {
                            failed_missing_inputs.fetch_add(1, Ordering::Relaxed);
                        }
                        GeneFailureKind::Processing => {
                            failed_processing.fetch_add(1, Ordering::Relaxed);
                        }
                        GeneFailureKind::Panic => {
                            failed_panics.fetch_add(1, Ordering::Relaxed);
                        }
                    }
                    error_lines
                        .lock()
                        .unwrap_or_else(|poisoned| poisoned.into_inner())
                        .push(format!(
                            "{gene}\t{}\t{}",
                            failure.kind.as_str(),
                            failure.message
                        ));
                    ("fail", false)
                }
            };
            if mergeable {
                mergeable_genes
                    .lock()
                    .unwrap_or_else(|poisoned| poisoned.into_inner())
                    .push(gene.clone());
            }
            resume_decisions
                .lock()
                .unwrap_or_else(|poisoned| poisoned.into_inner())
                .push(ResumeDecision {
                    gene: gene.clone(),
                    action: resume_action,
                    reason: result.resume_reason,
                });
            if let Some(record) = result.downsample {
                downsample_records
                    .lock()
                    .unwrap_or_else(|poisoned| poisoned.into_inner())
                    .push(record);
            }

            let done_now = done.fetch_add(1, Ordering::Relaxed) + 1;
            let every = args.runtime.progress_every.max(1);
            if done_now.is_multiple_of(every) || done_now == total {
                eprintln!(
                    "progress: {done_now}/{total} (processed={}, skipped={}, errors={})",
                    processed.load(Ordering::Relaxed),
                    skipped.load(Ordering::Relaxed),
                    errors.load(Ordering::Relaxed),
                );
            }
        }));
    }

    for (worker_idx, handle) in handles.into_iter().enumerate() {
        if handle.join().is_err() {
            errors.fetch_add(1, Ordering::Relaxed);
            failed_panics.fetch_add(1, Ordering::Relaxed);
            error_lines
                .lock()
                .unwrap_or_else(|poisoned| poisoned.into_inner())
                .push(format!(
                    "<worker-{worker_idx}>\tpanic\tworker thread terminated unexpectedly"
                ));
        }
    }
    if let Some(tx) = heartbeat_stop_tx {
        let _ = tx.send(());
    }
    if let Some(handle) = heartbeat_handle {
        if handle.join().is_err() {
            errors.fetch_add(1, Ordering::Relaxed);
            failed_panics.fetch_add(1, Ordering::Relaxed);
            error_lines
                .lock()
                .unwrap_or_else(|poisoned| poisoned.into_inner())
                .push("<heartbeat>\tpanic\theartbeat thread terminated unexpectedly".to_owned());
        }
    }

    let mut mergeable_genes = cloned_records(&mergeable_genes);
    mergeable_genes.sort();
    mergeable_genes.dedup();

    ExecutionReport {
        processed: processed.load(Ordering::Relaxed),
        skipped: skipped.load(Ordering::Relaxed),
        skipped_completed_outputs: skipped_completed_outputs.load(Ordering::Relaxed),
        skipped_empty_reads: skipped_empty_reads.load(Ordering::Relaxed),
        skipped_no_usable_reads: skipped_no_usable_reads.load(Ordering::Relaxed),
        rejected_read_tracks: rejected_read_tracks.load(Ordering::Relaxed),
        genes_with_rejected_reads: genes_with_rejected_reads.load(Ordering::Relaxed),
        prepare_rejected_read_tracks: 0,
        errors: errors.load(Ordering::Relaxed),
        failed_missing_inputs: failed_missing_inputs.load(Ordering::Relaxed),
        failed_processing: failed_processing.load(Ordering::Relaxed),
        failed_panics: failed_panics.load(Ordering::Relaxed),
        elapsed: started.elapsed(),
        worker_count,
        error_lines: cloned_records(&error_lines),
        downsample_records: cloned_records(&downsample_records),
        resume_decisions: cloned_records(&resume_decisions),
        mergeable_genes,
    }
}
