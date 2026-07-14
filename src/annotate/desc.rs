use std::collections::{HashMap, HashSet};

use crate::annotate::addgene::{add_gene, AddGeneOpts};
use crate::model::{Strand, Transcript};

const FIGURE2_UTR_MIN_PERCENT: u128 = 5;

#[derive(Clone, Copy, Debug)]
pub struct DescOpts {
    pub offset_bp: u32,
    pub end_shift_bp: u32,
    pub fusion_fraction_read: f64,
    pub fusion_fraction_ref: f64,
}

impl Default for DescOpts {
    fn default() -> Self {
        Self {
            offset_bp: 10,
            end_shift_bp: 0,
            fusion_fraction_read: 0.1,
            fusion_fraction_ref: 0.1,
        }
    }
}

impl DescOpts {
    /// Validate fusion overlap fractions and materialize typed offsets.
    pub fn validate(self) -> Result<(), crate::config::ParameterError> {
        let _ = crate::config::BasePairOffset::new(self.offset_bp);
        let _ = crate::config::BasePairOffset::new(self.end_shift_bp);
        crate::config::UnitFraction::new(
            "fusion read overlap fraction",
            self.fusion_fraction_read,
        )?;
        crate::config::UnitFraction::new(
            "fusion reference overlap fraction",
            self.fusion_fraction_ref,
        )?;
        Ok(())
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DescRow {
    pub isoform_id: String,
    pub ref_id: String,
    pub gene: String,
    pub miss: String,
    pub extra: String,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Class4Row {
    pub isoform_id: String,
    pub class: String,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FusionRow {
    pub isoform_id: String,
    pub genes: Vec<String>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Class12Row {
    pub isoform_id: String,
    pub class: String,
}

#[derive(Clone, Debug, Default)]
pub struct DescResult {
    pub desc_rows: Vec<DescRow>,
    pub class4_rows: Vec<Class4Row>,
    pub fusion_rows: Vec<FusionRow>,
    pub class12_rows: Vec<Class12Row>,
}

fn ttype(tx: &Transcript) -> Option<&str> {
    tx.metadata().transcript_type()
}

fn is_isoform_anno(tx: &Transcript) -> bool {
    matches!(ttype(tx), Some("isoform_anno"))
}

fn gene_name(tx: &Transcript) -> &str {
    tx.metadata().gene_id().unwrap_or("none")
}

fn genes(tx: &Transcript) -> Vec<&str> {
    let gene_field = gene_name(tx);
    if gene_field == "none" {
        return Vec::new();
    }
    gene_field
        .split("||")
        .map(str::trim)
        .filter(|gene| !gene.is_empty() && *gene != "none")
        .collect()
}

fn same_comparison_locus(tx: &Transcript, reference: &Transcript) -> bool {
    tx.chrom == reference.chrom
        && (tx.strand == Strand::Unknown
            || reference.strand == Strand::Unknown
            || tx.strand == reference.strand)
}

fn junction_positions(tx: &Transcript) -> Vec<u32> {
    let mut boundaries: Vec<u32> = Vec::with_capacity(tx.exons.len() * 2);
    for exon in &tx.exons {
        boundaries.push(exon.start.get());
        boundaries.push(exon.end.get());
    }

    match tx.strand {
        Strand::Plus | Strand::Unknown => {}
        Strand::Minus => boundaries.reverse(),
    }

    if boundaries.len() <= 2 {
        return Vec::new();
    }

    boundaries[1..boundaries.len() - 1].to_vec()
}

fn ordered_boundary_matches(a: &[u32], b: &[u32], offset: u32) -> Vec<(usize, usize)> {
    crate::matching::ordered_one_to_one_matches_by(a.len(), b.len(), |a_idx, b_idx| {
        let delta = a[a_idx].abs_diff(b[b_idx]);
        (delta <= offset).then_some(u64::from(delta))
    })
}

fn junctions_equal(a: &[u32], b: &[u32], offset: u32) -> bool {
    a.len() == b.len() && ordered_boundary_matches(a, b, offset).len() == a.len()
}

fn compare_ei_by_boundary(a: &[u32], reference: &[u32], offset: u32) -> (Vec<usize>, Vec<usize>) {
    let matches = ordered_boundary_matches(a, reference, offset);
    let mut matched_a = vec![false; a.len()];
    let mut matched_reference = vec![false; reference.len()];
    for (a_idx, reference_idx) in matches {
        matched_a[a_idx] = true;
        matched_reference[reference_idx] = true;
    }

    let missed_order = matched_reference
        .iter()
        .enumerate()
        .filter_map(|(idx, matched)| (!matched).then_some(idx))
        .collect();
    let extra_order = matched_a
        .iter()
        .enumerate()
        .filter_map(|(idx, matched)| (!matched).then_some(idx))
        .collect();

    (missed_order, extra_order)
}

fn group_site(indices: &[usize]) -> Vec<Vec<usize>> {
    let mut groups: Vec<Vec<usize>> = Vec::new();
    for &idx in indices {
        match groups.last_mut() {
            Some(group) if idx == group.last().copied().unwrap_or(idx).saturating_add(1) => {
                group.push(idx);
            }
            _ => groups.push(vec![idx]),
        }
    }
    groups
}

fn has_new_junction(tx: &Transcript, refs: &[Transcript], offset: u32) -> bool {
    let junctions = junction_positions(tx);
    if junctions.is_empty() {
        return false;
    }

    let mut ref_junctions: Vec<u32> = Vec::new();
    for reference in refs {
        ref_junctions.extend(junction_positions(reference));
    }
    ref_junctions.sort_unstable();
    ref_junctions.dedup();
    if tx.strand == Strand::Minus {
        ref_junctions.reverse();
    }

    ordered_boundary_matches(&junctions, &ref_junctions, offset).len() != junctions.len()
}

fn class4(tx: &Transcript, refs: &[Transcript], offset: u32) -> String {
    if has_new_junction(tx, refs, offset) {
        return "new_junction".to_owned();
    }

    let tx_j = junction_positions(tx);
    for reference in refs {
        let ref_j = junction_positions(reference);
        if junctions_equal(&tx_j, &ref_j, offset)
            && (!tx_j.is_empty() || crate::interval::exonic_overlap_bp(tx, reference) > 0)
        {
            let tx_len = tx.tx_end.get().abs_diff(tx.tx_start.get());
            let ref_len = reference.tx_end.get().abs_diff(reference.tx_start.get());
            if tx_len >= ref_len {
                return format!("all_matched>=_{}", reference.name);
            }
            return format!("all_matched_<_{}", reference.name);
        }
    }

    "new_combination".to_owned()
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Figure2UtrClass {
    Extra,
    Missing,
}

fn summed_exon_len(tx: &Transcript) -> u64 {
    tx.exons.iter().map(|exon| u64::from(exon.len())).sum()
}

/// Figure 2I contains two UTR buckets, not separate 5'/3' end buckets.
/// The paper requires a difference of at least 5% of the reference's
/// summed exon length.
fn figure2_utr_class(
    tx: &Transcript,
    reference: &Transcript,
    offset: u32,
) -> Option<Figure2UtrClass> {
    let tx_j = junction_positions(tx);
    let ref_j = junction_positions(reference);
    if !junctions_equal(&tx_j, &ref_j, offset)
        || (tx_j.is_empty() && crate::interval::exonic_overlap_bp(tx, reference) == 0)
    {
        return None;
    }

    let tx_len = summed_exon_len(tx);
    let ref_len = summed_exon_len(reference);
    let difference = tx_len.abs_diff(ref_len);
    if difference == 0
        || u128::from(difference) * 100 < u128::from(ref_len) * FIGURE2_UTR_MIN_PERCENT
    {
        return None;
    }

    if tx_len > ref_len {
        Some(Figure2UtrClass::Extra)
    } else {
        Some(Figure2UtrClass::Missing)
    }
}

/// Detect the Figure 2L shape: one or more splice boundaries are replaced at
/// the same ordinal positions, without adding or removing a complete exon.
fn is_alternative_splice_site(
    tx: &Transcript,
    reference: &Transcript,
    locus_refs: &[Transcript],
    offset: u32,
) -> bool {
    if !has_new_junction(tx, locus_refs, offset) {
        return false;
    }

    let tx_j = junction_positions(tx);
    let ref_j = junction_positions(reference);
    if tx_j.is_empty() || tx_j.len() != ref_j.len() {
        return false;
    }

    let (missed, extra) = compare_ei_by_boundary(&tx_j, &ref_j, offset);
    !missed.is_empty() && missed == extra
}

/// A retained intron is a reference intron spanned by one query exon.  This
/// geometric test also covers retention of the first or last intron, which a
/// boundary-index-only rule can otherwise mistake for a terminal missing exon.
fn retains_reference_intron(tx: &Transcript, reference: &Transcript) -> bool {
    reference.exons.windows(2).any(|pair| {
        let intron_start = pair[0].end;
        let intron_end = pair[1].start;
        tx.exons
            .iter()
            .any(|exon| exon.start < intron_start && exon.end > intron_end)
    })
}

fn find_nearest_ref<'a>(tx: &Transcript, refs: &'a [Transcript], offset: u32) -> &'a Transcript {
    let tx_j = junction_positions(tx);
    let mut best_idx: usize = 0;

    let (best_missed, best_extra) =
        compare_ei_by_boundary(&tx_j, &junction_positions(&refs[0]), offset);
    let mut best_metric = (
        group_site(&best_extra).len(),
        best_extra.len(),
        group_site(&best_missed).len(),
        best_missed.len(),
    );
    let mut best_end_delta_sum = u64::from(tx.tx_start.get().abs_diff(refs[0].tx_start.get()))
        + u64::from(tx.tx_end.get().abs_diff(refs[0].tx_end.get()));

    for (idx, reference) in refs.iter().enumerate().skip(1) {
        let (missed, extra) = compare_ei_by_boundary(&tx_j, &junction_positions(reference), offset);
        let metric = (
            group_site(&extra).len(),
            extra.len(),
            group_site(&missed).len(),
            missed.len(),
        );
        let end_delta_sum = u64::from(tx.tx_start.get().abs_diff(reference.tx_start.get()))
            + u64::from(tx.tx_end.get().abs_diff(reference.tx_end.get()));

        let better = metric < best_metric
            || (metric == best_metric
                && (end_delta_sum < best_end_delta_sum
                    || (end_delta_sum == best_end_delta_sum
                        && reference.name < refs[best_idx].name)));

        if better {
            best_idx = idx;
            best_metric = metric;
            best_end_delta_sum = end_delta_sum;
        }
    }

    &refs[best_idx]
}

fn idx_minus_two_or_last(group: &[usize]) -> usize {
    if group.len() >= 2 {
        group[group.len() - 2]
    } else {
        group[group.len() - 1]
    }
}

fn desc_ei_by_boundary(
    tx: &Transcript,
    reference: &Transcript,
    offset: u32,
    end_shift_bp: u32,
) -> (Vec<String>, Vec<String>) {
    let tx_j = junction_positions(tx);
    let ref_j = junction_positions(reference);

    let (missed, extra) = compare_ei_by_boundary(&tx_j, &ref_j, offset);
    let group_miss = group_site(&missed);
    let group_extra = group_site(&extra);
    let splice_equal = group_miss.is_empty() && group_extra.is_empty();

    let ref_last_exon_no = ref_j.len() / 2 + 1;

    let mut miss_desc: Vec<String> = Vec::new();
    let mut extra_desc: Vec<String> = Vec::new();

    if group_miss.is_empty() {
        miss_desc.push("No miss exon.".to_owned());
    }

    for junctions in group_miss {
        if junctions[0] == 0 {
            miss_desc.push(format!(
                "5 primer miss: exon 1 to {}",
                junctions[junctions.len() - 1] / 2 + 1
            ));
        }

        if junctions[junctions.len() - 1] == ref_j.len().saturating_sub(1) {
            miss_desc.push(format!(
                "3 primer miss: exon {} to {}",
                junctions[0] / 2 + 2,
                ref_last_exon_no
            ));
        }

        if junctions[0] != 0 && junctions[junctions.len() - 1] != ref_j.len().saturating_sub(1) {
            let end = idx_minus_two_or_last(&junctions);
            if junctions[0] % 2 == 0 {
                miss_desc.push(format!(
                    "Intron retention: intron {} to {}",
                    junctions[0] / 2 + 1,
                    end / 2 + 1
                ));
            } else {
                miss_desc.push(format!(
                    "Exon miss: exon {} to {}",
                    junctions[0] / 2 + 2,
                    end / 2 + 2
                ));
            }
        }
    }

    if group_extra.is_empty() {
        extra_desc.push("No extra exon.".to_owned());
    }

    for junctions in group_extra {
        if junctions[0] == 0 {
            extra_desc.push(format!(
                "5 primer extra: exon 1 to {}",
                junctions[junctions.len() - 1] / 2 + 1
            ));
        }

        if junctions[junctions.len() - 1] == tx_j.len().saturating_sub(1) {
            let end = idx_minus_two_or_last(&junctions);
            extra_desc.push(format!(
                "3 primer extra: exon {} to {}",
                junctions[0] / 2 + 2,
                end / 2 + 2
            ));
        }

        if junctions[0] != 0 && junctions[junctions.len() - 1] != tx_j.len().saturating_sub(1) {
            let end = idx_minus_two_or_last(&junctions);
            if junctions[0] % 2 == 0 {
                extra_desc.push(format!(
                    "Intron extra: intron {} to {}",
                    junctions[0] / 2 + 1,
                    end / 2 + 1
                ));
            } else {
                extra_desc.push(format!(
                    "Exon extra: exon {} to {}",
                    junctions[0] / 2 + 2,
                    end / 2 + 2
                ));
            }
        }
    }

    if splice_equal {
        push_end_shift_tags(tx, reference, end_shift_bp, &mut miss_desc, &mut extra_desc);
    }

    (miss_desc, extra_desc)
}

fn push_end_shift_tags(
    tx: &Transcript,
    reference: &Transcript,
    end_shift_bp: u32,
    miss_desc: &mut Vec<String>,
    extra_desc: &mut Vec<String>,
) {
    if end_shift_bp == 0 {
        return;
    }

    let threshold = end_shift_bp as i64;

    let strand = match tx.strand {
        Strand::Unknown => reference.strand,
        strand => strand,
    };

    match strand {
        Strand::Minus => {
            let diff_5p = tx.tx_end.get() as i64 - reference.tx_end.get() as i64;
            if diff_5p >= threshold {
                extra_desc.push(format!("5 end extension: {}bp", diff_5p as u32));
            } else if diff_5p <= -threshold {
                miss_desc.push(format!("5 end truncation: {}bp", (-diff_5p) as u32));
            }

            let diff_3p = tx.tx_start.get() as i64 - reference.tx_start.get() as i64;
            if diff_3p <= -threshold {
                extra_desc.push(format!("3 end extension: {}bp", (-diff_3p) as u32));
            } else if diff_3p >= threshold {
                miss_desc.push(format!("3 end truncation: {}bp", diff_3p as u32));
            }
        }
        Strand::Plus | Strand::Unknown => {
            let diff_5p = tx.tx_start.get() as i64 - reference.tx_start.get() as i64;
            if diff_5p <= -threshold {
                extra_desc.push(format!("5 end extension: {}bp", (-diff_5p) as u32));
            } else if diff_5p >= threshold {
                miss_desc.push(format!("5 end truncation: {}bp", diff_5p as u32));
            }

            let diff_3p = tx.tx_end.get() as i64 - reference.tx_end.get() as i64;
            if diff_3p >= threshold {
                extra_desc.push(format!("3 end extension: {}bp", diff_3p as u32));
            } else if diff_3p <= -threshold {
                miss_desc.push(format!("3 end truncation: {}bp", (-diff_3p) as u32));
            }
        }
    }
}

fn desc_to_text(desc: &[String]) -> String {
    desc.join(";")
}

fn flow_desc(tx: &Transcript, refs: &[Transcript], opts: DescOpts) -> DescRow {
    let reference = find_nearest_ref(tx, refs, opts.offset_bp);
    let (miss_desc, extra_desc) =
        desc_ei_by_boundary(tx, reference, opts.offset_bp, opts.end_shift_bp);

    DescRow {
        isoform_id: tx.name.clone(),
        ref_id: reference.name.clone(),
        gene: gene_name(reference).to_owned(),
        miss: desc_to_text(&miss_desc),
        extra: desc_to_text(&extra_desc),
    }
}

fn describe_impl(isoforms: &[Transcript], references: &[Transcript], opts: DescOpts) -> DescResult {
    let mut iso_by_gene: HashMap<String, Vec<Transcript>> = HashMap::new();
    for tx in isoforms {
        for gene in genes(tx) {
            iso_by_gene
                .entry(gene.to_owned())
                .or_default()
                .push(tx.clone());
        }
    }

    let mut ref_by_gene: HashMap<String, Vec<Transcript>> = HashMap::new();
    for tx in references {
        for gene in genes(tx) {
            ref_by_gene
                .entry(gene.to_owned())
                .or_default()
                .push(tx.clone());
        }
    }

    let mut class4_rows: Vec<Class4Row> = Vec::new();
    let mut desc_rows: Vec<DescRow> = Vec::new();
    let mut fullmatch_more: HashSet<String> = HashSet::new();
    let mut fullmatch_less: HashSet<String> = HashSet::new();
    let mut alternative_splice_site: HashSet<String> = HashSet::new();
    let mut geometric_intron_retention: HashSet<String> = HashSet::new();
    let mut genes: Vec<String> = iso_by_gene.keys().cloned().collect();
    genes.sort();

    for gene in genes {
        let Some(refs) = ref_by_gene.get(&gene) else {
            continue;
        };
        if refs.is_empty() {
            continue;
        }

        let Some(isoforms) = iso_by_gene.get(&gene) else {
            continue;
        };

        for tx in isoforms {
            let locus_refs: Vec<Transcript> = refs
                .iter()
                .filter(|reference| same_comparison_locus(tx, reference))
                .cloned()
                .collect();
            if locus_refs.is_empty() {
                continue;
            }

            let nearest_ref = find_nearest_ref(tx, &locus_refs, opts.offset_bp);
            match figure2_utr_class(tx, nearest_ref, opts.offset_bp) {
                Some(Figure2UtrClass::Extra) => {
                    fullmatch_more.insert(tx.name.clone());
                }
                Some(Figure2UtrClass::Missing) => {
                    fullmatch_less.insert(tx.name.clone());
                }
                None => {}
            }
            if is_alternative_splice_site(tx, nearest_ref, &locus_refs, opts.offset_bp) {
                alternative_splice_site.insert(tx.name.clone());
            }
            if retains_reference_intron(tx, nearest_ref) {
                geometric_intron_retention.insert(tx.name.clone());
            }

            class4_rows.push(Class4Row {
                isoform_id: tx.name.clone(),
                class: class4(tx, &locus_refs, opts.offset_bp),
            });
            desc_rows.push(flow_desc(tx, &locus_refs, opts));
        }
    }

    class4_rows.sort_by(|a, b| {
        a.isoform_id
            .cmp(&b.isoform_id)
            .then_with(|| a.class.cmp(&b.class))
    });
    desc_rows.sort_by(|a, b| {
        a.isoform_id
            .cmp(&b.isoform_id)
            .then_with(|| a.ref_id.cmp(&b.ref_id))
    });
    let mut fusion_rows: Vec<FusionRow> = Vec::new();
    let mut fusion: HashSet<String> = HashSet::new();
    let mut isoforms_for_fusion: Vec<Transcript> = isoforms.to_vec();
    for tx in &mut isoforms_for_fusion {
        tx.metadata_mut().set_gene_id("none");
    }
    let fusion_annotated = add_gene(
        &isoforms_for_fusion,
        references,
        AddGeneOpts {
            fraction_read: opts.fusion_fraction_read,
            fraction_ref: opts.fusion_fraction_ref,
        },
    );
    for tx in &fusion_annotated {
        let genes: Vec<String> = gene_name(tx)
            .split("||")
            .map(str::trim)
            .filter(|g| !g.is_empty() && *g != "none")
            .map(|g| g.to_owned())
            .collect();
        if genes.len() > 1 {
            fusion.insert(tx.name.clone());
            fusion_rows.push(FusionRow {
                isoform_id: tx.name.clone(),
                genes,
            });
        }
    }

    fusion_rows.sort_by(|a, b| a.isoform_id.cmp(&b.isoform_id));

    // A fusion is Figure 2M, not a UTR change against either participating
    // gene.  Keep the legacy priority order below while making those evidence
    // sets biologically disjoint.
    fullmatch_more.retain(|isoform_id| !fusion.contains(isoform_id));
    fullmatch_less.retain(|isoform_id| !fusion.contains(isoform_id));

    let new_junction: HashSet<String> = class4_rows
        .iter()
        .filter(|row| row.class.contains("new_junction"))
        .map(|row| row.isoform_id.clone())
        .collect();

    let mut extra3: HashSet<String> = desc_rows
        .iter()
        .filter(|row| row.extra.contains("3 primer extra"))
        .map(|row| row.isoform_id.clone())
        .collect();
    let mut extra5: HashSet<String> = desc_rows
        .iter()
        .filter(|row| row.extra.contains("5 primer extra"))
        .map(|row| row.isoform_id.clone())
        .collect();
    let mut intron_retention: HashSet<String> = desc_rows
        .iter()
        .filter(|row| row.miss.contains("Intron retention"))
        .map(|row| row.isoform_id.clone())
        .collect();
    let mut exon_miss: HashSet<String> = desc_rows
        .iter()
        .filter(|row| row.miss.contains("Exon miss"))
        .map(|row| row.isoform_id.clone())
        .collect();
    let mut exon_extra: HashSet<String> = desc_rows
        .iter()
        .filter(|row| row.extra.contains("Exon extra"))
        .map(|row| row.isoform_id.clone())
        .collect();
    let mut miss3: HashSet<String> = desc_rows
        .iter()
        .filter(|row| row.miss.contains("3 primer miss"))
        .map(|row| row.isoform_id.clone())
        .collect();
    let mut miss5: HashSet<String> = desc_rows
        .iter()
        .filter(|row| row.miss.contains("5 primer miss"))
        .map(|row| row.isoform_id.clone())
        .collect();

    intron_retention.extend(geometric_intron_retention);

    // A one-for-one boundary replacement is Figure 2L.  The legacy textual
    // parity descriptions can also resemble terminal/exon/retention events;
    // remove those secondary interpretations so `new_junction` survives the
    // documented later-wins bucketization.
    for set in [
        &mut miss5,
        &mut miss3,
        &mut extra5,
        &mut extra3,
        &mut intron_retention,
        &mut exon_miss,
        &mut exon_extra,
    ] {
        set.retain(|isoform_id| !alternative_splice_site.contains(isoform_id));
    }

    let reference: HashSet<String> = isoforms
        .iter()
        .filter(|tx| is_isoform_anno(tx))
        .map(|tx| tx.name.clone())
        .collect();

    let mut isoform_ids: Vec<String> = isoforms.iter().map(|tx| tx.name.clone()).collect();
    isoform_ids.sort();
    isoform_ids.dedup();

    let mut class12_rows: Vec<Class12Row> = Vec::new();
    for isoform_id in isoform_ids {
        let mut class: Option<&str> = None;
        // Figure 2 evidence is ordered from low to high priority.  This loop
        // deliberately does not break: when several rules match, each later
        // match overwrites the earlier one (legacy TrackCluster behavior).
        for (label, set) in [
            ("new_junction", &new_junction),
            ("5'missing", &miss5),
            ("3'missing", &miss3),
            ("5'extra", &extra5),
            ("3'extra", &extra3),
            ("intron_retention", &intron_retention),
            ("inner_miss_exon", &exon_miss),
            ("inner_extra_exon", &exon_extra),
            ("fusion_gene", &fusion),
            ("full_matched<", &fullmatch_less),
            ("full_matched>=", &fullmatch_more),
            ("reference", &reference),
        ] {
            if set.contains(&isoform_id) {
                class = Some(label);
            }
        }
        let Some(class) = class else {
            continue;
        };
        class12_rows.push(Class12Row {
            isoform_id,
            class: class.to_owned(),
        });
    }

    class12_rows.sort_by(|a, b| a.isoform_id.cmp(&b.isoform_id));

    DescResult {
        desc_rows,
        class4_rows,
        fusion_rows,
        class12_rows,
    }
}

/// Describe isoforms, returning invalid option errors at the library boundary.
pub fn try_describe(
    isoforms: &[Transcript],
    references: &[Transcript],
    opts: DescOpts,
) -> Result<DescResult, crate::config::ParameterError> {
    opts.validate()?;
    Ok(describe_impl(isoforms, references, opts))
}

pub fn describe(isoforms: &[Transcript], references: &[Transcript], opts: DescOpts) -> DescResult {
    try_describe(isoforms, references, opts)
        .unwrap_or_else(|error| panic!("invalid description options: {error}"))
}

#[cfg(test)]
mod tests {
    use crate::model::{Bed12Attrs, Coord, Interval, Strand, Transcript};

    use super::*;

    fn make_tx(name: &str, exons: &[(u32, u32)], ttype: &str, gene: &str) -> Transcript {
        make_tx_strand(name, Strand::Plus, exons, ttype, gene)
    }

    fn make_tx_strand(
        name: &str,
        strand: Strand,
        exons: &[(u32, u32)],
        ttype: &str,
        gene: &str,
    ) -> Transcript {
        let tx_start = exons.iter().map(|(s, _)| *s).min().unwrap_or(0);
        let tx_end = exons.iter().map(|(_, e)| *e).max().unwrap_or(0);
        let exons = exons
            .iter()
            .map(|(s, e)| Interval::new(Coord::new(*s), Coord::new(*e)).unwrap())
            .collect::<Vec<_>>();

        Transcript::new(
            "chr1".to_owned(),
            strand,
            Coord::new(tx_start),
            Coord::new(tx_end),
            name.to_owned(),
            exons,
            Bed12Attrs {
                score: 0,
                thick_start: Coord::new(tx_start),
                thick_end: Coord::new(tx_end),
                item_rgb: "0".to_owned(),
                extra_fields: vec![
                    "none".to_owned(),
                    "none".to_owned(),
                    "none".to_owned(),
                    "-1,".to_owned(),
                    ttype.to_owned(),
                    gene.to_owned(),
                    "none".to_owned(),
                    "none".to_owned(),
                ],
            },
        )
        .unwrap()
    }

    fn figure2_class(
        strand: Strand,
        query_exons: &[(u32, u32)],
        reference_exons: &[(u32, u32)],
    ) -> Option<String> {
        let query = make_tx_strand("query", strand, query_exons, "nanopore_read", "GENE1");
        let reference = make_tx_strand(
            "reference",
            strand,
            reference_exons,
            "isoform_anno",
            "GENE1",
        );
        describe(
            &[query],
            &[reference],
            DescOpts {
                offset_bp: 0,
                ..DescOpts::default()
            },
        )
        .class12_rows
        .into_iter()
        .find(|row| row.isoform_id == "query")
        .map(|row| row.class)
    }

    #[test]
    fn description_never_matches_a_reused_gene_id_across_chromosomes() {
        let local_reference = make_tx(
            "local_ref",
            &[(100, 110), (130, 140)],
            "isoform_anno",
            "SHARED_GENE",
        );
        let mut cross_chrom_reference = make_tx(
            "cross_chrom_exact",
            &[(100, 110), (120, 140)],
            "isoform_anno",
            "SHARED_GENE",
        );
        cross_chrom_reference.chrom = "chr2".to_owned();
        let query = make_tx(
            "query",
            &[(100, 110), (120, 140)],
            "nanopore_read",
            "SHARED_GENE",
        );

        let result = describe(
            &[query],
            &[cross_chrom_reference, local_reference],
            DescOpts {
                offset_bp: 0,
                ..DescOpts::default()
            },
        );
        assert_eq!(result.desc_rows[0].ref_id, "local_ref");
        assert_eq!(result.class4_rows[0].class, "new_junction");
    }

    #[test]
    fn known_strand_query_never_matches_opposite_strand_reference() {
        let wrong_strand_exact = make_tx_strand(
            "minus_exact",
            Strand::Minus,
            &[(100, 200)],
            "isoform_anno",
            "SHARED_GENE",
        );
        let right_strand_distant = make_tx_strand(
            "plus_local",
            Strand::Plus,
            &[(300, 400)],
            "isoform_anno",
            "SHARED_GENE",
        );
        let query = make_tx_strand(
            "plus_query",
            Strand::Plus,
            &[(100, 200)],
            "nanopore_read",
            "SHARED_GENE",
        );

        let result = describe(
            &[query],
            &[wrong_strand_exact, right_strand_distant],
            DescOpts {
                offset_bp: 0,
                ..DescOpts::default()
            },
        );
        assert_eq!(result.desc_rows[0].ref_id, "plus_local");
        assert_eq!(result.class4_rows[0].class, "new_combination");
    }

    #[test]
    fn class4_detects_new_junction() {
        let refs = vec![make_tx(
            "ref",
            &[(100, 110), (120, 130), (140, 150)],
            "isoform_anno",
            "GENE1",
        )];
        let novel = make_tx(
            "novel",
            &[(100, 110), (115, 130), (140, 150)],
            "nanopore_read",
            "GENE1",
        );

        let class = class4(&novel, &refs, 0);
        assert_eq!(class, "new_junction");
    }

    #[test]
    fn fuzzy_boundary_matching_is_one_to_one() {
        // Every boundary has a nearby counterpart, but no one-to-one alignment
        // can cover all three: the first two query boundaries both depend on
        // reference boundary 100, while query boundary 111 is the only match
        // for both remaining reference boundaries.
        let query = [99, 101, 111];
        let reference = [100, 110, 112];

        assert!(!junctions_equal(&query, &reference, 2));
        let (missed, extra) = compare_ei_by_boundary(&query, &reference, 2);
        assert_eq!(missed.len(), 1);
        assert_eq!(extra.len(), 1);
    }

    #[test]
    fn fuzzy_boundary_matching_preserves_minus_strand_order() {
        let reference = make_tx_strand(
            "ref_minus",
            Strand::Minus,
            &[(100, 110), (120, 130), (140, 150)],
            "isoform_anno",
            "GENE1",
        );
        let query = make_tx_strand(
            "query_minus",
            Strand::Minus,
            &[(100, 111), (121, 131), (141, 150)],
            "nanopore_read",
            "GENE1",
        );

        assert!(!has_new_junction(
            &query,
            std::slice::from_ref(&reference),
            1
        ));
        assert!(class4(&query, &[reference], 1).starts_with("all_matched"));
    }

    #[test]
    fn flow_desc_mentions_5prime_miss_for_truncation() {
        let refs = vec![make_tx(
            "ref",
            &[(100, 110), (120, 130), (140, 150)],
            "isoform_anno",
            "GENE1",
        )];
        let trunc = make_tx("trunc", &[(120, 130), (140, 150)], "nanopore_read", "GENE1");

        let row = flow_desc(
            &trunc,
            &refs,
            DescOpts {
                offset_bp: 0,
                end_shift_bp: 0,
                ..DescOpts::default()
            },
        );
        assert!(row.miss.contains("5 primer miss"));
    }

    #[test]
    fn flow_desc_adds_end_shift_tags_on_plus_strand_when_enabled() {
        let refs = vec![make_tx(
            "ref",
            &[(100, 110), (120, 130), (140, 150)],
            "isoform_anno",
            "GENE1",
        )];
        let longer = make_tx(
            "longer",
            &[(90, 110), (120, 130), (140, 160)],
            "nanopore_read",
            "GENE1",
        );

        let row = flow_desc(
            &longer,
            &refs,
            DescOpts {
                offset_bp: 0,
                end_shift_bp: 5,
                ..DescOpts::default()
            },
        );
        assert!(row.extra.contains("5 end extension"));
        assert!(row.extra.contains("3 end extension"));
    }

    #[test]
    fn flow_desc_adds_end_shift_tags_on_minus_strand_when_enabled() {
        let refs = vec![make_tx_strand(
            "ref",
            Strand::Minus,
            &[(100, 110), (120, 130), (140, 150)],
            "isoform_anno",
            "GENE1",
        )];
        let longer = make_tx_strand(
            "longer",
            Strand::Minus,
            &[(90, 110), (120, 130), (140, 160)],
            "nanopore_read",
            "GENE1",
        );

        let row = flow_desc(
            &longer,
            &refs,
            DescOpts {
                offset_bp: 0,
                end_shift_bp: 5,
                ..DescOpts::default()
            },
        );
        assert!(row.extra.contains("5 end extension"));
        assert!(row.extra.contains("3 end extension"));
    }

    #[test]
    fn end_shift_is_disabled_when_threshold_is_zero() {
        let refs = vec![make_tx(
            "ref",
            &[(100, 110), (120, 130), (140, 150)],
            "isoform_anno",
            "GENE1",
        )];
        let longer = make_tx(
            "longer",
            &[(90, 110), (120, 130), (140, 160)],
            "nanopore_read",
            "GENE1",
        );

        let row = flow_desc(
            &longer,
            &refs,
            DescOpts {
                offset_bp: 0,
                end_shift_bp: 0,
                ..DescOpts::default()
            },
        );
        assert!(!row.extra.contains("end extension"));
        assert!(!row.miss.contains("end truncation"));
    }

    #[test]
    fn class12_matches_figure2_terminal_exon_events_on_both_strands() {
        let reference = [(100, 140), (200, 240), (300, 340), (400, 440)];

        let cases = [
            (
                Strand::Plus,
                vec![(200, 240), (300, 340), (400, 440)],
                "5'missing",
            ),
            (
                Strand::Plus,
                vec![(100, 140), (200, 240), (300, 340)],
                "3'missing",
            ),
            (
                Strand::Plus,
                vec![(20, 60), (100, 140), (200, 240), (300, 340), (400, 440)],
                "5'extra",
            ),
            (
                Strand::Plus,
                vec![(100, 140), (200, 240), (300, 340), (400, 440), (480, 520)],
                "3'extra",
            ),
            (
                Strand::Minus,
                vec![(100, 140), (200, 240), (300, 340)],
                "5'missing",
            ),
            (
                Strand::Minus,
                vec![(200, 240), (300, 340), (400, 440)],
                "3'missing",
            ),
            (
                Strand::Minus,
                vec![(100, 140), (200, 240), (300, 340), (400, 440), (480, 520)],
                "5'extra",
            ),
            (
                Strand::Minus,
                vec![(20, 60), (100, 140), (200, 240), (300, 340), (400, 440)],
                "3'extra",
            ),
        ];

        for (strand, query, expected) in cases {
            assert_eq!(
                figure2_class(strand, &query, &reference).as_deref(),
                Some(expected),
                "strand={strand:?}, query={query:?}"
            );
        }
    }

    #[test]
    fn class12_matches_figure2_internal_events_on_both_strands() {
        let reference = [(100, 140), (200, 240), (300, 340), (400, 440)];
        let reference_without_b = [(100, 140), (300, 340), (400, 440)];

        for strand in [Strand::Plus, Strand::Minus] {
            assert_eq!(
                figure2_class(strand, &[(100, 140), (300, 340), (400, 440)], &reference,)
                    .as_deref(),
                Some("inner_miss_exon")
            );
            assert_eq!(
                figure2_class(strand, &reference, &reference_without_b).as_deref(),
                Some("inner_extra_exon")
            );
            assert_eq!(
                figure2_class(strand, &[(100, 140), (200, 340), (400, 440)], &reference,)
                    .as_deref(),
                Some("intron_retention")
            );

            // Figure 2K also applies to the first and last intron.  These
            // cases used to be mistaken for terminal missing exons.
            assert_eq!(
                figure2_class(strand, &[(100, 240), (300, 340), (400, 440)], &reference,)
                    .as_deref(),
                Some("intron_retention")
            );
            assert_eq!(
                figure2_class(strand, &[(100, 140), (200, 240), (300, 440)], &reference,)
                    .as_deref(),
                Some("intron_retention")
            );

            // Figure 2L is a boundary replacement, not a whole-exon event.
            assert_eq!(
                figure2_class(
                    strand,
                    &[(100, 140), (210, 240), (300, 340), (400, 440)],
                    &reference,
                )
                .as_deref(),
                Some("new_junction")
            );
            assert_eq!(
                figure2_class(
                    strand,
                    &[(100, 140), (200, 250), (300, 340), (400, 440)],
                    &reference,
                )
                .as_deref(),
                Some("new_junction")
            );
        }
    }

    #[test]
    fn figure2_utr_classes_require_five_percent_summed_exon_difference() {
        let reference = [(100, 140), (200, 240), (300, 340), (400, 440)];

        for strand in [Strand::Plus, Strand::Minus] {
            // Reference summed exon length is 160 bp, so 8 bp is exactly 5%.
            assert_eq!(
                figure2_class(
                    strand,
                    &[(92, 140), (200, 240), (300, 340), (400, 440)],
                    &reference,
                )
                .as_deref(),
                Some("full_matched>=")
            );
            assert_eq!(
                figure2_class(
                    strand,
                    &[(108, 140), (200, 240), (300, 340), (400, 440)],
                    &reference,
                )
                .as_deref(),
                Some("full_matched<")
            );

            assert_eq!(
                figure2_class(
                    strand,
                    &[(93, 140), (200, 240), (300, 340), (400, 440)],
                    &reference,
                ),
                None
            );
            assert_eq!(figure2_class(strand, &reference, &reference), None);
        }
    }

    #[test]
    fn optional_end_shift_tags_do_not_change_figure2_class12() {
        let reference = make_tx(
            "reference",
            &[(100, 140), (200, 240), (300, 340), (400, 440)],
            "isoform_anno",
            "GENE1",
        );
        let query = make_tx(
            "query",
            &[(92, 140), (200, 240), (300, 340), (400, 440)],
            "nanopore_read",
            "GENE1",
        );

        let result = describe(
            &[query],
            &[reference],
            DescOpts {
                offset_bp: 0,
                end_shift_bp: 5,
                ..DescOpts::default()
            },
        );

        assert!(result.desc_rows[0].extra.contains("5 end extension"));
        assert_eq!(result.class12_rows[0].class, "full_matched>=");
    }

    #[test]
    fn class12_later_matching_terminal_rule_overwrites_earlier_rule() {
        let reference = [(100, 140), (200, 240), (300, 340), (400, 440)];
        let missing_both_ends = [(200, 240), (300, 340)];

        assert_eq!(
            figure2_class(Strand::Plus, &missing_both_ends, &reference).as_deref(),
            Some("3'missing")
        );
    }

    #[test]
    fn find_nearest_ref_breaks_metric_ties_by_end_delta_then_name() {
        let tx = make_tx(
            "tx",
            &[(90, 110), (120, 130), (140, 150)],
            "nanopore_read",
            "GENE1",
        );
        let ref_far = make_tx(
            "ref_far",
            &[(50, 110), (120, 130), (140, 150)],
            "isoform_anno",
            "GENE1",
        );
        let ref_near = make_tx(
            "ref_near",
            &[(100, 110), (120, 130), (140, 150)],
            "isoform_anno",
            "GENE1",
        );

        let refs = vec![ref_far.clone(), ref_near.clone()];
        let chosen = find_nearest_ref(&tx, &refs, 0);
        assert_eq!(chosen.name, "ref_near");

        let refs = vec![ref_near, ref_far];
        let chosen_name = find_nearest_ref(&tx, &refs, 0);
        assert_eq!(chosen_name.name, "ref_near");
    }

    #[test]
    fn find_nearest_ref_end_distance_covers_the_full_u32_coordinate_domain() {
        let tx = make_tx(
            "tx",
            &[(4_294_967_000, 4_294_967_100)],
            "nanopore_read",
            "GENE1",
        );
        let ref_far = make_tx("ref_far", &[(100, 200)], "isoform_anno", "GENE1");
        let ref_near = make_tx(
            "ref_near",
            &[(4_294_966_900, 4_294_967_000)],
            "isoform_anno",
            "GENE1",
        );

        let refs = [ref_far, ref_near];
        let chosen = find_nearest_ref(&tx, &refs, 0);
        assert_eq!(chosen.name, "ref_near");
    }

    #[test]
    fn describe_marks_fusion_when_isoform_overlaps_multiple_genes() {
        let refs = vec![
            make_tx("ref1", &[(100, 200)], "isoform_anno", "GENE1"),
            make_tx("ref2", &[(210, 310)], "isoform_anno", "GENE2"),
        ];
        let iso = make_tx("read1", &[(150, 260)], "nanopore_read", "none");

        let res = describe(&[iso], &refs, DescOpts::default());
        assert_eq!(res.fusion_rows.len(), 1);
        assert_eq!(res.fusion_rows[0].isoform_id, "read1");
        assert_eq!(res.fusion_rows[0].genes, vec!["GENE1", "GENE2"]);

        assert_eq!(res.class12_rows.len(), 1);
        assert_eq!(res.class12_rows[0].isoform_id, "read1");
        assert_eq!(res.class12_rows[0].class, "fusion_gene");
    }

    #[test]
    fn figure2_fusion_is_not_reclassified_as_a_utr_event() {
        let refs = vec![
            make_tx("ref1", &[(100, 200)], "isoform_anno", "GENE1"),
            make_tx("ref2", &[(210, 310)], "isoform_anno", "GENE2"),
        ];
        let iso = make_tx("read1", &[(100, 310)], "nanopore_read", "GENE1||GENE2");

        let result = describe(&[iso], &refs, DescOpts::default());
        assert_eq!(result.fusion_rows[0].isoform_id, "read1");
        assert_eq!(result.class12_rows[0].class, "fusion_gene");
    }

    #[test]
    fn annotated_reference_overwrites_other_class12_evidence() {
        let reference = make_tx(
            "reference",
            &[(100, 140), (200, 240), (300, 340)],
            "isoform_anno",
            "GENE1",
        );
        let annotated_with_new_site = make_tx(
            "annotated",
            &[(100, 140), (210, 240), (300, 340)],
            "isoform_anno",
            "GENE1",
        );

        let result = describe(
            &[annotated_with_new_site],
            &[reference],
            DescOpts {
                offset_bp: 0,
                ..DescOpts::default()
            },
        );
        assert_eq!(result.class12_rows[0].class, "reference");
    }

    #[test]
    fn describe_does_not_mark_fusion_without_overlap_even_if_gene_field_is_multi() {
        let refs = vec![
            make_tx("ref1", &[(100, 200)], "isoform_anno", "GENE1"),
            make_tx("ref2", &[(210, 310)], "isoform_anno", "GENE2"),
        ];
        let iso = make_tx("read1", &[(120, 180)], "nanopore_read", "GENE1||GENE2");

        let res = describe(&[iso], &refs, DescOpts::default());
        assert!(res.fusion_rows.is_empty());
        assert!(res
            .class12_rows
            .iter()
            .all(|row| row.class != "fusion_gene"));
    }

    #[test]
    fn rejects_invalid_fusion_fractions_at_library_boundary() {
        let opts = DescOpts {
            fusion_fraction_read: f64::INFINITY,
            ..DescOpts::default()
        };
        assert!(try_describe(&[], &[], opts).is_err());

        let opts = DescOpts {
            fusion_fraction_ref: -0.1,
            ..DescOpts::default()
        };
        assert!(try_describe(&[], &[], opts).is_err());
    }
}
