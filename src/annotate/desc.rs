use std::collections::{HashMap, HashSet};

use crate::annotate::addgene::{add_gene, AddGeneOpts};
use crate::model::{Strand, Transcript};

const TTYPE_COL: usize = 4;
const GENE_NAME_COL: usize = 5;

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
    tx.extra_fields.get(TTYPE_COL).map(|value| value.as_str())
}

fn is_isoform_anno(tx: &Transcript) -> bool {
    matches!(ttype(tx), Some("isoform_anno"))
}

fn gene_name(tx: &Transcript) -> &str {
    tx.extra_fields
        .get(GENE_NAME_COL)
        .map(|value| value.as_str())
        .unwrap_or("none")
}

fn set_extra(tx: &mut Transcript, idx: usize, value: String) {
    if tx.extra_fields.len() <= idx {
        tx.extra_fields.resize(idx + 1, "none".to_owned());
    }
    tx.extra_fields[idx] = value;
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

fn junctions_equal(a: &[u32], b: &[u32], offset: u32) -> bool {
    let mut matched_a: HashSet<u32> = HashSet::new();
    let mut matched_b: HashSet<u32> = HashSet::new();

    for &i in a {
        for &j in b {
            if i.abs_diff(j) <= offset {
                matched_a.insert(i);
                matched_b.insert(j);
            }
        }
    }

    a.iter().copied().collect::<HashSet<u32>>() == matched_a
        && b.iter().copied().collect::<HashSet<u32>>() == matched_b
}

fn fuzzy_intersection(a: &[u32], b: &[u32], offset: u32) -> HashMap<u32, u32> {
    let mut match_dic: HashMap<u32, u32> = HashMap::new();
    for &i in a {
        for &j in b {
            if i.abs_diff(j) <= offset {
                match_dic.insert(i, j);
            }
        }
    }
    match_dic
}

fn compare_ei_by_boundary(a: &[u32], reference: &[u32], offset: u32) -> (Vec<usize>, Vec<usize>) {
    let match_dic = fuzzy_intersection(a, reference, offset);
    let junction_new: Vec<u32> = a
        .iter()
        .copied()
        .map(|pos| match_dic.get(&pos).copied().unwrap_or(pos))
        .collect();

    let mut posdic_a: HashMap<u32, usize> = HashMap::new();
    for (idx, pos) in junction_new.iter().copied().enumerate() {
        posdic_a.insert(pos, idx);
    }

    let mut posdic_ref: HashMap<u32, usize> = HashMap::new();
    for (idx, pos) in reference.iter().copied().enumerate() {
        posdic_ref.insert(pos, idx);
    }

    let reference_set: HashSet<u32> = reference.iter().copied().collect();
    let junction_set: HashSet<u32> = junction_new.iter().copied().collect();

    let mut missed_order: Vec<usize> = reference_set
        .difference(&junction_set)
        .filter_map(|pos| posdic_ref.get(pos).copied())
        .collect();
    missed_order.sort_unstable();

    let mut extra_order: Vec<usize> = junction_set
        .difference(&reference_set)
        .filter_map(|pos| posdic_a.get(pos).copied())
        .collect();
    extra_order.sort_unstable();

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

    let match_dic = fuzzy_intersection(&junctions, &ref_junctions, offset);
    let matched: HashSet<u32> = match_dic.keys().copied().collect();
    junctions.iter().any(|pos| !matched.contains(pos))
}

fn class4(tx: &Transcript, refs: &[Transcript], offset: u32) -> String {
    if has_new_junction(tx, refs, offset) {
        return "new_junction".to_owned();
    }

    let tx_j = junction_positions(tx);
    for reference in refs {
        let ref_j = junction_positions(reference);
        if junctions_equal(&tx_j, &ref_j, offset) {
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
    let mut best_end_delta_sum = tx.tx_start.get().abs_diff(refs[0].tx_start.get())
        + tx.tx_end.get().abs_diff(refs[0].tx_end.get());

    for (idx, reference) in refs.iter().enumerate().skip(1) {
        let (missed, extra) = compare_ei_by_boundary(&tx_j, &junction_positions(reference), offset);
        let metric = (
            group_site(&extra).len(),
            extra.len(),
            group_site(&missed).len(),
            missed.len(),
        );
        let end_delta_sum = tx.tx_start.get().abs_diff(reference.tx_start.get())
            + tx.tx_end.get().abs_diff(reference.tx_end.get());

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

pub fn describe(isoforms: &[Transcript], references: &[Transcript], opts: DescOpts) -> DescResult {
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
            class4_rows.push(Class4Row {
                isoform_id: tx.name.clone(),
                class: class4(tx, refs, opts.offset_bp),
            });
            desc_rows.push(flow_desc(tx, refs, opts));
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
        set_extra(tx, GENE_NAME_COL, "none".to_owned());
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

    let end_ext3: HashSet<String> = desc_rows
        .iter()
        .filter(|row| row.extra.contains("3 end extension"))
        .map(|row| row.isoform_id.clone())
        .collect();
    let end_ext5: HashSet<String> = desc_rows
        .iter()
        .filter(|row| row.extra.contains("5 end extension"))
        .map(|row| row.isoform_id.clone())
        .collect();
    let end_trunc3: HashSet<String> = desc_rows
        .iter()
        .filter(|row| row.miss.contains("3 end truncation"))
        .map(|row| row.isoform_id.clone())
        .collect();
    let end_trunc5: HashSet<String> = desc_rows
        .iter()
        .filter(|row| row.miss.contains("5 end truncation"))
        .map(|row| row.isoform_id.clone())
        .collect();

    let mut end_shift_any: HashSet<String> = HashSet::new();
    end_shift_any.extend(end_ext3.iter().cloned());
    end_shift_any.extend(end_ext5.iter().cloned());
    end_shift_any.extend(end_trunc3.iter().cloned());
    end_shift_any.extend(end_trunc5.iter().cloned());

    let fullmatch_more_all: HashSet<String> = class4_rows
        .iter()
        .filter(|row| row.class.contains("all_matched>"))
        .map(|row| row.isoform_id.clone())
        .collect();
    let fullmatch_less_all: HashSet<String> = class4_rows
        .iter()
        .filter(|row| row.class.contains("all_matched_<"))
        .map(|row| row.isoform_id.clone())
        .collect();

    let fullmatch_more: HashSet<String> = fullmatch_more_all
        .difference(&end_shift_any)
        .cloned()
        .collect();
    let fullmatch_less: HashSet<String> = fullmatch_less_all
        .difference(&end_shift_any)
        .cloned()
        .collect();

    let new_junction: HashSet<String> = class4_rows
        .iter()
        .filter(|row| row.class.contains("new_junction"))
        .map(|row| row.isoform_id.clone())
        .collect();

    let extra3: HashSet<String> = desc_rows
        .iter()
        .filter(|row| row.extra.contains("3 primer extra"))
        .map(|row| row.isoform_id.clone())
        .collect();
    let extra5: HashSet<String> = desc_rows
        .iter()
        .filter(|row| row.extra.contains("5 primer extra"))
        .map(|row| row.isoform_id.clone())
        .collect();
    let intron_retention: HashSet<String> = desc_rows
        .iter()
        .filter(|row| row.miss.contains("Intron retention"))
        .map(|row| row.isoform_id.clone())
        .collect();
    let exon_miss: HashSet<String> = desc_rows
        .iter()
        .filter(|row| row.miss.contains("Exon miss"))
        .map(|row| row.isoform_id.clone())
        .collect();
    let exon_extra: HashSet<String> = desc_rows
        .iter()
        .filter(|row| row.extra.contains("Exon extra"))
        .map(|row| row.isoform_id.clone())
        .collect();
    let miss3: HashSet<String> = desc_rows
        .iter()
        .filter(|row| row.miss.contains("3 primer miss"))
        .map(|row| row.isoform_id.clone())
        .collect();
    let miss5: HashSet<String> = desc_rows
        .iter()
        .filter(|row| row.miss.contains("5 primer miss"))
        .map(|row| row.isoform_id.clone())
        .collect();

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
        for (label, set) in [
            ("new_junction", &new_junction),
            ("5'missing", &miss5),
            ("3'missing", &miss3),
            ("5'extra", &extra5),
            ("3'extra", &extra3),
            ("intron_retention", &intron_retention),
            ("inner_miss_exon", &exon_miss),
            ("inner_extra_exon", &exon_extra),
            ("5'end_truncation", &end_trunc5),
            ("3'end_truncation", &end_trunc3),
            ("5'end_extension", &end_ext5),
            ("3'end_extension", &end_ext3),
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
}
