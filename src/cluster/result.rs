use crate::model::Transcript;

#[derive(Debug)]
pub struct ClusterResult {
    pub isoforms: Vec<Transcript>,
    pub read_to_isoform: Vec<(String, String)>,
    pub unused: Vec<Transcript>,
}
