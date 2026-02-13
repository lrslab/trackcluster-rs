pub const SAMPLE_DELIM: &str = "::";

pub fn tagged_read_name(sample: &str, read_name: &str) -> String {
    format!("{sample}{SAMPLE_DELIM}{read_name}")
}

pub fn split_tagged_read_name(name: &str) -> Option<(&str, &str)> {
    let (sample, read_name) = name.split_once(SAMPLE_DELIM)?;
    if sample.is_empty() || read_name.is_empty() {
        return None;
    }
    Some((sample, read_name))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tagged_name_roundtrip() {
        let tagged = tagged_read_name("S1", "readA");
        assert_eq!(tagged, "S1::readA");
        assert_eq!(split_tagged_read_name(&tagged), Some(("S1", "readA")));
    }

    #[test]
    fn split_rejects_missing_parts() {
        assert_eq!(split_tagged_read_name("S1"), None);
        assert_eq!(split_tagged_read_name("::readA"), None);
        assert_eq!(split_tagged_read_name("S1::"), None);
    }
}
