version 1.0

struct RegionParams {
  Int?   filter
  Int?   trim_ends
  Int?   pile_size
  Float? min_var_frequency
  Int?   max_alignments_per_read

  Int? max_reads_per_guide
  Int? iterations
  Int? seed

  Int? max_consensus_reads

  Int?   max_amplicon_size
  Float? min_read_qv
  String? off_target_groups
  Float? min_cluster_frequency
  Int?   min_cluster_read_count
  Float? max_uchime_score

  Boolean? skip_chimera_detection
  Boolean? skip_consensus
}

struct PbaaParamsFile {
  String version

  Map[String, RegionParams] regions
}