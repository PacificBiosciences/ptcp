# Reference data information

## TRE resources

### STRchive-disease-loci.hg38.TRGT.strchive-\<strchive-version\>.ptcp-\<ptcp-version\>-FXS_FMR1_AGG.bed

A TRGT BED for use as `trgt_bed`. It is STRchive's TRGT BED as published, plus a second FMR1 row using PTCP's definition from `meta/trgt`:

```
chrX	147912049	147912111	ID=FXS_FMR1;MOTIFS=CGG;STRUC=<TR>
chrX	147912049	147912111	ID=FXS_FMR1_AGG;MOTIFS=CGG,AGG;STRUC=<TR>
```

STRchive models FMR1 with the single motif `CGG`, so AGG interruptions are counted as part of the CGG repeat. PTCP's `CGG,AGG` definition reports the AGG count per allele in `MC`. Both rows cover the same span, so FMR1 is called with and without interruptions. The PTCP row is renamed `FXS_FMR1_AGG` because TRIDs must be unique.

```bash
STRCHIVE_VERSION=v2.29.0
PTCP_VERSION=3.2.1
STRCHIVE_BED=STRchive-disease-loci.hg38.TRGT.bed
PTCP_BED=meta/trgt/PureTarget_repeat_expansion_panel_2.0.repeat_definition.GRCh38.bed
OUTPUT_BED=STRchive-disease-loci.hg38.TRGT.strchive-${STRCHIVE_VERSION}.ptcp-${PTCP_VERSION}-FXS_FMR1_AGG.bed

wget "https://raw.githubusercontent.com/dashnowlab/STRchive/${STRCHIVE_VERSION}/data/catalogs/${STRCHIVE_BED}"

{
  cat "${STRCHIVE_BED}"
  grep 'ID=FXS_FMR1;' "${PTCP_BED}" \
  | sed 's~ID=FXS_FMR1;~ID=FXS_FMR1_AGG;~'
} \
| sort -k1,1V -k2,2n \
> "${OUTPUT_BED}"

# Confirm all TRIDs are unique; should be no output
grep -oE 'ID=[^;]*;' "${OUTPUT_BED}" | sort | uniq -d
```

With `v2.29.0` the result has 83 loci: STRchive's 82 and `FXS_FMR1_AGG`. STRchive releases are listed at <https://github.com/dashnowlab/STRchive/releases>.
