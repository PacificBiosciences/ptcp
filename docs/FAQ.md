# Frequently asked questions

## How does PTCP handle PureTarget-specific configurations?

PTCP is specifically designed for PureTarget data and automatically configures all tools with the optimal settings. 

**TRGT configuration**: PTCP automatically:
- Uses the `--preset targeted` option for TRGT analysis
- Includes fail reads when available to improve coverage of expanded alleles
- Disables quality filtering that might exclude important repeat expansion reads
- Uses cluster-based genotyping for better allele assignment

**Paraphase configuration**: PTCP automatically configures Paraphase with settings optimized for targeted sequencing data, including:
- `--targeted`: Enables targeted sequencing mode for improved analysis of specific genomic regions
- `--write-nocalls-in-vcf`: Paraphase will write no-call sites in the VCFs, marked with LowQual filter

## Why are fail reads important for analyzing tandem repeats with PureTarget data?

The PureTarget protocol produces insert sizes of about 5 kb, but large expansions of loci like *FXN*, *C9orf72*, *DMPK* and *CNBP* produce much larger molecules that may not produce reads reaching HiFi quality thresholds at typical movie times. PTCP automatically includes fail reads when available, which can significantly increase the coverage of expanded alleles and prevent allelic dropouts.

## Where can I learn more about the underlying tool configurations?

For detailed information about how TRGT analyzes PureTarget data, see the [TRGT PureTarget documentation](https://github.com/PacificBiosciences/trgt/blob/main/docs/puretarget.md).

For information about Paraphase configuration with targeted sequencing data, see the [Paraphase targeted data documentation](https://github.com/PacificBiosciences/paraphase/blob/main/docs/targeted_data.md).

## Can I use PTCP to analyze PureTarget repeat expansion panel 2.0 data?

Yes, because TRGT is part of PTCP, and the regions file contains Carrier Screening and repeat expansion panel 2.0 regions, you can use PTCP to process repeat expansion panel 2.0 data.

## Can I run PTCP with the GRCh37 reference?

Yes. Be sure to use the correct configuration file for the corresponding reference genome. Refer to the [Reference genome section](./Input_files.md#13-reference-genome) of the Input files page for more details.