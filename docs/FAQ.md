# Frequently asked questions

## How does PTCP handle PureTarget-specific configurations?

PTCP is specifically designed for PureTarget data and automatically configures all tools with the optimal settings. 

**TRGT configuration**: PTCP automatically:
- Uses the `--preset targeted` option for TRGT analysis
- Includes fail reads when provided to improve coverage of expanded alleles
- Disables quality filtering that might exclude important repeat expansion reads
- Uses cluster-based genotyping for better allele assignment

**Paraphase configuration**: PTCP automatically configures Paraphase with settings optimized for targeted sequencing data, including:
- `--targeted`: Enables targeted sequencing mode for improved analysis of specific genomic regions
- `--write-nocalls-in-vcf`: Paraphase will write no-call sites in the VCFs, marked with LowQual filter

## Why is *F8* calling handled outside of Paraphase?

Paraphase has explicit support for *F8* in WGS, but in PureTarget we only use Paraphase to remap reads onto a single target for the *F8* introns. The actual *F8* inversion calling is performed by a method using clipping patterns tied to the exact PureTarget guide RNA design, which is too PureTarget-specific to live inside Paraphase. Note that we cannot do meaningful small variant calling for *F8* here either, as Paraphase is only remapping reads to a single *F8* target and the haplotypes it produces are not informative for small variants, so annotation is out of scope for *F8*.

## Why is there an *HBA*-specific workflow if Paraphase supports *HBA1/2*?

While Paraphase can classify common *HBA1/2* events (such as 3p7/4p2 del/dup) using consistent soft-clip start/end positions on haplotypes and label them as *HBA1*, *HBA2*, or the hybrid patterns. It does not scan arbitrary split reads or discover novel breakpoints, so larger or non-canonical structural variants need a dedicated SV caller like Sawfish that uses split/supplementary alignments and depth signals. That functionality is provided via the additional *HBA* workflow.

## Why are fail reads important for analyzing tandem repeats with PureTarget data?

The PureTarget protocol produces insert sizes of about 5 kb, but large expansions of loci like *FXN*, *C9orf72*, *DMPK* and *CNBP* produce much larger molecules that may not produce reads reaching HiFi quality thresholds at typical movie times. Including fail reads (when available) can significantly increase the coverage of expanded alleles and prevent allelic dropouts.

## Where can I learn more about the underlying tool configurations?

For detailed information about how TRGT analyzes PureTarget data, see the [TRGT PureTarget documentation](https://github.com/PacificBiosciences/trgt/blob/main/docs/puretarget.md).

For information about Paraphase configuration with targeted sequencing data, see the [Paraphase targeted data documentation](https://github.com/PacificBiosciences/paraphase/blob/main/docs/targeted_data.md).

## Should I validate ptcp-qc CSV/JSON calls against the underlying outputs?

Yes. The ptcp-qc JSON/CSV reports are summaries intended for QC triage and reporting, and they should be treated as pointers to supporting evidence. For any surprising or relevant result, validate it in the primary tool outputs (for example: `trgt_vcf` + TRGT plots + `trgt_spanning_bam` for tandem repeat calls; `paraphase_json`/BAM/VCFs for Paraphase; and `sawfish_vcf` for larger SVs such as HBA deletions).

## Can I use PTCP to analyze PureTarget repeat expansion panel 2.0 data?

Yes, because TRGT is part of PTCP, and the regions file contains Carrier Screening and repeat expansion panel 2.0 regions, you can use PTCP to process repeat expansion panel 2.0 data.

## Can I run PTCP with the GRCh37 reference?

Yes. Be sure to use the correct configuration file for the corresponding reference genome. Refer to the [Reference genome section](./Input_files.md#13-reference-genome) of the Input files page for more details.
