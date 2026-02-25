# Frequently Asked Questions

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

Paraphase has explicit support for *F8* in WGS, but in PureTarget we only use Paraphase to remap reads onto a single target for the *F8* introns. The actual *F8* inversion calling is performed by a method that uses clipping patterns tied to the exact PureTarget guide RNA design, which is too PureTarget-specific to live inside Paraphase. Note that we cannot do meaningful small variant calling for *F8* here either, as Paraphase is only remapping reads to a single *F8* target and the haplotypes it produces are not informative for small variants, so annotation is out of scope for *F8*.

## Why is there an *HBA*-specific workflow if Paraphase supports *HBA1/2*?

Paraphase can classify common *HBA1/2* events (for example 3p7/4p2 del/dup) because those variants tend to produce consistent haplotype-level patterns, such as reproducible soft-clip start/end positions and characteristic hybrid configurations that can be labeled as *HBA1*, *HBA2*, or a known hybrid. Paraphase is not designed to scan arbitrary split reads or to discover novel breakpoints, so larger or non-canonical structural variants are better handled by a dedicated SV caller such as Sawfish, which uses split and supplementary alignments plus depth signals. The additional *HBA* workflow provides that SV-calling capability.

To make split-read based SV calling work reliably in the *HBA* region, we remap the locus reads before running Sawfish. In a highly homologous locus like *HBA1/2*, the aligner's choices about when to emit a supplementary alignment (versus soft-clipping), how to place alignment breakpoints, and how to assign primary/secondary/supplementary status can vary depending on parameters. Sawfish depends on those supplementary and split alignments as breakpoint evidence. If reads are represented as soft-clipped or the split placement shifts, the caller can lose support or move the breakpoint. We evaluated using pbmm2 for this remapping, but we could not find parameters that produced consistently stable split and supplementary alignments across our *HBA1/2* test cases. Minimap2 produced more consistent supplementary alignment behavior in this context, so it is used for the remapping step. We may revisit pbmm2 in a future release.

## Why are fail reads important for analyzing tandem repeats with PureTarget data?

The PureTarget protocol produces insert sizes of about 5 kb, but large expansions of loci like *FXN*, *C9orf72*, *DMPK* and *CNBP* produce much larger molecules that may not produce reads reaching HiFi quality thresholds at typical movie times. Including fail reads (when available) can significantly increase the coverage of expanded alleles and prevent allelic dropouts.

## Where can I learn more about the underlying tool configurations?

For detailed information about how TRGT analyzes PureTarget data, see the [TRGT PureTarget documentation](https://github.com/PacificBiosciences/trgt/blob/main/docs/puretarget.md).

For information about Paraphase configuration with targeted sequencing data, see the [Paraphase targeted data documentation](https://github.com/PacificBiosciences/paraphase/blob/main/docs/targeted_data.md).

## Should I validate ptcp-qc CSV/JSON calls against the underlying outputs?

Yes. The ptcp-qc JSON/CSV reports are summaries intended for QC triage and reporting, and they should be treated as pointers to supporting evidence. For any surprising or relevant result, validate it in the primary tool outputs (for example: `trgt_vcf` + TRGT plots + `trgt_spanning_bam` for tandem repeat calls; `paraphase_json`/BAM/VCFs for Paraphase; and `sawfish_vcf` for larger SVs such as *HBA* deletions).

## Can I use PTCP to analyze PureTarget repeat expansion panel 2.0 data?

Yes, because TRGT is part of PTCP, and the regions file contains Carrier Screening and repeat expansion panel 2.0 regions, you can use PTCP to process repeat expansion panel 2.0 data.

## Can I run PTCP with the GRCh37 reference?

Yes. Be sure to use the correct configuration file for the corresponding reference genome. Refer to the [Reference genome section](./Input_files.md#13-reference-genome) of the Input files page for more details.

## How is sample sex used in PTCP?

PTCP does not infer sex from the data. It uses the `sex` column in your sample sheet exactly as provided. Downstream logic treats **only** the exact value `M` as male; every other value (including `F`, `f`, `Male`, or missing) is treated as female. For **autosomal** targets, `sex` has no effect.

- **TRGT**: PTCP maps `sex` to a TRGT karyotype and passes it via `--karyotype` (`XY` if `sex == M`, otherwise `XX`). This sets the expected chrX ploidy: `XY` samples are modeled as haploid on X (one allele), while `XX` samples are modeled as diploid on X (two alleles).
- **Paraphase and *F8* inversion calling**: PTCP passes the same `sex` value into the Paraphase subworkflow and into the *F8* inversion caller (`f8_inversion.py --sex`). More broadly, correct sex is required to interpret X-linked Paraphase outputs because copy-number expectations differ between XX and XY.

If `sex` is missing or incorrect, chrX calls can look wrong, for example apparent two-allele behavior on X in an XY sample, or unexpected copy numbers. PTCP currently models only XX vs XY; atypical karyotypes (for example XXY) and mosaicism are not handled and should be interpreted with extra care.

## Why is coverage unexpectedly low, or mapped read counts look wrong?

If a locus (for example *CYP21A2*) shows unexpectedly low coverage or almost no reads in PTCP's `*.mapped.bam`, a likely explanation is a **reference genome mismatch**, where the FASTA you used for mapping includes ALT contigs and/or decoys. In an ALT/decoy-rich GRCh38, reads from duplicated sequence, close paralogs, or segmental duplications can align preferentially to ALT/decoy contigs instead of the primary chromosome. The reads are still in the BAM, but they are distributed across different reference sequences, so depth at the primary locus looks artificially low.

A quick check is to count the number of reference sequences in the BAM header:

```bash
samtools view -H sample.mapped.bam | awk '$1=="@SQ"{c++} END{print c}'
```

For the PTCP-supported GRCh38 "no ALT" reference genome, you should see 195 contigs, and PTCP supports only the reference genomes listed in our reference genome repository (see the [Reference genome section](./Input_files.md#13-reference-genome)). If the contig count in your BAM header does not match the supported reference, re-run PTCP with the supported FASTA and its matching index files. If you want PTCP and another workflow/tool to produce comparable locus-level mapping counts, make sure both use the exact same reference build and contig set, including whether ALT/decoy contigs are present. This also applies to hg19.
