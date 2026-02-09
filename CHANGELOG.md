# Changelog

Any changes to PTCP are noted below:

## [3.2.0] - 2026-02-09

### Added

- QC CSV reporting: ptcp-qc JSON reports are now also summarized into CSV outputs written alongside existing QC outputs (qc.coverage.csv, qc.trgt.csv, qc.paraphase.csv, qc.aggregate.csv).

### Changed
- **Tool updates:** 
    - TRGT -> 5.0.0 [(Release notes)](https://github.com/PacificBiosciences/trgt/releases/tag/v5.0.0)
        - Includes improved consensus generation for long repeat expansions changes introduced in; see also changes introduced in [4.1.0](https://github.com/PacificBiosciences/trgt/releases/tag/v4.1.0)
    - Paraphase -> 3.4.0 [(Release notes)](https://github.com/PacificBiosciences/paraphase/releases/tag/v3.4.0)
    - Sawfish -> 2.2.1 [(Release notes)](https://github.com/PacificBiosciences/sawfish/releases/tag/v2.2.1)
    - ptcp-qc -> 1.0.0 [(Release notes)](https://github.com/PacificBiosciences/ptcp-qc/releases/tag/1.0.0)
        - Adds improved annotation for known HBA deletions, including the most likely deletion name when detected (e.g., if the observed breakpoints and deletion length match the 3.7 kb deletion (`3p7del`), it will be annotated as `3p7del`)

- Variant annotation with Havanno is now deterministic to improve reproducibility across runs
- SMN homology correction supports batch mode for higher-throughput processing
- Repeat annotation VCF updated with additional HBB variants
- Sample sheet format simplified to a two-column mapping (sample to sex)
- TRGT TRID extraction for plotting now uses bcftools instead of awk parsing
- BAM collection improved across storage backends to reduce low-probability sample mis-association risk

### Fixed

- *SMN* copy number correction: fixed an edge case where the correction logic could report a single copy
- miniwdl-slurm updated to address rare Slurm-related failures observed in production environments
- *GBA* and *RPGR* target intervals adjusted to prevent 1 to 2 bp edge clipping artifacts during Paraphase realignment
- Fixed a rare naming truncation issue during BAM normalization. Outputs were not affected, and sample names/labels are now preserved correctly throughout normalization and reporting

## [3.1.2] - 2025-10-03

### Added
- Add paraphase_annotation_vcf link to the template JSON

### Fixed
- Update to the current PTCP quay.io URL in the main.wdl and within documentation
- Update dead link in Havanno documentation

## [3.1.1] - 2025-09-17

### Added
- Extensive updates to documentation relating to input and output files, FAQ, and running PTCP on HPC systems

## [3.1.0] - 2025-09-05

### Changed
- **Tool updates:** 
    - Paraphase -> 3.3.4 [(Release notes)](https://github.com/PacificBiosciences/paraphase/releases/tag/v3.3.4)
        - Haplotype assignments are now deterministic
    - ptcp-qc -> 0.8.0
        - Adds a new `genotype_adjusted` field to the `paraphase_results` of **HBA**, this will adjust the Paraphase genotyping if large deletions are detected by Sawfish in **HBA**
    - F8 calling script: renamed the output "breakpoints" to "spans"
- **Tool configurations:**
    - Updated Paraphase gene configurations for **GBA**, **f8inv22**, and **f8inv1**
    - Updated the variant list VCF

### Fixed
- Update requirements to make PTCP python 3.13 compatible

___

## [3.0.0] - 2025-07-31

### Added
 - Sawfish for calling structural variants, initially focused on the HBA locus
 - The F8 inversion script now reports breakpoint regions
 - The documentation now includes a section on running the pipeline on DNAnexus

### Changed
- **Tool updates:** 
    - TRGT -> 4.0.0 [(Release notes)](https://github.com/PacificBiosciences/trgt/releases/tag/v4.0.0)
    - Paraphase -> 3.3.2 [(Release notes)](https://github.com/PacificBiosciences/paraphase/releases/tag/v3.3.2)
    - ptcp-qc -> 0.7.0
    - Removed the -N 1 flag from the pbmm2 alignment step

___

## [2.2.1] - 2025-06-24

### Added
 - **Reference genome support:** 
    - Variant list VCF for hg37 (beyond just hg38)

### Fixed
 - **Workflow:** 
    - Corrected the `havanno_json` definition in the Paraphase task which could otherwise lead to errors
 - **Tool updates:** 
    - Correctly annotate overlapping pathogenic variants during small-variant annotation

___

## [2.2.0] - 2025-06-10

### Added
- **Reference‐genome support:**
    - Support for hg37 (in addition to hg38), corresponding configuration files under `meta/{paraphase, trgt, ptcp-qc}` for hg37 and hg38
- **Workflow:**
    - Per-haplotype small-variation annotation of Paraphase VCF output ([see: Havanno readme](docker/ptcp/scripts/havanno/README.md)); optionally triggered if an annotation VCF is provided. An example variant list VCF is provided in `meta/variant_list`
    - Annotation step for Paraphase SMN1/2 output to scale copy number based on coverage at other targets ([see: SMN_homology readme](docker/ptcp/scripts/smn-homology/README.md)); optionally triggered if a regression-weights JSON is provided
- **Utilities:** Scripts to generate sample sheets and PTCP workflow input JSON files

### Changed
- **Tool updates:**
    - ptcp-qc -> v0.6.0
        - Expand JSON output to include:
            - Genome version, target-BED name, ptcp-qc version, timestamps
            - Per-sample X∕Autosomal and Y∕Autosomal ratio calculations (off-target reads)
            - Paraphase reporting now includes absolute read counts, unique read counts, and fractional counts for reads mapping to multiple haplotypes
- **Tool configurations:**
    - Update Paraphase gene definition for **HBB** to prevent haplotype truncation
___

## [2.1.0] - 2025-04-30

### Changed
- **Tool updates:**
    - TRGT -> 3.0.0 [(Release notes)](https://github.com/PacificBiosciences/trgt/releases/tag/v3.0.0)
        - Modifies repeat-motif detection: segmentation now matches only perfect STR motifs, while still tolerating imperfections in VNTR motifs
        - Introduces flexible motif matching for improved detection and visualization of repeat interruptions
        - Includes all upstream improvements from TRGT [2.1.0](https://github.com/PacificBiosciences/trgt/releases/tag/v2.1.0) for faster genotyping and plotting
    - Paraphase → v3.3.0 (**Full release**) [(Release notes)](https://github.com/PacificBiosciences/paraphase/releases/tag/v3.3.0)
        - Release consolidating beta fixes and targeted analysis improvements
        - Also incorporates the latest[ 3.3.1 feature-branch](https://github.com/PacificBiosciences/paraphase/tree/v3.3.1) updates
    - F8 calling module: refactored and streamlined logic for clarity, added comprehensive unit tests
- **Tool configurations:**
    - Updated Paraphase gene definitions for **F8**, **GBA**, and added the new `gene1_cn2` flag to **GBA** to adjust copy-number calculation when only one haplotype is detected

___

## [2.0.0] - 2025-04-10

### Added
- **Breaking changes:**
    - Renamed workflow to **PureTarget Carrier Pipeline (PTCP)** 
    - Renamed QC tool from **manatee-qc** to **ptcp-qc**

### Changed
- **Tool updates:**
    - TRGT -> 2.0.0 [(Release notes)](https://github.com/PacificBiosciences/trgt/releases/tag/v2.0.0)
        - Up to 200x genotyping speed-up in targeted datasets
    - Upgraded QC tool (now ptcp-qc) from manatee-qc v0.4.0 -> ptcp-qc v0.5.0
- **Tool configurations:**
    - Updated Paraphase gene definitions for **CYP21** and **GBA**
    - Updated TRGT tandem repeat BED file
    - Updated ptcp-qc target BED file

___

## [1.6.0] - 2025-04-02

### Added
- **manatee-qc** task which analyzes generated outputs, producing both per-sample and aggregate QC reports
- **tar_outputs** task which consolidates upstream per sample output files into a single .tar archive for input into the manatee-qc task
- Configuration file for manatee-qc at `meta/manatee-qc/manatee_qc.GRCh38.bed`
- manatee-qc binary of version v0.4.0 in the docker/manatee/binaries directory

### Changed
- **Tool updates:**
    - Paraphase -> latest commit on 3.3.0 branch

### Removed
- Redundant Jasmine task from the workflow.

___

## [1.5.0] - 2025-03-14

### Changed
- **Tool updates:**
    - Paraphase -> 3.3.0 [(Release notes)](https://github.com/PacificBiosciences/paraphase/tree/v3.3)
        - Now reports the breakpoints of structural variants detected in HBA within the JSON output
        - Improved handling of coverage variability:
            - Parameters now accept frequencies instead of absolute read counts for setting minimal support thresholds
                - Renamed `--min-read-haplotype` → `--min-haplotype-frequency`
                - Renamed `--min-read-variant` → `--min-variant-frequency`
- **Tool configurations:**
    - Updated Paraphase gene definitions for SMN1: align with new guide definitions, HBA: now accounts for the shorter homology haplotypes that were reported as unknown haplotypes

### Removed
- Sawfish has been removed from the workflow as Paraphase now reports HBA structural variant breakpoints

___

## [1.4.0] - 2025-02-27

### Changed
- **Tool updates:**
    - Paraphase -> 3.2.1 [(Release notes)](https://github.com/PacificBiosciences/paraphase/releases/tag/v3.2.1)
        - Adds the RN (region name) tag to output BAM files, allowing any read to be matched to a specific analyzed region
    - TRGT -> 1.5.1 [(Release notes)](https://github.com/PacificBiosciences/trgt/releases/tag/v1.5.1)
    - Sawfish -> 0.12.10 [(Release notes)](https://github.com/PacificBiosciences/sawfish/releases/tag/v0.12.10)
    - F8 inversion script updates
        - Simplified target configuration
        - Now reports unclassified reads
        - Standardized naming to align with Paraphase F8 targets
- **Tool configurations:**
    - Updated the SMN1 gene definition to align with new guides

___

## [1.3.0] - 2025-02-07

### Changed
- **Tool updates:**
    - Smrtlink -> v25.1.0.257715
    - Paraphase -> 3.2.0 [(Release notes)](https://github.com/PacificBiosciences/paraphase/releases/tag/v3.2.0)
        - Enabled new flag `--write-nocalls-in-vcf`, Paraphase will now write no-call sites in the VCFs, marked with LowQual filter
    - Sawfish -> 0.12.9 [(Release notes)](https://github.com/PacificBiosciences/sawfish/releases/tag/v0.12.9)
- **Tool configurations:**
    - Update HBA gene definition

___

## [1.2.0] - 2025-01-23

### Added
- **Breaking changes**: Update BAM file collection
    - Sample matching is now based on the filename (`movie.barcode`) instead of sample ids from BAM header
        - Example: `movie.hifi_reads.barcode.bam` and `movie.fail_reads.barcode.bam` will now be merged based on their filenames
        - Assumes a single HiFi BAM and at most one fail reads BAM for each sample
        - Changed the sample sheet formatting from `movie,barcode,sample_id,sex` to `bam_name,bam_id,sex`, where `bam_id` can be derived from `bam_name` by removing `.hifi_reads.`, `.fail_reads.`, or `.reads`
- Introduced basic BAM validation using samtools quickcheck to detect malformed headers or truncated files
- Added metadata to all workflows and tasks. Defined disk space requirements in the runtime of each task

### Changed
- **Tool updates:**
    - TRGT -> 1.5.0 [(Release notes)](https://github.com/PacificBiosciences/trgt/releases/tag/v1.5.0)

### Removed
- `infer_chry` task
- `pbmm2_index` task (handled automatically by pbmm2, which builds an index if required)

___

## [1.1.0] - 2025-01-23

### Added
- Sawfish subworkflow for breakpoint detection of structural variations (e.g., HBA deletions)

### Changed
- **Breaking change**: Revert the change to collect_inputs (can specify HiFi and fail reads again in the inputs.json)
- **Tool updates:**
    - TRGT -> 1.4.1 [(Release notes)](https://github.com/PacificBiosciences/trgt/releases/tag/v1.4.1)

___

## [1.0.0] - 2025-01-23

### Added
- **Breaking changes**: Major restructuring of the workflow, key updates:
    - Replace the PBAA-based subworkflow with a Paraphase-based subworkflow for hard gene analysis
    - Minimize alignments, all subworkflows now use the same PBMM2 aligned BAM file

### Changed
- **Tool updates:**
- TRGT -> 1.4.0 [(Release notes)](https://github.com/PacificBiosciences/trgt/releases/tag/v1.4.0)
