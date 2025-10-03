# Havanno
HAplotype Vriant ANNOtation is a program that annotates paraphase output for
variants of interest. It produces a JSON reporting, for each haplotype

1. the variants found in a query VCF, stratified by haplotype and pathogenic and
   pseudogene-predictive variants,
2. the total insertion and deletion size (as some gene-pseudogene pairs have
   different sizes).

## Installation

This tool only uses standard python libraries so it should work with any python
version >=3.0.0 (tested on 3.12).


## Running
An example run
```
havanno --variant-vcf PURETARGET_VARIANTS_VCF --paraphase-dir DIR
```

## Output

The output looks like this, with JSON entry per paraphase locus
```
{
  "cyp21": {
    "cyp21_hap1": {
      "num_pathogenic_variants": 1,
      "num_pseudogene_variants": 2,
      "total_insertion_size": 0,
      "total_deletion_size": 0,
      "variants": {
        "pathogenic": [
          "CYP21A2_c518T>A]_830_841del_TCCTGGAAGGGC_M"
        ],
        "pseudogene": [
          "CYP21A2_c293-13_C>G_Pseudo_M",
          "CYP21A2_c518T>A]_830_841del_TCCTGGAAGGGC_M"
        ]
      }
    },
    "cyp21_hap2": {
      "num_pathogenic_variants": 0,
      "num_pseudogene_variants": 2,
      "total_insertion_size": 0,
      "total_deletion_size": 0,
      "variants": {
        "pathogenic": [],
        "pseudogene": [
          "CYP21A2_c.118C>T",
          "CYP21A2_c.138C>T"
        ]
      }
    },
    "cyp21_hap3": {
      "num_pathogenic_variants": 6,
      "num_pseudogene_variants": 18,
      "total_insertion_size": 0,
      "total_deletion_size": 128,
      "variants": {
        "pathogenic": [
          "CYP21A2_p.Pro31Leu_p.Pro31Leu_c.92C>T_M",
          "CYP21A2_c332_339del_p.Gly111fs",
          "CYP21A2_c518T>A]_830_841del_TCCTGGAAGGGC_M",
          "CYP21A2_p.Ile237Asn_p.Ile237Asn_c.710T>A_Val238Cluster_M",
          "CYP21A2_p.Val238Glu_p.Val238Glu_c.713T>A_Val238Cluster_M",
          "CYP21A2_p.Met240Lys_p.Met240Lys_c.719T>A_Val238Cluster_M"
        ],
        "pseudogene": [
          "CYP21A2_Pseudo4",
          "CYP21A2_Pseudo5",
          "CYP21A2_Pseudo6",
          "CYP21A2_Pseudo7",
          "CYP21A2_Pseudo8",
          "CYP21A2_Pseudo9",
          "CYP21A2_Pseudo10",
          "CYP21A2_Pseudo11",
          "CYP21A2_Pseudo12",
          "CYP21A2_Pseudo13",
          "CYP21A2_Pseudo14",
          "CYP21A2_Pseudo15",
          "CYP21A2_Pseudo16",
          "CYP21A2_c332_339del_p.Gly111fs",
          "CYP21A2_Pseudo17",
          "CYP21A2_c518T>A]_830_841del_TCCTGGAAGGGC_M",
          "CYP21A2_Pseudo18",
          "CYP21A2_Pseudo19"
        ]
      }
    },
    "cyp21_hap4": {
      "num_pathogenic_variants": 8,
      "num_pseudogene_variants": 18,
      "total_insertion_size": 0,
      "total_deletion_size": 128,
      "variants": {
        "pathogenic": [
          "CYP21A2_p.Pro31Leu_p.Pro31Leu_c.92C>T_M",
          "CYP21A2_c332_339del_p.Gly111fs",
          "CYP21A2_c518T>A]_830_841del_TCCTGGAAGGGC_M",
          "CYP21A2_p.Ile237Asn_p.Ile237Asn_c.710T>A_Val238Cluster_M",
          "CYP21A2_p.Val238Glu_p.Val238Glu_c.713T>A_Val238Cluster_M",
          "CYP21A2_p.Met240Lys_p.Met240Lys_c.719T>A_Val238Cluster_M",
          "CYP21A2_p.Gln319Ter_p.Gln319Ter_c.955C>T_M",
          "CYP21A2_p.Arg357Trp_p.Arg357Trp_c.1069C>T_M"
        ],
        "pseudogene": [
          "CYP21A2_Pseudo4",
          "CYP21A2_Pseudo5",
          "CYP21A2_Pseudo6",
          "CYP21A2_Pseudo7",
          "CYP21A2_Pseudo8",
          "CYP21A2_Pseudo9",
          "CYP21A2_Pseudo10",
          "CYP21A2_Pseudo11",
          "CYP21A2_Pseudo12",
          "CYP21A2_Pseudo13",
          "CYP21A2_Pseudo14",
          "CYP21A2_Pseudo15",
          "CYP21A2_Pseudo16",
          "CYP21A2_c332_339del_p.Gly111fs",
          "CYP21A2_Pseudo17",
          "CYP21A2_c518T>A]_830_841del_TCCTGGAAGGGC_M",
          "CYP21A2_Pseudo18",
          "CYP21A2_Pseudo19"
        ]
      }
    }
  }
}
```

Here paraphase reported 4 haplotypes. Hap1 contains a pathogenic mutation. Hap2
is a gene but contains two pseudogene-specific mutations. Haps 3 and 4 are
clearly pseudogenes given the number of pseudogene-predictive SNPs, the total
deletion size and the high number of both pseudogene-specific and pathogenic
mutations.

more documentation:


```
usage: havanno [-h] --variant-vcf VCF --paraphase-dir DIR

annotate paraphase outputs based on a variant list

options:
  -h, --help           show this help message and exit
  --variant-vcf VCF    variant list VCF file with boolean flags PATHO and PSEUDO
  --paraphase-dir DIR  Path to the paraphase output directory
```

## Where to find the --variants-vcf parameter?

It is up  to the user to define which variants they want to query. The
PureTarget variant list used to assess the accuracy of the protocol is
maintained [here](https://github.com/PacificBiosciences/ptcp/tree/main/meta/variant_list).
