#!/usr/bin/env python3
"Create dark genes IGV session files from a JSON file containing data to render."

import argparse
import json
from jinja2 import Environment, FileSystemLoader


TEMPLATE = """\
<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<Session genome="hg38" locus="chrX:154864801-155055733" version="8">
    <Resources> {% for resource in resources -%} <Resource path="{{ resource.path }}"
            type="{{ resource.type }}" /> {% endfor -%} </Resources>
    <Panel height="200" name="DataPanel" width="1600">
        <Track attributeKey="{{ sample_id }} F8 inversion VCF"
            clazz="org.broad.igv.variant.VariantTrack"
            colorScale="ContinuousColorScale;0.0;0.0;255,255,255;0,0,178" displayMode="EXPANDED"
            fontSize="10" groupByStrand="false" id="{{ f8inversionvcf }}"
            name="{{ sample_id }} F8 inversion VCF" siteColorMode="ALLELE_FREQUENCY" visible="true" />
        <Track attributeKey="{{ sample_id }} consensus VCF"
            clazz="org.broad.igv.variant.VariantTrack"
            colorScale="ContinuousColorScale;0.0;0.0;255,255,255;0,0,178" displayMode="EXPANDED"
            fontSize="10" groupByStrand="false" id="{{ consensusvcf }}"
            name="{{ sample_id }} consensus VCF" siteColorMode="ALLELE_FREQUENCY" visible="true" />
    </Panel>
    <Panel height="80" name="Panel1718672712439" width="1600">
        <Track attributeKey="{{ sample_id }} consensus alignments"
            clazz="org.broad.igv.sam.AlignmentTrack" color="185,185,185" experimentType="THIRD_GEN"
            fontSize="10" id="{{ consensusbam }}" name="{{ sample_id }} consensus alignments"
            visible="true">
            <RenderOptions />
        </Track>
    </Panel>
    <Panel height="800" name="Panel1718672723602" width="1600">
        <Track attributeKey="{{ sample_id }} HiFi alignments"
            clazz="org.broad.igv.sam.AlignmentTrack" color="185,185,185" displayMode="SQUISHED"
            experimentType="THIRD_GEN" fontSize="10" id="{{ hifi_readsbam }}"
            name="{{ sample_id }} HiFi alignments" visible="true">
            <RenderOptions colorOption="YC_TAG" groupByOption="PHASE" />
        </Track>
    </Panel>
    <Panel height="120" name="FeaturePanel" width="1600">
        <Track attributeKey="Reference sequence" clazz="org.broad.igv.track.SequenceTrack"
            fontSize="10" id="Reference sequence" name="Reference sequence"
            sequenceTranslationStrandValue="+" shouldShowTranslation="false" visible="true" />
        <Track attributeKey="Refseq Genes" clazz="org.broad.igv.track.FeatureTrack"
            colorScale="ContinuousColorScale;0.0;836.0;255,255,255;0,0,178" fontSize="10"
            groupByStrand="false"
            id="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/ncbiRefSeq.txt.gz"
            name="Refseq Genes" visible="true" />
        <Track attributeKey="pbaa guide regions" clazz="org.broad.igv.track.FeatureTrack"
            colorScale="ContinuousColorScale;0.0;0.0;255,255,255;0,0,178" fontSize="10"
            groupByStrand="false" id="{{ pbaa_guides }}" name="pbaa guide regions" visible="true" />
        <Track attributeKey="masked regions" clazz="org.broad.igv.track.FeatureTrack"
            colorScale="ContinuousColorScale;0.0;0.0;255,255,255;0,0,178" fontSize="10"
            groupByStrand="false" id="{{ pbaa_mask }}" name="masked regions" visible="true" />
    </Panel>
    <PanelLayout dividerFractions="0.17794486215538846,0.24561403508771928,0.8922305764411027" />
    <HiddenAttributes>
        <Attribute name="DATA FILE" />
        <Attribute name="DATA TYPE" />
        <Attribute name="NAME" />
    </HiddenAttributes>
</Session>
"""


def replace_path_prefix(path, old_path_prefix, new_path_prefix):
    "Replace the old path prefix with the new path prefix."
    if old_path_prefix is None or new_path_prefix is None:
        return path
    return path.replace(old_path_prefix, new_path_prefix)


def extract_relevant_data(data, old_path_prefix, new_path_prefix):
    "Extract the relevant data from the JSON file."
    pbaa_guides = replace_path_prefix(
        data["puretarget_darkgenes.pbaa_guides"], old_path_prefix, new_path_prefix
    )
    pbaa_mask = replace_path_prefix(
        data["puretarget_darkgenes.pbaa_mask"], old_path_prefix, new_path_prefix
    )
    for i, sample in enumerate(data["puretarget_darkgenes.sample_names"]):
        resources = []
        sample_id = sample
        f8inversionvcf = replace_path_prefix(
            data["puretarget_darkgenes.f8_pbaa_vcf"][i],
            old_path_prefix,
            new_path_prefix,
        )
        consensusvcf = replace_path_prefix(
            data["puretarget_darkgenes.minipileup_vcfs"][i],
            old_path_prefix,
            new_path_prefix,
        )
        consensusbam = replace_path_prefix(
            data["puretarget_darkgenes.consensus_bams"][i],
            old_path_prefix,
            new_path_prefix,
        )
        hifi_readsbam = replace_path_prefix(
            data["puretarget_darkgenes.painted_bams"][i],
            old_path_prefix,
            new_path_prefix,
        )
        resources.append(
            {
                "path": replace_path_prefix(
                    f8inversionvcf, old_path_prefix, new_path_prefix
                ),
                "type": "vcf",
            }
        )
        resources.append(
            {
                "path": replace_path_prefix(
                    consensusvcf, old_path_prefix, new_path_prefix
                ),
                "type": "vcf",
            }
        )
        resources.append(
            {
                "path": replace_path_prefix(
                    consensusbam, old_path_prefix, new_path_prefix
                ),
                "type": "bam",
            }
        )
        resources.append(
            {
                "path": replace_path_prefix(
                    hifi_readsbam, old_path_prefix, new_path_prefix
                ),
                "type": "bam",
            }
        )
        resources.append(
            {
                "path": replace_path_prefix(
                    pbaa_guides, old_path_prefix, new_path_prefix
                ),
                "type": "bed",
            }
        )
        resources.append(
            {
                "path": replace_path_prefix(
                    pbaa_mask, old_path_prefix, new_path_prefix
                ),
                "type": "bed",
            }
        )
        yield {
            "resources": resources,
            "sample_id": sample_id,
            "f8inversionvcf": f8inversionvcf,
            "consensusvcf": consensusvcf,
            "consensusbam": consensusbam,
            "hifi_readsbam": hifi_readsbam,
            "pbaa_guides": pbaa_guides,
            "pbaa_mask": pbaa_mask,
        }


def render_dark_genes_igv_session(resources):
    "Render the dark genes IGV session template with resources."
    env = Environment(loader=FileSystemLoader("."))
    template = env.from_string(TEMPLATE)
    return template.render(resources)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", help="JSON file containing data to render.")
    parser.add_argument("prefix", help="Output file prefix.")
    parser.add_argument("--old-path-prefix", "-o", help="Old path prefix.")
    parser.add_argument("--new-path-prefix", "-n", help="New path prefix.")
    args = parser.parse_args()
    with open(args.input, "r") as infile:
        data = json.load(infile)
        for sample in extract_relevant_data(
            data, args.old_path_prefix, args.new_path_prefix
        ):
            with open(
                ".".join([args.prefix, sample["sample_id"], "dark_genes.igv.xml"]), "w"
            ) as outfile:
                outfile.write(render_dark_genes_igv_session(sample))


if __name__ == "__main__":
    main()
