#!/usr/bin/env python

configfile: "config.yaml"

rule target:
    input:
        ".fonts_registered.txt",
        "figures/rnaseq_summary.pdf",
        "figures/zf_chipseq_coverage-affinity-dependent-peaks.pdf"

rule register_fonts:
    input:
        fonts_path = config["fonts_path"],
    output:
        output_path = ".fonts_registered.txt"
    conda:
        "envs/plot.yaml"
    script:
        "scripts/register_fonts.R"

rule rnaseq_figures:
    input:
        high_affinity = config["rnaseq"]["high-affinity"],
        low_affinity = config["rnaseq"]["low-affinity"],
        low_affinity_w_clamp = config["rnaseq"]["low-affinity-w-clamp"],
        fonts = ".fonts_registered.txt"
    output:
        volcano = "figures/rnaseq_volcano.pdf",
        maplot = "figures/rnaseq_maplot.pdf",
        summary = "figures/rnaseq_summary.pdf",
    conda:
        "envs/plot.yaml"
    params:
        fdr = config["rnaseq"]["fdr"]
    script:
        "scripts/rnaseq_figures.R"

rule chipseq_coverage:
    input:
        coverage = config["chipseq_coverage"]["coverage"],
        summit_annotation = config["chipseq_coverage"]["peak_annotation"],
        transcript_annotation = config["chipseq_coverage"]["transcripts"],
        orf_annotation = config["chipseq_coverage"]["orfs"],
    output:
        coverage = "figures/zf_chipseq_coverage-affinity-dependent-peaks.pdf"
    conda:
        "envs/plot.yaml"
    script:
        "scripts/chipseq_coverage.R"


