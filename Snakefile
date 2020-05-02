#!/usr/bin/env python

configfile: "config.yaml"

rule target:
    input:
        ".fonts_registered.txt",
        "figures/rnaseq_summary.pdf",
        "figures/zf_chipseq_coverage-affinity-dependent-peaks.pdf",
        "figures/zf_venus_reporter_datavis.pdf",
        "figures/rnaseq_summary/rnaseq_summary.pdf",
        "figures/rnaseq_summary_alternate/rnaseq_summary_scatter_highlight_motifs.pdf",
        "figures/zf_chipseq_global_coverage/zf_chipseq_global_coverage.pdf",

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
        output_path = ".fonts_registered.txt",
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

rule zf_venus_reporter_datavis:
    input:
        theme = config["theme_path"],
        fonts = ".fonts_registered.txt",
        target_annotation = config["zf_venus_reporter_datavis"]["target_annotation"],
        transcript_annotation = config["zf_venus_reporter_datavis"]["transcript_annotation"],
        orf_annotation = config["zf_venus_reporter_datavis"]["orf_annotation"],
        motif_annotation = config["zf_venus_reporter_datavis"]["motif_annotation"],
        data = config["zf_venus_reporter_datavis"]["data"],
    output:
        plot = "figures/zf_venus_reporter_datavis.pdf"
    conda:
        "envs/plot.yaml"
    script:
        "scripts/zf_venus_reporter_datavis.R"

rule rnaseq_summary:
    input:
        theme = config["theme_path"],
        fonts = ".fonts_registered.txt",
        high_affinity = config["rnaseq_summary"]["high-affinity"],
        low_affinity = config["rnaseq_summary"]["low-affinity"],
        low_affinity_w_clamp = config["rnaseq_summary"]["low-affinity-w-clamp"],
        motif_results = config["rnaseq_summary"]["motif_results"],
    output:
        pdf = "figures/rnaseq_summary/rnaseq_summary.pdf",
    params:
        fdr = config["rnaseq_summary"]["rnaseq_fdr"],
    conda:
        "envs/ggforce.yaml"
    script:
        "scripts/rnaseq_summary.R"

rule rnaseq_summary_alternate:
    input:
        theme = config["theme_path"],
        fonts = ".fonts_registered.txt",
        high_affinity = config["rnaseq_summary"]["high-affinity"],
        low_affinity = config["rnaseq_summary"]["low-affinity"],
        low_affinity_w_clamp = config["rnaseq_summary"]["low-affinity-w-clamp"],
        motif_results = config["rnaseq_summary"]["motif_results"],
    output:
        scatter_motifs = "figures/rnaseq_summary_alternate/rnaseq_summary_scatter_highlight_motifs.pdf",
        scatter_chip = "figures/rnaseq_summary_alternate/rnaseq_summary_scatter_highlight_chip.pdf",
    params:
        fdr = config["rnaseq_summary"]["rnaseq_fdr"],
    conda:
        "envs/plot.yaml"
    script:
        "scripts/rnaseq_summary_alternate.R"


rule zf_chipseq_global_coverage:
    input:
        theme = config["theme_path"],
        fonts = ".fonts_registered.txt",
        data = config["zf_chipseq_global_coverage"]["data"],
    output:
        pdf = "figures/zf_chipseq_global_coverage/zf_chipseq_global_coverage.pdf",
    conda:
        "envs/plot.yaml"
    script:
        "scripts/zf_chipseq_global_coverage.R"


