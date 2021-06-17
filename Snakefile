#!/usr/bin/env python

configfile: "config.yaml"

rule target:
    input:
        ".fonts_registered.txt",
        "figures/rnaseq_summary.pdf",
        expand("figures/zf_chipseq_coverage/zf_chipseq_coverage-affinity-dependent-peaks-{dataset}.pdf", dataset=config['chipseq_coverage']),
        expand("figures/zf_chipseq_coverage_ratio_motifs/zf_chipseq_coverage-ratio-affinity-dependent-peaks-{dataset}.pdf", dataset=config['chipseq_coverage_ratio']),
        expand("figures/zf_chipseq_coverage_ratio_motifs/zf_chipseq_coverage-ratio-affinity-dependent-peaks-{dataset}-freescale.pdf", dataset=config['chipseq_coverage_ratio']),
        "figures/zf_venus_reporter_datavis.pdf",
        "figures/rnaseq_summary/rnaseq_summary.pdf",
        "figures/rnaseq_summary_alternate/rnaseq_summary_scatter_highlight_motifs.pdf",
        "figures/zf_chipseq_global_coverage/zf_chipseq_global_coverage.pdf",
        "figures/expression_vs_chip_enrichment/expression_vs_chip_enrichment.pdf",
        "figures/motif_distance_vs_rnaseq/motif_distance_vs_rnaseq.pdf",
        "figures/rna_vs_chip_windows/rna_vs_chip_windows.pdf",
        expand("figures/chipseq_global_abundance/chipseq_global_abundance_{dataset}.pdf", dataset=config["chipseq_global_abundance"]),
        expand("figures/reporter_chipseq_qc/reporter_chipseq_qc_{dataset}.pdf", dataset=config['reporter_chipseq_qc']),
        "figures/rnaseq_summary/rnaseq_summary_mutants.pdf",
        "figures/rnaseq_summary/rnaseq_summary_mutants_all.pdf",
        'figures/chip_hits_w_rnaseq_info/chip_joined.tsv',
        "figures/rnaseq_volcano_custom.pdf",
        "figures/rnaseq_summary/rnaseq_summary_13_6.pdf",
        "figures/rnaseq_volcano_custom_13_6.pdf",
        "figures/rnaseq_maplot_13_6.pdf",

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
        fonts = ".fonts_registered.txt",
        coverage = lambda wc: config["chipseq_coverage"][wc.dataset]["coverage"],
        summit_annotation = lambda wc: config["chipseq_coverage"][wc.dataset]["peak_annotation"],
        transcript_annotation = lambda wc: config["chipseq_coverage"][wc.dataset]["transcripts"],
        orf_annotation = lambda wc: config["chipseq_coverage"][wc.dataset]["orfs"],
    output:
        coverage = "figures/zf_chipseq_coverage/zf_chipseq_coverage-affinity-dependent-peaks-{dataset}.pdf"
    conda:
        "envs/plot.yaml"
    script:
        "scripts/chipseq_coverage.R"

rule chipseq_coverage_ratio:
    input:
        fonts = ".fonts_registered.txt",
        coverage = lambda wc: config["chipseq_coverage_ratio"][wc.dataset]["coverage"],
        summit_annotation = lambda wc: config["chipseq_coverage_ratio"][wc.dataset]["peak_annotation"],
        transcript_annotation = lambda wc: config["chipseq_coverage_ratio"][wc.dataset]["transcripts"],
        orf_annotation = lambda wc: config["chipseq_coverage_ratio"][wc.dataset]["orfs"],
        motif_annotation = lambda wc: config["chipseq_coverage_ratio"][wc.dataset]["motifs"],
    output:
        coverage = "figures/zf_chipseq_coverage_ratio_motifs/zf_chipseq_coverage-ratio-affinity-dependent-peaks-{dataset}.pdf"
    params:
        filter_groups = lambda wc: config["chipseq_coverage_ratio"][wc.dataset]["filter_groups"],
    conda:
        "envs/plot.yaml"
    script:
        "scripts/chipseq_coverage_ratio_motifs.R"

rule chipseq_coverage_ratio_freescale:
    input:
        fonts = ".fonts_registered.txt",
        coverage = lambda wc: config["chipseq_coverage_ratio"][wc.dataset]["coverage"],
        summit_annotation = lambda wc: config["chipseq_coverage_ratio"][wc.dataset]["peak_annotation"],
        transcript_annotation = lambda wc: config["chipseq_coverage_ratio"][wc.dataset]["transcripts"],
        orf_annotation = lambda wc: config["chipseq_coverage_ratio"][wc.dataset]["orfs"],
        motif_annotation = lambda wc: config["chipseq_coverage_ratio"][wc.dataset]["motifs"],
    output:
        coverage = "figures/zf_chipseq_coverage_ratio_motifs/zf_chipseq_coverage-ratio-affinity-dependent-peaks-{dataset}-freescale.pdf"
    params:
        filter_groups = lambda wc: config["chipseq_coverage_ratio"][wc.dataset]["filter_groups"],
    conda:
        "envs/plot.yaml"
    script:
        "scripts/chipseq_coverage_ratio_motifs_freescale.R"

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

rule expression_vs_chip_enrichment:
    input:
        theme = config["theme_path"],
        fonts = ".fonts_registered.txt",
        high_rna = config["expression_vs_chip_enrichment"]["high_rna"],
        high_chip = config["expression_vs_chip_enrichment"]["high_chip"],
        low_rna = config["expression_vs_chip_enrichment"]["low_rna"],
        low_chip = config["expression_vs_chip_enrichment"]["low_chip"],
        clamp_rna = config["expression_vs_chip_enrichment"]["clamp_rna"],
        clamp_chip = config["expression_vs_chip_enrichment"]["clamp_chip"],
    output:
        pdf = "figures/expression_vs_chip_enrichment/expression_vs_chip_enrichment.pdf",
    conda:
        "envs/plot.yaml"
    script:
        "scripts/expression_vs_chip_enrichment.R"

rule motif_distance_vs_rnaseq:
    input:
        theme = config["theme_path"],
        fonts = ".fonts_registered.txt",
        rnaseq = config["motif_distance_vs_rnaseq"]["rnaseq"],
        motifs = config["motif_distance_vs_rnaseq"]["motifs"],
    output:
        pdf = "figures/motif_distance_vs_rnaseq/motif_distance_vs_rnaseq.pdf",
    conda:
        "envs/plot.yaml"
    script:
        "scripts/motif_distance_vs_rnaseq.R"

rule rna_vs_chip_windows:
    input:
        theme = config["theme_path"],
        fonts = ".fonts_registered.txt",
        rna = config["rna_vs_chip_windows"]["rna"],
        chip = config["rna_vs_chip_windows"]["chip"],
    output:
        pdf = "figures/rna_vs_chip_windows/rna_vs_chip_windows.pdf",
    conda:
        "envs/plot.yaml"
    script:
        "scripts/rna_vs_chip_windows.R"

rule chipseq_global_abundance:
    input:
        theme = config["theme_path"],
        fonts = ".fonts_registered.txt",
        data = lambda wc: config["chipseq_global_abundance"][wc.dataset],
    output:
        pdf = "figures/chipseq_global_abundance/chipseq_global_abundance_{dataset}.pdf",
    conda:
        "envs/plot.yaml"
    script:
        "scripts/chipseq_abundance_barplot.R"

rule reporter_chipseq_qc:
    input:
        theme = config["theme_path"],
        fonts = ".fonts_registered.txt",
        data = lambda wc: config["reporter_chipseq_qc"][wc.dataset],
    output:
        pdf = "figures/reporter_chipseq_qc/reporter_chipseq_qc_{dataset}.pdf",
    conda:
        "envs/plot.yaml"
    script:
        "scripts/reporter_chipseq_qc.R"

rule rnaseq_summary_mutants:
    input:
        theme = config["theme_path"],
        fonts = ".fonts_registered.txt",
        high_affinity = config["rnaseq_summary"]["high-affinity"],
        low_affinity = config["rnaseq_summary"]["low-affinity"],
        low_affinity_w_clamp = config["rnaseq_summary"]["low-affinity-w-clamp"],
        high_affinity_zev_mutant = config["rnaseq_summary_mutants"]["high_affinity_zev_mutant"],
        high_affinity_zftf_mutant = config["rnaseq_summary_mutants"]["high_affinity_zftf_mutant"],
    output:
        pdf = "figures/rnaseq_summary/rnaseq_summary_mutants.pdf",
    params:
        fdr = config["rnaseq_summary"]["rnaseq_fdr"],
    conda:
        "envs/ggforce.yaml"
    script:
        "scripts/rnaseq_summary_mutants_new.R"

rule rnaseq_summary_mutants_all:
    input:
        theme = config["theme_path"],
        fonts = ".fonts_registered.txt",
        high_affinity = config["rnaseq_summary"]["high-affinity"],
        low_affinity = config["rnaseq_summary"]["low-affinity"],
        low_affinity_w_clamp = config["rnaseq_summary"]["low-affinity-w-clamp"],
        high_affinity_zev_mutant = config["rnaseq_summary_mutants"]["high_affinity_zev_mutant"],
        high_affinity_zftf_mutant = config["rnaseq_summary_mutants"]["high_affinity_zftf_mutant"],
    output:
        pdf = "figures/rnaseq_summary/rnaseq_summary_mutants_all.pdf",
    params:
        fdr = config["rnaseq_summary"]["rnaseq_fdr"],
    conda:
        "envs/ggforce.yaml"
    script:
        "scripts/rnaseq_summary_mutants_all.R"

rule join_chip:
    input:
        affinity_dependent = config['chip_hits_w_rnaseq_info']['affinity_dependent_peaks'],
        high = config['chip_hits_w_rnaseq_info']['chip_results']['high'],
        low = config['chip_hits_w_rnaseq_info']['chip_results']['low'],
        low_clamp = config['chip_hits_w_rnaseq_info']['chip_results']['low_clamp'],
    output:
        tsv = 'figures/chip_hits_w_rnaseq_info/chip_joined.tsv'
    conda:
        "envs/plot.yaml"
    script:
        "scripts/join_chip.R"

rule join_rna:
    input:
        high = config['chip_hits_w_rnaseq_info']['rna_results']['high'],
        low = config['chip_hits_w_rnaseq_info']['rna_results']['low'],
        low_clamp = config['chip_hits_w_rnaseq_info']['rna_results']['low_clamp'],
    output:
        tsv = 'figures/chip_hits_w_rnaseq_info/rna_joined.tsv'
    conda:
        "envs/plot.yaml"
    script:
        "scripts/join_rna.R"

rule chip_hits_w_rnaseq_info:
    input:
        chip = 'figures/chip_hits_w_rnaseq_info/chip_joined.tsv',
        rna = 'figures/chip_hits_w_rnaseq_info/rna_joined.tsv'
    output:
        'figures/chip_hits_w_rnaseq_info/chip_hits_w_rnaseq_info.tsv'
    shell: """
        bedtools closest \
            -a <(tail -n +2 {input.chip} | \
                 awk 'BEGIN{{OFS="\t"}} {{print $0, NR}}' | \
                 sort -k1,1 -k2,2n) \
            -b <(tail -n +2 {input.rna} | \
                 sort -k1,1 -k2,2n) \
            -k 2 | \
        sort -k8,8n | \
        cut -f8 --complement | \
        cat <(paste <(head -n 1 {input.chip}) <(head -n 1 {input.rna})) - > \
        {output}
        """

rule rnaseq_volcano_custom:
    input:
        theme = config["theme_path"],
        low_v_high = config["rnaseq_volcano_custom"]["low_v_high"],
        low_clamp_v_high = config["rnaseq_volcano_custom"]["low_clamp_v_high"],
        low_clamp_v_low = config["rnaseq_volcano_custom"]["low_clamp_v_low"],
        fonts = ".fonts_registered.txt"
    output:
        volcano = "figures/rnaseq_volcano_custom.pdf",
    conda:
        "envs/plot.yaml"
    params:
        fdr = config["rnaseq"]["fdr"]
    script:
        "scripts/rnaseq_volcano_custom.R"

rule rnaseq_summary_13_6:
    input:
        theme = config["theme_path"],
        fonts = ".fonts_registered.txt",
        high_affinity = config["rnaseq_summary_13_6"]["high-affinity"],
        low_affinity = config["rnaseq_summary_13_6"]["low-affinity"],
        low_affinity_w_clamp = config["rnaseq_summary_13_6"]["low-affinity-w-clamp"],
        # motif_results = config["rnaseq_summary"]["motif_results"],
    output:
        pdf = "figures/rnaseq_summary/rnaseq_summary_13_6.pdf",
    params:
        fdr = config["rnaseq_summary"]["rnaseq_fdr"],
    conda:
        "envs/ggforce.yaml"
    script:
        "scripts/rnaseq_summary_13-6.R"

rule rnaseq_volcano_custom_13_6:
    input:
        theme = config["theme_path"],
        low_v_high = config["rnaseq_volcano_custom_13_6"]["low_v_high"],
        low_clamp_v_high = config["rnaseq_volcano_custom_13_6"]["low_clamp_v_high"],
        low_clamp_v_low = config["rnaseq_volcano_custom_13_6"]["low_clamp_v_low"],
        fonts = ".fonts_registered.txt"
    output:
        volcano = "figures/rnaseq_volcano_custom_13_6.pdf",
    conda:
        "envs/plot.yaml"
    params:
        fdr = config["rnaseq"]["fdr"]
    script:
        "scripts/rnaseq_volcano_custom.R"

rule rnaseq_figures_13_6:
    input:
        high_affinity = config["rnaseq_13_6"]["high-affinity"],
        low_affinity = config["rnaseq_13_6"]["low-affinity"],
        low_affinity_w_clamp = config["rnaseq_13_6"]["low-affinity-w-clamp"],
        fonts = ".fonts_registered.txt"
    output:
        volcano = "figures/rnaseq_volcano_13_6.pdf",
        maplot = "figures/rnaseq_maplot_13_6.pdf",
        summary = "figures/rnaseq_summary_13_6.pdf",
    conda:
        "envs/plot.yaml"
    params:
        fdr = config["rnaseq"]["fdr"]
    script:
        "scripts/rnaseq_figures.R"
