
import = function(path,
                  fdr,
                  strain_id){
    read_tsv(path) %>%
        mutate(significant=log10_padj > -log10(fdr),
               strain=strain_id) %>%
        return()
}

main = function(theme_path = "custom_theme.R",
                input_high_affinity = "high-affinity-ZF-v-reporter-only_rnaseq-libsizenorm-verified_genes_plus_Venus-diffexp-results-all.tsv",
                input_low_affinity = "low-affinity-ZF-v-reporter-only_rnaseq-libsizenorm-verified_genes_plus_Venus-diffexp-results-all.tsv",
                input_low_affinity_w_clamp = "low-affinity-ZF-with-clamp-v-reporter-only_rnaseq-libsizenorm-verified_genes_plus_Venus-diffexp-results-all.tsv",
                motif_results = "high-affinity-upregulated-transcripts-v-all-transcripts_control_unmergedFIMOresults.tsv.gz",
                fdr=0.1,
                # panel_letter = "A",
                # fig_width = 11.4,
                # fig_height= 6,
                scatter_motifs_out,
                scatter_chip_out){

    source(theme_path)

    motif_hits = read_tsv(motif_results) %>%
        filter(motif_start >= 0) %>%
        pull(region_id) %>%
        unique()

    chip_hits = c("VRP1",
                  "Venus",
                  "PHO8",
                  "FTR1",
                  "PMA2",
                  "PRY3",
                  "ENO1",
                  "AIM3")

    high = import(input_high_affinity,
                  fdr,
                  "high-affinity")
    low = import(input_low_affinity,
                 fdr,
                 "low-affinity")
    clamp = import(input_low_affinity_w_clamp,
                   fdr,
                   "low-affinity + clamp")
    gene_list = bind_rows(high,
                          low,
                          clamp) %>%
        filter(significant) %>%
        distinct(name) %>%
        pull(name)

    df = bind_rows(high,
                   clamp) %>%
        left_join(low,
                  by=c("chrom", "start", "end", "name", "strand"),
                  suffix=c("_y", "_x")) %>%
        filter(name %in% gene_list) %>%
        mutate(motif_hit = name %in% motif_hits,
               chip_hit = name %in% chip_hits)

    rnaseq_summmary_alt_motifs = ggplot(data=df,
           aes(x=log2_foldchange_x,
               y=log2_foldchange_y,
               color=strain_y)) +
        geom_vline(xintercept=0,
                   size=0.2,
                   color="gray70") +
        geom_hline(yintercept=0,
                   size=0.2,
                   color="gray70") +
        geom_abline(slope=1,
                    intercept=0,
                    size=0.2,
                    color="gray70") +
        geom_point(aes(shape=motif_hit),
                   alpha=0.6,
                   size=1,
                   stroke=0.4) +
        scale_x_continuous(name=quote("log"[2] ~ textstyle(frac("low-affinity",
                                                                "reporter-only"))),
                           breaks=scales::pretty_breaks(4)) +
        scale_y_continuous(name=quote("log"[2] ~ textstyle(frac("strain",
                                                                "reporter-only"))),
                           breaks=scales::pretty_breaks(6)) +
        scale_color_viridis_d(end=0.8,
                              name=NULL) +
        scale_shape_manual(values=c(1, 16),
                           name="motif hit") +
        theme_default +
        theme(panel.grid=element_blank(),
              axis.title.y=element_text(angle=0,
                                        vjust=0.5),
              legend.position="left")

    rnaseq_summmary_alt_chip = ggplot(data=df,
           aes(x=log2_foldchange_x,
               y=log2_foldchange_y,
               color=strain_y)) +
        geom_vline(xintercept=0,
                   size=0.2,
                   color="gray70") +
        geom_hline(yintercept=0,
                   size=0.2,
                   color="gray70") +
        geom_abline(slope=1,
                    intercept=0,
                    size=0.2,
                    color="gray70") +
        geom_point(aes(shape=chip_hit),
                   alpha=0.6,
                   size=1,
                   stroke=0.4) +
        scale_x_continuous(name=quote("log"[2] ~ textstyle(frac("low-affinity",
                                                                "reporter-only"))),
                           breaks=scales::pretty_breaks(4)) +
        scale_y_continuous(name=quote("log"[2] ~ textstyle(frac("strain",
                                                                "reporter-only"))),
                           breaks=scales::pretty_breaks(6)) +
        scale_color_viridis_d(end=0.8,
                              name=NULL) +
        scale_shape_manual(values=c(1, 16),
                           name="ChIP hit") +
        theme_default +
        theme(panel.grid=element_blank(),
              axis.title.y=element_text(angle=0,
                                        vjust=0.5),
              legend.position="left")

    ggsave(scatter_motifs_out,
           plot=rnaseq_summmary_alt_motifs,
           width=16,
           height=9,
           units="cm",
           device=cairo_pdf)
    ggsave(scatter_chip_out,
           plot=rnaseq_summmary_alt_chip,
           width=16,
           height=9,
           units="cm",
           device=cairo_pdf)

}

main(theme_path = snakemake@input[["theme"]],
     input_high_affinity=snakemake@input[["high_affinity"]],
     input_low_affinity=snakemake@input[["low_affinity"]],
     input_low_affinity_w_clamp=snakemake@input[["low_affinity_w_clamp"]],
     motif_results=snakemake@input[["motif_results"]],
     fdr=snakemake@params[["fdr"]],
     # panel_letter=snakemake@params[["panel_letter"]],
     # fig_width=snakemake@params[["fig_width"]],
     # fig_height=snakemake@params[["fig_height"]],
     scatter_motifs_out= snakemake@output[["scatter_motifs"]],
     scatter_chip_out =snakemake@output[["scatter_chip"]])

