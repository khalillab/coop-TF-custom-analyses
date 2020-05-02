main = function(theme_path = "custom_theme.R",
     rnaseq_path = "high-affinity-ZF-v-reporter-only_rnaseq-libsizenorm-verified_genes_plus_Venus-diffexp-results-all.tsv",
     motif_path = "high-affinity-upregulated-transcripts-v-all-transcripts_control_unmergedFIMOresults.tsv.gz",
     pdf_out="test.pdf"){
    source(theme_path)

    df = read_tsv(motif_path) %>%
        filter(motif_start >= 0) %>%
        inner_join(read_tsv(rnaseq_path),
                   by=c("chrom", "region_id"="name", "region_strand"="strand")) %>%
        mutate(motif_position = (motif_start + motif_end) / 2,
               relative_motif_position = if_else(region_strand=="+",
                                                 motif_position - region_end,
                                                 region_start - motif_position)) %>%
        group_by(chrom, start, end, region_id) %>%
        mutate(n_motifs = n(),
               n_motifs_labeled = paste(n_motifs, "motifs"),
               category=case_when(log10_padj > -log10(0.1) & log2_foldchange > 0 ~ "upregulated",
                                  log10_padj > -log10(0.1) & log2_foldchange < 0 ~ "downregulated",
                                  TRUE ~ "not significantly changed"))

    plot = ggplot(data=df,
           aes(x=relative_motif_position,
               y=log2_foldchange)) +
        geom_hline(yintercept=0,
                   size=0.2,
                   color="gray70") +
        geom_vline(xintercept=0,
                   size=0.2,
                   color="gray70") +
        geom_point(aes(color=region_id),
                   size=0.5,
                   alpha=0.5,
                   shape=16) +
        geom_smooth(size=0.2) +
        scale_x_continuous(labels=function(x) if_else(x==0, "TSS", paste(x, "bp")),
                           breaks=scales::pretty_breaks(2),
                           name="motif position relative to TSS") +
        scale_y_continuous(quote("log"[2] ~ textstyle(frac("high-affinity",
                                                           "reporter-only")))) +
        facet_grid(n_motifs_labeled~category,
                   margins=c("category")) +
        theme_default +
        theme(legend.position="none",
              panel.grid=element_blank(),
              axis.title.y=element_text(angle=0,
                                        vjust=0.5),
              strip.text.y=element_text(angle=0))

    ggsave(pdf_out,
           plot=plot,
           width=16*1.5,
           height=9*1.5,
           units="cm",
           device=cairo_pdf)
}

main(theme_path = snakemake@input[["theme"]],
     rnaseq_path = snakemake@input[["rnaseq"]],
     motif_path = snakemake@input[["motifs"]],
     pdf_out= snakemake@output[["pdf"]])

