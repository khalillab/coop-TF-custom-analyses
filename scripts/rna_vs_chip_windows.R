main = function(theme_path = "custom_theme.R",
                rna_path = "high-affinity-ZF-v-reporter-only_rnaseq-libsizenorm-verified_genes_plus_Venus-diffexp-results-all.tsv",
                chip_path = "high-affinity-ZF-v-reporter-only_ZF-chipseq-libsizenorm-verified_genes_plus_Venus_upstream_windows-diffbind-results-all.tsv",
                pdf_out="test.pdf"){
    source(theme_path)

    df = read_tsv(rna_path) %>%
        left_join(read_tsv(chip_path),
                  by=c("chrom", "name", "strand"),
                  suffix=c("_rna", "_chip"))

    scatter = ggplot(data=df,
           aes(y=log2_foldchange,
               x=log2FC_enrichment)) +
        geom_vline(xintercept=0,
                   color="gray70",
                   size=0.2) +
        geom_hline(yintercept=0,
                   color="gray70",
                   size=0.2) +
        stat_binhex(geom="point",
                    aes(color=..count..),
                    bins=200,
                    shape=16,
                    size=0.3,
                    alpha=0.8) +
        scale_x_continuous(name=quote("ZF ChIP enrichment in [-300 bp, TSS]: log"[2] ~
                                          textstyle(frac("high-affinity",
                                                         "reporter-only"))),
                           breaks=scales::pretty_breaks(5)) +
        scale_y_continuous(name=quote(atop("RNAseq:",
                                           "log"[2] ~ textstyle(frac("high-affinity",
                                                                     "reporter-only")))),
                           breaks=scales::pretty_breaks(4)) +
        scale_color_viridis_c(option="cividis") +
        theme_default +
        theme(panel.grid=element_blank(),
              legend.position="none",
              axis.title.y=element_text(angle=0,
                                        vjust=0.5))

    ggsave(pdf_out,
           plot=scatter,
           width=8.5,
           height=8.5*9/16,
           unit="cm",
           device=cairo_pdf)
}

main(theme_path = snakemake@input[["theme"]],
     rna_path = snakemake@input[["rna"]],
     chip_path = snakemake@input[["chip"]],
     pdf_out = snakemake@output[["pdf"]])

