import = function(df,
                  rna_path,
                  chip_path,
                  strain_id){
    read_tsv(rna_path) %>%
        inner_join(read_tsv(chip_path),
                   by=c("chrom", "start", "end", "name", "strand"),
                   suffix=c("_rna", "_chip")) %>%
        mutate(strain=strain_id) %>%
        bind_rows(df, .) %>%
        return()
}

main = function(theme_path = "custom_theme.R",
                high_rna_path = "high-affinity-ZF-v-reporter-only_rnaseq-libsizenorm-verified_genes_plus_Venus-diffexp-results-all.tsv",
                high_chip_path = "high-affinity-ZF-v-reporter-only_ZF-chipseq-libsizenorm-verified_genes_plus_Venus-diffbind-results-all.tsv",
                low_rna_path = "low-affinity-ZF-v-reporter-only_rnaseq-libsizenorm-verified_genes_plus_Venus-diffexp-results-all.tsv",
                low_chip_path = "low-affinity-ZF-v-reporter-only_ZF-chipseq-libsizenorm-verified_genes_plus_Venus-diffbind-results-all.tsv",
                clamp_rna_path = "low-affinity-ZF-with-clamp-v-reporter-only_rnaseq-libsizenorm-verified_genes_plus_Venus-diffexp-results-all.tsv",
                clamp_chip_path = "low-affinity-ZF-with-clamp-v-reporter-only_ZF-chipseq-libsizenorm-verified_genes_plus_Venus-diffbind-results-all.tsv",
                pdf_out="test.pdf"){
    source(theme_path)
    library(cowplot)

    df = tibble() %>%
        import(high_rna_path,
               high_chip_path,
               "high-affinity") %>%
        import(low_rna_path,
               low_chip_path,
               "low-affinity") %>%
        import(clamp_rna_path,
               clamp_chip_path,
               "low-affinity + clamp")

    unnormalized = ggplot(data=df,
           aes(x=condition_expr,
               # y=log2FC_enrichment)) +
               y=condition_enrichment)) +
        stat_binhex(geom="point",
                    aes(color=..count..),
                    bins=150,
                    size=0.2,
                    shape=16,
                    alpha=0.8) +
        geom_smooth(size=0.2,
                    color="orange") +
        # facet_wrap(~strain,
                   # ncol=1) +
        facet_grid(.~strain) +
        scale_color_viridis_c() +
        scale_x_log10("transcript abundance") +
        scale_y_continuous(limits=c(-1.5, NA),
                           oob=scales::squish,
                           name=quote("ZF: log"[2] ~ frac("IP", "input"))) +
        theme_default +
        theme(legend.position="none",
              panel.grid=element_blank(),
              axis.title.y=element_text(angle=0,
                                        vjust=0.5))

    normalized = ggplot(data=df,
           aes(x=condition_expr,
               y=log2FC_enrichment)) +
               # y=condition_enrichment)) +
        stat_binhex(geom="point",
                    aes(color=..count..),
                    bins=150,
                    size=0.2,
                    shape=16,
                    alpha=0.8) +
        geom_smooth(size=0.2,
                    color="orange") +
        facet_grid(.~strain) +
        # facet_wrap(~strain,
                   # ncol=1) +
        scale_color_viridis_c() +
        scale_x_log10("transcript abundance") +
        scale_y_continuous(limits=c(-1, NA),
                           oob=scales::squish,
                           name=quote("ZF: log"[2] ~ frac("strain", "reporter-only"))) +
        theme_default +
        theme(legend.position="none",
              panel.grid=element_blank(),
              axis.title.y=element_text(angle=0,
                                        vjust=0.5))

    plot = plot_grid(unnormalized,
                     normalized,
                     ncol=1,
                     align="v",
                     axis="lr")

    ggsave(pdf_out,
           plot=plot,
           width=16,
           height=9,
           units="cm",
           device=cairo_pdf)
}

main(theme_path = snakemake@input[["theme"]],
     high_rna_path = snakemake@input[["high_rna"]],
     high_chip_path = snakemake@input[["high_chip"]],
     low_rna_path = snakemake@input[["low_rna"]],
     low_chip_path = snakemake@input[["low_chip"]],
     clamp_rna_path = snakemake@input[["clamp_rna"]],
     clamp_chip_path = snakemake@input[["clamp_chip"]],
     pdf_out= snakemake@output[["pdf"]])

