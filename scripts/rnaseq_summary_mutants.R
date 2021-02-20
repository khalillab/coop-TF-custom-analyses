
import = function(path, strain_id){
    read_tsv(path) %>%
        mutate(strain=strain_id) %>%
        return()
}

main = function(theme_path = "custom_theme.R",
                input_high_affinity = "high-affinity-ZF-v-reporter-only_rnaseq-libsizenorm-verified_genes_plus_Venus-diffexp-results-all.tsv",
                input_low_affinity = "low-affinity-ZF-v-reporter-only_rnaseq-libsizenorm-verified_genes_plus_Venus-diffexp-results-all.tsv",
                input_low_affinity_w_clamp = "low-affinity-ZF-with-clamp-v-reporter-only_rnaseq-libsizenorm-verified_genes_plus_Venus-diffexp-results-all.tsv",
                input_high_affinity_zev_mutant = "high-affinity-ZF-ZEV-mutant-v-reporter-only_rnaseq-libsizenorm-verified_genes_plus_Venus-diffexp-results-all.tsv",
                input_high_affinity_zftf_mutant = "high-affinity-ZF-ZFTF-mutant-v-reporter-only_rnaseq-libsizenorm-verified_genes_plus_Venus-diffexp-results-all.tsv",
                fdr=0.1,
                panel_letter = "A",
                fig_width = 11.4,
                fig_height= 6,
                pdf_out="test.pdf",
                grob_out="test.Rdata"){

    source(theme_path)

    df = import(input_high_affinity,
                "high-affinity") %>%
        # bind_rows(import(input_low_affinity,
        #                  "low-affinity")) %>%
        # bind_rows(import(input_low_affinity_w_clamp,
        #                  "low-affinity + clamp")) %>%
        bind_rows(import(input_high_affinity_zev_mutant,
                         "high-affinity ZEV mutant")) %>%
        bind_rows(import(input_high_affinity_zftf_mutant,
                         "high-affinity ZFTF mutant")) %>%
        mutate(significant = (log10_padj > -log10(fdr))) %>%
        group_by(chrom, start, end, name, strand) %>%
        filter(any(significant)) %>%
        ungroup()

    gene_order = df %>%
        filter(strain == "high-affinity") %>%
        arrange(log2_foldchange) %>%
        pull(name)

    df %<>%
        mutate(name = ordered(name, levels=gene_order))

    reporter_coord_x = df %>%
        filter(name=="Venus") %>%
        slice(1) %>%
        pull(name) %>%
        as.integer()

    rnaseq_summary = ggplot(data = df,
           aes(x=name,
               y=log2_foldchange,
               ymin=log2_foldchange - lfc_SE,
               ymax=log2_foldchange + lfc_SE,
               color=strain)) +
        geom_vline(xintercept=reporter_coord_x,
                   size=0.4,
                   color="gray70",
                   alpha=0.4) +
        geom_hline(yintercept = 0,
                   color="grey70",
                   size=0.2) +
        geom_linerange(alpha=0.4,
                       size=0.2) +
        geom_point(size=0.1,
                   alpha=0.8,
                   shape=16) +
        scale_color_viridis_d(end=0.85,
                              name=NULL,
                              guide=guide_legend(override.aes = list(size=0.5))) +
        # scale_color_manual(values=c("#2b398d",
        #                             "#91bfda",
        #                             "#d23c28",
        #                             "red",
        #                             "orange"),
        #                    name=NULL,
        #                    guide=guide_legend(override.aes = list(size=0.5))) +
        scale_x_discrete(expand=c(0, 2),
                         name="genes differentially expressed vs. reporter-only") +
        scale_y_continuous(name=expression("log"[2] ~ bgroup("[",
                                                             textstyle(atop("fold-change vs.",
                                                                  "reporter-only")),
                                                             "]")),
                           breaks=scales::pretty_breaks(6)) +
        annotate(geom="text",
                 x=reporter_coord_x-7,
                 y=9.375,
                 vjust=0.45,
                 hjust=1,
                 label="reporter",
                 size=7/72*25.4,
                 family="FreeSans") +
        annotate(geom="segment",
                 x=reporter_coord_x-7.7,
                 xend=reporter_coord_x-2.5,
                 y=9.375,
                 yend=9.375,
                 size=0.3,
                 arrow=arrow(type="open",
                             length=unit(2, "pt"))) +
        labs(tag=panel_letter) +
        theme_default +
        theme(panel.grid=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank(),
              axis.title.y=element_text(angle=0,
                                        vjust=0.5),
              legend.position=c(0.5,0.75),
              legend.spacing.x=unit(-4, "pt"))

    ggsave(pdf_out,
           plot=rnaseq_summary,
           width=fig_width,
           height=fig_height,
           units="cm",
           device=cairo_pdf)
    # save(rnaseq_summary,
    #      file=grob_out)
}

main(theme_path = snakemake@input[["theme"]],
     input_high_affinity=snakemake@input[["high_affinity"]],
     input_low_affinity=snakemake@input[["low_affinity"]],
     input_low_affinity_w_clamp=snakemake@input[["low_affinity_w_clamp"]],
     input_high_affinity_zev_mutant = snakemake@input[["high_affinity_zev_mutant"]],
     input_high_affinity_zftf_mutant = snakemake@input[["high_affinity_zftf_mutant"]],
     fdr=snakemake@params[["fdr"]],
     panel_letter=snakemake@params[["panel_letter"]],
     fig_width=11.4,
     fig_height=6,
     pdf_out= snakemake@output[["pdf"]],
     grob_out=snakemake@output[["grob"]])
