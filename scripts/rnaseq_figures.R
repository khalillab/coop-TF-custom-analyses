library(tidyverse)
library(magrittr)
# library(ggpmisc)
# library(ggrepel)

import = function(path, strain_id){
    read_tsv(path) %>%
        mutate(strain=strain_id) %>%
        return()
}

main = function(input_high_affinity,
                input_low_affinity,
                input_low_affinity_w_clamp,
                fdr,
                output_volcano,
                output_maplot,
                output_summary){

    df = import(input_high_affinity,
                "high-affinity ZF") %>%
        bind_rows(import(input_low_affinity,
                         "low-affinity ZF")) %>%
        bind_rows(import(input_low_affinity_w_clamp,
                         "low-affinity ZF with clamp")) %>%
        mutate(significant = (log10_padj > -log10(fdr)))

    volcano = ggplot() +
        geom_vline(xintercept=0,
                   color="grey80",
                   size=0.5) +
        geom_point(data = df,
                   aes(x=log2_foldchange,
                       y=log10_padj,
                       alpha=significant,
                       color=significant),
                   size=0.75,
                   shape=16) +
        facet_grid(.~strain) +
        scale_y_continuous(limits=c(0, NA),
                           expand=c(0, 2),
                           name=expression("-log"[10] * "FDR")) +
        scale_x_continuous(name=expression("log"[2] ~ bgroup("[",
                                                             atop("fold-change vs.",
                                                                  "reporter-only"),
                                                             "]"))) +
        scale_alpha_manual(values=c(0.2, 0.6)) +
        scale_color_manual(values=c("grey40", "#440154FF")) +
        theme_light() +
        theme(text=element_text(color="black",
                                size=10),
              strip.background=element_blank(),
              strip.text=element_text(color="black"),
              legend.position="none",
              axis.text=element_text(color="black"),
              axis.title=element_text(size=10),
              axis.title.y=element_text(angle=0, vjust=0.5),
              panel.grid=element_blank())

    ggsave(output_volcano,
           plot=volcano,
           width=16,
           height=9,
           units="cm")

    maplot = ggplot() +
        geom_hline(yintercept=0,
                   color="grey80") +
        geom_point(data = df,
                   aes(x=mean_expr,
                       y=log2_foldchange,
                       alpha=significant,
                       color=significant),
                   # size=0.4,
                   size=1,
                   shape=16) +
        facet_grid(.~strain) +
        scale_x_log10(#limits=c(0, NA),
                           #expand=c(0, 2),
                           name="mean expression",
                           breaks=c(1e1, 1e3, 1e5),
                           labels=c(expression(10^1),
                                    expression(10^3),
                                    expression(10^5))
                           # limits=c(0, 3e5),
                           # oob=scales::squish
                           ) +
        scale_y_continuous(name=expression("log"[2] ~ bgroup("[",
                                                             atop("fold-change vs.",
                                                                  "reporter-only"),
                                                             "]"))) +
        scale_alpha_manual(values=c(0.2, 0.6)) +
        scale_color_manual(values=c("grey40", "#440154FF")) +
        theme_light() +
        theme(text=element_text(color="black",
                                size=10),
              strip.background=element_blank(),
              strip.text=element_text(color="black"),
              legend.position="none",
              axis.text=element_text(color="black"),
              axis.title=element_text(size=10),
              axis.title.y=element_text(angle=0, vjust=0.5),
              panel.grid=element_blank())

    ggsave(output_maplot,
           plot=maplot,
           width=16,
           height=7,
           units="cm")

    # plot = ggplot() +
    #     geom_point(data = df %>%
    #                    filter(! significant),
    #                aes(x=log2_foldchange,
    #                    y=log10_padj),
    #                alpha=0.4,
    #                size=0.4,
    #                color="grey40") +
    #     geom_point(data = df %>%
    #                    filter(significant),
    #                aes(x=log2_foldchange,
    #                    y=log10_padj,
    #                    color=(log2_foldchange > 0)),
    #                size=0.4,
    #                alpha=0.8) +
    #     stat_dens2d_labels(data = df %>% filter(significant),
    #                        aes(x=log2_foldchange,
    #                            y=log10_padj,
    #                            label=name,
    #                            fill=(log2_foldchange > 0)),
    #                        geom="label_repel",
    #                        h=c(0.5, 20),
    #                        size=6/72*25.4,
    #                        box.padding=unit(0, "pt"),
    #                        label.r=unit(0.5, "pt"),
    #                        label.size=NA,
    #                        label.padding=unit(0.8, "pt"),
    #                        ylim=c(-log10(fdr), NA),
    #                        force=0.5,
    #                        segment.size=0.2,
    #                        segment.alpha=0.5) +
    #     facet_grid(.~strain) +
    #     scale_y_continuous(limits=c(0, NA),
    #                        expand=c(0, 2),
    #                        name=expression("-log"[10] * "FDR")) +
    #     scale_x_continuous(name=expression("log"[2] ~ bgroup("[",
    #                                                          atop("fold-change vs.",
    #                                                               "reporter-only"),
    #                                                          "]"))) +
    #     theme(legend.position = "none",
    #           axis.title.y = element_text(angle=0, vjust=0.5))

    df %<>%
        group_by(chrom, start, end, name, strand) %>%
        filter(any(significant)) %>%
        ungroup()

    gene_order = df %>%
        filter(strain == "high-affinity ZF") %>%
        arrange(log2_foldchange) %>%
        pull(name)

    df %<>%
        mutate(name = ordered(name, levels=gene_order))

    summary_plot = ggplot(data = df,
           aes(x=name,
               y=log2_foldchange,
               ymin=log2_foldchange - lfc_SE,
               ymax=log2_foldchange + lfc_SE,
               color=strain)) +
        geom_hline(yintercept = 0,
                   color="grey70",
                   size=0.2) +
        geom_linerange(alpha=0.4,
                       size=0.4) +
        geom_point(size=0.1,
                   alpha=0.95) +
        scale_color_viridis_d(end=0.85,
                              name=NULL) +
        scale_x_discrete(expand=c(0, 5),
                         name="genes differentially expressed vs. reporter-only") +
        scale_y_continuous(name=expression("log"[2] ~ bgroup("[",
                                                             atop("fold-change vs.",
                                                                  "reporter-only"),
                                                             "]")),
                           breaks=scales::pretty_breaks(6)) +
        annotate(geom="text",
                 x=300,
                 y=9.355,
                 vjust=0.5,
                 hjust=1,
                 label="reporter",
                 size=8/72*25.4,
                 family="FreeSans") +
        annotate(geom="segment",
                 x=302,
                 xend=317,
                 y=9.355,
                 yend=9.355,
                 size=0.3,
                 arrow=arrow(type="closed",
                             length=unit(3, "pt"))) +
        theme_light() +
        theme(text=element_text(size=8, family="FreeSans"),
              legend.position=c(0.5,0.75),
              axis.text.y=element_text(color="black", size=8),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank(),
              axis.title=element_text(size=8),
              axis.title.y=element_text(angle=0,
                                        vjust=0.5),
              legend.text=element_text(size=8),
              legend.spacing.x=unit(-2, "pt"),
              panel.grid=element_blank())

    ggsave(output_summary,
           plot=summary_plot,
           width=16,
           height=9,
           units="cm",
           device=cairo_pdf)
}

main(input_high_affinity=snakemake@input[["high_affinity"]],
     input_low_affinity=snakemake@input[["low_affinity"]],
     input_low_affinity_w_clamp=snakemake@input[["low_affinity_w_clamp"]],
     fdr=snakemake@params[["fdr"]],
     output_volcano=snakemake@output[["volcano"]],
     output_maplot=snakemake@output[["maplot"]],
     output_summary=snakemake@output[["summary"]])

