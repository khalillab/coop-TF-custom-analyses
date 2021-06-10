
import = function(path, strain_id){
    read_tsv(path) %>%
        mutate(strain=strain_id) %>%
        return()
}
sub_distance = function(a,b){
    return(mapply(function(x,y) sum(x!=y),strsplit(a,""),strsplit(b,"")))
}

main = function(theme_path = "custom_theme.R",
                input_high_affinity = "high-affinity-ZF-v-reporter-only_rnaseq-libsizenorm-verified_genes_plus_Venus-diffexp-results-all.tsv",
                input_low_affinity = "low-affinity-ZF-v-reporter-only_rnaseq-libsizenorm-verified_genes_plus_Venus-diffexp-results-all.tsv",
                input_low_affinity_w_clamp = "low-affinity-ZF-with-clamp-v-reporter-only_rnaseq-libsizenorm-verified_genes_plus_Venus-diffexp-results-all.tsv",
                motif_results = "high-affinity-upregulated-transcripts-v-all-transcripts_control_unmergedFIMOresults.tsv.gz",
                fdr=0.1,
                # panel_letter = "A",
                fig_width = 17.4,
                fig_height= 6,
                pdf_out="test.pdf"){

    source(theme_path)
    library(ggforce)
    library(cowplot)

    motif_hits = read_tsv(motif_results) %>%
        filter(motif_start >= 0) %>%
        mutate(n_mismatch = sub_distance(match_sequence, 'GACGCTGCT'))

    chip_hits = c("Venus",
                  "HMS1",
                  "AMF1",
                  "PRY1",
                  "REB1",
                  "VRP1",
                  "FTR1",
                  "DAN1",
                  "PAC11",
                  "OCH1",
                  "YDL144C",
                  "SUR2")

    df = import(input_high_affinity,
                "high-affinity") %>%
        bind_rows(import(input_low_affinity,
                         "low-affinity")) %>%
        bind_rows(import(input_low_affinity_w_clamp,
                         "low-affinity + clamp")) %>%
        mutate(significant = (log10_padj > -log10(fdr))) %>%
        group_by(chrom, start, end, name, strand) %>%
        filter(any(significant)) %>%
        ungroup()

    gene_order = df %>%
        filter(strain == "high-affinity") %>%
        arrange(log2_foldchange) %>%
        pull(name)

    df %<>%
        mutate(name = ordered(name, levels=gene_order),
               motif_hit = name %in% (motif_hits %>% filter(n_mismatch < 2) %>% pull(region_id)),
               chip_hit = name %in% chip_hits)

    reporter_coord_x = df %>%
        # filter(name=="Venus") %>%
        filter(name %in% chip_hits) %>%
        distinct(name, .keep_all=TRUE) #%>%
        # pull(name) %>%
        # as.integer()

    rnaseq_summary = ggplot(data = df,
                            aes(x=as.numeric(name),
                                y=log2_foldchange,
                                # ymin=log2_foldchange - lfc_SE,
                                # ymax=log2_foldchange + lfc_SE,
                                color=strain,
                                fill=strain)) +
        geom_vline(xintercept=as.integer(reporter_coord_x[['name']]),
                   size=0.4,
                   color="gray70",
                   alpha=0.4) +
        geom_text(data=reporter_coord_x,
                  aes(x=as.integer(name),
                      label=name),
                  color='black',
                  y=-2,
                  angle=90,
                  family='FreeSans',
                  size=5/72*25.4) +
        geom_hline(yintercept = 0,
                   color="grey70",
                   size=0.2) +
        # geom_linerange(alpha=0.4,
                       # size=0.2) +
        geom_point(#aes(shape=motif_hit),
                   size=0.5,
                   alpha=0.8,
                   stroke=0.1) +
        facet_zoom(xlim=c(114, 214),
                   zoom.size=1) +
        scale_fill_viridis_d(end=0.85,
                              name=NULL,
                              guide=guide_legend(override.aes = list(size=1))) +
        scale_color_viridis_d(end=0.85,
                              name=NULL,
                              guide=guide_legend(override.aes = list(size=1))) +
        # scale_shape_manual(values=c(1, 21),
                           # guide=FALSE) +
        scale_x_continuous(expand=c(0.01, 0),
                           name="genes differentially expressed vs. reporter-only") +
        scale_y_continuous(name=expression("log"[2] ~ bgroup("[",
                                                             textstyle(atop("fold-change vs.",
                                                                  "reporter-only")),
                                                             "]")),
                           breaks=scales::pretty_breaks(6)) +
        # annotate(geom="text",
        #          x=reporter_coord_x-7,
        #          y=9.375,
        #          vjust=0.45,
        #          hjust=1,
        #          label="reporter",
        #          size=7/72*25.4,
        #          family="FreeSans") +
        # annotate(geom="segment",
        #          x=reporter_coord_x-7.7,
        #          xend=reporter_coord_x-2.5,
        #          y=9.375,
        #          yend=9.375,
        #          size=0.3,
        #          arrow=arrow(type="open",
        #                      length=unit(2, "pt"))) +
        # labs(tag=panel_letter) +
        theme_default +
        theme(panel.grid=element_blank(),
              axis.text.x=element_blank(),
              axis.text.y=element_text(size=5),
              axis.ticks.x=element_blank(),
              axis.title.y=element_text(angle=0,
                                        vjust=0.5),
              legend.position=c(0.44,0.9),
              legend.spacing.x=unit(-5, "pt"),
              legend.background=element_blank(),
              strip.background = element_rect(fill="gray90",
                                              linetype="blank"),
              panel.spacing=unit(1, "pt"))

    density = ggplot(data=mutate(df, test="a"),
           aes(x=log2_foldchange,
               color=strain,
               fill=strain)) +
        coord_flip() +
        geom_vline(xintercept=0,
                   size=0.2,
                   color="gray70") +
        geom_density(size=0.2,
                     alpha=0.1) +
        facet_grid(test~.,
                   margins=TRUE) +
        scale_y_continuous(limits=function(x) c(0, x[2] * 1.05),
                           expand=c(0.01,0)) +
        scale_fill_viridis_d(end=0.85,
                              name=NULL,
                              guide=FALSE) +
        scale_color_viridis_d(end=0.85,
                              name=NULL,
                              guide=FALSE) +
        theme_default +
        theme(panel.grid=element_blank(),
              axis.title.y=element_blank(),
              axis.text=element_blank(),
              axis.ticks=element_blank(),
              strip.text=element_blank(),
              panel.spacing.y=unit(7.7, "pt"),
              plot.margin=margin(l=-6, r=2, unit="pt"))


    plots = plot_grid(rnaseq_summary,
                     density,
                     nrow=1,
                     align="h",
                     axis="tb",
                     rel_widths=c(1, 0.1))

    ggsave(pdf_out,
           plot=plots,
           width=fig_width,
           height=fig_height,
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
     pdf_out= snakemake@output[["pdf"]])
