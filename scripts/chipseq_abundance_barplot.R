
main = function(theme_path = "custom_theme.R",
                data_path = "ZF-chipseq_spikein-counts.tsv",
                pdf_out="test.pdf",
                # grob_out="test.Rdata",
                fig_width=4,
                fig_height=6,
                panel_letter="",
                plot_title="ZF ChIP-seq"){
    source(theme_path)

    df = read_tsv(data_path) %>%
        mutate(abundance = (experimental_counts_IP / spikein_counts_IP) *
                   (spikein_counts_input / experimental_counts_input),
               group=ordered(group,
                             levels=c('reporter-only',
                                      'high-affinity-ZF',
                                      'low-affinity-ZF',
                                      'low-affinity-ZF-with-clamp',
                                      'high-affinity-ZF-noVP16',
                                      'low-affinity-ZF-noVP16',
                                      'low-affinity-ZF-with-clamp-noVP16')))

    baseline_mean = df %>%
        filter(group == "reporter-only") %>%
        summarize(mean=mean(abundance)) %>%
        pull(mean)

    df %<>%
        mutate(scaled_abundance=scales::rescale(abundance,
                                                from=c(0, baseline_mean)))

    df_summary = df %>%
        group_by(group) %>%
        summarize(mean_scaled_abundance = mean(scaled_abundance),
               sd_scaled_abundance = sd(scaled_abundance)) %>%
        mutate(label_direction = if_else(mean_scaled_abundance > 0.5,
                                         -1,
                                         1))

    chipseq_abundance_barplot = ggplot() +
        geom_col(data=df_summary,
                 aes(x=group,
                     y=mean_scaled_abundance,
                     fill=group),
                 alpha=0.5,
                 width=0.65) +
        geom_errorbar(data=df_summary,
                      aes(x=group,
                          ymin=mean_scaled_abundance - sd_scaled_abundance,
                          ymax=mean_scaled_abundance + sd_scaled_abundance),
                      width=0.2,
                      size=0.3,
                      alpha=0.8) +
        geom_jitter(data=df,
                    aes(x=group,
                        y=scaled_abundance),
                    width=0.2,
                    height=0,
                    shape=16,
                    size=0.7,
                    alpha=0.8) +
        geom_text(data=df_summary,
                  aes(x=group,
                      y=mean_scaled_abundance + (sd_scaled_abundance + 0.35) * label_direction,
                      label=paste0(sprintf("%.2f", mean_scaled_abundance),
                                   "\n±",
                                  sprintf("%.2f", sd_scaled_abundance))),
                  size=5/72*25.4,
                  family="FreeSans") +
        scale_x_discrete(name=NULL,
                         expand=c(0,0.5)) +
        scale_y_continuous(limits=function(x) c(0, x[2] * 1.05),
                           breaks=scales::pretty_breaks(5),
                           expand=c(0,0),
                           name="total signal") +
        scale_fill_tableau(guide=FALSE) +
        labs(tag=panel_letter,
             title=plot_title) +
        theme_default +
        theme(panel.grid=element_blank(),
              axis.text.x=element_text(angle=20,
                                       hjust=0.9),
              axis.ticks.x=element_blank(),
              plot.title=element_text(hjust=0))

    ggsave(pdf_out,
           plot=chipseq_abundance_barplot,
           width=fig_width,
           height=fig_height,
           units="cm",
           device=cairo_pdf)
    # save(chipseq_abundance_barplot,
    #      file=grob_out)
}

main(theme_path=snakemake@input[["theme"]],
     data_path=snakemake@input[["data"]],
     pdf_out=snakemake@output[["pdf"]],
     # grob_out=snakemake@output[["grob"]],
     # fig_width=snakemake@params[["fig_width"]],
     fig_width=12,
     fig_height=12,
     panel_letter="",
     plot_title="ZF ChIP-seq")

