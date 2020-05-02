
main = function(theme_path = "custom_theme.R",
                target_annotation = "Venus.bed",
                transcript_annotation_path = "Scer_transcripts_w_verifiedORFs-withVenus-URA3_reporter_strain.bed",
                orf_annotation_path = "Scer_nondubious_ORFs-withVenus-URA3_reporter_strain.bed",
                motif_annotation_path = "42-10_9bp_fimo_7e4.bed",
                data_path = "Venus_all-assays.tsv.gz",
                plot_out = "test.pdf"){
    source(theme_path)
    library(ggridges)
    library(cowplot)

    df = read_tsv(data_path,
                  col_names=c("group", "sample", "annotation",
                              "assay", "index", "position", "signal")) %>%
        group_by(group, assay, position) %>%
        summarize(signal = mean(signal)) %>%
        ungroup %>%
        mutate(group=ordered(group,
                             levels=c("reporter-only",
                                      "high-affinity-ZF",
                                      "low-affinity-ZF",
                                      "low-affinity-ZF-with-clamp"),
                             labels=c("reporter only",
                                      "high-affinity",
                                      "low-affinity",
                                      "low-affinity + clamp")))

    df_chip = df %>%
        filter(assay=="ZF-ChIPseq")

    df_rnaseq = df %>%
        filter(assay %in% c("RNA-seq-sense-wholeread",
                            "RNA-seq-antisense-wholeread")) %>%
        ungroup() %>%
        mutate(assay = str_remove(assay, "RNA-seq-"),
               assay = str_remove(assay, "-wholeread")) %>%
        pivot_wider(names_from=assay,
                    values_from=signal)

    x_limits = range(df[["position"]])

    motif_annotations = read_tsv(target_annotation,
             col_names=c("chrom", "left", "right", "name", "score", "strand")) %>%
        select(chrom, left, right, strand) %>%
        left_join(read_tsv(motif_annotation_path,
                           col_names=c("chrom", "motif_left", "motif_right", "motif_name",
                                       "motif_score", "motif_strand")),
                  by=c("chrom")) %>%
        mutate_at(vars(motif_left, motif_right),
                  ~((. - ifelse(strand=="+", left, right)) * ifelse(strand=="+", 1, -1) / 1e3)) %>%
        mutate(m_left = ifelse(strand=="+", motif_left, motif_right),
               m_right = ifelse(strand=="+", motif_right, motif_left),
               motif_left = m_left,
               motif_right = m_right) %>%
        select(-c(m_left, m_right)) %>%
        filter(between(motif_left, x_limits[1], x_limits[2]) |
                   between(motif_right, x_limits[1], x_limits[2]) |
                   (motif_left <= x_limits[1] & motif_right >= x_limits[2])) %>%
        mutate_at(vars(motif_left, motif_right),
                  ~pmax(x_limits[1], pmin(x_limits[2], .)))

    transcript_annotations = read_tsv(target_annotation,
             col_names=c("chrom", "left", "right", "name", "score", "strand")) %>%
        select(chrom, left, right, strand) %>%
        left_join(read_tsv(transcript_annotation_path,
                           col_names=c("chrom", "transcript_left", "transcript_right", "transcript_name",
                                       "score", "transcript_strand")) %>%
                      select(-score),
                  by=c("chrom")) %>%
        mutate_at(vars(transcript_left, transcript_right),
                  ~((. - ifelse(strand=="+", left, right)) * ifelse(strand=="+", 1, -1) / 1e3)) %>%
        mutate(t_left = ifelse(strand=="+", transcript_left, transcript_right),
               t_right = ifelse(strand=="+", transcript_right, transcript_left),
               transcript_left = t_left,
               transcript_right = t_right) %>%
        select(-c(t_left, t_right)) %>%
        filter(between(transcript_left, x_limits[1], x_limits[2]) |
                   between(transcript_right, x_limits[1], x_limits[2]) |
                   (transcript_left <= x_limits[1] & transcript_right >= x_limits[2])) %>%
        left_join(read_tsv(orf_annotation_path,
                           col_names=c("chrom", "orf_left", "orf_right", "transcript_name",
                                       "score", "transcript_strand")) %>%
                      select(-score),
                  by=c("chrom", "transcript_name", "transcript_strand")) %>%
        mutate_at(vars(orf_left, orf_right),
                  ~((. - ifelse(strand=="+", left, right)) * ifelse(strand=="+", 1, -1) / 1e3)) %>%
        mutate(o_left = ifelse(strand=="+", orf_left, orf_right),
               o_right = ifelse(strand=="+", orf_right, orf_left),
               orf_left = o_left,
               orf_right = o_right) %>%
        select(-c(o_left, o_right)) %>%
        mutate(transcript_start = if_else(transcript_strand==strand,
                                          transcript_left,
                                          transcript_right),
               transcript_end = if_else(transcript_strand==strand,
                                        transcript_right,
                                        transcript_left),
               notch = ifelse(transcript_strand==strand,
                              (orf_right - orf_left) * 0.85 + orf_left,
                              orf_right - (orf_right - orf_left) * 0.85),
               orf_start = if_else(transcript_strand==strand,
                                   orf_left,
                                   orf_right),
               orf_end = if_else(transcript_strand==strand,
                                 orf_right,
                                 orf_left),
               label_xpos = if_else(is.na(notch),
                                    (transcript_start + transcript_end) / 2,
                                    (orf_start + orf_end) / 2)) %>%
        select(transcript_name,
               transcript_start,
               transcript_end,
               orf_start,
               orf_end,
               notch,
               label_xpos) %>%
        mutate_at(vars(transcript_start, transcript_end),
                  ~pmax(x_limits[1], pmin(x_limits[2], .))) %>%
        mutate(transcript_y = (row_number() %% 2) * 1.3)

    orf_annotations = transcript_annotations %>%
        filter(! is.na(notch)) %>%
        expand(nesting(transcript_name,
                       transcript_start,
                       transcript_end,
                       orf_start,
                       notch,
                       orf_end,
                       label_xpos),
               orf_y=c(-1, 1)) %>%
        pivot_longer(c(orf_start, notch, orf_end),
                     names_to="type",
                     values_to="orf_x") %>%
        mutate(orf_y = ifelse(type=="orf_end",
                              0,
                              orf_y)) %>%
        distinct() %>%
        group_by(transcript_name) %>%
        arrange(type,
                orf_y,
                .by_group=TRUE) %>%
        mutate(order = c(4,2,3,5,1)) %>%
        arrange(order,
                .by_group=TRUE) %>%
        mutate(orf_x = pmax(x_limits[1], pmin(x_limits[2], orf_x))) %>%
        left_join(transcript_annotations %>%
                      select(transcript_name, transcript_y),
                  by="transcript_name") %>%
        mutate(orf_y = orf_y + transcript_y)

    annotation_plot = ggplot() +
        geom_rect(data=motif_annotations,
                  aes(xmin=motif_left,
                      xmax=motif_right,
                      alpha=motif_score),
                  ymin=-Inf,
                  ymax=Inf,
                  fill="gray70") +
        scale_alpha_continuous(range=c(0.2, 0.6),
                               guide=FALSE) +
        geom_segment(data=transcript_annotations,
                     aes(x=transcript_start,
                         xend=transcript_end,
                         y=transcript_y,
                         yend=transcript_y),
                     color="gray50") +
        geom_polygon(data=orf_annotations,
                     aes(x=orf_x,
                         y=orf_y,
                         group=transcript_name),
                     fill="grey80") +
        geom_text(data=transcript_annotations,
                  aes(x=label_xpos,
                      y=transcript_y,
                      label=transcript_name,
                      vjust=if_else(is.na(notch),
                                    -0.2, 0.5)),
                  size=8/72*25.4,
                  family="FreeSans",
                  fontface="italic") +
        geom_vline(xintercept=x_limits[2],
                   color="gray70",
                   size=0.4) +
        scale_x_continuous(expand=c(0,0),
                           limits=x_limits) +
        scale_y_continuous(expand=c(0.2, 0)) +
        theme_void() +
        theme(axis.line.y=element_line(color="gray70",
                                       size=0.2))

    rnaseq_plot = ggplot() +
        geom_vline(data = transcript_annotations %>%
                       filter(transcript_start == 0) %>%
                       pivot_longer(cols=c(transcript_start, transcript_end),
                                    names_to="type",
                                    values_to="x"),
                   aes(xintercept=x),
                   color="gray70",
                   size=0.2) +
        # geom_rect(data=motif_annotations,
        #           aes(xmin=motif_left,
        #               xmax=motif_right,
        #               alpha=motif_score),
        #           ymin=-Inf,
        #           ymax=Inf,
        #           fill="gray70") +
        # scale_alpha_continuous(range=c(0.2, 0.4),
        #                        guide=FALSE) +
        geom_line(data = df_rnaseq,
                        aes(x=position,
                            color=group,
                            y=sense),
                  alpha=0.9) +
        geom_line(data = df_rnaseq,
                        aes(x=position,
                            color=group,
                            y=-antisense),
                  alpha=0.9) +
        geom_hline(yintercept=0,
                   color="gray70",
                   size=0.2) +
        scale_x_continuous(expand=c(0,0),
                           labels=function(x) ifelse(x==0, "TSS", paste(x, "kb")),
                           name=NULL) +
        scale_y_continuous(name="RNA-seq\n signal") +
        scale_color_ptol(name=NULL) +
        theme_default +
        theme(panel.grid=element_blank(),
              axis.title.y=element_text(angle=0,
                                        vjust=0.5),
              legend.justification=c(0,1),
              legend.position=c(0.03,0.80),
              legend.background=element_blank(),
              legend.key.height=unit(10, "pt"),
              legend.spacing.y=unit(0, "pt"))

    chip_plot = ggplot() +
        geom_vline(data = transcript_annotations %>%
                       filter(transcript_start == 0) %>%
                       pivot_longer(cols=c(transcript_start, transcript_end),
                                    names_to="type",
                                    values_to="x"),
                   aes(xintercept=x),
                   color="gray70",
                   size=0.2) +
        # geom_rect(data=motif_annotations,
        #           aes(xmin=motif_left,
        #               xmax=motif_right,
        #               alpha=motif_score),
        #           ymin=-Inf,
        #           ymax=Inf,
        #           fill="gray70") +
        # scale_alpha_continuous(range=c(0.2, 0.4),
        #                        guide=FALSE) +
        geom_line(data=df_chip,
                       aes(x=position,
                           y=signal,
                           color=group)) +
        geom_vline(xintercept=max(df_chip[["position"]]),
                   color="gray70",
                   size=0.2) +
        scale_x_continuous(expand=c(0,0),
                           name=NULL) +
        scale_y_continuous(breaks=scales::pretty_breaks(4),
                           name="ZFTF ChIP\nenrichment") +
        scale_color_ptol(guide=FALSE) +
        theme_default +
        theme(axis.text.x=element_blank(),
              axis.title.y=element_text(angle=0,
                                        vjust=0.5),
              panel.grid=element_blank())

    summary_plot = plot_grid(annotation_plot,
                     chip_plot,
                     rnaseq_plot,
                     align="v",
                     axis="lr",
                     ncol=1,
                     rel_heights=c(0.15,1,1))

    ggsave(plot_out,
           plot=summary_plot,
           width=16,
           height=9,
           units="cm",
           device=cairo_pdf)
}

main(theme_path = snakemake@input[["theme"]],
     target_annotation = snakemake@input[["target_annotation"]],
     transcript_annotation_path = snakemake@input[["transcript_annotation"]],
     orf_annotation_path = snakemake@input[["orf_annotation"]],
     motif_annotation_path = snakemake@input[["motif_annotation"]],
     data_path = snakemake@input[["data"]],
     plot_out = snakemake@output[["plot"]])

