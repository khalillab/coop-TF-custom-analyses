library(tidyverse)

main = function(input_coverage='affinity-dependent-peaks-filtered-allsamples-allannotations-ZF-chipseq-spikenorm-ratio.tsv.gz',
                input_summit_annotation='peaks_cluster-1.bed',
                input_transcript_annotation='Scer_transcripts_w_verifiedORFs-withVenus-URA3_reporter_strain.bed',
                input_orf_annotation='Scer_nondubious_ORFs-withVenus-URA3_reporter_strain.bed',
                motif_annotation='allmotifs.bed',
                filter_groups = c('high-affinity-ZF-noVP16',
                                  'low-affinity-ZF-noVP16',
                                  'low-affinity-ZF-with-clamp-noVP16'),
                output_path='test.pdf'){

    df = read_tsv(input_coverage,
                  col_names=c("group",
                              "sample",
                              "type",
                              "annotation",
                              "index",
                              "position",
                              "signal")) %>%
        filter(group %in% filter_groups) %>%
        group_by(group,
                 type,
                 annotation,
                 index,
                 position) %>%
        summarize(signal=mean(signal)) %>%
        ungroup() %>%
        mutate(group=ordered(group,
                             levels=c("high-affinity-ZF",
                                      "high-affinity-ZF-noVP16",
                                      "low-affinity-ZF",
                                      "low-affinity-ZF-noVP16",
                                      "low-affinity-ZF-with-clamp",
                                      "low-affinity-ZF-with-clamp-noVP16",
                                      "reporter-only"),
                             labels=c("high-affinity",
                                      "high-affinity, no VP16",
                                      "low-affinity",
                                      "low-affinity, no VP16",
                                      "coop. assembly",
                                      "coop. assembly, no VP16",
                                      "reporter-only")))

    # #%>%
    #     # mutate(signal=scales::rescale(signal)) %>%
    #     spread(key=group,
    #            value=signal) %>%
    #     mutate_at(c("high-affinity-ZF",
    #                 "low-affinity-ZF",
    #                 "low-affinity-ZF-with-clamp"),
    #               # ~(./(`reporter-only`+0.01))) %>%
    #               ~(.-`reporter-only`)) %>%
    #     select(-`reporter-only`) %>%
    #     gather(key=group,
    #            value=signal,
    #            -c(annotation,
    #               index,
    #               position)) %>%
    #     mutate(group=ordered(group,
    #                          levels=c("high-affinity-ZF",
    #                                   "low-affinity-ZF",
    #                                   "low-affinity-ZF-with-clamp"),
    #                          labels=c("high-affinity ZF",
    #                                   "low-affinity ZF",
    #                                   "low-affinity ZF\nwith clamp")))


    motifs = read_tsv(input_summit_annotation,
                           col_names=c("chrom",
                                       "summit_start",
                                       "summit_end",
                                       "summit_name",
                                       "score",
                                       "strand")) %>%
        select(-c(score,
                  strand)) %>%
        rowid_to_column("index") %>%
        left_join(read_tsv(motif_annotation,
                           col_names=c("chrom",
                                       "motif_left",
                                       "motif_right",
                                       "motif_name",
                                       "score",
                                       "strand",
                                       "junk",
                                       "sequence")) %>%
                      select(-junk),
                  by="chrom") %>%
        filter((motif_left >= (summit_start - 2010) & motif_left <= (summit_end + 2010)) |
                   (motif_right >= (summit_start - 2010) & motif_right <= (summit_end + 2010)) |
                   (motif_left <= (summit_start - 2010) & motif_right >= (summit_end + 2010))) %>%
        filter(score > 4) %>%
        mutate_at(vars(motif_left, motif_right),
                  ~((. - summit_start) / 1000)) %>%
        group_by(index) %>%
        mutate(y=(((row_number()-1) %% 4) * 0.1) - 0.5)

    annotations = read_tsv(input_summit_annotation,
                           col_names=c("chrom",
                                       "summit_start",
                                       "summit_end",
                                       "summit_name",
                                       "score",
                                       "strand")) %>%
        select(-c(score,
                  strand)) %>%
        rowid_to_column("index") %>%
        left_join(read_tsv(input_transcript_annotation,
                           col_names=c("chrom",
                                       "transcript_left",
                                       "transcript_right",
                                       "transcript_name",
                                       "score",
                                       "strand")) %>%
                      select(-score),
                  by="chrom") %>%
        filter((transcript_left >= (summit_start - 2010) & transcript_left <= (summit_end + 2010)) |
                   (transcript_right >= (summit_start - 2010) & transcript_right <= (summit_end + 2010)) |
                   (transcript_left <= (summit_start - 2010) & transcript_right >= (summit_end + 2010))) %>%
        left_join(read_tsv(input_orf_annotation,
                           col_names=c("chrom",
                                       "orf_left",
                                       "orf_right",
                                       "transcript_name",
                                       "score",
                                       "strand")) %>%
                      select(-score),
                 by=c("chrom",
                      "transcript_name",
                      "strand")) %>%
        mutate(transcript_start = ifelse(strand=="+",
                                         transcript_left,
                                         transcript_right),
               transcript_end = ifelse(strand=="+",
                                       transcript_right,
                                       transcript_left),
               notch = ifelse(strand=="+",
                              (orf_right - orf_left) * 0.85 + orf_left,
                              orf_right - (orf_right - orf_left) * 0.85),
               orf_start = ifelse(strand=="+",
                                  orf_left,
                                  orf_right),
               orf_end = ifelse(strand=="+",
                                orf_right,
                                orf_left)) %>%
        select(-c(transcript_left,
                  transcript_right,
                  orf_left,
                  orf_right,
                  strand)) %>%
        mutate_at(c("transcript_start",
                    "transcript_end",
                    "orf_start",
                    "orf_end",
                    "notch"),
                  ~((. - summit_start) / 1000)) %>%
        mutate(label_position = (orf_start + orf_end)/2) %>%
        expand(nesting(index,
                       transcript_name,
                       transcript_start,
                       transcript_end,
                       orf_start,
                       notch,
                       orf_end,
                       label_position),
               orf_y=c(-1, 1)) %>%
        gather(key=type,
               value=orf_x,
               c(orf_start,
                 notch,
                 orf_end)) %>%
        mutate(orf_y = ifelse(type=="orf_end",
                              0,
                              orf_y)) %>%
        distinct() %>%
        group_by(index,
                 transcript_name) %>%
        arrange(type,
                orf_y,
                .by_group=TRUE) %>%
        mutate(order = c(4,2,3,5,1)) %>%
        arrange(order,
                .by_group=TRUE) %>%
        mutate(ymin=min(df['signal']),
               range=max(df['signal']) - min(df['signal']),
               ymin=ymin-0.1*range)
        # left_join(df %>%
        #               group_by(index) %>%
        #               summarize(ymin=min(signal),
        #                         range=max(signal)-min(signal)) %>%
        #               mutate(ymin=ymin-0.1*range),
        #           by="index")

    plot = ggplot() +
        geom_segment(data = annotations %>%
                         distinct(transcript_name, .keep_all=TRUE),
                     aes(y=ymin,
                         yend=ymin,
                         x=pmax(-2.01, pmin(2.01, transcript_start)),
                         xend=pmax(-2.01, pmin(2.01, transcript_end))),
                     color="grey50") +
        geom_polygon(data = annotations,
                     aes(x=pmax(-2.01, pmin(2.01, orf_x)),
                         y=orf_y * 0.07 * range + ymin,
                         group=transcript_name),
                     fill="grey80") +
        geom_text(data = annotations %>%
                      filter((range(pmax(-2.01, pmin(2.01, orf_x)))[2] -
                                  range(pmax(-2.01, pmin(2.01, orf_x)))[1]) > 0.4) %>%
                      distinct(transcript_name, .keep_all=TRUE) %>%
                      filter(between(label_position, -2.01, 2.01)),
                  aes(x=label_position,
                      y=ymin,
                      label=transcript_name),
                  size=6/72*25.4,
                  fontface="italic",
                  family="FreeSans",
                  vjust=0.6) +
        geom_point(data=motifs %>%
                       filter(score>4),
                   aes(x=(motif_left + motif_right)/2,
                       y=y),
                   alpha=0.6,
                   size=0.4,
                   shape=16) +
        geom_line(data = df,
                  aes(x=position,
                      y=signal,
                      color=group),
                  size=0.4,
                  alpha=0.8,
                  position=position_identity()) +
        facet_wrap(~index,
                   # scales="free_y") +
                   scales="fixed") +
        scale_color_viridis_d(end=0.85,
        # scale_color_brewer(palette = 'Paired',
                              name="relative ChIP enrichment:",
                              guide=guide_legend(label.position="top",
                                                 label.vjust=0.5,
                                                 title.vjust=0.8,
                                                 override.aes=list(size=2))) +
        scale_linetype_discrete(name=NULL,
                                guide=guide_legend(label.position="top")) +
        scale_x_continuous(expand=c(0,0),
                           breaks=scales::pretty_breaks(n=2),
                           labels=function(x) case_when(x==0 ~ "peak summit",
                                                        x==2 ~ paste(x, "kb"),
                                                        TRUE ~ as.character(x)),
                           name=NULL) +
        scale_y_continuous(breaks=scales::pretty_breaks(n=3),
                           limits=c(NA, NA),
                           # limits=c(NA, 0.4),
                           name=NULL) +
        theme_light() +
        theme(text=element_text(color="black", size=8, family="FreeSans"),
              axis.text=element_text(color="black"),
              axis.text.y=element_text(size=4,
                                       margin=margin(r=0.5, unit="pt")),
              strip.text=element_blank(),
              panel.grid=element_blank(),
              panel.spacing.x=unit(10, "pt"),
              panel.spacing.y=unit(1, "pt"),
              legend.position="top",
              legend.direction="horizontal",
              legend.spacing.y=unit(0, "pt"),
              legend.text=element_text(margin=margin(b=-5, unit="pt")),
              legend.margin=margin(b=-11, unit="pt"),
              axis.title.y=element_text(angle=0, vjust=0.5))

    ggsave(output_path,
           plot=plot,
           width=17.4,
           height=9,
           units="cm",
           device=cairo_pdf)
}

main(input_coverage=snakemake@input[["coverage"]],
     input_summit_annotation=snakemake@input[["summit_annotation"]],
     input_transcript_annotation=snakemake@input[["transcript_annotation"]],
     input_orf_annotation=snakemake@input[["orf_annotation"]],
     motif_annotation=snakemake@input[['motif_annotation']],
     filter_groups=snakemake@params[['filter_groups']],
     output_path=snakemake@output[["coverage"]])

