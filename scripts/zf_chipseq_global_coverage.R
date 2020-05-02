main = function(theme_path = "custom_theme.R",
                pdf_out = "test.pdf",
                data_path = "ZF_chipseq_union-bedgraph-libsizenorm-midpoint-window-1000-allsamples.tsv.gz"){
    source(theme_path)

    df = read_tsv(data_path) %>%
        select_at(vars(-contains("input"))) %>%
        pivot_longer(-name,
                     names_to="sample",
                     values_to="signal") %>%
        mutate(group=str_remove(sample, "-[0-9]")) %>%
        group_by(name, group) %>%
        summarize(signal=mean(signal)) %>%
        separate(name,
                 into=c("chrom", "start", "end"),
                 sep="-",
                 convert=TRUE) %>%
        filter(chrom != "chrM") %>%
        mutate(chrom=ordered(chrom,
                             levels=c("chrI", "chrII", "chrIII", "chrIV", "chrV", "chrVI", "chrVII", "chrVIII", "chrIX", "chrX",
                                      "chrXI", "chrXII", "chrXIII", "chrXIV", "chrXV", "chrXVI")),
               group=ordered(group,
                             levels=c("reporter-only-ZF-IP",
                                      "high-affinity-ZF-ZF-IP",
                                      "low-affinity-ZF-ZF-IP",
                                      "low-affinity-ZF-with-clamp-ZF-IP"),
                             labels=c("reporter only",
                                      "high-affinity ZF",
                                      "low-affinity ZF",
                                      "low-affinity ZF + clamp")))

    chrom_sizes = df %>%
        group_by(chrom) %>%
        summarize(size=max(end)) %>%
        mutate(cum_start = cumsum(size) - size)

    n_chroms = length(chrom_sizes[["chrom"]])

    df %<>%
        left_join(chrom_sizes,
                  by="chrom") %>%
        mutate(position = (start + end) / 2 + cum_start)

    plot = ggplot(data=df,
           aes(x=position,
               y=signal,
               color=chrom,
               fill=chrom)) +
        geom_area(alpha=0.5,
                  size=0.1) +
        facet_grid(group ~ .) +
        scale_x_continuous(breaks=chrom_sizes[["cum_start"]] + (chrom_sizes[["size"]]/ 2),
                           minor_breaks=chrom_sizes[["cum_start"]],
                           labels=chrom_sizes[["chrom"]],
                           expand=c(0.01,0)#,
                           # sec.axis=dup_axis()
                           ) +
        scale_y_continuous(limits=quantile(df[["signal"]],
                                           probs=c(0, 0.999)),
                           oob=scales::squish,
                           expand=c(0,0),
                           name="IP ChIP signal") +
        scale_color_manual(values=rep(viridisLite::viridis(2, end=0.3),
                                      ceiling(n_chroms/2)),
                           guide=FALSE) +
        scale_fill_manual(values=rep(viridisLite::viridis(2, end=0.3),
                                     ceiling(n_chroms/2)),
                          guide=FALSE) +
        theme_default +
        theme(strip.text.y=element_text(angle=0,
                                        hjust=0),
              axis.text=element_text(color="black"),
              axis.text.x=element_text(angle=30,
                                       hjust=1,
                                       vjust=1),
              axis.text.x.top=element_text(angle=30,
                                           hjust=0,
                                           vjust=0),
              axis.title.x=element_blank(),
              legend.position="none",
              panel.grid=element_blank())

    ggsave(pdf_out,
           plot=plot,
           width=16*2,
           height=9*2,
           units="cm",
           device=cairo_pdf)
}

main(theme_path = snakemake@input[["theme"]],
     pdf_out = snakemake@output[["pdf"]],
     data_path = snakemake@input[["data"]])

