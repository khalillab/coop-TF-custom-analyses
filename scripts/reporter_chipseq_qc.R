main = function(theme_path = "custom_theme.R",
                data_path = "reporter-binding-site-allsamples-allannotations-ZF-chipseq-spikenorm-protection.tsv.gz",
                pdf_out='test.pdf'){

    source(theme_path)

    df = read_tsv(data_path,
                  col_names=c('group', 'sample', 'source', 'annotation', 'index', 'position', 'signal')) %>%
        mutate(zf_type=if_else(str_detect(group, 'noVP16'), 'no VP16', 'with VP16'),
               group=str_remove(group, '-noVP16'))

    df %<>%
        bind_rows(df %>%
                      filter(group=='reporter-only') %>%
                      mutate(zf_type='no VP16')) %>%
        mutate(zf_type=ordered(zf_type, levels=c('with VP16', 'no VP16')))

    plot = ggplot(data=df,
           aes(x=position,
               y=signal,
               group=interaction(sample, source),
               color=group,
               linetype=source)) +
        geom_vline(xintercept=0,
                   color='gray70',
                   size=0.2) +
        geom_line(alpha=0.7) +
        facet_grid(group~zf_type) +
        scale_x_continuous(expand=c(0,0),
                           labels=function(x) if_else(x==0, '42-10\nbinding sites', paste(as.character(x), 'kb')),
                           name=NULL) +
        scale_y_continuous(name='spike-in normalized signal') +
        scale_color_brewer(palette='Set1',
                           guide=FALSE) +
        theme_default +
        theme(strip.text.y=element_text(angle=0,
                                        hjust=0),
              panel.grid=element_blank(),
              legend.title=element_blank(),
              legend.position=c(0.99,0.99),
              legend.justification=c(1,1))
    ggsave(pdf_out,
           plot=plot,
           width=12,
           height=12,
           units='cm',
           device=cairo_pdf)
}

main(theme_path = snakemake@input[["theme"]],
     data_path = snakemake@input[["data"]],
     pdf_out= snakemake@output[["pdf"]])
