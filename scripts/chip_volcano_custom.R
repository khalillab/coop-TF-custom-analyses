import = function(path, strain_id){
    read_tsv(path) %>%
        mutate(strain=strain_id) %>%
        return()
}

main = function(theme_path='custom_theme.R',
                low_v_high='low-affinity-ZF-v-high-affinity-ZF_ZF-chipseq-libsizenorm-peaks-diffbind-results-all.tsv',
                low_clamp_v_high='low-affinity-ZF-with-clamp-v-high-affinity-ZF_ZF-chipseq-libsizenorm-peaks-diffbind-results-all.tsv',
                #low_clamp_v_low='low-affinity-ZF-with-clamp-v-low-affinity-ZF_rnaseq-libsizenorm-verified_genes_plus_Venus-diffexp-results-all.tsv',
                fdr=0.1,
                output_volcano){

    source(theme_path)

    df = import(low_v_high,
                "\"log\"[2] ~ frac(\"low-affinity ZF\", \"high-affinity ZF\")") %>%
        bind_rows(import(low_clamp_v_high,
                         "\"log\"[2] ~ frac(\"low-affinity ZF + clamp\", \"high-affinity ZF\")")) %>%
        mutate(significant = (log10_padj > -log10(fdr)),
               strain = fct_inorder(strain, ordered=TRUE))

    volcano = ggplot() +
        geom_vline(xintercept=0,
                   color="grey80",
                   size=0.5) +
        geom_point(data = df,
                   aes(x=log2FC_enrichment,
                       y=log10_padj,
                       alpha=significant,
                       color=significant),
                   size=0.75,
                   shape=16) +
        facet_grid(.~strain,
                   labeller = label_parsed,
                   switch='x') +
        scale_y_continuous(limits=c(0, NA),
                           expand=c(0, 2),
                           name=expression("-log"[10] * "FDR")) +
        scale_x_continuous(name=NULL) +
        scale_alpha_manual(values=c(0.2, 0.6)) +
        scale_color_manual(values=c("grey40", "#440154FF")) +
        theme_light() +
        theme(text=element_text(color="black",
                                size=10),
              strip.placement = 'outside',
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

}

main(theme_path=snakemake@input[["theme"]],
     low_v_high=snakemake@input[["low_v_high"]],
     low_clamp_v_high=snakemake@input[["low_clamp_v_high"]],
     fdr=snakemake@params[["fdr"]],
     output_volcano=snakemake@output[["volcano"]])

