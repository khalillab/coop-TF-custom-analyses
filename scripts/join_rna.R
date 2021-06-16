library(tidyverse)

high = snakemake@input[["high"]]
low = snakemake@input[["low"]]
low_clamp = snakemake@input[["low_clamp"]]
output = snakemake@output[["tsv"]]

import = function(path){
    read_tsv(path) %>%
        select(chrom, start, end, name, rna_log10_padj_high_v_reporter=log10_padj, rna_enrichment=log2_foldchange) %>%
        return()
}

import(high) %>%
    left_join(import(low) %>% rename(rna_log10_padj_low_v_reporter=rna_log10_padj_high_v_reporter),
              by=c("chrom", "start", "end", "name"),
              suffix=c("_high", "_low")) %>%
    left_join(import(low_clamp) %>% rename(rna_log10_padj_low_clamp_v_reporter=rna_log10_padj_high_v_reporter)) %>%
    rename(rna_enrichment_low_clamp = rna_enrichment) %>%
    write_tsv(output)

