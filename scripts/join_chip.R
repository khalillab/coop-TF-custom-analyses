library(tidyverse)

filter = snakemake@input[["affinity_dependent"]]
high = snakemake@input[["high"]]
low = snakemake@input[["low"]]
low_clamp = snakemake@input[["low_clamp"]]
output = snakemake@output[["tsv"]]

import = function(path){
    read_tsv(path) %>%
        # select(chrom, start, end, name, chip_enrichment=condition_enrichment) %>%
        select(chrom, start, end, name, chip_enrichment=log2FC_enrichment) %>%
        return()
}

import(high) %>%
    left_join(import(low),
              by=c("chrom", "start", "end", "name"),
              suffix=c("_high", "_low")) %>%
    left_join(import(low_clamp)) %>%
    rename(chip_enrichment_low_clamp = chip_enrichment) %>%
    semi_join(import(filter) %>% select(-chip_enrichment)) %>%
    write_tsv(output)

