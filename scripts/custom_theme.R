library(tidyverse)
library(magrittr)
library(extrafont)
library(ggthemes)

theme_default = theme_light() +
    theme(text=element_text(size=8,
                            color="black",
                            family="FreeSans"),
          axis.text=element_text(size=7,
                                 color="black"),
          axis.text.x=element_text(margin=margin(1, 0, 0, 0, "pt")),
          axis.text.y=element_text(margin=margin(0, 0.5, 0, 0.5, "pt")),
          axis.title.x=element_text(size=7,
                                    margin=margin(1, 0, 0, 0, "pt")),
          axis.title.y=element_text(size=7,
                                    margin=margin(0, 1, 0, 0, "pt")),
          plot.title=element_text(size=8,
                                  margin=margin(0, 0, 2, 0, "pt")),
          plot.subtitle=element_text(size=7,
                                     margin=margin(0, 0, 0, 0, "pt")),
          plot.margin=margin(11/2, 11/2, 11/2, 11/2, "pt"),
          legend.text=element_text(size=7),
          legend.key.height=unit(10, "pt"),
          legend.box.margin=margin(0, 0, 0, 0, "pt"),
          legend.margin=margin(0, 0, 0, 0, "pt"),
          strip.background=element_blank(),
          strip.text=element_text(color="black"),
          plot.tag=element_text(size=9,
                                face="bold"))
