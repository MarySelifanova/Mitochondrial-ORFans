library(tidyverse)
library(ggpubr)
library(ggrepel)

#Set working directory
setwd("/home/dima/Desktop/R")


# Read table - CHECK SEPARATORS!
head = TRUE  #set HEADER

inn <- read.table("nucl_dist.tsv", header=head, sep ="\t")
cai <-  read.table("cai_last_corr.tsv", header=head, sep ="\t")
seqstats <-  read.table("seqstats.tsv", header=head, sep ="\t")

my_theme <- theme_bw() +
  theme(axis.text.x = element_text(size=14, angle=90,hjust=0.95,vjust=0.5), 
        axis.text.y = element_text(size=14), 
        strip.text.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_blank(),
        legend.position = "none")   

inn <- inn %>% mutate(gene = factor(gene, levels = c(
  "cox1",
  "ORFan3",
  "nad6",
  "atp8",
  "ORFan1",
  "nad3",
  "cytb",
  "nad4",
  "cox3",
  "atp6",
  "ORFan2",
  "nad5",
  "nad4l",
  "cox2",
  "nad1",
  "nad2"
)))

cai <- cai %>% mutate(gene = factor(gene, levels = c(
  "cox1",
  "ORFan3",
  "nad6",
  "atp8",
  "ORFan1",
  "nad3",
  "cytb",
  "nad4",
  "cox3",
  "atp6",
  "ORFan2",
  "nad5",
  "nad4l",
  "cox2",
  "nad1",
  "nad2"
)))

seqstats <- seqstats %>% mutate(gene = factor(gene, levels = c(
  "cox1",
  "ORFan3",
  "nad6",
  "atp8",
  "ORFan1",
  "nad3",
  "CR",
  "cytb",
  "nad4",
  "cox3",
  "atp6",
  "ORFan2",
  "nad5",
  "nad4l",
  "cox2",
  "nad1",
  "nad2"
)))


unique(inn$gene)

Fig1C <- inn %>% ggplot(aes(x = ds, 
                   y = dn, 
                   fill = group,
                   color = group)) + 
  geom_point(shape = 21, alpha = 0.6, size = 3) + 
  xlim(0,1) +
  ylim(0,1) +
  geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), 
               lty = 2, 
               color = 'red') +
  geom_text_repel(aes(label = gene),
                      size = 6, 
                  max.overlaps = Inf, 
                  # alpha = 0.6,
                  segment.alpha = 0.6,
                  min.segment.length = 0,
                  size = 3) +
  my_theme +
  theme(axis.title.x = element_text(size = 14)) 
  

Fig1A <- inn %>% ggplot(aes(x = gene, 
                   y = p_nucl, 
                   fill = group)) +
  geom_point(size = 2, shape = 21) + 
  geom_segment(aes(x = gene, y = 0, xend = gene, yend = p_nucl), 
               lty = 3) +
  geom_errorbar(aes(ymin=p_nucl + se, ymax=p_nucl - se), width=.2) +
  ylim(0,0.6) +
  my_theme +
  ylab("Nucleotide diversity") +
  theme(axis.text.x = element_text(angle=90))

Fig1B <- inn %>% ggplot(aes(x = gene, 
                            y = p_aa, 
                            fill = group)) +
  geom_point(size = 2, shape = 21) + 
  geom_segment(aes(x = gene, y = 0, xend = gene, yend = p_aa), 
               lty = 3) +
  geom_errorbar(aes(ymin=p_aa + se.1, ymax=p_aa - se.1), width=.2) +
  ylim(0,0.6) +
  ylab("A.a. diversity") +
  my_theme +
  theme(axis.text.x = element_text(angle=90))

#### CAI
glimpse(cai)

Fig1D <- cai %>% ggplot(aes(x = gene, y = cai)) + 
  geom_boxplot() + 
  geom_point(aes(color = oragnism)) +
  ylim(0.5,1) +
  ylab("CAI") +
  my_theme + 
  theme(axis.text.x = element_text(angle=90),
        legend.position = "right")

Fig1AB <- ggarrange(Fig1A, Fig1B, ncol = 1)
Fig1ABC <- ggarrange(Fig1AB, Fig1C, nrow = 1)
Fig1ABCD <- ggarrange(Fig1ABC, Fig1D, ncol = 1)

#### Genome content stats 
seqstats_long <- seqstats %>%
  pivot_longer(names_to = "genome_param", values_to = "value", gc_content:gc_skew)

FigureS5B <- seqstats_long %>% ggplot(aes(x = gene, y = value, fill = species)) +
  geom_point(size = 3, shape = 21) + 
  geom_line(aes(group = species, color = species), alpha = 0.2) + 
  my_theme + 
  theme(legend.position = "right") + 
  facet_wrap( ~ genome_param, ncol = 1, scales = "free_y", strip.position = "right")
