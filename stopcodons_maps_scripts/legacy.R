#!/usr/bin/env Rscript

source("utils.R")
source("shannon.R")

library(Biostrings)

library(RColorBrewer)

library(IRanges)
library(rtracklayer)
library(circlize)



## end-to-end function to plot coverage


plot_cov <- function(cov_path, color, maximum) {
  coverage <- read_delim(cov_path, delim = "\t", col_names = F)
  
  colnames(coverage) <- c("name", "position", "coverage")
  
  breaks <- (0:ceiling(GenomeSize/step))*step
  
  coverage %>%
    mutate(chunk=cut(position, breaks)) %>%
    group_by(chunk, .drop=F) %>%
    summarize(mean_cov=mean(coverage), z_cov=(sd(coverage)/mean(coverage))) %>%
    replace_na(list(mean_cov=0, z_cov=0)) %>%
    mutate(chunk_center=(head(breaks, -1)+step/2))-> avg_cov
  
  
  circos.track(ylim = c(0, log10(maximum)))
  
  cov_palette <- colorRamp(c(color, "black"))
  max_z_cov <- max(avg_cov$z_cov[!is.na(avg_cov$z_cov)])
  cov_cols <- cov_palette(avg_cov$z_cov/max_z_cov)
  cov_cols_hex <- apply(cov_cols, 1, function(x) {
    if (!is.na(x[1])) {
      paste0("#", format(as.hexmode(round(x[1])*256*256 + round(x[2])*256 + round(x[3])), width=6, upper.case=T))
    } else {
      "#000000"
    }
  })
  
  
  circos.lines(x = avg_cov$chunk_center, y = log10(avg_cov$mean_cov),
               area = TRUE, type = 'h', col = cov_cols_hex)
  par(cex = 0.8)
  circos.yaxis(side = "right", at=seq(0, maximum, length.out = 3))
  
}



#xxxxxxxx hardcode zone - here you choose which assembly to use for further analysis

annelidGenome <- readDNAStringSet("spades_tr1.fasta")
annelidGenome <- readDNAStringSet("mitocontig_solov.fasta")

annelidGenome <- readDNAStringSet("assemblies/final_assembly_wsbs.fasta")
annelidGenome <- readDNAStringSet("assembly_brevipalpa.fasta")
annelidGenome <- readDNAStringSet("assembly_hoplura.fasta")
annelidGenome <- readDNAStringSet("assembly_websteri.fasta")
annelidGenome <- readDNAStringSet("bocham.fasta")
annelidGenome <- readDNAStringSet("lindaspio.fasta")
annelidGenome <- readDNAStringSet("marneg.fasta")

#xxxxxxx end of hardcode zone

######### Here you count sequence statistics for chunks of genome


seqstat_window <- 200 # set the window for GC content
information_window <- 100
step <- 100
GenomeSize <- length(annelidGenome[[1]])
steps <- GenomeSize %/% step

stops <- count_stops(annelidGenome[[1]])

seqstat_chunks <- get_chunks(annelidGenome[[1]], step, seqstat_window)
seq_stats <- decorate_pos(seqstat_chunks, count_seq_stats)

info_chunks <- get_chunks(annelidGenome[[1]], information_window, step)
info_stats <- decorate_pos(info_chunks, count_2_mers)

left_join(seq_stats, info_stats, by="pos")


seq_stats <- seq_stats %>% mutate(cumulativeATskew =
                                   cumsum(seq_stats$ATskew)) %>%
  mutate(cumulativeGCskew = cumsum(seq_stats$GCskew))


################ Circular maps
### BioCircus
png(out_file)
init_without_sectors(GenomeSize)

annot_name <- "mitocontig_solovki.bed"
annot_name <- "genes_fastas/wsbs_with_orfans.bed"

annot_name <- "genes_fastas/wsbs.bed"
annot_name <- "genes_fastas/polbre.bed"
annot_name <- "genes_fastas/polhop.bed"
annot_name <- "genes_fastas/polweb.bed"
annot_name <- "genes_fastas/bocham.bed"


annot_name <- "lindaspio.bed"
annot_name <- "marneg.bed"

annot <- read_annot(annot_name)


init_with_sectors(annot, GenomeSize)
# circos.track(ylim = c(0,1),bg.border = "red")
# circos.text(mtDNAfeatures$coordinates,0.7,mtDNAfeatures$X_Gene,
# niceFacing = TRUE, facing = "clockwise", cex = 0.5)
par(cex=0.8)

# annotation tracks

plot_annot(annot)

# stop codons tracks
plot_stops(stops, c("TAA", "TAG"), "red")
plot_stops(stops, c("TTA", "CTA"), "black")

# plot potential orfs
plot_annot(search_potential_orfs(stops, c("TAA", "TAG")), labels = "inside")


trnas <- annot[annot$type == "trna",]
pot_orf <- search_potential_orfs(stops, c("TAA", "TAG"))

pot_orf <- search_potential_orfs(stops, c("TTA", "CTA"))

save_bed(pot_orf, "wsbs_corrected_no_stops.bed")
save_bed(pot_orf, "polbre_no_stops.bed")
save_bed(pot_orf, "polhop_no_stops.bed")
save_bed(pot_orf, "polweb_no_stops.bed")
save_bed(pot_orf, "solov_no_stops.bed")

pot_orf %>%
  group_by(name) %>%
  group_modify(function(x, y) {
    row <- x[1,]
    trnas_overlapping <- (((trnas$end - 1) %% GenomeSize) + 1 > row$start) & (trnas$start < ((row$end - 1) %% GenomeSize) + 1)
    trna_subset <- trnas[trnas_overlapping, ]
    breaks <- sort(c(trna_subset$start, trna_subset$end, row$start, row$end))
    if(breaks[1] < row$start) {
      breaks <- breaks[3:length(breaks)]
    }
    break_mat <- matrix(breaks, ncol=2, byrow = T)
    colnames(break_mat) <- c("start", "end")
    break_mat <- as_tibble(break_mat)
    
    break_mat %>%
      mutate(i=1:n()) %>%
      filter(end-start > 210) -> filtered_break_mat
    
    orf_end <- row$end
    residuals <- (orf_end - filtered_break_mat$start) %% 3
    res_add <- c(1, 2, 0)
    filtered_break_mat$start <- filtered_break_mat$start + ((residuals + 1) %% 3)
    
    
    filtered_break_mat %>%
      mutate(chr=row$chr, strand=row$strand, type=row$type, new_name=paste0(y$name, ".", seq(1, n(), length.out = n()))) %>%
      select(chr, start, end, new_name, strand, type) -> cutted_orf
      
    cutted_orf
  }) %>%
  mutate(name=new_name) %>%
  select(chr, start, end, name, strand, type)-> cutted_pot_orfs

plot_annot(cutted_pot_orfs, labels = "inside")
save_bed(cutted_pot_orfs, "wsbs_corrected_trn_punct.bed")
save_bed(cutted_pot_orfs, "polbre_trn_punct.bed")
save_bed(cutted_pot_orfs, "polhop_trn_punct.bed")
save_bed(cutted_pot_orfs, "polweb_trn_punct.bed")
save_bed(cutted_pot_orfs, "solov_trn_punct.bed")
save_bed(cutted_pot_orfs, "bocham_trn_punct.bed")


renamed_fn <- "polwsbs_trn_punct_renamed.bed.tmp"
genome_name <- "final_assembly_wsbs.fasta"

file_groups <- list()
file_groups$species <- c("polwsbs", "polweb", "polbre", "polhop", "bocham")
file_groups$renamed_fn <- c("polwsbs_trn_punct_renamed.bed.tmp", "polweb_trn_punct_renamed.bed.tmp", "polbre_trn_punct_renamed.bed.tmp", "polhop_trn_punct_renamed.bed.tmp", "bocham_trn_punct_renamed.bed.tmp")
file_groups$assembly <- c("final_assembly_wsbs.fasta", "assembly_websteri.fasta", "assembly_brevipalpa.fasta", "assembly_hoplura.fasta", "bocham.fasta")
file_groups <- as_tibble(file_groups)

library(reshape2)

file_groups %>%
  group_by(species) %>%
  group_modify(function(.x, .y) {
    print(.y$species)
    annelidGenome <- readDNAStringSet(.x$assembly)
    cutted_pot_orfs <- as_tibble(read_annot(.x$renamed_fn))
    cutted_pot_orfs$start <- as.double(cutted_pot_orfs$start)
    cutted_pot_orfs$end <- as.double(cutted_pot_orfs$end)
    cutted_pot_orfs <- get_nsr_seq(annelidGenome[[1]], cutted_pot_orfs)
    cutted_pot_orfs <- count_nsr_seq_stats(cutted_pot_orfs)
    
    cutted_pot_orfs %>%
      select(name, gc_content, at_skew, gc_skew) %>%
      melt()
  }) -> for_plot

for_plot %>%
  pivot_wider(names_from = variable, values_from = value) -> wider_table

write_tsv(wider_table, "seqstats.tsv")

protlevels <- c("nad6", "orfan1", "nad3", "CR", "cytb", "nad4", "cox3", "atp6", "orfan2", "nad5", "nad4L", "cox2", "nad1", "nad2", "cox1", "orfan3")
for_plot$name <- factor(for_plot$name, levels=protlevels)

for_background <- unique(for_plot %>% ungroup() %>% select(name, variable))
for_background$value <- 1

orfan_coords <- grepl("orfan", for_background$name)
for_background$is_orfan <- orfan_coords

for_plot %>%
  group_by(name, variable) %>%
  summarise(avg=mean(value)) -> for_avg

ggplot(for_plot) +
  geom_rect(data=for_background, mapping=aes(xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=is_orfan), alpha=0.1) +
  geom_hline(data=for_avg, mapping=aes(yintercept=avg), show.legend = F, size=0.25) +
  scale_fill_manual(values=c("gray", "red")) +
  facet_grid(cols=vars(name), rows=vars(variable), scales = "free_y") +
  geom_point(aes(y=value, x=0, color=species), stat = "identity") + 
  ggtitle("5 species genes' sequence stats") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) -> p

gp <- ggplotGrob(p)
facet.columns <- gp$layout$l[grepl("panel", gp$layout$name)]
gp$widths[facet.columns] <- gp$widths[facet.columns] * x.var

diststats <- read_tsv("nucl_dist.txt")
diststats %>% 
  select("gene", "p_nucl", "p_aa", "dn/ds") -> protstats

orfan_coords <- grepl("orfan", protstats$gene)
protstats$is_orfan <- orfan_coords
protstats$gene <- factor(protstats$gene, levels = protlevels)
colnames(protstats) <- c("gene", "p_nucl", "p_aa", "dnds", "is_orfan")
protstats %>%
  select(gene, p_nucl, p_aa, dnds, is_orfan) %>%
  melt() -> for_prot_plot

ggplot(for_prot_plot) +
  geom_rect(mapping=aes(xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=is_orfan), alpha=0.1) +
  scale_fill_manual(values=c("gray", "red")) +
  facet_grid(cols=vars(gene), rows=vars(variable), scales = "free_y", labeller = labeller(variable = c("p_nucl"="nucleotide p-distance", "p_aa"="aminoacid p-distance", "dnds"="dn/ds"))) +
  geom_point(aes(y=value, x=0), stat = "identity") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

# # GC track
circos.track(ylim = c(0,max(seq_stats$GC)))
circos.lines(x = seq_stats$pos, y = seq_stats$GC, area = TRUE, type =
               'h', col = 'darkblue', lwd = 2)
par(cex = 0.8)
circos.yaxis(side = "right", at=seq(0, 0.4, length.out = 3))

# chi track
circos.par(track.height=0.15)
circos.track(ylim = c(0,max(seq_stats$chi)))
circos.lines(x=seq_stats$pos, y=seq_stats$chi, area=T, type="l", col = "cyan")
par(cex = 0.8)
circos.yaxis(side = "left", at=seq(0, floor(max(seq_stats$chi)), length.out = 3))
# entropy loss track
circos.track(ylim = c(0, max(seq_stats$entropy_loss)))
circos.lines(x=seq_stats$pos, y=seq_stats$entropy_loss, area=T, type="l", col = "gray")
par(cex = 0.8)

circos.yaxis(side = "right", at=c(0, floor(max(seq_stats$entropy_loss)*100)/100))
# cg dimers track
circos.par(track.height=0.15)
log_ratio <- log2(seq_stats$cg_rat)
log_ratio[log_ratio == -Inf] = -3
circos.track(ylim = c(-3, 3))
circos.lines(x=seq_stats$pos, y=log_ratio, area=F, type="l", col = "darkgreen")
circos.segments(0, 0, max(seq_stats$pos),0, type = "l", lwd = 1, lty =
                  2, col  = "black")
par(cex = 0.8)

circos.yaxis(side = "left", at=c(-3, 0, 3))

## to plot Coverage from real data
plot_cov("spades_tr1.tsv", "red", 150)
plot_cov("final_assembly_coverage.tsv", "green", 150)
plot_cov("solomap11.tsv", "green", 75)
plot_cov("solomap12.tsv", "green", 75)
plot_cov("solomap21.tsv", "green")
plot_cov("solomap22.tsv", "green")
plot_cov("masha_coverage.tsv", "green", 50)
plot_cov("masha_coverage_soft.tsv", "green", 50)
plot_cov("masha_coverage_softer.tsv", "green", 50)
plot_cov("masha_coverage_softest.tsv", "green", 50)
plot_cov("solov_self_coverage.tsv", "green", 50)

plot_cov("new_coverage/solovki_dora_coverage.tsv", "red", 10)
plot_cov("new_coverage/solS1_coverage.tsv", "red", 10)
plot_cov("new_coverage/solS2_coverage.tsv", "red", 10)
plot_cov("new_coverage/solS3_coverage.tsv", "red", 1000)
plot_cov("new_coverage/spio_coverage.tsv", "blue", 1000)

text(0, 0, paste0(round(GenomeSize/1000), "kb"), cex = 1)

dev.off()


# brevipalpa
title("Mitobim on brevipalpa reference trim1")
# hopula
title("Mitobim on hopula reference trim1")
# websteri
title("Mitobim on websteri reference trim1")
# for solov polydora coverage data
title("Coverage of reads from solovian polydora on our assembly.")

title("Solovian polydora assembly stats")

par(cex=1.25)
title("Stop codons density in different reading frames in corrected wsbs polydora")

# cov paths
cov_path <- "solomap11.tsv"
cov_path <- "solomap12.tsv"
cov_path <- "solomap21.tsv"
cov_path <- "solomap22.tsv"


# other plots

library(ggplot2)
stops$codon_distances %>%
  group_by(codon) %>%
  summarize(delta_sd=sd(delta)) -> dist_sd

stops$codon_distances %>%
  group_by(codon) %>%
  filter((max(delta) > 2000) || is_stop) -> filtered_dists

ggplot(filtered_dists, aes(x=codon, y=delta)) + geom_boxplot()

dist_sd %>%
  filter(delta_sd > 500) -> dist_sd_filtered

ggplot(dist_sd_filtered, aes(x=codon, y=delta_sd)) + geom_bar(stat="identity")

stops$codon_distances %>%
  filter(!is_stop) -> cod_dists

ggplot(cod_dists, aes(y=delta)) + geom_histogram()
