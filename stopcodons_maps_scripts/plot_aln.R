library(tidyverse)

library(RColorBrewer)
library(Biostrings)
library(IRanges)
library(rtracklayer)
library(circlize)

aln_path <- "~/knorre/align_polydora.fa"
primer_starts <- c(5, 5510, 100043, 14634)
annot_path <- "~/knorre/spades_tr1_moved.bed"

alignment <- readDNAStringSet(aln_path)
aln_len <- length(alignment[[1]])
ngenomes <- length(alignment)

alignment <- lapply(alignment, as.vector)
alignment <- do.call(bind_cols, alignment)

nuc_diversity <- function(col) {
  sim_pairs <- sum((table(col) - 1)*table(col)/2)
  all_pairs <- length(col)*(length(col) - 1)/2
  return(1 - sim_pairs/all_pairs)
}

apply(alignment, 1, function(x) {
  has_gap <- "-" %in% x
  
  sum(sapply(table(x[x != "-"]), function(count) {
    p <- count/ngenomes
    return(-p*log2(p))
  })) -> enthropy
  return(tibble(has_gap=has_gap, enthropy=enthropy, nuc_diversity=nuc_diversity(x[x != "-"])))
}) %>%
  bind_rows()  -> aln_stats

step=100

aln_stats %>%
  mutate(group=(1:aln_len %/% step)) %>%
  group_by(group) %>%
  summarise(mean_ent=mean(enthropy), gap_ratio=sum(has_gap)/n()) -> slice_summaries

struct_path <- "~/knorre/dora_struct.tsv"
struct <- read.csv(struct_path, header=F)

slice_summaries$struct_scor <- struct$V1

no_gap2_gap <- function(gap_seq, coords)
{ 
  new_coords<-which(gap_seq!='-')
  return(new_coords[coords])
}

mtDNAfeatures <- read.table(annot_path, sep = "\t")

colnames(mtDNAfeatures) <- c("node", "start", "end", "name", "xxx",
                             "strand")

mtDNAfeatures$start <- no_gap2_gap(alignment$`NODE_1_length_17700_cov_18.721508/1-17700`, mtDNAfeatures$start)
mtDNAfeatures$end <- no_gap2_gap(alignment$`NODE_1_length_17700_cov_18.721508/1-17700`, mtDNAfeatures$end)

mtDNAfeatures$chr <- "a"
mtDNAfeatures$end[mtDNAfeatures$end < mtDNAfeatures$start] <- mtDNAfeatures$end[mtDNAfeatures$end < mtDNAfeatures$start] + aln_len
mtDNAfeatures$col <- "pink"
mtDNAfeatures$col[3] <- "green"

affilations <- rep("intergenic", aln_len)

get_range <- function(x) {
  as.integer(x[2]):(as.integer(x[3]) - 1)
}

intergenic <- rep(T, aln_len)

trn_lines <- mtDNAfeatures[startsWith(mtDNAfeatures$name, "trn"),]
trn_pos <- unlist(apply(trn_lines, 1, get_range))
trn_pos[trn_pos > aln_len] <- trn_pos[trn_pos > aln_len] - aln_len
intergenic[trn_pos] <- F

prot_lines <- mtDNAfeatures[!startsWith(mtDNAfeatures$name, "trn"),]
prot_pos <- unlist(apply(prot_lines, 1, get_range))
prot_pos[prot_pos > aln_len] <- prot_pos[prot_pos > aln_len] - aln_len
intergenic[prot_pos] <- F

inter_pos <- which(intergenic)

mean_not_na <- function(x) {
  mean(x[!is.na(x)])
}

avg_prot_div <- mean_not_na(aln_stats$nuc_diversity[prot_pos])
avg_trna_div <- mean_not_na(aln_stats$nuc_diversity[trn_pos])
avg_inter_div <- mean_not_na(aln_stats$nuc_diversity[inter_pos])

hole1_pos <- 2295:3476
hole2_pos <- 4635:5490
hole3_pos <- 9827:10382
hole4_pos <- 17433:18000

hole1_div <- mean_not_na(aln_stats$nuc_diversity[hole1_pos])
hole2_div <- mean_not_na(aln_stats$nuc_diversity[hole2_pos])
hole3_div <- mean_not_na(aln_stats$nuc_diversity[hole3_pos])
hole4_div <- mean_not_na(aln_stats$nuc_diversity[hole4_pos])

circos.clear()
circos.par("gap.degree" = 0,
           "cell.padding" = c(0, 0, 0, 0),
           "start.degree" = 90,
           "track.height" = 0.1)

circos.initialize("a", xlim = c(0, aln_len))

circos.genomicLabels(mtDNAfeatures[, c("chr", "start", "end", "name")],
                     labels.column = 4,
                     side = "outside",
                     # col = FeatureColor[mtDNAfeatures$type],
                     cex = 1)

circos.par("track.height" = 0.05)
par(cex=0.8)
circos.track(ylim = c(0,1), bg.border = "white")
for (row in 1:nrow(mtDNAfeatures)) {
  circos.arrow(mtDNAfeatures$start[row],
               mtDNAfeatures$end[row],
               y= 0.5,
               width = 0.8,
               arrow.head.width = 0.8,
               arrow.head.length = 100,
               arrow.position =
                 ifelse(mtDNAfeatures$strand[row] == "+", "end", "start"),
               # # col =
               # FeatureColor[mtDNAfeatures$type[row]],
               col = mtDNAfeatures$col[row]
  )
}

circos.par("track.height" = 0.1)
circos.track(ylim=c(0, log2(ngenomes)))
circos.lines((slice_summaries$group + 0.5)*step, slice_summaries$mean_ent, col = "black", area = T)

circos.track(ylim=c(0, 1))
circos.lines((slice_summaries$group + 0.5)*step, slice_summaries$gap_ratio, col="red")

circos.par("track.height" = 0.05)
circos.track(ylim=c(0, 1))
circos.points(primer_starts, rep(c(0.5), length(primer_starts)), col="magenta")

circos.par("track.height" = 0.1)
circos.track(ylim=c(-1, 1))
circos.lines((slice_summaries$group + 0.5)*step, slice_summaries$struct_scor, col="blue", area=T)

text(0, 0, paste0(round(aln_len/1000), "kb"), cex = 3)
