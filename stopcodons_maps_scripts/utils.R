library(tidyverse)
library(Biostrings)

compl <- list(A="T", C="G", G="C", "T"="A")

read_annot <- function(annot_path, GenomeSize) {
  annotation <- read.table(annot_path, sep = "\t")
  
  colnames(annotation) <- c("node", "start", "end", "name", "xxx",
                            "strand")
  annot_bed <- annotation %>%
    mutate(chr = "a") %>%
    select(chr, start, end, name, strand)
  
  trna <- sapply(annot_bed$name, function(x) {substr(x, 1, 3) == "trn"})
  rrna <- sapply(annot_bed$name, function(x) {substr(x, 1, 3) == "rrn"})
  origin <- sapply(annot_bed$name, function(x) {substr(x, 1, 2) == "OH"})
  CR <- sapply(annot_bed$name, function(x) {x == "CR"})
  orfan <- sapply(annot_bed$name, function(x) {substr(x, 1, 3) == "orf"})
  annot_bed$type <- "protein"
  annot_bed$type[trna] <- "trna"
  annot_bed$type[rrna] <- "rrna"
  annot_bed$type[origin] <- "origin"
  annot_bed$type[CR] <- "control"
  annot_bed$type[orfan] <- "orfan"
  annot_bed$name[trna] <- substring(annot_bed$name[trna], 4)
  annot_bed %>%
    as_tibble() %>%
    mutate(end = ifelse(end < start, end + GenomeSize, end)) %>%
    mutate(start = ifelse(end >= GenomeSize, start - GenomeSize, start), end = ifelse(end >= GenomeSize, end - GenomeSize, end)) %>%
    arrange(start) %>%
    mutate(start = start + 1)
}


save_bed <- function(rranges, filename) {
  for_bed <- rranges
  for_bed$start <- for_bed$start - 1
  write_tsv(for_bed, filename, col_names = F)
}

circ_slice <- function(genome, start, end) {
  if (end <= length(genome)) {
    as.character(genome[start:end])
  } else {
    paste0(as.character(genome[start:length(genome)]), genome[1:(end %% length(genome))])
  }
}

get_nsr_seq <- function(genome, nsr) { # nsr - no-stop region
  nsr %>%
    group_by(name) %>%
    mutate(sequence=circ_slice(genome, start, end))
}

count_seq_stats <- function(seq) {
  let_table <- table(as.vector(seq))
  a <- let_table["A"]
  t <- let_table["T"]
  g <- let_table["G"]
  c <- let_table["C"]
  tibble(gc_content=(g + c)/sum(let_table), at_skew=(a - t)/(a + t), gc_skew=(g - c)/(g + c))
}

get_chunks <- function(genome, step, window) {
  genome_size <- length(genome)
  steps <- genome_size %/% step
  
  lapply(seq(1, steps), function(i) {
    position <- i*(step)
    circ_range <- (position-(window/2)):(position+(window/2))
    list(pos=position, seq=genome[((circ_range - 1) %% genome_size) + 1])
  })
}


decorate_pos <- function(chunklist, f) {
  lapply(chunklist, function(chunk) {
    as_tibble(f(chunk$seq)) %>%
      mutate(pos=chunk$pos)
  }) %>% bind_rows
}