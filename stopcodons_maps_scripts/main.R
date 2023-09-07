source("utils.R")
source("shannon.R")
source("codon_utils.R")
source("plot_utils.R")

seqstat_window <- 200 # set the window for GC content
information_window <- 100
step <- 100

read_data <- function(x, y) {
  x <- as.list(x)
  genome <- readDNAStringSet(x$sequence)[[1]]
  c(list(seq=genome, annot=read_annot(x$annotation, length(genome))), as.list(y))
    
}

count_stats <- function(sample_data, reverse=F, codons=NULL) {
  if (!reverse) {
    stop_codons <- c("TAA", "TAG")
  } else { stop_codons <- c("TTA", "CTA")}
  
  if (!is.null(codons)) {
    stop_codons <- codons
  }
  
  sample_data$genome_len <- length(sample_data$seq)
  stops <- count_stops(sample_data$seq)
  stops$codons %>%
    filter(codon %in% stop_codons) %>%
    mutate(chunk = cut(i, seq(1, sample_data$genome_len, step), include.lowest = T)) %>%
    group_by(chunk, frame) %>%
    summarise(nstop=n()) %>%
    mutate(pos=as.integer(chunk)*step) %>%
    ungroup() %>%
    select(pos, frame, nstop) %>% 
    pivot_wider(names_from = frame, names_prefix="frame_", values_from = nstop, values_fill = 0) -> stop_data
    
  
  seqstat_chunks <- get_chunks(sample_data$seq, step, seqstat_window)
  seq_stats <- decorate_pos(seqstat_chunks, count_seq_stats)
  
  info_chunks <- get_chunks(sample_data$seq, information_window, step)
  info_stats <- decorate_pos(info_chunks, count_2_mers)
  
  all_stats <- left_join(seq_stats, info_stats, by="pos") %>% left_join(stop_data, by="pos")
  all_stats %>% mutate(across(starts_with("frame"), ~ replace_na(.x, 0))) -> all_stats
  sample_data$seq <- NULL
  sample_data$stats <- all_stats
  sample_data
}

draw_graphs <- function(sample_data) {
  svglite(file.path("output", "stop_plots", paste0(sample_data$short_name, ".svg")))
  plot_stops(sample_data)
  dev.off()
  
  svglite(file.path("output", "coverage_plots", paste0(sample_data$short_name, ".svg")))
  plot_seqstats(sample_data)
  dev.off()
}

data_cor <- read_tsv("data_correspondance.tsv")
dir.create("output", showWarnings = F)
dir.create("output/coverage_plots", showWarnings = F)
dir.create("output/stop_plots", showWarnings = F)
dir.create("output/trn_punct_data", showWarnings = F)
dir.create("output/trn_punct_plots", showWarnings = F)
dir.create("output/counted_data", showWarnings = F)

data_cor %>%
  group_by(name, short_name) %>%
  group_map(read_data) %>%
  map(count_stats, codons=c("TTA", "TTG", "CTA", "CTT", "CTG", "CTC")) -> data

data %>%
  lapply(function(sample_data) {
    write_tsv(sample_data$stats, file.path("output", "counted_data", paste0(sample_data$short_name, "_reverse", ".tsv")))
  })

data %>%
  lapply(draw_graphs)
