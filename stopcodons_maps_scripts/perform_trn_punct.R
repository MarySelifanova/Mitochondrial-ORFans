source("sample_data.R")

data_cor <- read_tsv("data_correspondance.tsv")

data_cor

count_stop_codons <- function(chunk, reverse=FALSE) {
  stops <- c("TAA", "TAG")
  if (reverse) stops <- c("TTA", "CTA")
  chunk %>%
    mutate(is_stop=codon %in% stops) %>%
    group_by(frame) %>%
    summarise(nstop=sum(is_stop)) -> res
  
  ret <- res$nstop %>% as.list()
  names(ret) <- as.character(res$frame)
  ret
}

data_cor %>%
  rowwise() %>%
  group_split() %>%
  map(function(drow) {
    WGSData(drow$sequence) %>%
      consume_annot(base_annot=drow$annotation) %>%
      compute_codon_table() }) -> codon_frames

codon_frames %>%
  map(function(x) {
    x %>%
      nsr_annot(threshold = 300) %>%
      trn_punct(base_annot, nsr_annot, threshold = 300, adjust_start = TRUE)
  }) -> frames

frames %>%
  map(function(x) {
    x%>%
      init_summary() %>%
      pass_annots(base_annot, base_annot.t.nsr_annot) %>% 
      window_codon(nstops=count_stop_codons, .step = 50, .window = 50, .flatten = TRUE) %>%
      pop_summary()
  }) -> summaries


source("plot_utils.R")

single_plot <- function(short_name, name, summ) {
  print(summ)
  file.path("output", "trn_punct_plots", paste0(short_name, ".svg")) %>% svglite()
  plot_stops(summ)
  plot_annot(annotations(summ)$base_annot.t.nsr_annot, labels = "inside")
  title(name)
  dev.off()
}

data_cor %>%
  rowwise() %>%
  group_split() %>%
  map2(summaries, ~ single_plot(.x$short_name, .x$name, .y))

data_cor %>%
  rowwise() %>%
  group_split() %>%
  map2(summaries, function(.x, .y) {
    save_bed(annotations(.y)$base_annot.t.nsr_annot, paste0("output/trn_punct_data/", .x$short_name, ".bed"))
  })
