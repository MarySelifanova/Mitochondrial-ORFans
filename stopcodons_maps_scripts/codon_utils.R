library(tidyverse)

function(i) {
  lin_rg <- i:(i+2)
  genome[((lin_rg - 1) %% length(genome)) + 1]
}

count_stops <- function(genome) {
  map(1:length(genome), \(i) i:(i+2)) %>%
    map(\(i) ((i - 1) %% length(genome)) + 1) %>%
    map(\(i) genome[i]) %>%
    map_chr(as.character) %>%
    imap(\(x, idx) tibble(codon=x, i=idx)) %>%
    bind_rows() %>%
    mutate(frame = as_integer((i - 1) %% 3))
}

search_potential_orfs <- function(stops, stop_codons) {
  stops$codons %>%
    filter(codon %in% stop_codons) %>%
    group_by(frame) -> stop_codons
  
  stop_codons %>%
    summarise(coord=first(i)) -> first_stops
  
  stop_codons %>%
    group_modify(function(x, y) {
      first <- (x$i)[1]
      x %>%
        mutate(start=i+3, end=lead(i, default = first_stops$coord[((y$frame - GenomeSize) %% 3) + 1] + GenomeSize) + 2, delta=end-start)
    }) %>%
    filter(delta >= 330) %>%
    ungroup() %>%
    arrange(start) %>%
    filter(end > lag(end, default = 0)) %>%
    mutate(chr="a", strand="+", name=1:n(), type="protein") %>%
    select(chr, start, end, name, strand, type)
}