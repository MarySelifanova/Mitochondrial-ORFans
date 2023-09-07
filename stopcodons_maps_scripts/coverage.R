library(tidyverse)

get_cov_data <- function(cov_path, genome_size, step, window) {
  coverage <- read_delim(cov_path, delim = "\t", col_names = F)
  
  colnames(coverage) <- c("name", "position", "coverage")
  
  breaks <- (0:ceiling(GenomeSize/step))*step
  
  coverage %>%
    mutate(chunk=cut(position, breaks)) %>%
    group_by(chunk, .drop=F) %>%
    summarize(mean_cov=mean(coverage), z_cov=(sd(coverage)/mean(coverage))) %>%
    replace_na(list(mean_cov=0, z_cov=0)) %>%
    mutate(chunk_center=(head(breaks, -1)+step/2))-> avg_cov
  
  
}