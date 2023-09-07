library(tidyverse)

exp_freqs <- function(kmer, symuse) {
  exp_freq <- 1
  for (i in 1:2) {
    exp_freq <- exp_freq*symuse[kmer[i]]
  }
  
  return(exp_freq)
}

shen <- function(freqs) {
  sum(sapply((freqs[freqs != 0])/sum(freqs), function(x) {
    -x*log2(x)
  }))
}

all_dimers <- c()
bases <- c("A", "T", "G", "C")
for (i in 1:4) {
  for (j in 1:4) {
    di <- paste0(bases[i], bases[j])
    all_dimers <- c(all_dimers, di)
  }
} 

count_2_mers <- function(seq) {
  symuse <- table(as.vector(seq))
  symuse <- symuse/sum(symuse)
  k_mers <- c()
  for (i in 1:(length(seq) - 1)) {
    k_mers <- c(k_mers, do.call(paste0, as.list(seq[i:(i + 1)])))
  }
  k_mers <- factor(k_mers, levels = all_dimers)
  counts <- table(k_mers)
  back_freqs <- sapply(names(counts), function(x) {
    exp_freqs(strsplit(x, "")[[1]], symuse)
  })
  names(back_freqs) <- names(counts)
  expected_counts <- (length(seq) - 1)*back_freqs
  tibble(ent_loss=(shen(back_freqs) - shen(counts/sum(counts))), 
              chi=sum((counts - expected_counts)^2), 
              cg_rat=counts["CG"]/expected_counts["CG"])
}