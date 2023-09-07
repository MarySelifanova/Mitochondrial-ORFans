orfan_seqs <- readDNAStringSet("orfan_seqs/orfan1cut.fasta")

orfan_files <- c("orfan_seqs/orfan1cut.fasta", "orfan_seqs/orfan2cut.fasta", "orfan_seqs/orfan3cut.fasta")

stops_f <- c("TAA", "TAG")
stops_r <- c("TTA", "CTA")
lapply(1:length(orfan_files), function(j) {
  orfan_seqs <- readDNAStringSet(orfan_files[j])
  lapply(1:length(orfan_seqs), function(i) {
    stops <- count_stops(orfan_seqs[[i]])
    stops$codons %>%
      mutate(is_stop= codon %in% stops_f, is_stop_r= codon %in% stops_r) %>%
      group_by(frame) %>%
      summarise(nstops_f=100*sum(is_stop)/n(), nstops_r=100*sum(is_stop_r)/n()) %>%
      mutate(species=names(orfan_seqs)[[i]], orfan=paste0("orfan", j)) 
  }) -> stop_table
  
  do.call(rbind, stop_table)
}) -> stop_table_list

stop_table <- do.call(rbind, stop_table_list)
write_tsv(stop_table, "orfan_seqs/stopcount.tsv")
