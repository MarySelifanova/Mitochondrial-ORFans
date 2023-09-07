library(Biostrings)
library(tidyverse)
library(rlang)
library(reshape2)
source("utils.R")
source("codon_utils.R")
source("sample_summary.R")

# This is definition of S4 classes for storing genomes, it's annotations, coverage, codon tables (which are required for some types of
# analysis and take some time to compute) and performing analysis of them.

#========================

setClassUnion("MaybeTibble", members = c("tbl_df", "NULL"))
setClassUnion("MaybeSummary", members=c("SampleSummary", "NULL"))

setClass("SeqData", slots = c(data = "DNAStringSet", summary="MaybeSummary"), prototype = list(summary=NULL))
setClass("SSData", contains = "SeqData", slots = c(coverages="list", 
                                                   codon_table="MaybeTibble"), 
         prototype = list(coverages=list(), codon_table=NULL))
setClass("WGSData", contains = "SSData", slots = c(annotations="list"), prototype = list(annotations=list()))

#====================

setGeneric("seq", function(x) standardGeneric("seq"))
setGeneric("coverages", function(x) standardGeneric("coverages"))
setGeneric("has_codon_table", function(x) standardGeneric("has_codon_table"))
setGeneric("codon_table", function(x) standardGeneric("codon_table"))


setMethod("seq", "SSData", function(x) (x@data)[[1]])
setMethod("coverages", "SSData", function(x) {x@coverages})
setMethod("has_codon_table", "SSData", function(x) {!is.null(x@codon_table)})
setMethod("codon_table", "SSData", function(x) {x@codon_table})

setGeneric("annotations", function(x) standardGeneric("annotations"))
setMethod("annotations", "WGSData", function(x) x@annotations)
setMethod("annotations", "SampleSummary", function(x) x@annotations)

setMethod("genlen", "SSData", function(x) seq(x) %>% length())
#====================

to_crange <- function(start, end, len) {
  ((start:end - 1) %% len) + 1
}

setGeneric("subset_seq", function(x, start, end) standardGeneric("subset_seq"))
setMethod("subset_seq", "SSData", function(x, start, end) {
  crange <- to_crange(start, end, genlen(x))
  seq(x)[crange]
})

setGeneric("subset_codon_table", function(x, start, end) standardGeneric("subset_codon_table"))
setMethod("subset_codon_table", "SSData", function(x, start, end) {
  lrange <- start:end
  ret <- codon_table(x)[to_crange(start, end, genlen(x)),]
  
  ret$frame <- (lrange - 1) %% 3
  ret
})

setMethod("show", "SSData", function(x) {
  cat("This is object of a custom dora class for storing sequence data, it's metadata and computations results\n")
  cat(paste0("Seq length: ", length(seq(x)), "\n"))
  cat(paste0("Coverages: ", names(coverages(x)), "\n"))
  cat(paste0("Codon table status: ", if (has_codon_table(x)) {"computed"} else {"uncomputed"}, "\n"))
  invisible(x)
})

setMethod("show", "WGSData", function(x) {
  callNextMethod()
  ann_names <- annotations(x) %>% names() %>% as.list()
  paste0("Annotations: ", do.call(paste, ann_names)) %>% cat()
  invisible(x)
})

#====================

setGeneric("consume_annot", function(x, ...) standardGeneric("consume_annot"))

setMethod("consume_annot", "WGSData", function(x, ...) {
  x@annotations <- modifyList(x@annotations, 
                              map(list(...), function(fpath) {
                       read_annot(fpath, length(seq(x)))
                     }))
  x
})

#====================

setGeneric("consume_coverage", function(x, ...) standardGeneric("consume_coverage"))
setMethod("consume_coverage", "SSData", function(x, ...) {
  x@coverages <- modifyList(x@coverages, 
                            map(list(...), function(fpath) {
                              read_delim(fpath, delim = "\t", col_names = F)
                            }))
  x
})

#====================

setGeneric("compute_codon_table", function(x) standardGeneric("compute_codon_table"))
setMethod("compute_codon_table", "SSData", function(x) {
  x@codon_table <- count_stops(seq(x))
  x
})

#====================

setGeneric("nsr_annot", function(x, stop_codons=c("TAA", "TAG"), threshold=330, name="nsr_annot", circular=TRUE) standardGeneric("nsr_annot"))
setMethod("nsr_annot", "WGSData", function(x, stop_codons=c("TAA", "TAG"), threshold=330, name="nsr_annot", circular=TRUE) {
  if (!has_codon_table(x)) {
    stop("Cannot perform algorithm without computed codon table")
  }
  
  codon_table(x) %>%
    filter(codon %in% stop_codons) -> stops
  
  if (circular) {
    genlen <- genlen(x)
    stops %>%
      group_by(frame) %>%
      summarise(codon=last(codon), i=last(i)) %>%
      mutate(i = i - genlen) %>%
      mutate(frame = (i - 1) %% 3) %>%
      arrange(i) -> virt_stops
    
    stops <- bind_rows(virt_stops, stops)
  }
  stops %>%
    group_by(frame) %>%
    group_modify(~ {
      .x %>%
        mutate(start=lag(i), end=i, length = end - start + 1) %>%
        filter(!is.na(start)) %>%
        select(-i)
    }) %>%
    ungroup() %>%
    mutate(quan=order(length, decreasing=TRUE)/n()) %>%
    filter(length >= threshold) %>%
    arrange(start) %>%
    mutate(chr="a", type="protein", strand="+", name=as.character(row_number())) %>%
    select(chr, start, end, name, strand, type, quan) -> nsr_annot
  
  adlist <- list(nsr_annot)
  names(adlist) <- name
  x@annotations <- modifyList(x@annotations, adlist)
  x
  
})

#=============================


#====================

setGeneric("init_summary", function(x) standardGeneric("init_summary"))
setMethod("init_summary", "SSData", function(x) {
  x@summary <- new("SampleSummary", genlen = seq(x) %>% length())
  x
})

#====================

setGeneric("pass_annots", function(x, ...) standardGeneric("pass_annots"))
setMethod("pass_annots", "WGSData", function(x, ...) {
  expressions <- enquos(...) %>% eval_tidy(annotations(x))
  names <- map_chr(expressions, as_name)
  annots <- annotations(x)[names]
  x@summary@annotations <- modifyList(x@summary@annotations, annots)
  x
})

#====================

setGeneric("get_windows", function(x, step, window) standardGeneric("get_windows"))
setMethod("get_windows", "SSData", function(x, step=100, window=200) {
  sequence(length(seq(x)) %/% step, from=1, by=step) %>%
    map(\(i) c((i - window/2) + 1, i + window/2))
})

window_func <- function(funclist, flatten=FALSE) {
  function(.seq, i) {
    map(funclist, \(f) f(.seq)) -> listed_res
    
    if (flatten) {
      listed_res %>% unlist(recursive = FALSE) %>% as.list() -> listed_res
    }
    
    listed_res %>%
      modifyList(list(.i=i)) %>%
      as_tibble()
  }
}

window_func_factory <- function(f) {
  function(x, ..., .step=100, .window=200, .flatten=FALSE) {
    arglist <- list(...)
    get_windows(x, .step, .window) %>%
      map(\(range) f(range, x)) %>%
      imap(window_func(arglist, flatten=.flatten)) %>%
      bind_rows() %>%
      mutate(pos=(.i - 1)*.step + 1) %>%
      select(!.i) -> short_res
    
    short_res %>%
      pivot_longer(!pos, names_to = "stat") %>%
      select(stat, pos, value) -> stat_results
    
    new_columns <- short_res %>% select(-pos) %>% colnames()
    
    x@summary@seq_stats %>% bind_rows(stat_results) -> x@summary@seq_stats
    win_appendix <- rep(.window, length(new_columns)) %>% as.integer()
    names(win_appendix) <- new_columns
    
    step_appendix <- rep(.step, length(new_columns)) %>% as.integer()
    names(step_appendix) <- new_columns
    
    x@summary@seqstat_windows %>% c(win_appendix) -> x@summary@seqstat_windows
    x@summary@seqstat_steps %>% c(step_appendix) -> x@summary@seqstat_steps
    x
  }
}

setGeneric("window_seq", function(x, ..., .step=100, .window=200, .flatten=FALSE) standardGeneric("window_seq"))
setMethod("window_seq", "SSData", window_func_factory(\(range, x) x %>% subset_seq(range[1], range[2])))

#=========================
setGeneric("window_codon", function(x, ..., .step=100, .window=200, .flatten=FALSE) standardGeneric("window_codon"))
setMethod("window_codon", "SSData", window_func_factory(\(range, x) x %>% subset_codon_table(range[1], range[2])))

#====================

setGeneric("pop_summary", function(x) standardGeneric("pop_summary"))
setMethod("pop_summary", "SeqData", function(x) {
  x@summary
})

#====================

setGeneric("trn_punct", function(x, punctor, punctee, threshold=330, adjust_start=FALSE) standardGeneric("trn_punct"))

get_overlapping_names <- function(prot_annot, trn_annot) {
  prot_annot %>%
    rowwise() %>%
    group_split() %>%
    map(function(protrow) {
      trn_annot %>%
        filter(start < protrow$end & end > protrow$start) -> filt
      c(filt$name, protrow$name)
    })
}

counter_fac <- function(counter_type) {
  force(counter_type)
  function(acc, nxt) {
    if (nxt$type == counter_type & nxt$coord_type == "start") {
      acc <- acc + 1
    }
    
    if (nxt$type == counter_type & nxt$coord_type == "end") {
      acc <- acc - 1
    }
    
    acc
  }
}

find_protname <- function(annot) {
  annot %>%
    filter(type == "protein") %>%
    summarise(name = dplyr::first(name)) %>%
    unlist()
}

setMethod("trn_punct", "WGSData", function(x, punctor, punctee, threshold=330, adjust_start=FALSE) {
  punctor_frame <- enquo(punctor) %>% eval_tidy(annotations(x))
  punctee_frame <- enquo(punctee) %>% eval_tidy(annotations(x))
  
  gen_len <- genlen(x)
  
  trn_frame <- punctor_frame %>% filter(type=="trna")
  nsr_frame <- punctee_frame %>% filter(type=="protein")
  
  double_trn <- bind_rows(trn_frame %>% mutate(start=start-gen_len, end=end-gen_len, name=paste0(name, "-")), trn_frame)
  mixed_frame <- bind_rows(double_trn, nsr_frame)
  
  mixed_frame %>%
    select(start, end, name, type) %>%
    pivot_longer(start:end, names_to = "coord_type") %>%
    arrange(value) -> longer_frame
  
  groups_names <- get_overlapping_names(nsr_frame, double_trn)
  groups_names %>%
    map(\(names) longer_frame %>% filter(name %in% names)) %>%
    map(function(x) {
      x %>%
        rowwise() %>%
        group_split() -> chunk_list
      
      chunk_list %>%
        accumulate(counter_fac("trna"),
                   .init = 0, .dir = "forward") -> betw_trn
      
      chunk_list %>%
        accumulate(counter_fac("protein"),
                   .init = 0, .dir = "forward") -> betw_prot
      
      x$betw_trn <- betw_trn[-1]
      x$betw_prot <- betw_prot[-1]
      
      x %>%
        mutate(next_coord=lead(value)) %>%
        filter(betw_trn == 0 & !is.na(next_coord) & betw_prot == 1) %>%
        select(value, next_coord) %>%
        mutate(name=paste0(find_protname(x), ".", row_number())) -> ret
      
      protstart <- x %>% filter(type == "protein" & coord_type == "start") %>% select(value) %>% as.integer()
      
      ret %>%
        mutate(value = value + ((protstart - value) %% 3))
      
    }) %>%
    bind_rows() %>%
    filter(next_coord - value >= threshold) %>%
    dplyr::rename(start=value, end=next_coord) %>%
    mutate(chr="a", strand="+", type="protein") %>%
    select(chr, start, end, name, strand, type) -> annot_res
  
  if (adjust_start) {
    start_codons <- c("TTG", "ATA", "ATG", "ATC", "ATT", "GTG")
    annot_res %>%
      rowwise() %>%
      group_split() %>%
      map_int(function(annot_row) {
        subset_codon_table(x, annot_row$start, annot_row$end) %>%
        filter(frame == (annot_row$start - 1) %% 3 & codon %in% start_codons) %>%
        select(i) %>% dplyr::first() %>% as.integer()
      }) -> new_starts
    
    
    annot_res$start <- ((new_starts - annot_res$start) %% gen_len) + annot_res$start
    
    annot_res <- annot_res %>% filter((end - start) + 1 >= threshold)
  }
  
  uplist <- list(annot_res)
  names(uplist) <- paste0(as.character(enquo(punctor))[2], ".t.", as.character(enquo(punctee))[2])
  x@annotations <- modifyList(x@annotations, uplist)
  x
})

WGSData <- function(seqpath) {
  new("WGSData", data=readDNAStringSet(seqpath))
}

