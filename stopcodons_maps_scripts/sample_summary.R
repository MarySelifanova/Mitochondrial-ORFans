

# This is definition of S4 class for storing summary of sample_data. It is not supposed to be used in analysis, it's only purpose is being
# a convenient storage for annotations, sequence (or coverage) by-window statistics.

setClass("SampleSummary", slots=c(genlen="integer", 
                                  annotations="list", 
                                  seq_stats="tbl_df", 
                                  seqstat_windows="integer", 
                                  seqstat_steps="integer"),
         prototype = list(annotations=list(), 
                          seq_stats=tibble(stat=character(), pos=integer(), value=double()), 
                          seqstat_windows=integer(),
                          seqstat_steps=integer()))

setMethod("show", "SampleSummary", function(object) {
  cat("A custom sample_data summary object for dora project\n")
  cat("Annotations: ")
  cat(do.call(paste, object@annotations %>% names() %>% as.list()))
  cat("\n")
  cat("Window stats: ")
  cat(do.call(paste, object@seqstat_windows %>% names() %>% as.list()))
  cat("\n")
})

setGeneric("get_stats", function(x, ...) standardGeneric("get_stats"))
setMethod("get_stats", "SampleSummary", function(x, ...) {
  expressions <- enquos(...) %>% eval_tidy(x@seqstat_steps)
  stat_names <- map_chr(expressions, as_name)
  
  steps <- x@seqstat_steps[stat_names]
  is_same <- all(steps == steps[1])
  if (!is_same) stop("all stats should have the same step")
  
  x@seq_stats %>%
    filter(stat %in% stat_names) %>%
    pivot_wider(id_cols = pos, names_from = stat)
})

setGeneric("genlen", function(x) standardGeneric("genlen"))
setMethod("genlen", "SampleSummary", function(x) x@genlen)
