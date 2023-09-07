library(circlize)
library(tidyverse)
library(Cairo)
library(svglite)

init_without_sectors <- function(genome_len) {
  circos.clear()
  circos.par("gap.degree" = 0,
             "cell.padding" = c(0, 0, 0, 0),
             "start.degree" = 90,
             "track.height" = 0.1)
  circos.initialize("a", xlim = c(0, genome_len))
}

init_with_sectors <- function(annotation, genome_len) {
  annotation %>%
    filter(type %in% c("protein", "orfan")) %>%
    select(start, end) %>%
    gather() %>%
    select(value) %>%
    arrange(value) %>%
    filter(value - lag(value, default = 1) > 50) %>%
    mutate(prev_val = lag(value, default = 1), len=value - prev_val) -> breaks
  
  sizes <- breaks$len
  lim_mat <- breaks %>% 
    select(prev_val, value) %>% 
    as.matrix()
  
  alphabet <- strsplit("abcdefghijkmlnopqrstuvwxyz123456789~.,'/[]{}", "")[[1]]
  sector_names <- alphabet[1:length(sizes)]
  names(sizes) <- sector_names
  rownames(lim_mat) <- sector_names
  circos.clear()
  circos.par("gap.degree" = 5,
             "cell.padding" = c(0, 0, 0, 0),
             "start.degree" = 90,
             "track.height" = 0.1)
  sectors <- factor(sector_names, levels = sector_names)
  circos.initialize(sectors=sectors, x=sizes)
}

plot_annot <- function(annotation, genlen, labels="outside") {
  par(cex=1)
  colormap <- c(protein="pink", trna="yellow", rrna="gray", origin="green", control="blue", orfan="magenta")
  
  if (labels == "outside") {
    circos.genomicLabels(annotation,
                         labels.column = 4,
                         side = "outside")
  }
  
  circos.info()
  
  circos.par("track.height" = 0.05)
  
  circos.track(ylim = c(0,1), bg.border = "white")
  for (row in 1:nrow(annotation)) {
    start <- annotation$start[row]
    end <- annotation$end[row]
    
    circos.arrow(start,
                 end,
                 y= 0.5,
                 width = 0.8,
                 arrow.head.width = 0.8,
                 arrow.head.length = 100,
                 arrow.position =
                   ifelse(annotation$strand[row] == "+", "end", "start"),
                 col = colormap[annotation$type[row]]
    )
  }
  
  if (labels == "inside") {
    circos.genomicLabels(annotation,
                         labels.column = 4,
                         side = "inside")
  }
  
  
  circos.par("track.height" = 0.1)
}

get_plotting_function <- function(f, setup = function() {}) {
  function(samp_summary) {
    genlen(samp_summary) %>% init_without_sectors()
    plot_annot(annotations(samp_summary)$base_annot, genlen(samp_summary))
    setup()
    f(samp_summary)
  }
}

#============================

plot_stop_track <- function(pos, nstops, ylim, track_idx.r, gen_len, color="red") {
  circos.track(ylim=c(0, ylim))
  circos.lines(x=pos, y=ylim*(nstops > 0), area=T, type = 's', col="lightgray", lwd=0.5)
  circos.lines(x = pos, y = nstops, area = TRUE, type =
                 's', col = color, lwd = 0.5)
  
  circos.text(100, ylim/2, as.character(track_idx.r))
  circos.text(gen_len - 100, ylim/2, (track_idx.r - gen_len) %% 3 %>% as.character())
  circos.yaxis(side = "right", at=sequence(0, ylim, length.out = 3))
}

plot_stops <- get_plotting_function(function(samp_summary) {
  stop_summary <- samp_summary %>% get_stats(nstops.0, nstops.1, nstops.2)
  
  stop_summary %>%
    select(pos, starts_with("nstops")) %>%
    pivot_longer(!pos) %>%
    summarise(max = max(value)) -> summary
  
  max_stop <- summary$max[1]
  stop_summary %>%
    mutate(across(starts_with("nstops"), function(x) {
      plot_stop_track(pos, x, max_stop, cur_column() %>% 
                                          strsplit("\\.") %>%
                                          unlist() %>%
                                          last() %>%
                                          as.integer(), genlen(samp_summary))
      
    }))
}, setup = function() {
  circos.par(track.height=0.1, track.margin=c(0,0))
})

#==========================

plot_seqstats <- get_plotting_function(function(sample_data) {
  seqstats <- sample_data$stats
  # GC track
  circos.track(ylim = c(0,max(seqstats$gc_content)))
  circos.lines(x = seqstats$pos, y = seqstats$gc_content, area = TRUE, type =
                 'h', col = 'darkblue', lwd = 2)
  
  par(cex = 0.8)
  circos.yaxis(side = "right", at=seq(0, 0.4, length.out = 3))
  
  # pyrine-pyrimidine skew track
  ylim <- c(seqstats$gc_skew, seqstats$at_skew) %>% abs() %>% max()
  
  circos.track(ylim = c(-ylim, ylim))
  circos.lines(x = seqstats$pos, y = seqstats$at_skew, type = 'l', col
               = 'orange', lwd = 2)
  circos.segments(0, 0, max(seqstats$pos), 0, type = "l", lwd = 1, lty =
                    2, col  = "red")
  circos.lines(x = seqstats$pos, y = seqstats$gc_skew, type = 'l', col
               = 'red', lwd = 2)
  par(cex = 0.8)
  circos.yaxis(side = "left", at=seq(-0.4, 0.4, length.out = 3))
})


