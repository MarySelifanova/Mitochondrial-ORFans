
library(tidyverse)
# library(Rbowtie2)
# library(Rsamtools)
library(RColorBrewer)
library(Biostrings)
library(IRanges)
library(rtracklayer)
library(circlize)
library(genbankr)
library(zoo)
library(egg)


my_theme <- theme_bw() +
  theme(axis.text.x = element_text(size=14, angle=00,hjust=0.95,vjust=0.5), 
        axis.text.y = element_text(size=14), 
        strip.text.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_blank(),
        legend.position = "none")   


  #Set working directory & import files
  setwd("/media/dima/KINGSTON/Papers_submition/0000_Polydora/Data/")
  dir()
  mtDNA <- readGenBank("submission1.gb", ret.seq = TRUE)

  
  genes <- genes(mtDNA)
  features <- otherFeatures(mtDNA)
  sequence <- getSeq(mtDNA)
  GenomeSize = length(sequence[[1]])
  getSeq(mtDNA)
  ############# import and analyse coverage
  coverage <- read.csv("depth.txt", sep = "\t") 
  names(coverage) <- c("seq", "pos", "count")
  
  coverage <- left_join(cbind(data.frame(pos = seq(1,GenomeSize))), 
                        coverage,  
                        by = "pos") %>% 
    mutate(count = ifelse(is.na(count), 0, count)) %>%
    mutate(low_cov = ifelse(count > 2, NA, 1))
  
  coverage_smooth <- tibble(coverage) %>% 
    mutate(smooth_cov = rollmean(count, 25, na.pad = TRUE)) %>%
    filter(pos %% 25 == 0) %>% 
    mutate(smooth_cov = ifelse(smooth_cov < 1, 1, smooth_cov))
  
   ############ #Analyze complete genome GCcontent v0.1
   ########### to fix: handling genome tail, overlapping window, ...
   window = 100 #set the window for GC content
   step = 100 # 

   steps <- GenomeSize %/% step
  
   GCcontent =  data.frame(pos = as.numeric(), 
                           GC = as.numeric(), 
                           ATskew = as.numeric(),
                           GCskew = as.numeric(),
                           coverage = as.numeric()
                           )
   for (i in seq(1, steps))
   {
   position <- i*(step)
   chunk <- sequence[[1]][(position-(window/2)):
                            min((position+(window/2)),
                                length(sequence[[1]]))]
   GC <- letterFrequency(chunk, "GC", as.prob = TRUE)
   Aa <- letterFrequency(chunk, "A")
   Tt <- letterFrequency(chunk, "T")
   Gg <- letterFrequency(chunk, "G")
   Cc <- letterFrequency(chunk, "C")
   ATskew <- (Aa-Tt)/(Aa+Tt) # AT skew = (A − T)/(A + T)
   GCskew <- (Gg-Cc)/(Gg+Cc)  # GC skew = (G − C)/(G + C)
   chunk_cov <- mean(coverage$count[(position-(window/2)):
                                     min((position+(window/2)),
                                         length(sequence[[1]]))])
   
   GCcontent <- rbind(GCcontent,
                      data.frame(pos = position, 
                                 GC = GC, 
                                 ATskew = ATskew, 
                                 GCskew = GCskew,
                                 coverage = max(chunk_cov, 1)))
   }
   
   GCcontent <-GCcontent %>% mutate(cumulativeATskew = cumsum(GCcontent$ATskew)) %>%
      mutate(cumulativeGCskew = cumsum(GCcontent$GCskew))
   
   
   #####parse MtDNA genomic features

   
   genes_t <- as.data.frame(genes) %>%
     mutate(product = gene_id) %>%
     select(start, end, width, strand, type, product) 
     
   mtDNAfeatures <- as.data.frame(features) %>%
     mutate(product = ifelse(is.na(product), allele, product)) %>% 
     select(start, end, width, strand, type, product) %>%
     rbind(genes_t) %>%
     mutate(type = factor(type))
     
  FeatureColor <- brewer.pal(length(levels(mtDNAfeatures$type)),"Set2") %>%
      set_names(levels(mtDNAfeatures$type))

  
############################## 
################ Circular maps
############################## 

### BioCircus
 circos.clear()
 circos.par("gap.degree" = 0, 
             "cell.padding" = c(0, 0, 0, 0), 
             "start.degree" = 90,
             "track.height" = 0.1)
 circos.initialize("a", xlim = c(0, GenomeSize))
  # circos.track(ylim = c(0,1),bg.border = "red")
    # circos.text(mtDNAfeatures$coordinates,0.7,mtDNAfeatures$X_Gene, niceFacing = TRUE, facing = "clockwise", cex = 0.5)

  BEDmtDNAfeatures <- mtDNAfeatures %>% 
     # filter(type != "tRNA_gene") %>% 
     # mutate(X_Gene = ifelse(type == "tRNA_gene", str_sub(X_Gene, 1,2),X_Gene)) %>% 
     mutate(chr = "a") %>% 
     select(chr, start, end, product)
  
  circos.genomicLabels(BEDmtDNAfeatures, 
                       labels.column = 4, 
                       side = "outside", 
                       # col = FeatureColor[mtDNAfeatures$type],
                       cex = 1)
  
  circos.info()
  
  circos.par("track.height" = 0.05)
  
  circos.track(ylim = c(0,1),bg.border = "white")
  for (row in 1:nrow(mtDNAfeatures)) 
    {arrowsize = ifelse(mtDNAfeatures$type[row] ==  "tRNA", 35, 100)
          circos.arrow(mtDNAfeatures$start[row],
          mtDNAfeatures$end[row], 
          y= 0.5,
          width = 0.8,
          arrow.head.width = 0.8,
          arrow.head.length =  arrowsize,
          arrow.position = ifelse(mtDNAfeatures$strand[row] == "+", "end", "start"),
          col = FeatureColor[mtDNAfeatures$type[row]]
          )
  }

  
  circos.par("gap.degree" = 0, 
             "cell.padding" = c(0, 0, 0, 0), 
             "start.degree" = 90,
             "track.height" = 0.15)

  ### Coverage track
  # # 
  circos.track(ylim = c(0,max(log2(coverage_smooth$smooth_cov))))
  circos.lines(x = coverage_smooth$pos, 
               y = log2(coverage_smooth$smooth_cov), 
               area = TRUE, 
               type = 'h', 
               col = 'darkgreen', lwd = 2)
  
  # circos.lines(x = coverage$pos,
  #              y = coverage$low_cov,
  #              area = TRUE,
  #              type = 'h',
  #              col = 'red',
  #              lwd = 1)
  par(cex = 0.8)
  circos.yaxis(side = "right", labels.cex = 0.5)


   # # GC track
  circos.par("gap.degree" = 0, 
             "cell.padding" = c(0, 0, 0, 0), 
             "start.degree" = 90,
             "track.height" = 0.1)
  circos.track(ylim = c(0,max(GCcontent$GC)))
  circos.lines(x = GCcontent$pos, y = GCcontent$GC, area = TRUE, type = 'h', col = 'darkblue', lwd = 2)
  par(cex = 0.8)
  circos.yaxis(side = "right", labels.cex = 0.5)


  # #AT skew track
 circos.track(ylim = c(min(GCcontent$GCskew, na.rm = TRUE),-min(GCcontent$GCskew, na.rm = TRUE)))
 circos.lines(x = GCcontent$pos, y = GCcontent$ATskew, type = 'l', col = 'orange', lwd = 2)
 circos.segments(0, 0, max(GCcontent$pos),0, type = "l", lwd = 1, lty = 2, col  = "red")
 circos.lines(x = GCcontent$pos, y = GCcontent$GCskew, type = 'l', col = 'red', lwd = 2)
 par(cex = 0.8)
 circos.yaxis(side = "left", labels.cex = 0.5)

################################
################# GENE EXPRESSION
################################

    
expression <- mtDNAfeatures[1:42,] %>%
    filter(type != "tRNA") %>%
    filter(!is.na(product)) %>% 
    group_by(product) %>%
    mutate(sum_count = sum(coverage$count[start:end])) %>% 
    mutate(expression = sum_count/width)
  
Figure5B <- expression %>% 
  ggplot(aes(x =  log2(expression), y = product)) + 
  geom_col() +
  my_theme

# expression %>% 
#   ggplot(aes(x = width, y =  log2(expression), label = product)) + 
#   geom_text() +
#   my_theme

FigureSXX1 <- coverage %>% 
  ggplot(aes(x= pos, y = count)) +
  geom_line() +
  xlim(1774-300, 2317+100) +
  geom_vline(xintercept = c(1774,2317), size = 1.3, lty = 3,color = "coral4") +
  ylim(0, 275) +
  my_theme

FigureSXX2 <- coverage %>% 
  ggplot(aes(x= pos, y = count)) +
  geom_line() +
  xlim(4671-100, 5730+100) +
  geom_vline(xintercept = c(4671,5274, 5730), size = 1.3, lty = 3,color = "coral4") +
  ylim(0, 275) +
  my_theme

FigureSXX3 <- coverage %>% 
  ggplot(aes(x= pos, y = count)) +
  geom_line() +
  xlim(11968-100, 12394+100) +
  geom_vline(xintercept = c(11968,12394), size = 1.3, lty = 3,color = "coral4") +
  ylim(0, 275) +
  my_theme

FigureNad4l <- coverage %>% 
  ggplot(aes(x= pos, y = count)) +
  geom_line() +
  xlim(14256-100, 14538+100) +
  geom_vline(xintercept = c(14256,14538), size = 1.3, lty = 3,color = "coral4") +
  ylim(0, 275) +
  my_theme

FigureATP6 <- coverage %>% 
  ggplot(aes(x= pos, y = count)) +
  geom_line() +
  scale_x_continuous(breaks = seq(11000, 12100,  by = 50), 
                     limits = c(11170-100, 11967+100)) + 
  geom_vline(xintercept = c(11170,11967), size = 1.3, lty = 3,color = "coral4") +
  ylim(0, 1100) +
  my_theme

 ggarrange(FigureSXX1, FigureSXX2, FigureSXX3, FigureNad4l)