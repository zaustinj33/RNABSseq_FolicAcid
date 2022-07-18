## Histogram of m5C sites in overlap and union between replicates
library(ggplot2)
# set directories #
# parent directory
# location of 5mC files (filtered, overlap, and union files saves here)

# desired save directory for figures
save_dir <- paste0(getwd(),"/figures/")
dir.exists(save_dir)

overlap_files <- list.files(pattern = "_allOverlap.csv$")
print(overlap_files[1])

# plot histogram
histogram_plot <- function (overlap_file) {
  overlap_tmp <- read.csv(overlap_file)
  print(nrow(overlap_tmp))
  overlap_name <- gsub("_allOverlap.csv","",overlap_file)
  # select only ML columns
  overlap_tmp <- overlap_tmp[, grepl(overlap_name, names(overlap_tmp))]
  overlap_tmp <- overlap_tmp[, grepl("meth", names(overlap_tmp))]
  # create ML average
  overlap_tmp$avg_ML <- round(((overlap_tmp[,1]+overlap_tmp[,2])/2), 2)
  # just in case any snuck through the 0.1 filter
  overlap_tmp <- overlap_tmp[overlap_tmp$avg_ML >= 0.1,]
  
  histo <- ggplot(overlap_tmp) + geom_histogram(aes(x=(avg_ML), y = (..density..)),
      binwidth = .1, fill="lightblue", color="steelblue", boundary=0, closed="right", 
      position=position_nudge(x=0, y=0)) + 
    scale_x_continuous(breaks=c(seq(0,1,by=.1)), limits = c(0.1,1),expand = c(0.1,0)) +
    scale_y_continuous(breaks=c(seq(0,8)), limits=c(0,8), expand = c(0,0)) + 
    theme_bw() +
    theme(panel.grid =element_blank()) + 
    theme(panel.border = element_blank()) + 
    theme(axis.line = element_line(size=0.5, colour = "black")) + 
    xlab(paste0("methylation level of ", paste(overlap_name, "overlap"))) + 
    ylab("Percentage(%)")
  print(histo)
  #dev.off()
  ggsave(path = save_dir,filename =  paste0("Histogram_", overlap_name, "_overlap.png"),
         width = 8, height = 8,units = "cm", plot = histo)
}

histogram_plot(overlap_files[4])
lapply(overlap_files, histogram_plot)

histogram_plot_union <- function (union_file) {
  union_tmp <- read.csv(paste0("../Final_Data/", union_file))
  union_name <- gsub("_union.csv$","",union_file)
  colnames(union_tmp) <- c("X","group","cov_1","count_1","ML_1","cov_2","count_2","ML_2")
  
  union_tmp$ML_1 <- as.numeric(as.character(union_tmp$ML_1))
  union_tmp$ML_2 <- as.numeric(as.character(union_tmp$ML_2))
  # hashout if not including ML=0 transcripts
  union_tmp[is.na(union_tmp)] <- 0
  
  union_tmp$ML_1 <- round(union_tmp$ML_1, digits = 2)
  union_tmp$ML_2 <- round(union_tmp$ML_2, digits = 2)
  
  df <- data.frame(c(union_tmp$group, union_tmp$group), c(union_tmp$ML_1, union_tmp$ML_2))
  names(df) <- c("SeqID", "methLevel")
  union_final <- df[complete.cases(df),]
  
  #union_tmp$methLevel <- round(union_tmp$C_count/union_tmp$cov, digits=2)
  
  histo_union <- ggplot(union_final) + geom_histogram(aes(x=(methLevel), y = (..density..)),binwidth = .1, fill="lightblue", color="steelblue", boundary=0, closed="right", position=position_nudge(x=0, y=0)) + 
    scale_x_continuous(breaks=c(seq(0.1,1,by=.1)), expand = c(0.1,0)) + scale_y_continuous(breaks=c(seq(0,7)), limits=c(0,7), expand = c(0,0)) + theme_bw() +
    theme(panel.grid =element_blank()) + theme(panel.border = element_blank()) + theme(axis.line = element_line(size=0.5, colour = "black")) + 
    xlab(paste0("methylation level of ", union_name, " union")) + ylab("Precentage(%)")
  print(histo_union)
}

histogram_plot_union(union_3)

