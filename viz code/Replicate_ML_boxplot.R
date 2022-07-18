library(ggpubr)
library(ggplot2)
library(ggsci)
library(dplyr)

overlap_Depth_files <- list.files(pattern = "_allOverlapDepth.csv", full.names = T)
overlap_files <- list.files(pattern = "Overlap.csv", full.names = T)

conditions <- c("tLF","tMF","tHF", "pLF","pMF","pHF")

get_ML_info <- function(input_table, name){
  raw_info <- read.csv(input_table)
  ML_info <- raw_info[,grepl(name, colnames(raw_info))]
  ML_table <- data.frame(group = raw_info$group,
                         ML_1 = ML_info[[3]],
                         ML_2 = ML_info[[6]])
  return(ML_table)
}

merge_lists <- function(L) {
  pad.na <- function(x,len) {
    c(x,rep(NA,len-length(x)))
  }
  maxlen <- max(sapply(L,length))
  df <- do.call(data.frame,lapply(L,pad.na,len=maxlen))
  colnames(df) <- c("rep1_OL","rep2_OL","rep1_nonOL", "rep2_nonOL")
  #df$rep2_nonOL[!is.na(df$rep2_OL)] <- NA
  #df$rep1_nonOL[!is.na(df$rep1_OL)] <- NA
  return(df)
}

plot_replicateOL_ML  <- function(name){
  ## Get methylation level info of seperate groups (rep1 OL, rep2 OL, rep1 nonOL, rep2 nonOL)
  OL_name <- overlap_files[grepl(name, overlap_files)]
  OLdepth_name <- overlap_Depth_files[grepl(name, overlap_Depth_files)]
  OL_file <- get_ML_info(OL_name, name)
  OL_depth_file <- get_ML_info(OLdepth_name, name)
  non_overlap <- OL_depth_file[!OL_depth_file$group %in% OL_file$group,]

  # Add null column in case equal length?
  #non_overlap[nrow(non_overlap) + 1,] <- c('None', 0, 0)

  ## Merge lists to plottable DF
  plot_list <- list(OL_file$ML_1, OL_file$ML_2, non_overlap$ML_1, non_overlap$ML_2)
  plot_df <- stack(merge_lists(plot_list))

  ## set color
  mypal <- pal_npg("nrc", alpha = 0.7)(9)[c(1,3,2,4)]

  g <- ggplot(plot_df, aes(x=ind, y=values, fill=ind)) + geom_boxplot() + theme_classic() +
    scale_fill_manual(values=mypal)  +
    scale_y_continuous(limits = c(0,1,1), breaks = c(0,0.5,1)) +
    theme(axis.ticks.length=unit(.25, "cm")) +
    #scale_x_discrete(labels=c("total_ML" = paste('total',name), "poly_ML" = paste('poly',name))) +
    xlab("Fraction") +
    ylab("Methylation level")
  #stat_compare_means(comparisons = list(c('total_ML', 'poly_ML')), label.y = 0.9)
  g

  ggsave(filename = paste0(name,'all_depthOL_boxplot.png'),width = 4,height = 3,dpi = 200)


}

plot_replicateOL_ML('pMF')

lapply(conditions, plot_replicateOL_ML)



