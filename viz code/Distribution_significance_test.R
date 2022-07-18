library(ggplot2)
library(tidyr)
library(dplyr)

# Input data
OLdepth_files <- list.files(pattern = "_compare.*.csv$", full.names = T)  # Overlap with depth

# Get input distribution
# Plot distribution of transcript distribution's correlations
record_pvalues <- function(category, stats_df) {
  print(category)

  og_dist <- table(stats_df$status) / sum(table(stats_df$status))  # probability distribution
  test_dist <- table(stats_df$status[stats_df$status == category])  # number of successes

  # test 5UTR, exon, and 3UTR for significance
  p_table <- mapply(function(X, Y, dist_sum) {
    sapply(seq_along(X), function (row) binom.test(x = X[row], n = dist_sum, p = Y[row])$p.val)
  }, X=test_dist, Y=og_dist, dist_sum=sum(test_dist))
  return(list(p_table))
}

# save table of pvalues with counts, plot results
save_pvalues <- function(input_df){
  name <- gsub("LFC_ML_compare_","",basename(input_df))
  name <- gsub("005.csv", "", name)

  test <- read.csv(input_df)
  test <- test[!is.na(test$status),]  # dont change when switching type and status

  ## Counts
  dist_all <- test %>%
    group_by(type) %>%
    summarise(counts = n())
  dist_all$status <- 'all_sites'

  dist_summary <- test %>%
    group_by(status, type) %>%
    summarise(counts = n())
  dist_plot <- rbind(dist_all, dist_summary)
  # pivot by status - c("all_sites", "DMSonly","pos_corr", "neg_corr")
  dist_plot$type <- factor(dist_plot$type, levels = c("5UTR", "exon", "3UTR"))
  dist_plot$status <- factor(dist_plot$status, levels = c("all_sites", "pos_corr", "neg_corr", "DMSonly"))

  # Plot
  library(ggplot2)
  p <- ggplot(data = dist_plot, aes(x=status, y=counts, fill = type, label = counts)) +
    geom_bar(position = 'fill', stat = 'identity',  color = 'black') +
    theme_classic() +
    theme(axis.text.x=element_text(angle = -90)) +
    scale_fill_manual(values = c('yellow', 'orange', 'red')) + #c("#00FF00","#FF0000","#0000FF")+
    scale_y_continuous(expand = c(0, 0))
  p
  ggsave(paste0('distribution_plot_',name,'.png'), width = 5,height = 4,dpi = 300)

  ### Change type <-> status
  correlations <- unique(test$status)

  ## Binmonial test
  test_df <- as.data.frame(lapply(correlations, record_pvalues, stats_df = test))
  colnames(test_df) <-  correlations
  test_df$status <- row.names(test_df)

  # long form to merge with Counts
  save_pvals <- as.data.frame(pivot_longer(test_df, !status, names_to = "type", values_to = "pvalues"))
  save_stats <- left_join(dist_plot, save_pvals, by=c("status" = "status", "type" = "type"))

  write.csv(save_stats, paste0('distribution_pvals_',name,'.csv'), row.names = F)
  return(test_df)
}

save_stats <- lapply(OLdepth_files, save_pvalues)




