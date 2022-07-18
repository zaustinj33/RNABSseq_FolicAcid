## plot RNABS v RNAseq comparisons
library(ggplot2)
library(ggpubr)
library(ggsci)
library(dplyr)
library(DESeq2)
library(tximportData)
library(tximport)
library(data.table)

## 1 RNAseq
mypal = pal_npg("nrc", alpha = 0.7)(9)[5:8]
## Load counts
RNAseq_files <- list.files(pattern = "genes.results", path='RNAseq', full.names = T)
sample_table <- data.frame('sample' = c('pHF_rep1','pHF_rep2','pLF_rep1','pLF_rep2', 'pMF_rep1','pMF_rep2',
                                        'tHF_rep1','tHF_rep2','tMF_rep1','tMF_rep2', 'tLF_rep1','tLF_rep2'),
                           'treatment' = rep(c(rep('HF',2), rep('MF',2), rep('LF',2)),2),
                           'fraction' = c(rep('poly', 6), rep('total', 6)))
sample_table$group <- paste0(sample_table$treatment, '_', sample_table$fraction)

txi.rsem <- tximport(RNAseq_files, type = "rsem", txIn = FALSE, txOut = FALSE)
head(txi.rsem$counts)
txi.rsem$length[txi.rsem$length <= 0] <- 1

## Perform normalization
dds <- DESeqDataSetFromTximport(txi.rsem, sample_table, ~group)
dds <- dds[rowSums(counts(dds)) > 0,]
dds <- DESeq(dds)
resultsNames(dds)
# Create normalized counts table
dds.counts <- counts(dds, normalized = T)
colnames(dds.counts) <- sample_table$sample
dds.counts <- dds.counts[,c(1:6,9:12,7:8)]

get_DEG_results <- function (dds_results_object) {
  DEG_table <- data.frame(as.data.table(dds_results_object))
  rownames(DEG_table) <- rownames(dds)
  return(DEG_table)
}

LF_res <- get_DEG_results(results(dds, contrast = c("group","LF_poly", "LF_total")))
MF_res <- get_DEG_results(results(dds, contrast = c("group","MF_poly", "MF_total")))
HF_res <- get_DEG_results(results(dds, contrast = c("group","HF_poly", "HF_total")))

plot_DEGs <- function(data) {
  data$change <-  as.factor(ifelse(data$padj < 0.05 &
                                     abs(data$log2FoldChange) > log2(1.2) ,
                                   ifelse(data$log2FoldChange > log2(1.2) ,'Up','Down'),'No Change'))

  print(paste0("Number DownReg: ", length(data$change[data$change  == 'Down'])))
  print(paste0("Number UpReg: ", length(data$change[data$change  == 'Up'])))

  p <- ggplot(data = data, aes(x = log2FoldChange, y = -log10(padj), color = change)) +
    geom_point(alpha=0.8, size = 2) + labs(title="") +
    theme_bw(base_size = 15) +
    scale_x_continuous(limits = c(-30,30),breaks = c(-30,-15,0,15,30)) +
    scale_y_continuous(limits = c(0,170)) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          panel.border = element_rect(size = 2),
          plot.title = element_text(hjust = 0.5),
          #axis.text = element_text(size = 34)
    ) + scale_color_manual(name = "", values = c("green", "red", "black"), limits = c("Up", "Down", "No Change")) +
    geom_hline(yintercept=1.30103,linetype=3)+geom_vline(xintercept=c(-log2(1.2),log2(1.2)),linetype=3) #+
  #ylab('-log(p.adj)')+
  #xlab('rlog(FC)')
  p
}
plot_DEGs(HF_res)
ggsave(filename = 'DEGs_HF_volcano.png',width = 6,height = 5,dpi = 200)

#write.csv(HF_res, "HF_DEGs.csv")


## 2 Methylation level
ML_files <- list.files(pattern = "Compare.*_TrueOLtable005.csv", full.names = T, recursive = T)
library(tibble)
reduce_ML_files <- function(ML_file, treatment){
  ML_table <- read.csv(ML_file)
  row.names(ML_table) <- ML_table$group
  ML_reduced <- ML_table %>% select(contains(treatment))
  ML_reduced <- left_join(rownames_to_column(ML_reduced),
                          rownames_to_column(as.data.frame(ML_table)), by = 'rowname')
  ML_reduced$ML_1 <- (ML_reduced[[3]] + ML_reduced[[6]])/ (ML_reduced[[2]] + ML_reduced[[5]])
  ML_reduced$ML_2 <- (ML_reduced[[9]] + ML_reduced[[12]])/ (ML_reduced[[8]] + ML_reduced[[11]])
  return(ML_reduced)
}

LF_ML <- (reduce_ML_files(ML_file = 'Compare_tLFpLF_TrueOLtable005.csv', 'LF'))
MF_ML <- (reduce_ML_files(ML_file = 'Compare_tMFpMF_TrueOLtable005.csv', 'MF'))
HF_ML <- (reduce_ML_files(ML_file = 'Compare_tHFpHF_TrueOLtable005.csv', 'HF'))

comparisons <- c("LF", "MF", "HF")

LF_compare <- left_join(LF_ML, read.csv('RNAseq/LF_DEGs.csv'),
                        by = c("gene" = "X"))

MF_compare <- left_join(MF_ML, read.csv('RNAseq/MF_DEGs.csv'),
                        by = c("gene" = "X"))

HF_compare <- left_join(HF_ML, read.csv('RNAseq/HF_DEGs.csv'),
                        by = c("gene" = "X"))

## 3: Annotate and plot

# annotate quadrants
create_plot_table <- function(all_data, ML_cutoff, LFC_cutoff, padj_cutoff){
  # ML filter 0.1
  #all_data <- all_data[all_data$ML_1 >= 0.1 | all_data$ML_1 >= 0.1,]
  all_data$ML_diff <- all_data$ML_2 - all_data$ML_1 # diff [poly - total]
  #all_data <- all_data[complete.cases(all_data),]

  all_data$DEGsig <- ifelse(all_data$padj < padj_cutoff, "DEGsig","DEGnonsig")
  all_data$allSig <- ifelse(all_data$shape == 'sig' & all_data$DEGsig == 'DEGsig', 'both',
                            ifelse(all_data$shape == 'non-sig' & all_data$DEGsig == 'DEGsig', 'DEGonly',
                                   ifelse(all_data$shape == 'sig' & all_data$DEGsig == 'DEGnonsig', 'DMSOnly',
                                          'not signficant or DEG')))

  all_data$status <- ifelse(all_data$allSig == 'both' & all_data$ML_diff > ML_cutoff & all_data$log2FoldChange > LFC_cutoff, 'pos_corr',
                            ifelse(all_data$allSig == 'both' & all_data$ML_diff < ML_cutoff & all_data$log2FoldChange < LFC_cutoff, 'pos_corr',
                            ifelse(all_data$allSig == 'both' & all_data$ML_diff > ML_cutoff & all_data$log2FoldChange < LFC_cutoff, 'neg_corr',
                            ifelse(all_data$allSig == 'both' & all_data$ML_diff < ML_cutoff & all_data$log2FoldChange > LFC_cutoff, 'neg_corr',
                                   ifelse(all_data$allSig == 'DMSOnly', 'DMSonly','no_corr')))))

  all_data$status <- factor(all_data$status,
                            levels = c('pos_corr', 'neg_corr', 'DMSonly'))

  return(all_data)
}
plot_LF <- create_plot_table(LF_compare, 0, 0, 0.05)
write.csv(plot_LF, 'LFC_ML_compare_LF005.csv')

plot_HF <- create_plot_table(HF_compare, 0, 0, 0.05)
write.csv(plot_HF, 'LFC_ML_compare_HF005.csv')

plot_MF <- create_plot_table(MF_compare, 0, 0, 0.05)
write.csv(plot_MF, 'LFC_ML_compare_MF005.csv')


## 4 Plot
plot_LF <- read.csv('LFC_ML_compare_LF005.csv')
plot_MF <- read.csv('LFC_ML_compare_MF005.csv')
plot_HF <- read.csv('LFC_ML_compare_HF005.csv')

plot_DMS_DEG <- function(plot_table, name){
  p <- ggplot(plot_table[!is.na(plot_table$status) & (plot_table$status != 'no_corr'),]) +
    geom_point(aes(
    x = ML_diff, y = log2FoldChange, color = status)) + #, shape = allSig
    geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +theme_bw() +
    scale_x_continuous(limits = c(-1,1)) + scale_y_continuous(limits = c(-3,3)) +
    scale_color_manual(values = c('#00FF00','#FF0000','#0000FF','#808080')) +
    scale_shape_manual(values = c(8,4,16,1)) +
    theme(
      panel.grid = element_blank(),
      axis.text = element_text(size=14),
      legend.text = element_text(size = 14)
    ) +
    ggtitle(paste0('LFC x ML; t', name, 'v p', name))
  p
  #ggsave(filename = paste0('Point_compareLFCvML_', name, '005.png'),height = 10,width = 15,units = 'cm',dpi = 400)

}
plot_DMS_DEG(plot_HF, 'HF')

library(tidyr)
plot_stacked_bar <- function(){
  plot_table<-as.data.frame(sapply(list(plot_LF$status[!is.na(plot_LF$status) & (plot_LF$status != 'no_corr')],
                                        plot_MF$status[!is.na(plot_MF$status) & (plot_MF$status != 'no_corr')],
                                        plot_HF$status[!is.na(plot_HF$status) & (plot_HF$status != 'no_corr')]), table))
  colnames(plot_table) <- c('LF','MF','HF')
  plot_table$feature <- row.names(plot_table)

  df <- plot_table %>%
    tidyr::pivot_longer(!feature, names_to='sample', values_to='count') %>%
    dplyr::group_by(sample) %>%
    mutate(Percent = count / sum(count)*100) %>%
    ungroup()
  df$sample <- factor(df$sample, levels = c('LF', 'MF', 'HF'))
  df$feature <- factor(df$feature, levels = c('pos_corr', 'neg_corr', 'DMSonly'))

  ggplot(df, aes(x=sample,y=Percent,fill=feature)) + geom_bar(position = 'stack',stat='identity') +
    scale_fill_manual(values = c('#00FF00','#FF0000','#0000FF','#808080')) + theme_bw() +
    #scale_x_discrete(expand = c(0,0))+
    scale_y_continuous(expand = c(0,0))+#,breaks = c(0,200,400),limits = c(0,400)) +
    geom_text(aes(label = paste(round(count,2))),position = position_stack(vjust=0.5)) +
    theme(panel.grid = element_blank(),
          panel.border = element_blank(), axis.line = element_line(size = 1.25),
          axis.text = element_text(size=16))

}
plot_stacked_bar()
ggsave(filename = 'Bar_compareLFCxML_sig005.png',height = 10,width = 12,units = 'cm',dpi = 400)

