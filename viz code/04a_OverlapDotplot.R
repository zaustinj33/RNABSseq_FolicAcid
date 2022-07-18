library(ggpubr)
library(ggplot2)
library(MASS)
library(ggsci)
library(dplyr)

## Dotplot showing similarity of methylation level between two replicates for 
# all transcripts
# set directories if necessary #
# parent directory
save_dir <- paste0(getwd(),"/figures")

#Overlap from coprehensive union file
all_m5C <- read.csv("allm5C_libraries_0MLfilteredDepthAnno.csv")

## MANUALLY CHANGE SAVE FILE NAME AND MT FILTER IN FUNCTION ##
conditions <- c("tLF","tMF","tHF", "pLF","pMF","pHF")

# color density function
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

## dotplot ##
overlap_files <- list.files(pattern = "_allOverlap.csv", full.names = T)
#View(read.csv(overlap_files[[1]]))

## Comparing replicates
replicate_dotplot <- function (overlap) {
  name <- gsub("_allOverlap.csv","",basename(overlap))
  all_overlap <- read.csv(overlap)

  overlap_df <- all_overlap[,grepl(paste0(name), colnames(all_overlap))]
  # only grab methylLevel columns
  methylDF <- data.frame(group=all_overlap$group,
                         ML_1=as.numeric(overlap_df[[3]]),
                         ML_2=as.numeric(overlap_df[[6]])
  )
  methylDF <- methylDF[methylDF[,2] >= 0.1 | methylDF[,3] >= 0.1,]
  methylDF <- methylDF[complete.cases(methylDF),]
  print(nrow(methylDF))
  methylDF$density <- get_density(methylDF$ML_1, methylDF$ML_2, n = nrow(methylDF))
  methyl_anno <- left_join(all_overlap, methylDF, by='group')
  write.table(methyl_anno, paste0("Dotplot_", name, "_table.csv"), sep=',', quote = F, row.names = F)


  #methylDF$shape <- ifelse((methylDF$ML_1 > 0.1 & methylDF$ML_2 > 0.1), 1,0)

  p<-ggplot(methylDF) +
    geom_point(aes(x = ML_1, y = ML_2, color = density),size = 3) + theme_bw() +
    geom_abline(slope = 1, linetype = 'longdash', size = 1.5) +
    scale_color_gradientn(colors = c('blue','yellow','#DC0000FF'), values = c(0,0.5,1)) +
    theme(legend.title = element_blank(), legend.position = 'none', axis.text = element_text(size = 18),
          text = element_text(size = 18, colour = 'black'), axis.line = element_line(size = 1),
          panel.grid=element_blank(), panel.border = element_rect(size = 2), axis.ticks.length=unit(.25, "cm"),
          plot.margin = margin(1,1,1,1,'cm')) +
    scale_x_continuous(expand = c(0,0), breaks = c(0,0.5,1),limits = c(0:1)) +
    scale_y_continuous(expand = c(0,0), breaks = c(0,0.5,1),limits = c(0:1)) +
    #ggtitle(name) +
    xlab("m5c level of rep1") + 
    ylab("m5c level of rep2") +
    stat_cor(method="pearson",aes(x = ML_1, y = ML_2),label.x.npc = 0.05)
  p
  ggsave(path = save_dir,filename = paste0("/Dotplot_", name, "_true.png")
    ,plot = p,height = 5,width = 5,dpi = 400)

  return(methyl_anno)
}
plot <- replicate_dotplot(overlap_files[1])

lapply(overlap_files, replicate_dotplot)

## DMS of conditions -> volcano plot

## Comparing conditions
condition_dotplot <- function (overlap1,overlap2) {
  name1 <- gsub("overlap.*$","",basename(overlap1))
  name2 <- gsub("overlap.*$","",basename(overlap2))
  
  overlap1 <- read.csv(overlap1)
  overlap2 <- read.csv(overlap2)
  
  # only grab methylLevel columns
  methylDF <- full_join(overlap1[,c(7,8)],overlap2[,c(7,8)],by='group')
  #methylDF <- methylDF[methylDF[,2] > 0.1 | methylDF[,3] > 0.1,]
  methylDF <- methylDF[complete.cases(methylDF),]
  methylDF <- full_join(methylDF,info_0v2[,c(1,6)],by='group')
  methylDF$status[is.na(methylDF$status)] <- "None"
  
  methylDF$status <- factor(methylDF$status, 
                            levels = c('upReg_hyperMeth','upReg_hypoMeth',
                                       'downReg_hyperMeth','downReg_hypoMeth',
                                       'hypoMeth_only',"hyperMeth_only","None"))
  
  #methylDF$density <- get_density(methylDF[[2]], methylDF[[3]], n = nrow(methylDF))
  pearson <- cor.test(methylDF[[2]], methylDF[[3]], method="pearson") 
  methylDF$alpha <- ifelse(methylDF$status == 'None',1,1)
  
  p<-ggplot(methylDF) + 
    geom_point(aes(x = methylDF[[2]], y = methylDF[[3]], color = status, alpha=alpha),size = 2) + theme_bw() + 
    #geom_abline(slope = 1, linetype = 'longdash', size = 1.5) +
    geom_smooth(aes(x = methylDF[[2]], y = methylDF[[3]]), method = 'lm') +
    #scale_color_gradientn(colors = c('blue','yellow','#DC0000FF'), values = c(0,0.5,1)) +
    scale_color_manual(values = c('red','purple','green','blue','yellow','cyan', '#e0e0e0')) +
    theme(legend.title = element_blank(), legend.position = 'none', axis.text = element_text(size = 18),
          text = element_text(size = 18, colour = 'black'), axis.line = element_line(size = 1),
          panel.grid=element_blank(), panel.border = element_rect(size = 2), axis.ticks.length=unit(.25, "cm"),
          plot.margin = margin(1,1,1,1,'cm')) +
    scale_x_continuous(expand = c(0,0), breaks = c(0,0.5,1),limits = c(0:1)) +
    scale_y_continuous(expand = c(0,0), breaks = c(0,0.5,1),limits = c(0:1)) +
    ggtitle(paste(name1,'v',name2)) +
    xlab(paste("m5c level of",name1)) + 
    ylab(paste("m5c level of",name2)) +
    stat_cor(method="pearson",aes(x = methylDF[[2]], y = methylDF[[3]]),label.x.npc = 0.05)
  p
  ggsave(path = save_dir,filename =  paste0("/Dotplot_",
                                            name1,'v',name2, "_depthDMSlabels.png"),plot = p,height = 5,width = 5,dpi = 500)
}
plot <- condition_dotplot(overlap_files[1],overlap_files[2])
plot <- condition_dotplot(overlap_files[1],overlap_files[3])
plot <- condition_dotplot(overlap_files[3],overlap_files[2])

