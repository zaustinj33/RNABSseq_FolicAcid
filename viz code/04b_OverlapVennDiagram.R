library(grid)
library(VennDiagram)
library(scales)
library(gplots)

## Venn diagram to compare m5C replicates and overlap

# set directories #

all_m5C <- read.csv("allm5C_libraries_6Cfiltered.csv")
all_m5C_depth <- read.csv("allm5C_libraries_filteredDepth.csv")
conditions <- c("tLF","tMF","tHF", "pLF","pMF","pHF")

write_overlap_files <- function(condition){
  total_df <- all_m5C
  rownames(total_df) <- total_df$group
  print(condition)
  # only sites present in both replicates, or only in one replicate
  subset_df <- total_df[,grepl(paste0("C_count_",condition), colnames(total_df))]
  colnames(subset_df) <- gsub('C_count_','',colnames(subset_df))

  # remove rows with all NA
  subset_df <- subset_df[rowSums(is.na(subset_df)) != ncol(subset_df), ]
  input <- lapply(subset_df, function(i) (row.names(subset_df)[ !is.na(i) ]))

  png(file=paste0("figures/",condition,"_6Cvenn.png"))
  venn_table <- venn(input, intersections = TRUE)
  dev.off()

  intersections <- attr(venn_table,"intersections")
  union <- unlist(intersections)

  # all overlap sites depth. change intersections[[X]] to match desired overlap
  # first subset all annotated depth sites by sites only present in condition
  rownames(all_m5C_depth) <- all_m5C_depth$group
  subset_depth <- all_m5C_depth[,grepl(paste0("C_count_",condition), colnames(all_m5C_depth))]
  subset_depth <- subset_depth[rowSums(is.na(subset_depth)) != ncol(subset_depth), ]
  subset_depth$group <- rownames(subset_depth)
  # only sites present in both replicates, or only in one replicate

  OL <- intersections[[1]]
  all_int <- subset(total_df, group %in% OL)
  print(paste0("Overlap: ", nrow(all_int)))
  print(paste0("Rep1 OL rate: ", length(OL)/(length(intersections[[2]])+length(OL))))
  print(paste0("Rep2 OL rate: ", length(OL)/(length(intersections[[3]])+length(OL))))
  write.csv(all_int, paste0(condition,"_alloverlap.csv"), row.names = F)

  # Write union site file
  all_union <- subset(total_df, group %in% union)
  print(paste0("Union: ", nrow(all_union)))
  #write.csv(all_union, paste0(condition,"_allUnion.csv"), row.names = F)

  # only sites passing filters in at least one replicate
  passing_sites <- subset(all_m5C_depth, group %in% union)
  passing_sites <- passing_sites[,grepl(paste0("C_count_",condition), colnames(passing_sites))]
  passing_sites[passing_sites == 0] <- NA
  passing_sites <- passing_sites[rowSums(is.na(passing_sites)) != ncol(passing_sites),]

  # Write venn of depth sites
  input <- lapply(passing_sites, function(i) (row.names(passing_sites)[ !is.na(i) ]))

  png(file=paste0(condition,"_venn_depth.png"))
  venn_table <- venn(input, intersections = TRUE)
  dev.off()

  # Write union with depth
  all_unionDepth <- subset(all_m5C_depth, group %in% row.names(passing_sites))
  write.csv(all_unionDepth, paste0(condition,"_allUnionDepth.csv"), row.names = F)

  passing_sites <- passing_sites[complete.cases(passing_sites),]

  # Write all intersection with depth
  all_depth <- subset(all_m5C_depth, group %in% row.names(passing_sites))
  print(paste0("Overlap depth: ", nrow(all_depth)))
  write.csv(all_depth, paste0(condition,"_allOverlapDepth.csv"), row.names = F)

  return(intersections)
}
draw_venns <- lapply(conditions, write_overlap_files)

test <- write_overlap_files("tHF")
