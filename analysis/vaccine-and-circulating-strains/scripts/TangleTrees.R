#!/usr/bin/env Rscript

library(ape)
library(dplyr)
library(phytools)
library(dendextend)
library(magrittr)
library(lubridate)
library(rprojroot)
library(gridExtra)

# exports
source(sprintf("%s/utils.R", "scripts"))


# read metadata
prepare_metadata <- function(metadata_fn, lineages_fn) {
  
  # read metadata
  meta <- read.csv(metadata_fn, sep=",", header = TRUE)
  meta$taxa <- paste(meta$accession, meta$host,  meta$country, meta$date, sep="|")
  meta$date <- as.Date(as.character(meta$date), format = "%Y-%m-%d")
  meta$year <- sapply(strsplit(as.character(meta$date), "-", fixed = TRUE), `[`, 1)
  meta$sample_date <- decimal_date(meta$date)
  
  # fix strain labelled 13 as Clone 13
  meta$strain[meta$strain == '13'] <- 'Clone 13'
  
  # read lineages file
  lineages <- read.csv(lineages_fn, sep=",", header = TRUE)
  
  # merge
  metadata <- left_join(meta, lineages, by=c("accession"="Query"))

  
  return(metadata)
}


metadata_large <- prepare_metadata(
  metadata_fn = file.path(segmentsDir[1], "shared-strains/data/RVFV-L.csv"),
  lineages_fn = file.path(segmentsDir[1], "assignment/output-dir/report/lineages.csv")
  )

metadata_medium <- prepare_metadata(
  metadata_fn = file.path(segmentsDir[2], "shared-strains/data/RVFV-M.csv"),
  lineages_fn = file.path(segmentsDir[2], "assignment/output-dir/report/lineages.csv")
  )

metadata_np <- prepare_metadata(
  metadata_fn = file.path(segmentsDir[3], "shared-strains/data/RVFV-S.csv"),
  lineages_fn = file.path(segmentsDir[3], "assignment/output-dir/report/lineages.csv")
  )

metadata_nss <- prepare_metadata(
  metadata_fn = file.path(segmentsDir[3], "shared-strains/data/RVFV-S.csv"),
  lineages_fn = file.path(segmentsDir[3], "assignment/output-dir/report/lineages.csv")
  )



# function to read tree
read_tree <- function(tree_fn, metadata){
  tree <- read.tree(tree_fn)
  tree <- phytools::midpoint.root(tree)
  
  # match labels
  meta <- metadata[match(tree$tip.label, metadata$taxa), ]
  print(all(tree$tip.label == meta$taxa))
  tree$tip.label <- meta$strain
  # tree <- chronoMPL(tree)
  print(ape::is.ultrametric(tree))
  result <- list(tree, meta)
  return(result)
}

# read trees
tree_large <- read_tree(
  # tree = file.path(segmentsDir[1], "output-dir/riftL/iqtree/riftL.treefile"),
  tree = file.path(segmentsDir[1], "shared-strains/iqtree/RVFV-L.treefile"),
  metadata = metadata_large
  )
treeL <- tree_large[[1]]


tree_medium <- read_tree(
  # tree = file.path(segmentsDir[2], "output-dir/riftM/iqtree/riftM.treefile"),
  tree = file.path(segmentsDir[2], "shared-strains/iqtree/RVFV-M.treefile"),
  metadata = metadata_medium
)
treeM <- tree_medium[[1]]


tree_np <- read_tree(
  # tree = file.path(segmentsDir[3], "output-dir/riftS-NP/iqtree/riftS-NP.treefile"),
  tree = file.path(segmentsDir[3], "shared-strains/iqtree/RVFV-S-NP.treefile"),
  metadata = metadata_np
)
treeNP <- tree_np[[1]]


tree_nss <- read_tree(
  # tree = file.path(segmentsDir[3], "output-dir/riftS-NSS/iqtree/riftS-NSS.treefile"),
  tree = file.path(segmentsDir[3], "shared-strains/iqtree/RVFV-S-NSS.treefile"),
  metadata = metadata_nss
)
treeNSS <- tree_nss[[1]]

# get common tip labels
common_labels <- Reduce(intersect, list(treeL$tip.label, 
                                        treeM$tip.label,
                                        treeNP$tip.label,
                                        treeNSS$tip.label
                                        ))

if(length(common_labels) < 2) stop("Error: Fewer than 2 common taxa found across all trees.")

trees <- list(treeL, treeM, treeNP, treeNSS)


# finds tips containing vaccine names and colors their ancestral branches
color_vaccine_branches <- function(dend, 
                                   vaccine_pattern = "Smithburn|MP-12|Clone 13|vaccine", 
                                   color = "darkred", lwd = 3) {
  # Identify indices of labels that match the vaccine names
  labels_dend <- labels(dend)
  vaccine_indices <- grep(vaccine_pattern, labels_dend, ignore.case = TRUE)
  
  if(length(vaccine_indices) > 0) {
    # Color the specific branches leading to these tips
    dend <- color_branches(dend, labels = labels_dend[vaccine_indices], col = color, lwd = lwd)
  }
  return(dend)
}


# subset on common tips, normalize and color-code vaccine strains
process_tree <- function(tr, tips) {
  tr_sub <- keep.tip(tr, tips)
  dend <- as.dendrogram(as.hclust(compute.brlen(tr_sub, method = "Grafen")))
  
  # Set base branch width
  dend <- set(dend, "branches_lwd", 1.5)
  
  # Apply vaccine color coding
  dend <- color_vaccine_branches(dend)
  
  return(dend)
}


# convert to dendrograms for comparison
dends <- lapply(trees, process_tree, tips = common_labels)



# Calculate Global Entanglement
calculate_global_entanglement <- function(dend_list) {
  # Calculate pairwise entanglement for the sequence (L-M, M-SNSs, SNSs-SNP)
  e_lm   <- entanglement(dendlist(dend_list[[1]], dend_list[[2]]))
  e_ms   <- entanglement(dendlist(dend_list[[2]], dend_list[[3]]))
  e_ls   <- entanglement(dendlist(dend_list[[1]], dend_list[[3]]))
  
  # Global score is the mean of the chain
  global_score <- mean(c(e_lm, e_ms, e_ls))
  print(global_score)
  print(e_lm)
  print(e_ms)
  print(e_ls)
  
  return(list(
    Global_Score = global_score,
    Pairwise_Scores = c(L_vs_M = e_lm, M_vs_SNSs = e_ms, L_vs_SNSs = e_ls)
  ))
}

# untangle the chain to get the most conservatuve (lowest) entanglement scores

dl12 <- dendlist(dends[[1]], dends[[2]]) %>% dendextend::untangle(method = "step1side")
dl23 <- dendlist(dl12[[2]], dends[[3]]) %>% dendextend::untangle(method = "step1side")
dl34 <- dendlist(dl23[[2]], dends[[4]]) %>% dendextend::untangle(method = "step1side")

final_dends <- list(dl12[[1]], dl12[[2]], dl23[[2]], dl34[[2]])
results <- calculate_global_entanglement(final_dends)


# calculate individual Entanglement Scores
e_lm <- entanglement(dendlist(final_dends[[1]], final_dends[[2]]))
e_ms <- entanglement(dendlist(final_dends[[2]], final_dends[[3]]))
e_ls <- entanglement(dendlist(final_dends[[1]], final_dends[[3]]))

print(round(results$Pairwise_Scores, 4))



# plot tanglegrams
plot_panel <- function(d1, d2, main_l, main_r, score) {
  tanglegram(d1, d2, 
             main_left = main_l, 
             main_right = main_r,
             sub = paste("Entanglement Score:", round(score, 4)),
             common_subtrees_color_lines = TRUE, 
             highlight_distinct_edges = TRUE,
             highlight_branches_lwd = FALSE,
             lab.cex = 0.5, 
             margin_inner = 14, 
             columns_width = c(4, 2, 4),
             lwd = 1.2,
             cex_sub = 1.3)
}


pdf(file.path(outDir, "RVFV_Tanglegram.pdf"), width = 28, height = 12)
par(mfrow = c(1, 3), mar = c(10, 2, 10, 2))


plot_panel(final_dends[[1]], final_dends[[2]], "L", "M", e_lm)
plot_panel(final_dends[[2]], final_dends[[3]], "M", "S-NSs", e_ms)
plot_panel(final_dends[[3]], final_dends[[4]], "L", "S-NSs", e_ls)

dev.off()



