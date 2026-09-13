#!/usr/bin/env Rscript

library(tidyverse)
library(ggplot2)
library(pheatmap)


# exports
source(sprintf("%s/utils.R", "scripts"))


# INPUT FILES
files <- list(
  Smithburn = sprintf(
    "%s/output-dir/riftM/strain-types/riftM.strain_type.DQ380193.mutations.per.strain.singleton.csv",
    segmentsDir[[2]]),
  MP12 = sprintf(
    "%s/output-dir/riftM/strain-types/riftM.strain_type.DQ380208.mutations.per.strain.singleton.csv",
    segmentsDir[[2]]),
  Clone13 = sprintf(
    "%s/output-dir/riftM/strain-types/riftM.strain_type.DQ380213.mutations.per.strain.singleton.csv",
    segmentsDir[[2]])
)

aa_size <- 1197
step <- 200

# READ AND COMBINE DATA
mutation_df <- lapply(names(files), function(v){
  df <- read.csv(files[[v]])
  df$vaccine <- v
  return(df)
}) %>% bind_rows()

# Extract amino acid position
mutation_df$position <- as.numeric(gsub("\\D+", "", mutation_df$snps))

# ----------------------------------------------------
# Define structural regions
# ----------------------------------------------------

regions <- data.frame(
  start = c(1,155,561,691,1121),
  end   = c(154,560,690,1120,1197),
  region = c("NSm","Gn","Gn-Gc junction","Gc","Gc tail")
)

# assign region
mutation_df$region <- NA

for(i in 1:nrow(regions)){
  mutation_df$region[mutation_df$position >= regions$start[i] &
                       mutation_df$position <= regions$end[i]] <- regions$region[i]
}

# ----------------------------------------------------
# Define neutralizing epitope clusters
# ----------------------------------------------------

epitopes <- data.frame(
  start = c(180,301,421),
  end   = c(300,420,500),
  epitope = c("Epitope1","Epitope2","Epitope3")
)

mutation_df$epitope <- NA

for(i in 1:nrow(epitopes)){
  mutation_df$epitope[mutation_df$position >= epitopes$start[i] &
                        mutation_df$position <= epitopes$end[i]] <- epitopes$epitope[i]
}

# =====================================================
# 1 Mutation density plot (hotspots)
# =====================================================

density_df <- mutation_df %>%
  group_by(position) %>%
  summarise(count=n())

p_density <- ggplot(density_df, aes(position,count)) +
  geom_line() +
  geom_point() +
  scale_x_continuous(limits=c(0,aa_size),
                     breaks=seq(0,aa_size,by=step)) +
  labs(x="Amino acid position",
       y="Mutation density",
       title="RVFV M-segment mutation hotspots") +
  theme_bw()

print(p_density)

# =====================================================
# 2 Epitope mapped mutations
# =====================================================

epitope_df <- mutation_df %>%
  filter(!is.na(epitope)) %>%
  group_by(epitope) %>%
  summarise(count=n())

p_epitope <- ggplot(epitope_df,aes(epitope,count,fill=epitope)) +
  geom_bar(stat="identity") +
  labs(x="Neutralizing epitope cluster",
       y="Mutation count",
       title="Mutations within Gn neutralizing epitopes") +
  theme_bw()

print(p_epitope)

# =====================================================
# 3 Host mutation heatmap
# =====================================================

heatmap_df <- mutation_df %>%
  group_by(host,position) %>%
  summarise(count=n()) %>%
  pivot_wider(names_from=position,values_from=count,values_fill=0)

heatmap_matrix <- as.matrix(heatmap_df[,-1])
rownames(heatmap_matrix) <- heatmap_df$host

pheatmap(heatmap_matrix,
         cluster_rows=TRUE,
         cluster_cols=TRUE,
         main="Host-specific mutation heatmap")

# =====================================================
# 4 Structural domain mutation distribution
# =====================================================

domain_df <- mutation_df %>%
  group_by(region) %>%
  summarise(count=n())

p_domain <- ggplot(domain_df,aes(region,count,fill=region)) +
  geom_bar(stat="identity") +
  labs(x="M segment structural domain",
       y="Mutation count",
       title="Mutation distribution across RVFV glycoprotein domains") +
  theme_bw()

print(p_domain)

# =====================================================
# Save plots
# =====================================================

# ggsave("mutation_hotspots.png",p_density,width=10,height=4)
# ggsave("epitope_mutations.png",p_epitope,width=6,height=4)
# ggsave("domain_mutations.png",p_domain,width=6,height=4)

print(mutation_df[mutation_df$vaccine=='Smithburn',] |> 
        group_by(positions,snps) |> summarise(n=n()) |> 
        filter(n>=10), n=100)
