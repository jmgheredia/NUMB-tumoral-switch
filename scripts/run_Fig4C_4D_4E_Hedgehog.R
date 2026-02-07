# ============================================================
# run_Figure4C_Hedgehog_PCA_UMAP.R
#
# PCA, UMAP and clustering of Hedgehog pathway genes in PDMR
#
# Inputs:
# - data/processed/pdmr_isoform_expression_log2.csv
# - data/metadata/metadata_PDMR.csv
# - data/metadata/genes_by_pathway.csv
#
# Outputs (Figure 4C and related panels):
# - results/Figures/Figure_4/Hedgehog_pathway/*.tiff
#
# Intermediate:
# - results/processed/Figure_4/UMAP_coordinates_Hedgehog.csv
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(uwot)
  library(ggalluvial)
  library(cluster)
  library(RColorBrewer)
})

set.seed(123)

# ============================================================
# 0) Input / Output
# ============================================================

expr_file  <- file.path("data","processed","pdmr_isoform_expression_log2.csv")
meta_file  <- file.path("data","metadata","metadata_PDMR.csv")
genes_file <- file.path("data","metadata","genes_by_pathway.csv")

fig_dir  <- file.path("results","Figures","Figure_4","Hedgehog_pathway")
proc_dir <- file.path("results","processed","Figure_4")

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(proc_dir, recursive = TRUE, showWarnings = FALSE)

umap_cache_file <- file.path(proc_dir, "UMAP_coordinates_Hedgehog.csv")

stopifnot(
  file.exists(expr_file),
  file.exists(meta_file),
  file.exists(genes_file)
)

# ============================================
# 1. Load data
# ============================================

genes_by_pathway <- read.csv(genes_file, row.names = 1)
entire_dataset   <- read.csv(expr_file, check.names = FALSE)
metadata         <- read.csv(meta_file, sep = ";")

# ============================================
# 2. Clean IDs
# ============================================

entire_dataset$transcript_id <- sub("\\..*", "", entire_dataset$transcript_id)
entire_dataset$gene_id <- as.character(entire_dataset$gene_id)
genes_by_pathway$gene_symbol <- as.character(genes_by_pathway$gene_symbol)

metadata$PatientID.SpecimenID.SampleID <-
  trimws(as.character(metadata$PatientID.SpecimenID.SampleID))

colnames(entire_dataset) <- trimws(colnames(entire_dataset))

# ============================================
# 3. Merge + filter Hedgehog pathway genes
# ============================================

genes_transcripts_by_pathway <- merge(
  genes_by_pathway,
  entire_dataset,
  by.x = "gene_symbol",
  by.y = "gene_id"
)

Hedgehog_genes <- unique(subset(genes_transcripts_by_pathway,
                                pathway == "Hedgehog")$gene_symbol)

expression_Hedgehog <- entire_dataset %>%
  filter(gene_id %in% Hedgehog_genes)

# ============================================
# 4. Prepare data for PCA
# ============================================

rownames(expression_Hedgehog) <- expression_Hedgehog$transcript_id
expression_Hedgehog_clean <- expression_Hedgehog[, !(names(expression_Hedgehog) %in% c("gene_id", "transcript_id"))]
expression_Hedgehog_t <- t(expression_Hedgehog_clean)

expression_Hedgehog_t <- expression_Hedgehog_t[order(rownames(expression_Hedgehog_t)), ]

# ============================================
# 5. PCA
# ============================================

pca_Hedgehog <- prcomp(expression_Hedgehog_t, scale. = TRUE)

metadata_ordered <- metadata[
  match(rownames(expression_Hedgehog_t), metadata$PatientID.SpecimenID.SampleID),
]

pca_df <- as.data.frame(pca_Hedgehog$x)

pca_df$Passage <- metadata_ordered$Passage_of_this_sample
pca_df$Passage <- as.character(pca_df$Passage)
pca_df$Passage[pca_df$Passage %in% c("5","6","8","9")] <- "+5"
pca_df$Passage <- factor(pca_df$Passage,
                         levels=c("Resection","0","1","2","3","4","+5"))

var_exp <- (pca_Hedgehog$sdev)^2
var_exp_perc <- round(100 * var_exp / sum(var_exp), 1)

# ============================================
# 6. Plot PCA
# ============================================

p_pca <- ggplot(pca_df, aes(PC1, PC2, color = Passage)) +
  geom_point(size=3) +
  theme_minimal() +
  labs(
    title="PCA of Hedgehog pathway gene expression",
    x=paste0("PC1 (",var_exp_perc[1],"%)"),
    y=paste0("PC2 (",var_exp_perc[2],"%)")
  ) +
  scale_color_manual(values=c(
    "Resection"="black","0"="orange","1"="darkgoldenrod",
    "2"="green","3"="blue","4"="purple","+5"="red"
  )) +
  theme(plot.title=element_text(hjust=.5), panel.grid=element_blank())

ggsave(
  filename = file.path(fig_dir, "PCA_Hedgehog.tiff"),
  plot = p_pca,
  width = 6,
  height = 5,
  dpi = 300
)

# ============================================
# 7. UMAP (cached)
# ============================================

if (!file.exists(umap_cache_file)) {
  
  umap_result <- umap(
    expression_Hedgehog_t,
    n_neighbors = 15,
    min_dist = 0.1,
    metric = "euclidean",
    n_threads = 1,
    init = "spectral",
    ret_model = FALSE
  )
  
  umap_df <- as.data.frame(umap_result)
  colnames(umap_df) <- c("UMAP1","UMAP2")
  umap_df$SampleID <- rownames(expression_Hedgehog_t)
  
  write.csv(umap_df, umap_cache_file, row.names = FALSE)
  
} else {
  
  umap_df <- read.csv(umap_cache_file)
  
}

umap_df$Passage <- metadata_ordered$Passage_of_this_sample
umap_df$Passage <- as.character(umap_df$Passage)
umap_df$Passage[umap_df$Passage %in% c("5","6","8","9")] <- "+5"
umap_df$Passage <- factor(umap_df$Passage,
                          levels=c("Resection","0","1","2","3","4","+5"))

# ============================================
# 8. Plot UMAP
# ============================================

p_umap <- ggplot(umap_df, aes(UMAP1,UMAP2,color=Passage)) +
  geom_point(size=3) +
  theme_minimal() +
  labs(title="UMAP of Hedgehog pathway gene expression") +
  scale_color_manual(values=c(
    "Resection"="black","0"="orange","1"="darkgoldenrod",
    "2"="green","3"="blue","4"="purple","+5"="red"
  )) +
  theme(plot.title=element_text(hjust=.5), panel.grid=element_blank())

ggsave(
  filename = file.path(fig_dir, "UMAP_Hedgehog.tiff"),
  plot = p_umap,
  width = 6,
  height = 5,
  dpi = 300
)

# ============================================
# 9.1 Estimate number of clusters
# ============================================

wss <- sapply(1:10, function(k) {
  kmeans(umap_df[, 1:2], centers = k, nstart = 10)$tot.withinss
})

tiff(file.path(fig_dir, "Elbow_Kmeans_Hedgehog.tiff"),
     width = 6, height = 5, units = "in", res = 300)
plot(1:10, wss, type = "b", pch = 19, frame = FALSE,
     xlab = "Number of clusters K",
     ylab = "Total within-cluster sum of squares (WSS)")
dev.off()

silhouette_scores <- sapply(2:10, function(k) {
  km <- kmeans(umap_df[, 1:2], centers = k, nstart = 10)
  ss <- silhouette(km$cluster, dist(umap_df[, 1:2]))
  mean(ss[, 3])
})

tiff(file.path(fig_dir, "Silhouette_Kmeans_Hedgehog.tiff"),
     width = 6, height = 5, units = "in", res = 300)
plot(2:10, silhouette_scores, type = "b", pch = 19, frame = FALSE,
     xlab = "Number of clusters K",
     ylab = "Mean silhouette width",
     main = "K selection based on the silhouette index")
dev.off()

optimal_k <- which(diff(diff(wss)) == min(diff(diff(wss)))) - 2

# ============================================
# 9.2 k-means clustering
# ============================================

set.seed(123)
kmeans_result <- kmeans(umap_df[, 1:2], centers = optimal_k, nstart = 10)
umap_df$Cluster <- as.factor(kmeans_result$cluster)

if (optimal_k <= 8) {
  palette_name <- "Set2"
} else if (optimal_k <= 12) {
  palette_name <- "Set3"
} else {
  stop("Too many clusters for the predefined RColorBrewer palettes.")
}

cluster_colors_cb <- setNames(
  brewer.pal(n = optimal_k, name = palette_name)[1:optimal_k],
  levels(umap_df$Cluster)
)

p_clusters <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Cluster)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "k-means clusters on Hedgehog UMAP") +
  scale_color_manual(values = cluster_colors_cb) +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank())

ggsave(file.path(fig_dir, "UMAP_Clusters_Hedgehog.tiff"),
       p_clusters, width = 6, height = 5, dpi = 300)

p_clusters_ell <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Cluster)) +
  geom_point(size = 3) +
  stat_ellipse(type = "t", linetype = 2, size = 1) +
  theme_minimal() +
  labs(title = "k-means clusters on Hedgehog UMAP (with ellipses)") +
  scale_color_manual(values = cluster_colors_cb) +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank())

ggsave(file.path(fig_dir, "UMAP_Clusters_Hedgehog_Ellipses.tiff"),
       p_clusters_ell, width = 6, height = 5, dpi = 300)

# ============================================
# 9.3 Automatic cluster ordering
# ============================================

cluster_order <- umap_df %>%
  group_by(Cluster) %>%
  summarise(
    pct_resection = mean(Passage=="Resection"),
    mean_umap1 = mean(UMAP1),
    .groups="drop"
  ) %>%
  arrange(desc(pct_resection), mean_umap1)

umap_df$Cluster_ordered <- factor(
  umap_df$Cluster,
  levels = cluster_order$Cluster,
  labels = seq_len(optimal_k)
)

# ============================================
# 10. Composition barplot
# ============================================

cluster_table_norm <- umap_df %>%
  count(Cluster_ordered, Passage) %>%
  group_by(Cluster_ordered) %>%
  mutate(percentage = n/sum(n)*100)

p_bar <- ggplot(cluster_table_norm,
                aes(Cluster_ordered,percentage,fill=Passage)) +
  geom_bar(stat="identity") +
  theme_minimal() +
  scale_fill_manual(values=c(
    "Resection"="darkblue","0"="orange","1"="darkgoldenrod",
    "2"="green","3"="blue","4"="purple","+5"="red"
  )) +
  labs(title="Passage percentage distribution within each cluster",
       x="Cluster",y="Percentage of samples") +
  theme(plot.title=element_text(hjust=.5))

ggsave(file.path(fig_dir,"Distribution_Passages_Clusters_Hedgehog.tiff"),
       p_bar, width = 6, height = 5, dpi = 300)

# ============================================
# 11. Alluvial
# ============================================

umap_df$Passage <- recode(umap_df$Passage,"Resection"="R")

alluvial_data <- umap_df %>% count(Passage,Cluster_ordered)

custom_colors <- c(
  "R"="darkblue","0"="blue","1"="#9ecae1","2"="#fdbb84",
  "3"="#fc8d59","4"="#ef6548","+5"="#d7301f"
)

p_alluvial <- ggplot(alluvial_data,
                     aes(axis1=Passage,axis2=Cluster_ordered,y=n)) +
  geom_alluvium(aes(fill=Passage),width=1/12) +
  geom_stratum(width=1/12,fill="grey80",color="black") +
  geom_text(stat="stratum",aes(label=after_stat(stratum)),size=3) +
  scale_x_discrete(limits=c("Passage","Cluster"),expand=c(.05,.05)) +
  scale_fill_manual(values=custom_colors,drop=FALSE) +
  theme_minimal() +
  labs(title="Passages flowing into clusters (UMAP + k-means)",
       y="Number of samples") +
  theme(plot.title=element_text(hjust=.5,face="bold"),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

ggsave(file.path(fig_dir, "Alluvial_Clusters_Hedgehog.tiff"),
       p_alluvial, width = 6, height = 5, dpi = 300)

message("Figure 4C (Hedgehog PCA/UMAP) completed.")
