
library(tidyverse)
library(clusterCrit)
source("clone/runfactorization_custom.R")

# Load RNA-seq data
exp <- readRDS("data/single-cell/CellLines_RNAseqCounts.RDS", refhook = NULL) #ENS for genes and counts

# Apply log2 on RNA-seq data
exp <- log2(exp + 1)

# Load ATAC-seq data
atac_counts <- readRDS("data/single-cell/CellLines_ATACseqCounts.RDS", refhook = NULL) # peaks counts

# Load metadata
metadata <- readRDS("data/single-cell/CellLines_metadata.RDS", refhook = NULL)

# Rename columns from metadata
colnames(atac_counts) <- metadata[,1]

# Export RNA-seq data as tab-separated table
write.table(exp, "data/single-cell/CellLines_RNAseqCounts.txt", sep = "\t", col.names = TRUE, row.names = TRUE)

# Add a name ("probe") to the first column
system("sed -i '1s/^/probe\t/' data/single-cell/CellLines_RNAseqCounts.txt")

# Export ATAC-seq data as tab-separated table
write.table(atac_counts, "data/single-cell/CellLines_ATACseqCounts.txt", sep = "\t", col.names = TRUE, row.names = TRUE)

# Add a name ("probe") to the first column
system("sed -i '1s/^/probe\t/' data/single-cell/CellLines_ATACseqCounts.txt")

# Parameters for the plots
dot_size <- 1.5
dot_alpha <- 1.0
xlabel <- "Factor 1"
ylabel <- "Factor 2"

# Load annotations from the metadata
sample_annot <- metadata[, c("sample.rna", "celltype")]

# Folder for results
results_folder <- "clone/results_single_cell/"

# Create output folder
dir.create(results_folder, showWarnings = FALSE)

# Run factorization methods
out <- runfactorization(folder = "data/single-cell/",
                        file.names = c("CellLines_RNAseqCounts.txt", "CellLines_ATACseqCounts.txt"),
                        num.factors = 2, 
                        minPts = 10,
                        sep = "\t", 
                        filtering = "stringent")

# save(out, file = "clone/results_single_cell/factorization.RData")
# load("clone/results_single_cell/factorization.RData")

c_index <- numeric(0)

# For each factorization method
for(i in 1:length(out$factorizations)){
  
  # Get factorization result
  factors <- out$factorizations[[i]][[1]]
  
  # Delete NAs
  factors <- factors[!is.na(factors[,1]) & !is.na(factors[,2]) ,]
  sample_annotI <- sample_annot[!is.na(sample_annot[,1]) & !is.na(sample_annot[,2]), ]
  
  # Data to be plotted
  df <- data.frame(x =  factors[,1], y = factors[,2], color_by = sample_annotI[,2])
  # Plot results
  p <- ggplot(df, aes_string(x = "x", y = "y")) + 
    geom_point(aes_string(color = "color_by"), size=dot_size, alpha=dot_alpha) + 
    xlab(xlabel) + ylab(ylabel) +
    # scale_shape_manual(values=c(19,1,2:18)[seq_along(unique(shape_by))]) +
    theme(plot.margin = margin(20, 20, 10, 10), 
          axis.text = element_text(size = rel(1), color = "black"), 
          axis.title = element_text(size = 16), 
          axis.title.y = element_text(size = rel(1.1), margin = margin(0, 10, 0, 0)), 
          axis.title.x = element_text(size = rel(1.1), margin = margin(10, 0, 0, 0)), 
          axis.line = element_line(color = "black", size = 0.5), 
          axis.ticks = element_line(color = "black", size = 0.5),
          panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(),
          legend.key = element_rect(fill = "white"),
          legend.text = element_text(size = 16),
          legend.title = element_text(size =16)
    )
  p + scale_color_manual(values=c("#0072B2", "#D55E00", "#CC79A7"))
  # Export plot as JPEG image
  ggsave(paste0(results_folder, "plot_", out$method[i], ".jpg"))
  
  # Encode cell type annotations by numeric codes
  ann <- factor(sample_annotI[, 2], levels = c("HCT", "Hela", "K562"))
  ann <- as.integer(ann)
  # Compare factors and annotations
  c_index <- c(c_index, intCriteria(as.matrix(factors), as.integer(ann), crit = c("C_index"))$c_index)
  
}

# Build output table
report_cindex <- data.frame(method = out$method, cindex = c_index)

# Export results as one tab-separated table
write.table(report_cindex, file = paste0(results_folder, "singlecell_cindex.txt"), 
            sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

## Final plots
rgcca <- as.data.frame(out$factorizations[[1]][[1]])
rgcca$method <- "RGCCA"
colnames(rgcca)[1:2] <- c("factor1", "factor2")

mcia <- as.data.frame(out$factorizations[[2]][[1]])
mcia$method <- "MCIA"
colnames(mcia)[1:2] <- c("factor1", "factor2")

mofa <- as.data.frame(out$factorizations[[3]][[1]])
mofa$method <- "MOFA"
colnames(mofa)[1:2] <- c("factor1", "factor2")

intnmf <- as.data.frame(out$factorizations[[4]][[1]])
intnmf$method <- "intNMF"
colnames(intnmf)[1:2] <- c("factor1", "factor2")

tica <- as.data.frame(out$factorizations[[6]][[1]])
tica$method <- "tICA"
colnames(tica)[1:2] <- c("factor1", "factor2")

ubmi <- as.data.frame(out$factorizations[[7]][[1]])
ubmi$method <- "UBMI"
colnames(ubmi)[1:2] <- c("factor1", "factor2")

## Merge methods and plot
all_sc <- rbind(rgcca, mcia, mofa, intnmf, tica, ubmi) %>% 
  as.data.frame() %>% 
  rownames_to_column("id") %>% 
  mutate(type = case_when(grepl("HCT", id) ~ "HCT",
                          grepl("Hela", id) ~ "Hela",
                          grepl("K562", id) ~ "K562")) %>% 
  left_join(report_cindex, by = c("method" = "V1")) #%>% 
  # mutate(method = paste0(method, ": ", round(V2, 3)))

ggplot(all_sc, aes(x = factor1, y = factor2)) + 
  geom_point(aes(fill = type), color = "black", pch = 21, size = 3, alpha = 0.7) +
  labs(x = "Factor 1",
       y = "Factor 2") +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = "top",
        text = element_text(size = 15)) +
  facet_wrap(~ method, nrow = 2, scales = "free") +
  scale_fill_viridis_d()
  # POMA::scale_fill_poma_d()

ggsave(paste0(results_folder, "plot_overall.jpg"), width = 10)

# UBMI clusters
ubmi_clusters <- data.frame(sample = rownames(ubmi), ubmi, clust = paste0("cluster", out$UBMI.clusters))

ggplot(ubmi_clusters, aes(x = factor1, y = factor2)) + 
  geom_point(aes(fill = clust), color = "black", pch = 21, size = 3, alpha = 0.7) +
  labs(x = "Factor 1",
       y = "Factor 2") +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = "top",
        text = element_text(size = 15)) +
  scale_fill_viridis_d()

ggsave(paste0(results_folder, "plot_UBMI_clusters.jpg"), width = 5, height = 5)

# Cluster expression analysis
library(POMA)

ubmi_clusters_subset <- ubmi_clusters[ubmi_clusters$clust %in% c("cluster1", "cluster3"),]
exp_subset <- exp[, colnames(exp) %in% rownames(ubmi_clusters_subset)]

poma_obj <- PomaSummarizedExperiment(target = ubmi_clusters_subset[, c("sample", "clust")],
                                     features = t(exp_subset))

poma_obj %>% 
  PomaVolcano(method = "limma", xlim = 5) +
  theme(text = element_text(size = 15))

ggsave(paste0(results_folder, "volcanoplot_clusters_1_3_UBMI.jpg"), width = 8)

# Enrichment analysis
library(clusterProfiler)
library(org.Hs.eg.db)

limma_res <- poma_obj %>% 
  PomaLimma(contrast = "cluster1-cluster3")

limma_res %>% 
  filter(P.Value < 0.05) %>%
  count()

# saveRDS(limma_res, file = paste0(results_folder, "limma_res_cluster1_cluster3.Rds"))

## KEGG
feature_list <- limma_res %>%
  arrange(desc(logFC))

gene_list <- feature_list$logFC
names(gene_list) <- feature_list$feature

ids <- bitr(names(gene_list),
            fromType = "ENSEMBL",
            toType = "ENTREZID",
            OrgDb = org.Hs.eg.db)

feature_list_entrez <- feature_list %>%
  right_join(ids, by = c("feature" = "ENSEMBL"))

# Create a vector of the gene universe
kegg_gene_list <- feature_list_entrez$logFC

# Name vector with ENTREZ ids
names(kegg_gene_list) <- feature_list_entrez$ENTREZID

kegg_res <- gseKEGG(geneList      = kegg_gene_list,
                    organism      = "hsa",
                    nPerm         = 10000,
                    minGSSize     = 10,
                    maxGSSize     = 500,
                    pvalueCutoff  = 0.05,
                    pAdjustMethod = "fdr",
                    keyType       = "ncbi-geneid")

dotplot(kegg_res,
        showCategory = 20,
        title = "",
        split = ".sign") +
  facet_grid(. ~ .sign) +
  theme(legend.position = "top")

ggsave(paste0(results_folder, "GSEA_KEGG.jpg"), width = 10)

## Gene Ontology
gene_list <- limma_res$logFC
names(gene_list) <- limma_res$feature

# omit any NA values 
gene_list <- na.omit(gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list <- sort(gene_list, decreasing = TRUE)

gse <- gseGO(geneList = gene_list, 
             ont = "BP", 
             keyType = "ENSEMBL", 
             nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = org.Hs.eg.db, 
             pAdjustMethod = "fdr")

dotplot(gse, showCategory = 10, split = ".sign") + 
  facet_grid(. ~ .sign) +
  theme(legend.position = "top",
        text = element_text(size = 15))

ggsave(paste0(results_folder, "GSEA_GO.jpg"), width = 10)

