# Load Libraries
library(dplyr)
library(Seurat)
library(rhdf5)
library(SeuratDisk)

# Set paths
h5_file <- "/home/80023619/Moffitt_Internship_2023_02_13/H5_ISCVA_Pipeline/Input/GSE84465_h5friendly.h5seurat"
meta_file <- "/home/80023619/Moffitt_Internship_2023_02_13/H5_ISCVA_Pipeline/Input/GSE84465_metafile_with_annotation.txt"
customized_h5_file <- "/home/80023619/Moffitt_Internship_2023_02_13/H5_ISCVA_Pipeline/Input/customized_GSE84465_ISCVAM_H5.R"
output_h5_file <- "/home/80023619/Moffitt_Internship_2023_02_13/H5_ISCVA_Pipeline/Output/GSE84465_h5_2.h5"

# h5_file <- "J:\\Biostat_Interns\\Nobles_Gabrielle\\Post_processing_write_h5\\Input\\GSE84465_h5friendly.h5seurat"
# meta_file <- "J:\\Biostat_Interns\\Nobles_Gabrielle\\Post_processing_write_h5\\Input\\GSE84465_metafile_with_annotation.txt"
# customized_h5_file <- "J:\\Biostat_Interns\\Nobles_Gabrielle\\Post_processing_write_h5\\Input\\customized_GSE84465_ISCVAM_H5.R"
# output_h5_file <- 


# Load H5 Seurat Compliant File
seurat_obj <- LoadH5Seurat(h5_file,
                             assays = c("integrated", "RNA"),
                             reductions = c("pca", "tsne", "umap", "rpca"),
                             graphs = NULL,
                             images = NULL,
                             meta.data = FALSE,
                             commands = FALSE,
                             verbose = FALSE)

# Read Metadata File
metadata <- read.table(meta_file, header = TRUE, sep = "\t")
# colnames(metadata)[grep("seurat_clusters", colnames(metadata))] <- "Clusters"
colnames(metadata)[13] = "Clusters"


# Filter and Append Selected Columns to seurat_obj@meta.data
selected_columns <- colnames(metadata)[grepl("^pANN|^orig.ident|^S.Score|^G2M.Score|^percent.rp|^percent.mt|^tSNE|^RNA|^nCount_RNA|^nFeature_RNA", 
                                             colnames(metadata))]

seurat_obj@meta.data[selected_columns] <- metadata[selected_columns]

# Remove Selected Columns from Metadata
metadata[selected_columns] <- NULL
# Remove columns by column names from the metadata data frame
# metadata <- metadata[, -which(colnames(metadata) %in% c("orig.ident", "nCount_RNA", "nFeature_RNA"))]

# Create Seurat.list
seurat.list <- list(all = list(seurat = seurat_obj, covs = metadata))

# Load Functions from Customized ISCVAM H5 File
source(customized_h5_file)

# Convert Seurat.list to H5 File
write.h5(seurat.list = seurat.list, fn = output_h5_file)

# Check Dimensions of the H5 File
h5 <- H5Fopen(output_h5_file)
h5ls(h5)

# Check Covs in ISCVA Format
colnames(h5$'/artifacts/all/covs')

# Check DiscreteCovs and ContinuousCovs
h5$'/artifacts/all/discreteCovs'
h5$'/artifacts/all/continuousCovs'

# Close the H5 File
h5closeAll()
