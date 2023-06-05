# Load Library 
library(dplyr)
library(Seurat)
library(patchwork)
library(readr)
library(SingleR)
library(DoubletFinder)
library(Matrix)
library(fields)
library(DoubletFinder)
library(SeuratDisk)
library(SeuratData)




####---- Read in GSE141982_GSM4217362_SF12199, convert it to seurat object and read in the metadata----####

# File path

mtx_1 <-  "/home/80023619/Moffitt_Internship_2023_02_13/GSE141982/GSM4217362_SF12199_matrix.mtx.gz"
cells_1 <- "/home/80023619/Moffitt_Internship_2023_02_13/GSE141982/GSM4217362_SF12199_barcodes.tsv.gz"
features_1 <- "/home/80023619/Moffitt_Internship_2023_02_13/GSE141982/GSM4217362_SF12199_features.tsv.gz"

# ReadMtx loads sparse data 
counts_1 <- ReadMtx(mtx_1, cells_1, features_1)


# Create the Seurat object 
seurat_obj_1 <- CreateSeuratObject(counts = counts_1, project = "GSE141982_GSM4217362_SF12199", 
                                   min.cells = 3, min.features = 200)
# Preform data normalization 
seurat_obj_1 <- NormalizeData(seurat_obj_1, normalization.method = "LogNormalize",
                              scale.factor = 10000)
# Preform cell cycle scoring 
seurat_obj_1 <- CellCycleScoring(object = seurat_obj_1, g2m.features = cc.genes$g2m.genes,
                                 s.features = cc.genes$s.genes)
# Find highly  variable features 
seurat_obj_1 <- FindVariableFeatures(seurat_obj_1, selection.method = "vst", nfeatures = 2000)

# Scale the data 
all.genes <- rownames(seurat_obj_1)

seurat_obj_1 <- ScaleData(seurat_obj_1, verbose = FALSE, features = all.genes)

# Dimensional Reduction and Clustering 
seurat_obj_1 <- RunPCA(seurat_obj_1, npcs = 30, verbose = FALSE, seed.use = 42 )


# Preform doublet detection 
suppressMessages(require(DoubletFinder))

nExp <- round(ncol(seurat_obj_1) * 0.04)  # expect 4% doublets
seurat_obj_1 <- doubletFinder_v3(seurat_obj_1, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10)

# Read in the individualized metadata 
Meta_file_1 <- "/home/80023619/Moffitt_Internship_2023_02_13/GSE141982/GSE141982_GSM4217362_SF12199_metafile_with_annotations.txt"
metadata_1 <- read.table(Meta_file_1, header = T, sep = "\t")

#Add metafile columns to the seurat object
seurat_obj_1@meta.data$seurat_clusters_gabby_annotation <- metadata_1$seurat_clusters_gabby_annotation
seurat_obj_1@meta.data$Series_title <- "GSM4217362"
seurat_obj_1@meta.data$Patient_ID <- "SF12199"
seurat_obj_1@meta.data$Sample_type <- "SRA"
seurat_obj_1@meta.data$Organ <- "Brain"
seurat_obj_1@meta.data$Organism <- "homo_sapiens"
seurat_obj_1@meta.data$Age <- "74"
seurat_obj_1@meta.data$Gender <- "Male"
seurat_obj_1@meta.data$idh_status <- "Mutant"
seurat_obj_1@meta.data$Extracted_molecule <- "total_RNA"
seurat_obj_1@meta.data$Extraction_protocol <- "10x_genomics"
seurat_obj_1@meta.data$Library_strategy <- "RNAseq"
seurat_obj_1@meta.data$Library_source <- "transcriptomic"
seurat_obj_1@meta.data$Library_selection <- "cDNA"
seurat_obj_1@meta.data$Instrunment_model <- "Illumina_NovaSeq_600"
seurat_obj_1@meta.data$Platform_ID <- "GPL24676"


####---- Read in GSE141982_GSM4217363_Mixing, convert it to seurat object and read in the metadata----####

 # File path

mtx_2 <- "/home/80023619/Moffitt_Internship_2023_02_13/GSE141982/GSM4217363_Mixing_matrix.mtx.gz"
cells_2 <- "/home/80023619/Moffitt_Internship_2023_02_13/GSE141982/GSM4217363_Mixing_barcodes.tsv.gz"
features_2 <- "/home/80023619/Moffitt_Internship_2023_02_13/GSE141982/GSM4217363_Mixing_features.tsv.gz"

# ReadMtx loads in sparse data 
counts_2 <- ReadMtx(mtx_2, cells_2, features_2)


# Create the Seurat object 
seurat_obj_2 <- CreateSeuratObject(counts = counts_2, project = "GSE141982_GSM4217363_Mixing", 
                                   min.cells = 3, min.features = 200)

# Preform data normalization 
seurat_obj_2 <- NormalizeData(seurat_obj_2, normalization.method = "LogNormalize",
                              scale.factor = 10000)

# Preform cell cycle scoring 
seurat_obj_2 <- CellCycleScoring(object = seurat_obj_2, g2m.features = cc.genes$g2m.genes,
                                 s.features = cc.genes$s.genes)
#Find higly variable features 
seurat_obj_2 <- FindVariableFeatures(seurat_obj_2, selection.method = "vst", nfeatures = 2000)

# Scale the data 
all.genes <- rownames(seurat_obj_2)

seurat_obj_2 <- ScaleData(seurat_obj_2, verbose = FALSE, features = all.genes)

# Dimensional Reduction and Clustering 
seurat_obj_2 <- RunPCA(seurat_obj_2, npcs = 30, verbose = FALSE, seed.use = 42 )

# Run doublet detection 
suppressMessages(require(DoubletFinder))

nExp <- round(ncol(seurat_obj_2) * 0.04)  # expect 4% doublets
seurat_obj_2 <- doubletFinder_v3(seurat_obj_2, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10)


# Load in the metadata 
Meta_file_2 <- "/home/80023619/Moffitt_Internship_2023_02_13/GSE141982/GSE141982_GSM4217363_Mixing_metafile_with_cell_annotation_2023_01_04.txt"
metadata_2 <- read.table(Meta_file_2, header = T, sep = "\t")

#Add metafile columns to the seurat object
seurat_obj_2@meta.data$seurat_clusters_gabby_annotation <- metadata_2$seurat_clusters_gabby_annotation
seurat_obj_2@meta.data$Series_title <- "GSM4217363"
seurat_obj_2@meta.data$Patient_ID <- "Mixing"
seurat_obj_2@meta.data$Sample_type <- "SRA"
seurat_obj_2@meta.data$Organ <- "Brain"
seurat_obj_2@meta.data$Organism <- "homo_sapiens"
seurat_obj_2@meta.data$Age <- "34/44"
seurat_obj_2@meta.data$Gender <- "Male/Male"
seurat_obj_2@meta.data$idh_status <- "Mutant"
seurat_obj_2@meta.data$Extracted_molecule <- "total_RNA"
seurat_obj_2@meta.data$Extraction_protocol <- "10x_genomics"
seurat_obj_2@meta.data$Library_strategy <- "RNAseq"
seurat_obj_2@meta.data$Library_source <- "transcriptomic"
seurat_obj_2@meta.data$Library_selection <- "cDNA"
seurat_obj_2@meta.data$Instrunment_model <- "Illumina_NovaSeq_600"
seurat_obj_2@meta.data$Platform_ID <- "GPL24676"


####---- Read in GSE141982_GSM4217364_SF10432, convert it to seurat object and read in the metadata ----####

# File path


mtx_3 <-  "/home/80023619/Moffitt_Internship_2023_02_13/GSE141982/GSM4217364_SF10432_matrix.mtx.gz"
cells_3 <- "/home/80023619/Moffitt_Internship_2023_02_13/GSE141982/GSM4217364_SF10432_barcodes.tsv.gz"
features_3 <- "/home/80023619/Moffitt_Internship_2023_02_13/GSE141982/GSM4217364_SF10432_features.tsv.gz"

# ReadMtx loads in sparse data 
counts_3 <- ReadMtx(mtx_3, cells_3, features_3)


# Create the Seurat object 
seurat_obj_3 <- CreateSeuratObject(counts = counts_3, project = "GSE141982_GSM4217364_SF10432", 
                                   min.cells = 3, min.features = 200)


# Preform data normalization
seurat_obj_3 <- NormalizeData(seurat_obj_3, normalization.method = "LogNormalize",
                              scale.factor = 10000)


# Preform cell cycle scoring 
seurat_obj_3 <- CellCycleScoring(object = seurat_obj_3, g2m.features = cc.genes$g2m.genes,
                                 s.features = cc.genes$s.genes)
#Find higly variable features 
seurat_obj_3 <- FindVariableFeatures(seurat_obj_3, selection.method = "vst", nfeatures = 2000)

# Scale the data 
all.genes <- rownames(seurat_obj_3)

seurat_obj_3 <- ScaleData(seurat_obj_3, verbose = FALSE, features = all.genes)

# Dimensional Reduction and Clustering 
seurat_obj_3 <- RunPCA(seurat_obj_3, npcs = 30, verbose = FALSE, seed.use = 42 )

# Run doblet detection 
suppressMessages(require(DoubletFinder))

nExp <- round(ncol(seurat_obj_3) * 0.04)  # expect 4% doublets
seurat_obj_3 <- doubletFinder_v3(seurat_obj_3, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10)


# Load in the meta data 
Meta_file_3 <- "/home/80023619/Moffitt_Internship_2023_02_13/GSE141982/GSE141982_GSM4217364_SF10432_metafile_with_cell_annotation_2023_01_04.txt"
metadata_3 <- read.table(Meta_file_3, header = T, sep = "\t")


#Add metafile columns to the seurat object
seurat_obj_3@meta.data$seurat_clusters_gabby_annotation <- metadata_3$seurat_clusters_gabby_annotation
seurat_obj_3@meta.data$Series_title <- "GSM4217364"
seurat_obj_3@meta.data$Patient_ID <- "S10432"
seurat_obj_3@meta.data$Sample_type <- "SRA"
seurat_obj_3@meta.data$Organ <- "Brain"
seurat_obj_3@meta.data$Organism <- "homo_sapiens"
seurat_obj_3@meta.data$Age <- "50"
seurat_obj_3@meta.data$Gender <- "Female"
seurat_obj_3@meta.data$idh_status <- "Wildtype"
seurat_obj_3@meta.data$Extracted_molecule <- "total_RNA"
seurat_obj_3@meta.data$Extraction_protocol <- "10x_genomics"
seurat_obj_3@meta.data$Library_strategy <- "RNAseq"
seurat_obj_3@meta.data$Library_source <- "transcriptomic"
seurat_obj_3@meta.data$Library_selection <- "cDNA"
seurat_obj_3@meta.data$Instrunment_model <- "Illumina_NovaSeq_600"
seurat_obj_3@meta.data$Platform_ID <- "GPL24676"


####----- Merge the three  seurat objects together ----####

# To prevent over scaling or normaliZation of the data we comment it out 
combined.list <- lapply(X = c(seurat_obj_1, seurat_obj_2,seurat_obj_3), FUN = function(x) {
  # x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# Select integration features across all Seurat objects based on the variable features found in the previous step
features <- SelectIntegrationFeatures(object.list = combined.list)

# Perform data integration on each Seurat object in the combined.list
combined.list <- lapply(X = combined.list, FUN = function(x) {
  #x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE, seed.use = 42 )
  
})

# Find integration anchors across all Seurat objects using the integration features and the "rpca" reduction method with dimensions 1 to 50
anchors <- FindIntegrationAnchors(object.list = combined.list, 
                                  anchor.features = features, 
                                  reduction = "rpca", 
                                  dims = 1:50)

# this command creates an 'integrated' data assay
combined.integrated <- IntegrateData(anchorset = anchors, dims = 1:50)

####--- Run QC on combined data ----####
# Function to perform quality control on a Seurat object

# Input:
# - seurat: Seurat object
# - species: Species of the data (e.g., "Human" or "Mouse")
# - nFeature: Minimum number of features required for a cell to pass QC

qc.seurat <- function(seurat, species, nFeature) {
  
  # Define patterns for mitochondrial genes
  mt.pattern <- case_when(
    species == "Human" ~ "^MT-",
    species == "Mouse" ~ "^mt-",
    TRUE ~ "^MT-"
  )
  
  # Define patterns for ribosomal protein genes
  ribo.pattern <- case_when(
    species == "Human" ~ "^RP[LS]",
    species == "Mouse" ~ "^Rp[ls]",
    TRUE ~ "^RP[LS]"
  )
  
  # Calculate the percentage of mitochondrial genes per cell
  seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = mt.pattern, assay = 'RNA')
  
  # Calculate the percentage of ribosomal protein genes per cell
  seurat[["percent.rp"]] <- PercentageFeatureSet(seurat, pattern = ribo.pattern, assay = 'RNA')
  
  # Filter cells based on QC criteria (<= 20% mitochondrial genes and >= nFeature features)
  seurat <- seurat[, seurat[["percent.mt"]] <= 20 & seurat[["nFeature_RNA"]] >= nFeature]
  
  return(seurat)
}

# Perform QC on the combined integrated Seurat object
combined.integrated <- qc.seurat(combined.integrated, "Human", 500)

# Scale the data
combined.integrated <- ScaleData(combined.integrated, vars.to.regress = c("nFeature_RNA", "percent.mt"), verbose = FALSE)

# Retrieve the scaled counts
scaled_count <- combined.integrated@assays$integrated@scale.data

####---- Dimensional Reduction and Clustering ----####
# Perform PCA
combined.integrated <- RunPCA(combined.integrated, npcs = 30, verbose = FALSE, seed.use = 42)

# Perform UMAP
combined.integrated <- RunUMAP(combined.integrated, dims = 1:10, verbose = FALSE, seed.use = 42)

# Perform t-SNE
combined.integrated <- RunTSNE(combined.integrated, reduction = "pca", dims = 1:30, seed.use = 1)

# Find nearest neighbors
combined.integrated <- FindNeighbors(combined.integrated, reduction = "pca", dims = 1:30)

# Find clusters using different resolutions
combined.integrated <- FindClusters(combined.integrated, resolution = c(0.10, 0.15, 0.25, 0.50, 0.75))


####---- Single R cell annotations ----####
# Reference database from celldex
hpca.ref <- celldex::HumanPrimaryCellAtlasData()
dice.ref <- celldex::DatabaseImmuneCellExpressionData()
blueprint.ref <- celldex::BlueprintEncodeData()
monaco.ref <- celldex::MonacoImmuneData()
northern.ref <- celldex::NovershternHematopoieticData()
#cover seurat object to single cell expirment
sce <- as.SingleCellExperiment(DietSeurat(combined.integrated))
sce



# Using the refernce database from celldex auto-annotate the cell types
hpca.main <- SingleR(test = sce,assay.type.test = 1,ref = hpca.ref,labels = hpca.ref$label.main)
hpca.fine <- SingleR(test = sce,assay.type.test = 1,ref = hpca.ref,labels = hpca.ref$label.fine)
dice.main <- SingleR(test = sce,assay.type.test = 1,ref = dice.ref,labels = dice.ref$label.main)
dice.fine <- SingleR(test = sce,assay.type.test = 1,ref = dice.ref,labels = dice.ref$label.fine)
blue.main <- SingleR(test = sce,assay.type.test = 1,ref = blueprint.ref,labels = blueprint.ref$label.main)
blue.fine <- SingleR(test = sce,assay.type.test = 1,ref = blueprint.ref,labels = blueprint.ref$label.fine)
monaco.main <- SingleR(test = sce,assay.type.test = 1,ref = monaco.ref,labels = monaco.ref$label.main)
monaco.fine <- SingleR(test = sce,assay.type.test = 1,ref = monaco.ref,labels = monaco.ref$label.fine)
northern.main <- SingleR(test = sce,assay.type.test = 1,ref = northern.ref,labels = northern.ref$label.main)
northern.fine <- SingleR(test = sce,assay.type.test = 1,ref = northern.ref,labels = northern.ref$label.fine)


# Add the celldex annotations to the metadata
combined.integrated@meta.data$hpca.main   <- hpca.main$pruned.labels
combined.integrated@meta.data$hpca.fine   <- hpca.fine$pruned.labels
combined.integrated@meta.data$dice.main   <- dice.main$pruned.labels
combined.integrated@meta.data$dice.fine   <- dice.fine$pruned.labels
combined.integrated@meta.data$monaco.main   <- monaco.main$pruned.labels
combined.integrated@meta.data$monaco.fine   <- monaco.fine$pruned.labels
combined.integrated@meta.data$northern.main   <- northern.main$pruned.labels
combined.integrated@meta.data$northern.fine <- northern.fine$pruned.labels
combined.integrated@meta.data$blue.main   <- blue.main$pruned.labels
combined.integrated@meta.data$blue.fine <- blue.fine$pruned.labels


#embed the umap and tsne corrdinates in the meta.data 

UMAP <-  as.data.frame(Embeddings(object = combined.integrated[["umap"]]))
combined.integrated@meta.data$UMAP_1 <- UMAP$UMAP_1
combined.integrated@meta.data$UMAP_2 <- UMAP$UMAP_2

tsne <-  as.data.frame(Embeddings(object = combined.integrated[["tsne"]]))
combined.integrated@meta.data$tSNE_1 <- tsne$tSNE_1
combined.integrated@meta.data$tSNE_2 <- tsne$tSNE_2




####---- Output Files ----####

# Save the integrated Seurat object in h5 format
SaveH5Seurat(combined.integrated, filename = "/home/80023619/Moffitt_Internship_2023_02_13/GSE141982/GSE141982_merged_h5friendly",
             overwrite = TRUE, verbose = TRUE)

# Set the identity of cells in the integrated Seurat object to "seurat_clusters_gabby_annotation"
Idents(combined.integrated) <- "seurat_clusters_gabby_annotation"

# Subset the integrated Seurat object to include only cells identified as "Myeloid"
Myeloid_subset <- subset(x = combined.integrated, idents = "Myeloid")

# Write the metadata of the Myeloid subset to a tab-separated file
write.table(Myeloid_subset@meta.data, file= "/home/80023619/Moffitt_Internship_2023_02_13/GSE141982/GSE141982_merged_myeloid_metafile.txt",
            sep="\t")

# Save the Myeloid subset Seurat object in h5 format
SaveH5Seurat(Myeloid_subset, filename = "/home/80023619/Moffitt_Internship_2023_02_13/GSE141982_merged_myeloid_h5_friendly",
             overwrite = TRUE, verbose = TRUE)

# Subset the integrated Seurat object to include only cells identified as "Stroma"
Stroma_subset <- subset(x = combined.integrated, idents = "Stroma")

# Write the metadata of the Stroma subset to a tab-separated file
write.table(Stroma_subset@meta.data, file="/home/80023619/Moffitt_Internship_2023_02_13/GSE141982/GSE141982_merged_stroma_metafile.txt",
            sep="\t")

# Save the Stroma subset Seurat object in h5 format
SaveH5Seurat(Stroma_subset, filename = "/home/80023619/Moffitt_Internship_2023_02_13/GSE141982/GSE141982_merged_stroma_h5friendly",
             overwrite = TRUE, verbose = TRUE)

# Subset the integrated Seurat object to include only cells identified as "Tumor"
Tumor_subset <- subset(x = combined.integrated, idents = "Tumor")

# Write the metadata of the Tumor subset to a tab-separated file
write.table(Tumor_subset@meta.data, file="/home/80023619/Moffitt_Internship_2023_02_13/GSE141982/GSE141982_merged_tumor_metafile.txt",
            sep="\t")

# Save the Tumor subset Seurat object in h5 format
SaveH5Seurat(Tumor_subset, filename = "/home/80023619/Moffitt_Internship_2023_02_13/GSE141982/GSE141982_merged_tumor_h5friendly",
             overwrite = TRUE, verbose = TRUE)

# Write the metadata of the integrated Seurat object to a tab-separated file
write.table(combined.integrated@meta.data, file = "/home/80023619/Moffitt_Internship_2023_02_13/GSE141982/GSE141982_merged_metafile_with_annotations_v2.txt",
            sep = "\t")

# Write the scaled count matrix of the integrated Seurat object to a tab-separated file
write.table(scaled_count, file = "/home/80023619/Moffitt_Internship_2023_02_13/GSE141982/GSE141982_merged_scale_count_matrix",
            sep = "\t")
