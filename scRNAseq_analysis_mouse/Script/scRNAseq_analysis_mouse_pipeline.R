

# Load required libraries
library(dplyr)  # Data manipulation
library(Seurat)  # Single-cell RNA-seq analysis
library(patchwork)  # Plotting
library(readr)  # Reading data
library(SingleR)  # Cell type annotation
library(DoubletFinder)  # Doublet detection
library(Matrix)  # Sparse matrix operations
library(fields)  # Spatial data analysis
library(SeuratDisk)  # Saving Seurat objects
library(SeuratData)  # Example datasets
library(SingleCellExperiment)  # Working with single-cell data
library(celldex)  # Cell type reference databases

#####---- Parameter for loading in file ----####

# Define the path to the count file
# Specify the input folder path
# This dataset is GSE185420

# # Matrix 
# mtx <-  "J://Biostat_Interns//Nobles_Gabrielle//Github_folders//scRNAseq_analysis_mouse//Input//matrix.mtx.gz"
# # Barcodes 
# cells <- "J://Biostat_Interns//Nobles_Gabrielle//Github_folders//scRNAseq_analysis_mouse//Input//barcodes.tsv.gz"
# # Genes/Features 
# features <- "J://Biostat_Interns//Nobles_Gabrielle//Github_folders//scRNAseq_analysis_mouse//Input//features.tsv.gz"
# 

mtx <-  "/home/80023619/Moffitt_Internship_2023_02_13/scRNAseq_analysis_mouse/Input/matrix.mtx.gz"
# Barcodes 
cells <- "/home/80023619/Moffitt_Internship_2023_02_13/scRNAseq_analysis_mouse/Input/barcodes.tsv.gz"
# Genes/Features 
features <- "/home/80023619/Moffitt_Internship_2023_02_13/scRNAseq_analysis_mouse/Input/features.tsv.gz"

# ReadMtx loads sparse data 
counts <- ReadMtx(mtx, cells, features)

counts <- as.data.frame(counts)

# Define a function for quality control (QC) using Seurat
qc.seurat <- function(seurat, species, nFeature) {
  mt.pattern <- case_when(
    species == "Human" ~ "^MT-",
    species == "Mouse" ~ "^mt-",
    TRUE ~ "^MT-"
  )
  ribo.pattern <- case_when(
    species == "Human" ~ "^RP[LS]",
    species == "Mouse" ~ "^Rp[ls]",
    TRUE ~ "^RP[LS]"
  )
  
  # Calculate percentage of mitochondrial and ribosomal genes
  seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = mt.pattern, assay = "RNA")
  seurat[["percent.rp"]] <- PercentageFeatureSet(seurat, pattern = ribo.pattern, assay = "RNA")
  
  # Filter cells based on QC criteria
  seurat[, seurat[["percent.mt"]] <= 20 & seurat[["nFeature_RNA"]] >= nFeature]
}

# Create a Seurat object and perform QC
seurat_obj <- CreateSeuratObject(counts = counts,project = "GSE185420")

seurat_obj <- qc.seurat(seurat_obj, "Human", 500)

# Normalize data using LogNormalize method
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)

# Calculate cell cycle scores
seurat_obj <- CellCycleScoring(object = seurat_obj, g2m.features = cc.genes$g2m.genes, s.features = cc.genes$s.genes)

# Find highly variable features
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)

# Scale data by regressing out unwanted sources of variation
seurat_obj <- ScaleData(seurat_obj, vars.to.regress = c("nFeature_RNA", "percent.mt"), verbose = FALSE)
scaled_count <- seurat_obj@assays$RNA@scale.data

# Perform dimensionality reduction and clustering
seurat_obj <- RunPCA(seurat_obj, npcs = 30, verbose = FALSE, seed.use = 42)
seurat_obj <- RunTSNE(seurat_obj, reduction = "pca", dims = 1:30, seed.use = 1)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:10, verbose = FALSE, seed.use = 42)
seurat_obj <- FindNeighbors(seurat_obj, reduction = "pca", dims = 1:30)
seurat_obj <- FindClusters(seurat_obj, resolution = c(0.10, 0.15, 0.25,0.75))

# Find doublets using DoubletFinder
suppressMessages(require(DoubletFinder))
nExp <- round(ncol(seurat_obj) * 0.04)  # expect 4% doublets
seurat_obj <- doubletFinder_v3(seurat_obj, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10)

# Load mouse reference databases from celldex
Immgen.ref <- celldex::ImmGenData()
mouseRNAseq.ref <- celldex::MouseRNAseqData()
# Convert Seurat object to SingleCellExperiment
sce <- as.SingleCellExperiment(DietSeurat(seurat_obj))
sce

# Auto-annotate cell types using reference databases from celldex
ImmGen.main <- SingleR(test = sce, assay.type.test = 1, ref = Immgen.ref, labels = Immgen.ref$label.main)
ImmGen.fine <- SingleR(test = sce, assay.type.test = 1, ref = Immgen.ref, labels = Immgen.ref$label.fine)
mouseRNAseq.main <- SingleR(test = sce, assay.type.test = 1, ref = mouseRNAseq.ref, labels = mouseRNAseq.ref$label.main)
mouseRNAseq.fine <- SingleR(test = sce, assay.type.test = 1, ref = mouseRNAseq.ref, labels = mouseRNAseq.ref$label.fine)

# Add the celldex annotations to the metadata of the Seurat object
seurat_obj@meta.data$ImmGen.main <- ImmGen.main$pruned.labels
seurat_obj@meta.data$ImmGen.fine <- ImmGen.fine$pruned.labels
seurat_obj@meta.data$mouseRNAseq.main <- mouseRNAseq.main$pruned.labels
seurat_obj@meta.data$mouseRNAseq.fine <- mouseRNAseq.fine$pruned.labels


# Feature Plots check cell automated annaotations 

Microglia.plot <- FeaturePlot(seurat_obj,
                              reduction = "umap",
                              features = c("Itgam","P2ry12", "Tmem119", "Cx3cr1"),
                              label = TRUE,
                              order = TRUE,
                              label.size = 2

)






Astrocyte.plot <- FeaturePlot(seurat_obj,
                              reduction = "tsne",
                              features = c("S100b","Slc1a3", "Aldh1l1", "Gfap"),
                              label = T,
                              order = TRUE,
                              label.size = 3
)



Fibroblast.plot <- FeaturePlot(seurat_obj,
                               reduction = "umap",
                               features = c("Fn1","Tnc"),
                               label = T,
                               order = TRUE,
                               label.size = 3
)

Fibroblast.plot


Oligo.precur.plot <- FeaturePlot(seurat_obj,
                                    reduction = "umap",
                                    features = c("Lhfpl3","Megf11", "Pcdh15", "Cspg4"),
                                    label = T,
                                    order = TRUE,
                                    label.size = 3
)



oligo.plot <-  FeaturePlot(seurat_obj,
                                     reduction = "umap",
                                     features = c("Mog","Mbp","Mag","Cldn11"),
                                     label = T,
                                     order = TRUE,
                                     label.size = 3
)



oligo.plot2 <- FeaturePlot(seurat_obj,
                                     reduction = "tsne",
                                     features = c("Olig1","Olig2","Olig3"),
                                     label = T,
                                     order = T,
                                     label.size = 3)






# Neural Progintor Cells
npc.plot <-  FeaturePlot(seurat_obj,
                         reduction = "umap",
                         features = c("Cspg4","Rnf5","Sox2","Sox1"),
                         label = TRUE,
                         order = TRUE,
                         label.size = 3
)

# Neural Stem Cells
nsc.plot <-  FeaturePlot(seurat_obj,
                         reduction = "umap",
                         features = c("Sox9","Prom1","Nes"),
                         label = TRUE,
                         order = TRUE,
                         label.size = 3
)




macrophage.plot <-  FeaturePlot(seurat_obj,
                                reduction = "umap",
                                features = c("Cd68","Cd163","Cd14"),
                                label = TRUE,
                                order = TRUE,
                                label.size = 3
)


Tcell.plot  <-  FeaturePlot(seurat_obj,
                                reduction = "umap",
                                features = c("Cd2","Cd3d","Cd45"),
                                label = TRUE,
                                order = TRUE,
                                label.size = 3
)


Bcell.plot <-  FeaturePlot(seurat_obj,
                                reduction = "umap",
                                features = c("Cd19","Cd27","Cd24, Igd"),
                                label = TRUE,
                                order = TRUE,
                                label.size = 3
)

Endothelial.plot <-  FeaturePlot(seurat_obj,
                           reduction = "umap",
                           features = c("Pecam1","Vwf","A2m"),
                           label = TRUE,
                           order = TRUE,
                           label.size = 3
)





 Plot.lists <- list(Microglia.plot,Astrocyte.plot,Fibroblast.plot,Oligo.precur.plot,
                    oligo.plot,oligo.plot2,npc.plot,nsc.plot,macrophage.plot,
                    Tcell.plot, Bcell.plot, Endothelial.plot)


# Manual annotations
seurat_obj@meta.data$seurat_clusters_gabby_annotation <- "NA"

seurat_obj@meta.data[which(seurat_obj@meta.data$seurat_clusters %in% c(0:4, 7:8, 11:12, 15)),
                     "seurat_clusters_gabby_annotation"] <- "Tumor"
seurat_obj@meta.data[which(seurat_obj@meta.data$seurat_clusters %in% c(5,6,9)),
                     "seurat_clusters_gabby_annotation"] <- "Myeloid"
seurat_obj@meta.data[which(seurat_obj@meta.data$seurat_clusters %in% c(13)),
                     "seurat_clusters_gabby_annotation"] <- "Tcells"
seurat_obj@meta.data[which(seurat_obj@meta.data$seurat_clusters %in% c(14)),
                     "seurat_clusters_gabby_annotation"] <- "Bcells"
seurat_obj@meta.data[which(seurat_obj@meta.data$seurat_clusters %in% c(10)),
                     "seurat_clusters_gabby_annotation"] <- "Stroma"

# Add UMAP coordinates to the metadata
UMAP <- as.data.frame(Embeddings(object = seurat_obj[["umap"]]))
seurat_obj@meta.data$UMAP_1 <- UMAP$UMAP_1
seurat_obj@meta.data$UMAP_2 <- UMAP$UMAP_2

# Add tSNE coordinates to the metadata
tsne <- as.data.frame(Embeddings(object = seurat_obj[["tsne"]]))
seurat_obj@meta.data$tSNE_1 <- tsne$tSNE_1
seurat_obj@meta.data$tSNE_2 <- tsne$tSNE_1


# Specify the output folder path
output_folder <- "/home/80023619/Moffitt_Internship_2023_02_13/scRNAseq_analysis_mouse/Output/"

# Write metadata to a file
write.table(seurat_obj@meta.data, file = file.path(output_folder, "GSE185420_metadata"),
            sep = "\t")

# Write scaled count matrix to a file
# write.table(scaled_count, file = file.path(output_folder, "GSE185420_scaled_count"),
#             sep = "\t")

# Write seurat H5 file 
SaveH5Seurat(seurat_obj, filename = file.path(output_folder, "GSE185420_h5friendly"),
             overwrite = TRUE, verbose = TRUE)


# Subset and write metadata for myeloid cells
Idents(seurat_obj) <- "seurat_clusters_gabby_annotation"

myeloid_subset <- subset(x = seurat_obj, idents = "Myeloid")
write.table(myeloid_subset@meta.data, file = file.path(output_folder, "GSE185420_myeloid_metadata"),
            sep = "\t")

# Save myeloid subset as H5Seurat object
SaveH5Seurat(myeloid_subset, filename = file.path(output_folder, "GSE185420_myeloid_h5friendly"),
             overwrite = TRUE, verbose = TRUE)

# Subset and write metadata for tumor cells
tumor_subset <- subset(x = seurat_obj, idents = "Tumor")

write.table(tumor_subset@meta.data, file = file.path(output_folder, "GSE185420_tumor_metadata"),
            sep = "\t")

# Save tumor subset as H5Seurat object
SaveH5Seurat(tumor_subset, filename = file.path(output_folder, "GSE185420_tumor_h5friendly"),
             overwrite = TRUE, verbose = TRUE)

# Subset and write metadata for Tcells cells
tcell_subset <- subset(x = seurat_obj, idents = "Tcells")

write.table(tcell_subset@meta.data, file = file.path(output_folder, "GSE185420_tcell_metadata"),
            sep = "\t")

# Save tumor subset as H5Seurat object
SaveH5Seurat(tcell_subset, filename = file.path(output_folder, "GSE185420_tcell_h5friendly"),
             overwrite = TRUE, verbose = TRUE)



# Subset and write metadata for Tcells cells
bcell_subset <- subset(x = seurat_obj, idents = "Bcells")

write.table(bcell_subset@meta.data, file = file.path(output_folder, "GSE185420_bcell_metadata"),
            sep = "\t")

# Save tumor subset as H5Seurat object
SaveH5Seurat(bcell_subset, filename = file.path(output_folder, "GSE185420_bcell_h5friendly"),
             overwrite = TRUE, verbose = TRUE)




# Subset and write metadata for Tcells cells
stroma_subset <- subset(x = seurat_obj, idents = "Stroma")

write.table(stroma_subset@meta.data, file = file.path(output_folder, "GSE185420_stroma_metadata"),
            sep = "\t")

# Save tumor subset as H5Seurat object
SaveH5Seurat(stroma_subset, filename = file.path(output_folder, "GSE185420_stroma_h5friendly"),
             overwrite = TRUE, verbose = TRUE)


# 
# # Loop through each plot and generate a PDF file
# for (plot in Plot.lists) {
#   # Generate UMAP plot
#   umap_plot <- DimPlot(seurat_obj, reduction = "umap", label = TRUE, repel = FALSE,
#                        label.size = 3, seed = 42)
#   
#   # Set the PDF file name based on the plot name
#   plot_name <- deparse(substitute(plot))
#   pdf_file <- file.path(output_folder, paste0(plot_name, ".pdf"))
#   
#   # Create the PDF file and save the plot
#   pdf(file = pdf_file, width = 10, height = 8.5)
#   print(umap_plot)
#   dev.off()
# }
