# The H5 script is taking the H5 file and loading it into a seurat object and outputing: 
## Gene expression txt files for each cluster and excel file: containing unfiltered and filtered by p_val_adj < 0.05
## Gene expression txt file for cell types and excel file: containing unfiltered and filtered by p_val_adj < 0.05
## Subsetting 1000 matrix and outputing, raw, normalized, scaled and metadata 




#####---- R package Installation ----####
# Define the list of packages to install
packages <- c("dplyr", "Seurat", "patchwork", "readr", "SingleR", 
              "DoubletFinder", "Matrix", "fields", "SeuratDisk",
              "SeuratData", "openxlsx")

# Check if the packages are already installed
installed_packages <- row.names(installed.packages())
packages_to_install <- setdiff(packages, installed_packages)

# If any packages need to be installed, install them
if (length(packages_to_install) > 0) {
  install.packages(packages_to_install)
}

# Load in the necessary packages 
lapply(packages, library, character.only = TRUE)



####---- Load in the seurat object from the seurat H5 file ----####
# Set the working directory
setwd("/home/80023619/Moffitt_Internship_2023_02_13/Test2/Input")

# List all files in the working directory
files <- list.files()

# Find the path of the h5 file 
h5_file <- files[grep(".h5seurat", files)]

# Load the file as a Seurat object using the LoadH5Seurat function
seurat_obj <- LoadH5Seurat(h5_file,
                           assays = c("integrated", "RNA"),
                           reductions = c("pca", "tsne", "umap", "rpca"),
                           graphs = NULL,
                           images = NULL,
                           meta.data = T,
                           commands = F,
                           verbose = F)

####---- Create an output directory for all files to go to----####
setwd("/home/80023619/Moffitt_Internship_2023_02_13/Test2/")

#Define the output directory name
output_dir <- "output_files"

# Create an output directory if it doesn't already exist
if (!file.exists(output_dir)) {
  dir.create(output_dir)
}

# Assign the full output directory path to a variable
full_output_dir <- paste0(getwd(),"/", output_dir, sep = "/")

####---- Find markers and Gene Expression ----#### 

# Generate the Markers for each clusters 
# Set the Idents 
Idents(seurat_obj) <- "seurat_clusters"


# get the number of cells in each cluster
cluster_sizes <- table(Idents(seurat_obj))

# only include clusters with more than 2 cells (FindMaker will only work on clusters that have at least 3 cells )
clusters <- names(cluster_sizes)[cluster_sizes > 2]

wb <- createWorkbook()

# loop through each cluster and find marker genes
for (cluster_id in clusters) {
  
  print(paste("Finding marker genes for cluster", cluster_id))
  
  # use the FindMarker() function to get the gene markers
  unfiltered_markers <- FindMarkers(object = seurat_obj, ident.1 = cluster_id, min.pct = .25 )
  
  upreg_markers <- unfiltered_markers[unfiltered_markers$p_val_adj < 0.05, ]
  
  
  # write the marker genes to a separate file 
  unfilter_output_file<- file.path(full_output_dir, paste0("cluster_", cluster_id, "_unfilt_gene_expression.txt"))
  upreg_output_file <- file.path(full_output_dir, paste0("cluster_", cluster_id, "_upreg_gene_expression.txt"))
  
  
  write.table(unfiltered_markers, file = unfilter_output_file, sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
  write.table(upreg_markers, file = upreg_output_file, sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
  
  # add a new sheet to the workbook for unfiltered markers and write the marker genes to the sheet
  unfiltered_sheet_name <- paste0("Cluster_", cluster_id, "_unfilt")
  addWorksheet(wb, unfiltered_sheet_name)
  writeData(wb, unfiltered_sheet_name, unfiltered_markers, row.names = TRUE)  # include row.names in the output
  setColWidths(wb, unfiltered_sheet_name, cols = 1:ncol(unfiltered_markers), widths = "auto") # Autofit columns
  
  # add a new sheet to the workbook for upregulated markers and write the marker genes to the sheet
  upregulated_sheet_name <- paste0("Cluster_", cluster_id, "_upreg")
  addWorksheet(wb, upregulated_sheet_name)
  writeData(wb, upregulated_sheet_name, upreg_markers, row.names = TRUE)  # include row.names in the output
  setColWidths(wb, upregulated_sheet_name, cols = 1:ncol(upreg_markers), widths = "auto") # Autofit columns
  
}

# Save the workbook to an Excel file
xlsx_output_file <- file.path(full_output_dir, "marker_genes.xlsx")
saveWorkbook(wb, xlsx_output_file)




# Generate marker list  by cell types
wb2 <-  createWorkbook()


Idents(seurat_obj) <- "seurat_clusters_gabby_annotation"

cell_types <- unique(seurat_obj@meta.data$seurat_clusters_gabby_annotation)


for(cells in cell_types){
  print(paste("Finding marker genes for", cells))
  
  # use the FindMarker() function to get the top 10 marker genes
  unfiltered_markers <- FindMarkers(object = seurat_obj, ident.1 = cells, min.pct = .25 )
  
  upreg_markers <- unfiltered_markers[unfiltered_markers$p_val_adj < 0.05, ]
  
  
  # write the marker genes to a separate file 
  unfilter_output_file<- file.path(full_output_dir, paste0(cells, "_unfilt_gene_expression.txt"))
  upreg_output_file <- file.path(full_output_dir, paste0(cells, "_upreg_gene_expression.txt"))
  
  
  write.table(unfiltered_markers, file = unfilter_output_file, sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
  write.table(upreg_markers, file = upreg_output_file, sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
  
  # add a new sheet to the workbook for unfiltered markers and write the marker genes to the sheet
  unfiltered_sheet_name <- paste0(cells, "_unfilt")
  addWorksheet(wb2, unfiltered_sheet_name)
  writeData(wb2, unfiltered_sheet_name, unfiltered_markers, row.names = TRUE)  # include row.names in the output
  setColWidths(wb2, unfiltered_sheet_name, cols = 1:ncol(unfiltered_markers), widths = "auto") # Autofit columns
  
  # add a new sheet to the workbook for upregulated markers and write the marker genes to the sheet
  upregulated_sheet_name <- paste0(cells, "_upreg")
  addWorksheet(wb2, upregulated_sheet_name)
  writeData(wb2, upregulated_sheet_name, upreg_markers, row.names = TRUE)  # include row.names in the output
  setColWidths(wb2, upregulated_sheet_name, cols = 1:ncol(upreg_markers), widths = "auto") # Autofit columns
  
  
}

# Save the workbook to an Excel file
xlsx_output_cell_file <- file.path(full_output_dir, "cell_marker_genes.xlsx")
saveWorkbook(wb2, xlsx_output_cell_file)



####---- Random 1000 Cell matrix for Umap app ----####

set.seed(28)  # set the seed for reproducibility
output_dir <- full_output_dir
num_cells <- 1000

folder_path <- paste0(output_dir, "/", num_cells, "_matrix")
dir.create(folder_path, showWarnings = FALSE)

# Iterate over columns of the subset Seurat object
for (i in seq_len(ncol(seurat_obj))) {
  # Randomly sample cells
  random_cells <- sample(1:nrow(seurat_obj), size = num_cells, replace = FALSE)
  subset_seurat_obj <- seurat_obj[random_cells, ]
  
  # Extract data from the subset Seurat object
  if ("integrated" %in% colnames(subset_seurat_obj@assays[[i]])) {
    # If 'integrated' assay is available, retrieve the scaled data
    int_scaled_matrix <- as.data.frame(subset_seurat_obj@assays[[i]]@scale.data)
    #scaled_counts <- as.data.frame(subset_seurat_obj@assays[[i]]@scale.data)
  } else {
    # If 'integrated' assay is not available, retrieve the scaled RNA data
    int_scaled_matrix <- NULL
    scaled_counts <- as.data.frame(subset_seurat_obj@assays[[i]]@scale.data)
  }
  
  # Retrieve raw counts, normalized data, scaled counts, and metadata
  raw_counts <- as.data.frame(subset_seurat_obj@assays[[i]]@counts)
  normalized_data <- as.data.frame(subset_seurat_obj@assays[[i]]@data)
  scaled_counts <- as.data.frame(subset_seurat_obj@assays[[i]]@scale.data)
  metadata <- as.data.frame(subset_seurat_obj@meta.data)
  
  # Add row names as a column for Alyssa Umap app
  raw_counts <- tibble::rownames_to_column(raw_counts, var = "Genes")
  normalized_data <- tibble::rownames_to_column(normalized_data, var = "Genes")
  metadata <- tibble::rownames_to_column(metadata, var = "Barcodes")
  
  # Write int_scaled_matrix or scaled_counts to a file based on availability
  if (!is.null(int_scaled_matrix)) {
    int_scaled_matrix <- tibble::rownames_to_column(int_scaled_matrix, var = "Genes")
    write.table(int_scaled_matrix, file = paste0(folder_path, "/", num_cells, "_int_scaled_counts_", i, ".txt"),
                sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  } else {
    scaled_counts <- tibble::rownames_to_column(scaled_counts, var = "Genes")
    write.table(scaled_counts, file = paste0(folder_path, "/", num_cells, "_rna_ScaledCounts_", i, ".txt"),
                sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  }
  
  # Write raw_counts, normalized_data, and metadata to separate files
  write.table(raw_counts, file = paste0(folder_path, "/", num_cells, "_raw_counts_", i, ".txt"),
              sep = "\t", row.names = FALSE, col.names = TRUE)
  write.table(normalized_data, file = paste0(folder_path, "/", num_cells, "_normalized_counts_", i, ".txt"),
              sep = "\t", row.names = FALSE, col.names = TRUE)
  write.table(metadata, file = paste0(folder_path, "/", num_cells, "_metafile_", i, ".txt"),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}


#####---- Find All Markers ----####
# Set the identity of the cells in the Seurat object
Idents(seurat_obj) <- "seurat_clusters_gabby_annotation"

# Find markers for each cluster with a minimum percentage threshold
all.markers <- FindAllMarkers(seurat_obj, min.pct = 0.1)

# Calculate the difference between pct.1 and pct.2 and add a new column
all.markers$diff_pct <- all.markers$pct.1 - all.markers$pct.2

# Specify the file path for saving the all markers data
all_marker_file_path <- file.path(full_output_dir, "all_markers.txt")

# Write the all markers data to a tab-separated file
write.table(all.markers, file = all_marker_file_path, 
            sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

