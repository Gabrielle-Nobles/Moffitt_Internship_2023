# Moffitt_Internship_2023
![Single Cell RNA-seq Analysis Pipeline and Post-processing h5-ISCVA compliant ](https://github.com/Gabrielle-Nobles/Single_Cell_Brain_Tumor_Atlas/assets/97853225/1b6c6eaf-00f4-4afb-be27-4c7037deb263)
![Post Processing gene expression and thousand matrix subset](https://github.com/Gabrielle-Nobles/Single_Cell_Brain_Tumor_Atlas/assets/97853225/c1983c40-f786-4f0d-9107-5d56d04ade04)
![Combine Data scRNAseq Analysis Pipeline](https://github.com/Gabrielle-Nobles/Single_Cell_Brain_Tumor_Atlas/assets/97853225/4600ae97-c63d-476e-9925-b9f4c4006c72)


## Introduction 
The utilization of single-cell RNA sequencing (scRNA-seq) has emerged as a robust approach to examining the complexities of cellular heterogeneity and genetic expression with distinctive resolution. However, analyzing scRNAseq data can be complex and requires a standardized approach. This protocol outlines a comprehensive workflow for analyzing scRNAseq data using R and the Seurat package, along with other tools like SeuratDisk, pathwork, dplyr, Single R, and Celldex. The protocol covers essential steps such as quality control, normalization, data scaling, dimensionality reduction, clustering, and automated cell annotation. Following this protocol allows bioinformaticians and researchers to effectively analyze scRNAseq data, identify distinct cell populations, and gain valuable insights into cellular dynamics and functions. The output files generated by this protocol, including metadata, H5 Seurat files, cell subpopulation metadata, and ISCVA-compliant files, facilitate downstream analyses and enable integration with other analysis and visualization tools. This protocol provides a standardized and reproducible framework for scRNAseq analysis. 


## Installation 

### R dependencies 

### Input Files 
To achieve accurate results in single-cell RNA sequencing (scRNAseq) analysis, it's crucial to have high-quality and compatible input data. The input file should be a raw count matrix of single cells, with each row representing a gene and each column representing a specific cell. The numbers within the matrix indicate the raw expression counts or read counts of genes in each cell. Ensure that the input file is in a compatible format, such as a Comma-Separated Values (CSV) file, tab-delimited text file or a Tab-Seperated Values (TSV) file, as these formats are easily readable and processed by R and the Seurat package. If you have used the CellRanger pipeline to process your scRNAseq data, you should have three input files: matrix (matrix.mtx), barcode (barcode.tsv), and features (features.tsv).

**Parameter file to read in the barcode.tsv, matrix.mtx and genes.tsv from each individual patient**
### Single Cell RNA-seq Analysis 
1. Quality Control
-- Hao et al.,Cell 2021.[PMID: 34062119]
- Perform quality control (QC) steps to filter out low-quality cells and genes. QC calculates the percentages of mitochondrial genes and ribosomal protein genes for each cell, and then filters out cells that have high mitochondrial gene content and low detected features in the RNA assay.  
3. Normalization 
- Hao et al.,Cell 2021.[PMID: 34062119]
- After completing the QC steps, it's important to normalize the count matrix to adjust for differences in library sizes and sequencing depth. This is a crucial step in analyzing scRNAseq data because it accounts for variations in library sizes and sequencing depth between individual cells. Normalization ensures that expression values can be compared accurately and eliminates technical biases. 
- Our method employs the "LogNormalize" parameter, which performs a global-scaling normalization. It divides each gene's expression by the total expression of that gene across all cells, multiplies the data by a default scale factor of 10,000, and log-transforms it. This normalization method is valuable for preserving the relative differences between cells while normalizing gene expression across them.
4. Cell Cycle Scores
- Hao et al.,Cell 2021.[PMID: 34062119]
- To determine cell cycle activity in scRNAseq data, the expression levels of certain genes associated with the cell cycle are analyzed and cell cycle scores are assigned to individual cells.
5. Scaling Data 
- Hao et al.,Cell 2021.[PMID: 34062119]
- Scaling normalizes and transform gene expression values, which helps to remove unwanted technical variation and improve the accuracy of downstream analyses
- The scaling method used by default is the "LogNormalize" method, which performs a natural logarithm transformation followed by centering and scaling of the gene expression values.
- During the scaling process, the variables specified in 'vars.to.regress' (nFeature_RNA and percent.mt in this case) are regressed out. This means that any variation in the gene expression values that can be attributed to these variables is removed.
6. Prinicpal Component Analysis (PCA)
- Hao et al.,Cell 2021.[PMID: 34062119]
- Conduct PCA on the scaled data to reduce the dimensionality of the dataset while preserving the most significant sources of variation. This step helps identify major sources of heterogeneity within the dataset. 
7. Nearest Neighbor 
- Compute the nearest neighbors for each cell in the reduced PCA space. This step helps identify cells that are likely to be biologically similar based on their expression profiles
8. SNN clustering 
- Hao et al.,Cell 2021.[PMID: 34062119]
- Perform clustering of the cells using the shared nearest neighbor (SNN) optimization based clustering algorithm. This algorithm group cells into clusters based on their similarity in the PCA space.
9. Doublet Finder 
- McGinnis et al.,CellPress.2019.[PMID: 30954475]
10. Automated Cell Annotation using SingleR and Celldex
- Aran et al.,Nature Immunology.2019.[PMID: 30643263]
- Leverage external reference datasets and computational tools like SingleR and Celldex to automatically annotate cell types or states. SingleR compares the gene expression of each cell to a reference dataset, while Celldex predicts cell type annotations based on a cell type reference database. These annotations provide biological context to the identified cell clusters.
### Output files 
- **Metadata**
1. Cell Cycle Score (S.score,G2M.score, and Phase) 
2. QC percentage (percent.rp and perdent.mt)
3. Doublets (pANN and DF.classfication) 
4. Resoulutions (integrated_snn_res or RNA_snn_res) 
5. Umap and Tsne coordinates 
6. Single R annotations 
- hpca.main and hpca.fine 
- dice.main and dice.fine 
- monacco.main and moncacco.fine 
- nothern.main and northern.fine
- blue.main and blue.fine
7. Manually curated cell annotations (seurat_clusters_gabby_annotations) 
- **H5 Seurat-compliant file** 
-- Export the processed and analyzed data in the H5 Seurat file format. This file contains the expression values, dimensionality reduction results, clustering information, and metadata. It serves as a comprehensive representation of the analyzed single-cell RNAseq data.
- **Subset of cell populations**
--Create separate metadata and H5 Seurat files for each identified cell subpopulation or cluster. This division facilitates downstream analyses focused on specific cell populations of interest.
- **H5 ISCVA-compliant file** 
--The H5 ISCVA-compliant file is a specific file format designed to load and interact with the Interactive Single Cell Visual Analytics (ISCVA) application. 

## Post Processsing: H5-ISCVA compliant file 
### Input File 
- H5 seurat object 
- Metadata file 
- Customized ISCVAM H5 file
### Write H5-ISCVA Compliant Pipeline 
- Edit the metadata file and seurat_obj@metadata
1. Open the metadata file and locate the column named "seurat_clusters".
2. Change the column name from "seurat_clusters" to "Clusters".
3. Append the following columns from the metadata file to the **seurat_obj@meta.data:**
- pANN scores
- DF.classifications
- Cell Cycle Scores (S.Score, G2M.Score)
- t-SNE coordinates
- QC percentages (percent.rp, percent.mt)
4. Identify any columns in the metadata file that should be removed to prevent duplicate columns in the H5 file.
5. Remove the selected columns from the metadata file.
6. Create a Seurat.list object containing the seurat_obj and the metadata file.
7. Use the source() function to input the functions from the customized ISCVAM H5 file.
8. Use the write.h5() function to convert the seurat_list to an H5-ISCVA compliant file.
9. Check the metadata columns to ensure they are arranged in the correct format for ISCVA.

### Output File 
- H5 ISCVA-compliant file 

## Post Processing: Gene expression 
### Input File
The H5 Seurat file typically contains essential information such as the expression matrix, cell metadata, dimensionality reduction results, clustering information, and other annotations relevant to the dataset. Loading the H5 Seurat file ensures that all the necessary data and attributes are available for subsequent analysis steps.
## Gene expresssion 
### Output Files 
- Gene expression .txt files for each cluster within the default resoultion (also an excel file of all gene expressions) 
- Gene expression for cell type (also an excel file of all gene expressions)
- Find all marker .txt files 
## Post Processing: 1000 Matrix for Umap application 
### Input File 
- **H5 Seurat-compliant file** 
-- Load in the processed and analyzed data in the H5 Seurat file format. This file contains the expression values, dimensionality reduction results, clustering information, and metadata. It serves as a comprehensive representation of the analyzed single-cell RNAseq data.

### Steps 
1. Read in scRNA-seq data for each sample: 
- Create file paths fro the matrix, barcodes and features files for each sample 
- Use the 'ReadMTX' function to the load the sparse data into the 'counts" object 
- Create a Seurat object using the 'CreateSeuratObject' function 
2. Preform data normalization: 
- Nomralize the data using the 'NormalizeData' function specifying the normalizatiion method (ex. "LogNormalize") and scale factor.
3. Perform cell cycle scoring:
- Use the 'CellCycleScoring' function to calculate cell cycle socres based on cell cycle gene lists (ex. G2M and S phase gene lists) 
4. Find highly variable features: 
- Use the 'FindVariableFeatures'  function to identify highly variable features in the data, specifying the selection method (e.g., "vst") and number of features to select 
5. Scale the data:
- 'ScaleData' function to scale the data
6. Perform dimensional reduction and clustering:
- Use the 'RunPCA' function to perform PCA on the scaled data and set the number of principal components (npcs) 
7. Perform doblet detection 
- Use the 'doubletFinder_v3' function to calulate the expected number of doubletes based on the percentage of cells 
8. Read in the metadata for each sample: 
- Create a file path for the metadata file 
- Use the 'read.table' function to read the metadata into a dataframe.
9. Add columns to the Seurat object consisting of:
- Manually curated annotations, Patient ID , Sample type, Organ, Organism, Age, Gender, IDH status, Extracted molecule, Extraction protocol, Library strategy, Library source, Library selection, Instrument model and Platform ID.
10. Combine the 3 Seurat objects together:
- Create a list of the Seurat objects 
- Select integration features across all Seurat Objects based on the variable features found in step 4
- Combine the data using the 'FindIntegrationAnchors' and "IntegrateData" functions, specifying the integration features and dimensions
- Store the integrated data in 'combined.integrated' 
11. Preform quality control on the combined data 
- Inside the function, define patterns for mitochondrial and ribosomal protein genes based on the species.
- Calculate QC metrics such as total features, expressed features, and mitochondrial gene percentage.
- Filter out low-quality cells based on the defined criteria.
12. Scale the combined data to perform PCA 
13. Generate UMAP and TSNE coordinates using 'RunUMAP' and 'RunTSNE', embed them into the metadata 
14. Find the nearest neighbors and cluster the combined data 
15. Perform cell type annotation: 
- Use SingleR to annotate cell types beased on gene expression profiles 
- The reference database used are: 
-- Human Primary Cell Atlas Data, Database Immune Cell Expression Data, Blueprint Encode Data, Monaco Immune Data, Novershtern Hematopietic Data
16. Save the combined data as a h5 compliant seurat object using 'SaveH5Seurat' 
17. Set the identify to the manually curated annotations then subset each cell type

Note: The protocol assumes that the necessart input files(matrix, barcoed features and metadata) are acailable in the specificed filed paths. 

### Output File 
- Scaled count matrix
- Raw count matrix 
- Normalized count matrix 
- metadat file


## Combining multiple count matrices
### Input files 
- Matrix (matrix.mtx), barcode (barcode.tsv), and features (features.tsv)
- Metadata for each count matrix 
### Output Files 
- Combined data metadata
- Combined data H5 seurat object
- Subset of cell populations h5 seuarat objects an metadata  







## Links to ISCVA Brain Tumor Datasets
 

http://zc:3838/apps/iscva/brain/GSE70630/

http://zc:3838/apps/iscva/brain/GSE89567/

http://zc:3838/apps/iscva/brain/GSE135437/

http://zc:3838/apps/iscva/brain/GSE125969/

http://zc:3838/apps/iscva/brain/GSE174401/

http://zc:3838/apps/iscva/brain/GSE173278/

http://zc:3838/apps/iscva/brain/GSE84465/

http://zc:3838/apps/iscva/brain/GSE141982/

http://zc:3838/apps/iscva/brain/GSE139948/

http://zc:3838/apps/iscva/brain/GSE179373/

http://zc:3838/apps/iscva/brain/GSE162631/

http://zc:3838/apps/iscva/brain/GSE138794/

http://zc:3838/apps/iscva/brain/GSE186344/

http://zc:3838/apps/iscva/brain/GSE119926/

http://zc:3838/apps/iscva/brain/GSE103224/

http://zc:3838/apps/iscva/brain/GSE185420/

http://zc:3838/apps/iscva/brain/GSE122871/

http://zc:3838/apps/iscva/brain/GSE129730/
