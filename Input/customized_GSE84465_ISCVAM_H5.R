create.visualization.artifacts <- 
  function(seurat, covs.discrete, covs.continuous=NULL) {
    #pcs = data.frame(seurat@reductions$pca@cell.embeddings)
    #colnames(pcs) <- gsub("_", "", colnames(pcs))
    
    qc.features = c("nFeature_RNA", "nCount_RNA" ,"percent.mt", "percent.rp")
    qc.stats = seurat[[qc.features]]
    
    cell.cycle.features <- c("S.Score", "G2M.Score")
    cell.stats <- seurat[[cell.cycle.features]]
    
    # pann.features <- c("pANN_0.25_0.09_152","pANN_0.25_0.09_144","pANN_0.25_0.09_350",
    #               "pANN_0.25_0.09_233","pANN_0.25_0.09_46")
    pann.features = c("pANN_0.25_0.09_142")
    pann.stats <- seurat[[pann.features]]
    
    # doublet.features <- c("DF.classifications_0.25_0.09_152","DF.classifications_0.25_0.09_144",
    #                       "DF.classifications_0.25_0.09_350","DF.classifications_0.25_0.09_46",
    #                       "DF.classifications_0.25_0.09_233")
    # doublet.stats <- seurat[[doublet.features]]
    
    # resolution.features <- c("integrated_snn_res.0.1","integrated_snn_res.0.15",
    #                          "integrated_snn_res.0.25" ,"integrated_snn_res.0.5",
    #                          "integrated_snn_res.0.75")
    
    resolution.features <- c("RNA_snn_res.0.1", "RNA_snn_res.0.15", "RNA_snn_res.0.25", "RNA_snn_res.0.75" )
                             #"RNA_snn_res.0.5") 
    
    
    resolution.stats <- seurat[[resolution.features]]
    
    tsne.features <- c("tSNE_1","tSNE_2")
    tsne.cords <- seurat[[tsne.features]]
    
    
    clustering.features = grep("Clusters",colnames(covs.discrete), value=T)
    clustering = data.frame(as.character(covs.discrete[[clustering.features]]))
    colnames(clustering) = "Clusters"
    
    
    umap.cols = grep("UMAP_.*",colnames(covs.discrete), value=T)
    umap.cords = covs.discrete[umap.cols]
    colnames(umap.cords) = tolower(colnames(umap.cords))
    
    covs.discrete <- covs.discrete[,-1]
    celltypes.cols <- grep("seurat_clusters_gabby_annotation", colnames(covs.discrete), value=T)
    # celltypes.cols <- grep("Celltype.*", colnames(covs.discrete), value=T)
    celltypes <- covs.discrete[celltypes.cols]
    # colnames(celltypes) <- c("Celltype.malignancy","Celltype.major.lineage", "Celltype.minor.lineage")
    
    clinical.cols <- covs.discrete[, !colnames(covs.discrete) %in% c(clustering.features, umap.cols, celltypes.cols)]
    
    clinical.covs <- data.frame(lapply(clinical.cols, as.character))
    
    if(is.null(covs.continuous)) {
      covs = cbind(clinical.covs, clustering,
                   celltypes, 
                   id=names(seurat$orig.ident),
                   umap.cords,qc.stats,cell.stats,
                   pann.stats,resolution.stats,tsne.cords)
    } else {
      covs = cbind(covs.continuous,
                   covs.discrete, clustering, 
                   celltypes, 
                   id=names(seurat$orig.ident),
                   umap.cords,qc.stats,cell.stats,
                   pann.stats,resolution.stats,tsne.cords)
    }
    
    discreteCovs <- c(
      colnames(clinical.covs),
      colnames(celltypes),
      colnames(resolution.stats),
      colnames(clustering)
    )
    
    continuousCovs <- c(
      colnames(covs.continuous),
      #colnames(tsne.cords),
      colnames(umap.cords),
      colnames(cell.stats),
      colnames(pann.stats),
      #colnames(resolution.stats),
      colnames(tsne.cords),
      colnames(qc.stats))
    
    
    ret <- list(covs=covs, discreteCovs=discreteCovs, continuousCovs=continuousCovs)
  }


write.h5.artifact <- function(seurat, fn, covs, artifactName) {
  require(rhdf5)
  artifacts <- create.visualization.artifacts(seurat, covs)
  h5createGroup(fn,paste0("artifacts/",artifactName))
  # h5write(artifacts$pcs, fn, paste0("artifacts/",artifactName,"/pcs"))
  i <- sapply(artifacts$covs, is.factor)
  artifacts$covs[i] <- lapply(artifacts$covs[i], as.character)
  
  h5write(artifacts$covs, fn, paste0("artifacts/",artifactName, "/covs"))
  h5write(artifacts$discreteCovs, fn, paste0("artifacts/",artifactName, "/discreteCovs"))
  h5write(artifacts$continuousCovs, fn,paste0("artifacts/",artifactName, "/continuousCovs"))
}


write.h5 <- function(seurat.list, fn) {
  unlink(fn)
  seurat = seurat.list$all$seurat  
  require(rhdf5)
  h5createFile(fn)
  
  h5createGroup(fn,"matrix")
  h5write(colnames(seurat), fn, "matrix/barcodes" )
  h5write(seurat@assays$RNA@counts@i, fn, "matrix/indices")
  h5write(seurat@assays$RNA@counts@p, fn, "matrix/indptr")
  h5write(seurat@assays$RNA@counts@x, fn, "matrix/data")
  h5write(seurat@assays$RNA@counts@Dim, fn, "matrix/shape")
  h5write(seurat@assays$RNA@counts@Dimnames[[1]], fn, "matrix/gene_names")
  
  h5createGroup(fn,"artifacts")
  for(group in names(seurat.list)) {
    write.h5.artifact(seurat.list[[group]]$seurat, fn, seurat.list[[group]]$covs, group)
  }
}