# NOTE : De_Op integrated data in 'INMG_SingleCell/' have been processed only in version that
#       had not yet used QC infos to FILTER out 'bad cells'.
# Instead, the following is done by the present script:
#         1. open RAW data.
#         2. create 'raw' seurat object  (as usual, very first step)
#         3. as we already have QC tables into 'qcdoubl' folder, 
#             join by barcodes that info AND seu metadata.
#         4. Filter seu object (reject low QC barcodes) and run clustering..etc
#         5. This FILTERED seu is available to:
#             impute cell types  (by using 'impute)
#              and  do GO and LigPlot
# --
# JohaGL
library(tidyverse)
library(Seurat)
library(ggplot2)
library(sctransform)
library(cowplot)
library(RColorBrewer)
source(file="~/INMG_SingleCell/scripts/functions_stock.R", local=T)
SP = "~/LigRecPairs_GO/"  # sibling project, OUR WORKING DIR HERE
prloc = "~/INMG_SingleCell/"  
prefix = "FI-DeOp"

setwd(SP)  ## attention! this time all results will be saved into this sibling project

resu = "results_DeOpNEWseu/"
rdsdir = "rds_NEWseu/" # TODO: add .gitignore because heavy
nbDIM=15
resol=0.4
threshold = 0.5 #threshold log2foldChange for saving markers tables
permit = 15
newfeat.intgr = 6000

dir.create(resu,recursive=T)
dir.create(rdsdir, recursive=T)

# ================================================================================
print("initializing raw seurat objects")
# ================================================================================

dmich.f <- list.files(paste0(prloc,"data/DeMicheliD0/"), pattern="D0.txt", full.names = T)
oprescu.f <- list.files(paste0(prloc,"data/OprescuD0"), pattern="Noninjured_raw.txt$", full.names = T)

dmich = getMatrixFromtxt(dmich.f)
oprescu = getMatrixFromtxt(oprescu.f)
dim(dmich)
#[1] 19208  7028
dim(oprescu)
#[1] 19311  5670
rm(dmich,oprescu)

authors = c("DeMicheliD0", "OprescuD0") 
fi.name = c(dmich.f, oprescu.f)
names(fi.name) = authors

print("checking vector as you defined it, to run analysis")
print(fi.name)
setwd(prloc)

print("opening files, loading seu objects into list  and filtering")
muscle.list <- list()
muscle.list <- lapply(names(fi.name), function(auth){
  cgmat <- getMatrixFromtxt(fi.name[[auth]])
  seu <- CreateSeuratObject(cgmat, project=auth, min.cells=3,min.features = 150)
  return(seu)
})

# ================================================================================
print("adding QC infos to each seurat object meta.data")
# ================================================================================

print("retreiving QC and formatting")
# --------------------------------------------------------------------------------
txtallD0FINDER=paste0(prloc,"qcdoubl/spli_FINDER/TABLE_DOUBLFINDER_SPLITTED_4D0.txt")
txtallD0SCRAN=paste0(prloc, "qcdoubl/spli_SCRAN/TABLE_DOUBLETS_SCRAN_splitted_4D0.txt")
allD0finder = read.table(txtallD0FINDER, sep="\t",  header=T) 
allD0scran = read.table(txtallD0SCRAN, sep="\t",header=T)
scran.1 = allD0scran %>% filter((startsWith(as.character(barcode),"D0_")) |
                                    (startsWith(as.character(barcode),"Noninj")))
finder.1 = allD0finder %>% filter((startsWith(as.character(id),"D0_")) |
                                    (startsWith(as.character(id),"Noninj")))

finder.1$id = as.character(finder.1$id)
scran.1$barcode = as.character(scran.1$barcode)
joinedQC = left_join(scran.1, finder.1, by=c("barcode" ="id"))
print("if some NA introduced, the reason is: scran has ALL barcodes, whereas Finder imposed cutoff; fix the script that does DoubletsFinder step")


print("intersecting doublets; passing barcodes to rownames")
dim(joinedQC); dim(joinedQC %>% distinct()) # to check dims; to check no row duplication
joinedQC <- joinedQC %>% mutate(DOUBL_INTERSECT =  case_when(
  DFclass == "Doublet" & classific == "doublets" ~ "Doublet",
  TRUE ~ "Singlet"))
# selecting only useful columns for this analysis:
joinedQC <- joinedQC %>% select(barcode, sample, DOUBL_INTERSECT, is_cell, is_inf_outlier)
rownames(joinedQC) = joinedQC$barcode

opre.joinedQC = joinedQC %>% filter(startsWith(as.character(barcode),"Noninj"))
dmich.joinedQC = joinedQC %>% filter(startsWith(as.character(barcode),"D0_"))

QClist = list(dmich.joinedQC, opre.joinedQC)

rm(allD0finder, allD0scran, finder.1,scran.1, joinedQC, opre.joinedQC, dmich.joinedQC)
# --------------------------------------------------------------------------------


print("adding QC infos to meta.data (each object)")
# --------------------------------------------------------------------------------
for (i in 1:length(muscle.list)){
  muscle.list[[i]]@meta.data$barcode = rownames(muscle.list[[i]]@meta.data)
  tmp = left_join(muscle.list[[i]]@meta.data, QClist[[i]], by="barcode")
  rownames(tmp) = tmp$barcode
  if(all(rownames(tmp) == rownames(muscle.list[[i]]@meta.data))){
    muscle.list[[i]]@meta.data <- tmp
    rm(tmp)
  }else{print("something went wrong when assigning QC infos to meta.data")}
}
rm(QClist)
# --------------------------------------------------------------------------------

# ================================================================================
print("filtering out bad quality barcodes == yielding filtered seurat objects")
# ================================================================================
##### filtering (finally): 
muscle.listFI <- mapply( function(seu,auth){
  seu@meta.data[["orig.ident"]] <- rep(auth,length(seu@meta.data[["orig.ident"]]))
  seu[["percent.mt"]] <- PercentageFeatureSet(seu,pattern="^mt-")
  seu <- subset(seu, subset= DOUBL_INTERSECT != "Doublet" &  is_cell != F &
    is_inf_outlier != T & percent.mt < permit)
  return(seu)}, muscle.list, names(fi.name)
)
rm(muscle.list)
# ====================================================================================
# integration
print("integrating...")
# ====================================================================================
title <- ''; for (i in 1:length(authors)){title <- paste(title, authors[i])}
print(paste0("Integrating ",title," data. This is the main integration job, 
  reason: both on tibialis anterior m, similar experimental design"))
print(paste0("Will use *",nbDIM,"* princ components (dims), accordingly to  ", 
             "elbowplot and good balance identities-nbclusters seen in independent runs."))
print(paste0("   Using nbfeatures= ",newfeat.intgr))


muscle.listFI <- lapply(muscle.listFI, function(seu){
  seu <- NormalizeData(seu,verbose=FALSE)
  seu <- FindVariableFeatures(seu, selection.method="vst",
                              nfeatures= newfeat.intgr, verbose=FALSE)
  return(seu)
})

muscle.anchors <- FindIntegrationAnchors(object.list = muscle.listFI, dims=1:30)
muscle.integrated <- IntegrateData(anchorset = muscle.anchors, dims = 1:30)

# switch to integrated assay. The variable features of this assay are
# automatically set during IntegrateData
DefaultAssay(object = muscle.integrated) <- "integrated"

# =====================================================================================
# find clusters
# =====================================================================================
print("scaling, reducing dim (PCA), printing ElbowPlot into preprocessIntegrated.pdf")
muscle.integrated <- ScaleData(muscle.integrated, verbose = FALSE)
muscle.integrated <- RunPCA(muscle.integrated, npcs = 30, verbose = FALSE)
#visuals
pdf(paste0(resu,prefix,"_preprocessIntegrated.pdf"))
ElbowPlot(muscle.integrated)
DimHeatmap(muscle.integrated, dims=1:30, cells = 500, balanced = TRUE)
dev.off()

muscle.integrated <- KNNplusFITSNE(muscle.integrated, nbDIM, resol)
# =====================================================================================
# also run UMAP (slot FITSNE kept available)
# =====================================================================================
muscle.integrated <- RunUMAP(muscle.integrated, umap.method="uwot", dims=1:nbDIM)

saveRDS(muscle.integrated,file=paste0(rdsdir, prefix, "_integrated_seu_fitsne.rds"))
# =====================================================================================
# markers : both, positive and negative LFC
# =====================================================================================
muscle.integrated.markers <- FindAllMarkers(muscle.integrated, only.pos = FALSE, 
                                            min.pct = 0.25, logfc.threshold = 0.25) #time consuming
# save into two tables:

positive <- muscle.integrated.markers %>% filter(avg_logFC > threshold)
negative <- muscle.integrated.markers %>% filter(avg_logFC < -threshold)

write.table(positive, paste0(resu,prefix,"_ALLMARKERS_Pos_integratedD0.txt"), sep="\t")
write.table(negative, paste0(resu,prefix,"_ALLMARKERS_Neg_integratedD0.txt"), sep="\t")

# example getting 8 most relevant genes Pos and Neg :
most.neg <- negative %>% group_by(cluster) %>% top_n(n=4, wt= -avg_logFC)
most.pos <- positive %>% group_by(cluster) %>% top_n(n=4, wt= avg_logFC)

# =====================================================================================
# preliminary plots into results
# =====================================================================================

pdf(paste0(resu,prefix, "_FITSNEandUMAP.pdf"), width=13)
fit <- DimPlot(muscle.integrated, reduction="tsne", label=T,
               cols=definecolors(muscle.integrated@active.ident)) + 
  ggtitle(title) + theme(title=element_text(size=9))
uma <- DimPlot(muscle.integrated, reduction="umap", label=T,
               cols=definecolors(muscle.integrated@active.ident)) + 
  ggtitle(title) + theme(title=element_text(size=9))
plot_grid(fit,uma)
dev.off()

pos4top <- positive %>% group_by(cluster) %>% top_n(n=4, wt= avg_logFC)
pdf(paste0(resu,prefix,"_HEATMAP.pdf"), width=13)
DoHeatmap(muscle.integrated, features = pos4top$gene, 
          group.colors=definecolors(muscle.integrated@active.ident)) 
dev.off()
print("finished")

# ====================================================================================





 

