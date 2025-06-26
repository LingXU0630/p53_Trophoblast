*********************scRNA-seq*******************
*********************0.BGI

singularity exec dnbc4tools.sif dnbc4tools tools mkgtf --ingtf gencode.v32.primary_assembly.annotation.gtf --output gencode.v32.primary_assembly.annotation.genes.filter.gtf --type gene_type

dnbc4tools rna mkref --ingtf gencode.v32.primary_assembly.annotation.genes.filter.gtf --fasta GRCh38.primary_assembly.genome.fa --thread 10  --species Homo_sapiens 

dnbc4tools rna run \
	--cDNAfastq1 V350183467_L01_15_1.fq.gz,V350183467_L02_15_1.fq.gz,V350183467_L03_15_1.fq.gz,V350183467_L04_15_1.fq.gz \
	--cDNAfastq2 V350183467_L01_15_2.fq.gz,V350183467_L02_15_2.fq.gz,V350183467_L03_15_2.fq.gz,V350183467_L04_15_2.fq.gz \
	--oligofastq1 V350183467_L01_11_1.fq.gz,V350183467_L02_11_1.fq.gz,V350183467_L03_11_1.fq.gz,V350183467_L04_11_1.fq.gz \
	--oligofastq2 V350183467_L01_11_2.fq.gz,V350183467_L02_11_2.fq.gz,V350183467_L03_11_2.fq.gz,V350183467_L04_11_2.fq.gz \
	--genomeDir star_dir  \
	--name test --threads 10



*********************1.intergration

library(monocle3)
library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(magrittr)
library(Matrix)
library(CytoTRACE)
library(velocyto.R)
library(SeuratDisk)
library(reticulate)
library(tidydr)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(glmGamPoi)
library(sceasy)

EPSC <- readRDS("EPSC_umap.rds")
EPSCKO <- readRDS("EPSCKO_umap.rds")
TSC <- readRDS("TSC_umap.rds")
TSCKO <- readRDS("TSCKO_umap.rds")
d2 <- readRDS("d2_umap.rds")
kd2 <- readRDS("kd2_umap.rds")
d7s <- readRDS("d7s_umap.rds")
kd7s <- readRDS("kd7s_umap.rds")

seurat_list_wt <-list(EPSC,d2,d7s,TSC)
names(seurat_list_wt)<- c("EPSC","d2","d7s","TSC")

merged_seurat_wt <- merge(x = seurat_list_wt[[1]],
		       y = seurat_list_wt[2:length(seurat_list_wt)],
		       merge.data = TRUE)

seurat_list_ko <-list(EPSCKO,kd2,kd7s,TSCKO)
names(seurat_list_ko)<- c("EPSCKO","kd2","kd7s","TSCKO")

merged_seurat_ko <- merge(x = seurat_list_ko[[1]],
		       y = seurat_list_ko[2:length(seurat_list_ko)],
		       merge.data = TRUE)

seurat_list <-list(merged_seurat_wt,merged_seurat_ko)

seurat_list <- lapply(X = seurat_list, FUN = SCTransform)
SCT.features <- SelectIntegrationFeatures(object.list = seurat_list, nfeatures = 2000)
seurat_list <- PrepSCTIntegration(object.list = seurat_list, anchor.features = SCT.features)

SCT.anchors <- FindIntegrationAnchors(object.list = seurat_list, normalization.method = "SCT",anchor.features = SCT.features)
combined.sct <- IntegrateData(anchorset = SCT.anchors, normalization.method = "SCT")

combined.sct <- RunPCA(combined.sct,npcs = 20, verbose = FALSE)
combined.sct <- JackStraw(combined.sct, num.replicate = 100);combined.sct <- ScoreJackStraw(combined.sct, dims = 1:20)
p<-ElbowPlot(combined.sct);ggsave("pca.pdf",p,width=12,height=8)

combined.sct <- FindNeighbors(combined.sct, reduction = "pca", dims = 1:17)
combined.sct <- FindClusters(combined.sct, resolution = 0.1)
combined.sct <- RunUMAP(combined.sct, reduction = "pca", dims = 1:17)
combined.sct <- RunTSNE(combined.sct, reduction = "pca", dims = 1:17)


*********************2. rename cluster and Find marker genes
combined.sct <- FindClusters(combined.sct, resolution = 1)
xlcluster <- factor(
  as.character(combined.sct@meta.data$integrated_snn_res.1),
  levels = as.character(0:17),
  labels = c(
    "C6-TSC",  # 0 
    "C0-EPSC", # 1
    "C0-EPSC", # 2
    "C3-D7",   # 3
    "C6-TSC",  # 4
    "C0-EPSC", # 5
    "C2-D2",   # 6
    "C0-EPSC", # 7
    "C6-TSC",  # 8
    "C2-D2",   # 9
    "C2-D2",   # 10
    "C1-D2",   # 11
    "C0-EPSC", # 12
    "C2-D2",   # 13
    "C5-TSC",  # 14
    "C4-D7",   # 15
    "Other",   # 16
    "C7-TSC"   # 17
  )
)

combined.sct@meta.data$xlcluster<-xlcluster
Idents(combined.sct)<-xlcluster

DefaultAssay(combined.sct) <- "RNA"
for (i in unique(Idents(combined.sct))) {
   markers <- FindConservedMarkers(intergateTSC, ident.1 = i, grouping.var = "condition", verbose = TRUE)
   write.csv(markers,paste0("combined.sct_conserve_markers_cluster",i,"csv"))  
  }

saveRDS(combined.sct, file = "wtko_intergrate.rds")

intergateTSC<- subset(combined.sct,idents=c("C5-TSC","C6-TSC","C7-TSC"))
****TSC c6 VS c7
combined.sct.TSCWT.markers67 <- FindMarkers(intergateTSCWT,assay="RNA",min.pct = 0, logfc.threshold = 0,ident.1="C6-TSC",ident.2="C7-TSC")
write.csv(combined.sct.TSCWT.markers07,"combined.sct2000_CR5.markersTSCWT_07_RNA.csv")

****TSC c5 VS c6
combined.sct.TSCWT.markers67 <- FindMarkers(intergateTSCWT,assay="RNA",min.pct = 0, logfc.threshold = 0,ident.1="C5-TSC",ident.2="C6-TSC")
write.csv(combined.sct.TSCWT.markers67,"combined.sct2000_CR5.markersTSCWT_60_RNA.csv")

**** DEGs between conditions

combined.sct$condition_xl <- paste(combined.sct$condition, combined.sct$xlcluster, sep = "_")
Idents(combined.sct)      <- "condition_xl"
DefaultAssay(ccombined.sct) <- "RNA"

clusters <- c(
  "C0-EPSC","C1-D2","C2-D2","C3-D7",
  "C4-D7","C5-TSC","C6-TSC","C7-TSC"
)

out_dir <- "./scRNAseq-TSC/scEPSC"
for (cl in clusters) {
  markers <- FindMarkers(
    combined.sct,
    ident.1 = paste0("WT_", cl),
    ident.2 = paste0("KO_", cl),
    verbose = FALSE
  )
  fname <- file.path(
    out_dir,
    paste0("wt_vsko_", gsub("-", "_", cl), "_marker_genes.csv")
  )
  write.csv(markers, file = fname, row.names = FALSE)
}


*********************3.monocle3
Idents(combined.sct) <-"xlcluster"
sling<- subset(combined.sct,idents=c("C5-TSC","C6-TSC","C7-TSC")) 

cds <- as.cell_data_set(sling)
cds <- cluster_cells(cds)
p1 <- plot_cells(cds, show_trajectory_graph = FALSE)
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
p=p1|p2
ggsave("reduction_compare.pdf",plot=p,width=10,height=5)

cds <- learn_graph(cds)
get_earliest_principal_node <- function(cds, time_bin="C6-TSC"){
  cell_ids <- which(colData(cds)[, "xlcluster"] == time_bin)
  
  closest_vertex <-
  cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
  igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
  (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}
cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))
pdf("preprocess_method_PCA_umap_pseudotime.pdf")
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=TRUE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=1.5)
dev.off()


*********************4. CytoTRACE

Idents(combined.sct) <-"orig.ident"
ko.sling1<- subset(combined.sct,idents=c("KO"))
Idents(ko.sling1) <-"xlcluster"
ko.sling<- subset(ko.sling1,idents=c("C5-TSC","C6-TSC","C7-TSC")) 
int.embed.ko<-Embeddings(ko.sling,reduction="umap")
exp<- GetAssayData(ko.sling,assay="RNA");exp=as.data.frame(exp)
results <- CytoTRACE(exp,ncores = 8, subsamplesize = 3000)

source("plotCytoTRACE1.R")
pdf("TSCKO_remove4cells_color_reverse.pdf",width=8,height=5)
plotCytoTRACE1(results,emb=int.embed.ko)
dev.off()
#color=c("#4583b3","#f5b375","#c63696")

Idents(combined.sct) <-"orig.ident"
wt.sling1<- subset(combined.sct,idents=c("WT"))
Idents(wt.sling1) <-"xlcluster"
wt.sling<- subset(wt.sling1,idents=c("C5-TSC","C6-TSC","C7-TSC")) 
int.embed.wt<-Embeddings(wt.sling,reduction="umap")
exp<- GetAssayData(wt.sling,assay="RNA");exp=as.data.frame(exp)
results <- CytoTRACE(exp,ncores = 8, subsamplesize = 3000)
source("plotCytoTRACE1.R")
pdf("TSC_remove4cells_color_reverse.pdf",width=8,height=5)
plotCytoTRACE1(results,emb=int.embed.wt)
dev.off()

Idents(combined.sct) <-"xlcluster"
sling<- subset(combined.sct,idents=c("C5-TSC","C6-TSC","C7-TSC")) 
int.embed.wt<-Embeddings(sling,reduction="umap")
exp<- GetAssayData(sling,assay="RNA");exp=as.data.frame(exp)
results <- CytoTRACE(exp,ncores = 8, subsamplesize = 3000)
source("plotCytoTRACE1.R")
pdf("TSCall_remove4cells_color_reverse.pdf",width=8,height=5)
plotCytoTRACE1(results,emb=int.embed.wt)
dev.off()


*********************5. velocity

***velocyto
samtools sort -t CB -O BAM \
-o WT_cellsorted_possorted_genome_bam.bam \
WT_possorted_anno_decon_sorted.bam

velocyto run WT_possorted_genome_bam.bam gencode.v32.primary_assembly.annotation.gtf  \
-b filtered_feature_bc_matrix/barcodes.tsv.gz  -o velocyte  

samtools sort -t CB -O BAM \
-o KO_cellsorted_possorted_genome_bam.bam \
KO_possorted_anno_decon_sorted.bam

velocyto run KO_possorted_genome_bam.bam gencode.v32.primary_assembly.annotation.gtf  \
-b filtered_feature_bc_matrix/barcodes.tsv.gz  -o velocyte  

**** In python
conda activate dnbc4tools-1
import anndata
import scvelo as scv
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
%load_ext rpy2.ipython

TSC_WT = anndata.read_loom("TSC-WT.loom")
TSC_KO = anndata.read_loom("TSC-KO.loom")

wt_obs = pd.read_csv("TSCWT_cellID_obs.csv")
ko_obs = pd.read_csv("TSCKO_cellID_obs.csv")

TSC_WT_F = TSC_WT[np.isin(TSC_WT.obs.index,wt_obs["x"])]
TSC_KO_F = TSC_KO[np.isin(TSC_KO.obs.index,ko_obs["x"])]

TSC_WT_F.var_names_make_unique()
TSC_KO_F.var_names_make_unique()

TSC = TSC_WT_F.concatenate(TSC_KO_F)
TSC_index = pd.DataFrame(TSC.obs.index)
TSC_index.to_csv("TSC_index-loom.csv",index=False,sep=',')

umap = pd.read_csv("TSC_loomOrder_embedding_obs.csv")
clusters = pd.read_csv("TSC_loomOrder_cluster_obs.csv")
TSC.obsm['X_umap'] = umap.values
TSC.obs['clusters']=clusters.values

scv.pp.filter_and_normalize(TSC,min_shared_counts = 20, n_top_genes = 2000)
scv.pp.moments(TSC, n_pcs=30, n_neighbors=30)
scv.tl.velocity(TSC, mode = "stochastic")
scv.tl.velocity_graph(TSC,sqrt_transform = True)
scv.tl.velocity_pseudotime(TSC)
scv.pl.scatter(TSC, color = 'velocity_pseudotime', cmap = 'gnuplot', show=False,palette=["#4583b3","#c63696"])
plt.savefig('TSC_RNA-velo-pseudotime_color.png')
scv.pl.proportions(TSC)
scv.pl.velocity_embedding(TSC, basis = 'umap')
scv.pl.velocity_embedding_stream(TSC, basis = "umap", color = "clusters", legend_loc = "best", dpi = 150, show=False,palette=["steelblue","#a38cbd","#f7afb9","#c63596"])
plt.savefig('TSC_RNA-velo-color.svg',dpi=100,format="svg",bbox_inches="tight")



*********************6. cell cycle

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
DefaultAssay(combined.sct)<-"RNA"
combined.sct <- CellCycleScoring(combined.sct, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

metrics <-  c("nFeature_RNA", "S.Score", "G2M.Score", "percent.mt")

p<- FeaturePlot(combined.sct, 
            reduction = "umap", 
            features = metrics,
            pt.size = 0.4, 
            min.cutoff = 'q1',
            label = TRUE)

ggsave("wtko-e27t_SCT3000_intergrate_umap_RC5_xlcluster_metrics.pdf",p,width=10,height=10)



*********************7. projection of wtko_intergrate.rds based on in vivo data

sceasy::convertFormat(
  "reference_fixed.h5ad",
  from = "anndata",
  to   = "seurat",
  outFile = "reference_fixed.h5ad.rds"
)

invio <- readRDS("reference_fixed.h5ad.rds")
Idents(invio)="final_annot_all_troph_corrected"
DefaultAssay(invio)<-"RNA"
p1 <- DimPlot(invio, reduction = "umap_model", label = FALSE,alpha=0.7)+
	tidydr::theme_dr(xlength=0.2,ylength=0.2,arrow=arrow(length=unit(0.1,"inches"),type="closed"))+
	theme(panel.grid=element_blank(),axis.title=element_text(face=1,hjust=0.03,size=7))+NoLegend()

ggsave("reference_fixed.h5ad_reUMAP.pdf",plot=p1,width=5,height=5)


invio <- NormalizeData(invio)
invio <- RunUMAP(invio, reduction="pca",dims = 1:8, reduction.name = "umap_model",  reduction.key = "UMM_", return.model = TRUE)
VariableFeatures(invio)=rownames(invio)

combined.sct <- NormalizeData(combined.sct, verbose = FALSE)
combined.sct <- FindVariableFeatures(combined.sct, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
combined.sct <- ScaleData(combined.sct, features = VariableFeatures(invio), verbose = FALSE)
combined.sct <- RunPCA(combined.sct, dims=1:8, verbose = FALSE)

features=intersect(rownames(invio),rownames(combined.sct))
anchors <- FindTransferAnchors(
  reference            = invio,
  query                = combined.sct,
  reference.reduction  = "pca",
  dims                 = 1:8,
  features      = features
)
qry_mapped <- MapQuery(
  anchorset           = anchors,
  reference           = invio,
  query               = combined.sct,
  refdata             = list(celltype = "final_annot_all_troph_corrected"),
  reference.reduction = "pca",
  reduction.model     = "umap_model",
  new.reduction.name  ="ref.umap"
)
ref_emb <- Embeddings(invio, "umapxl") %>% as.data.frame() %>% mutate(ds="invio")
qry_emb <- Embeddings(qry_mapped, "ref.umap") %>% as.data.frame() %>% mutate(ds="combined.sct")
colnames(qry_emb)=c("UMAP_1","UMAP_2","ds")
coords <- rbind(ref_emb, qry_emb)
p <-ggplot(coords, aes(UMAP_1, UMAP_2, color = ds)) +
  geom_point(alpha = 0.5) +
  theme_classic()
ggsave("combined.sct_paper_invivoP13_umap.pdf",plot=p,width=18,height=6)



*********************8. invivo cell cycle

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
DefaultAssay(invio)<-"RNA"
invio <- CellCycleScoring(invio, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

metrics <-  c("S.Score", "G2M.Score")
Idents(invio)="Phase"
p <- DimPlot(invio, reduction = "umap", label = FALSE,repel = TRUE)
ggsave("invio_TEST_phase2.pdf",p,width=9,height=8)

p<- FeaturePlot(invio, 
            reduction = "umap", 
            features = metrics,
            pt.size = 0.4, 
            min.cutoff = 'q1',
            label = TRUE)

ggsave("invio_cycle.pdf",p,width=10,height=10)


library(dplyr)
df=data.frame(invio@meta.data)
data=df %>%
  count(final_annot_all_troph_corrected, Phase) 

write.csv(data,"invio_cell_cycle_summary.csv")
************************************below is figure plotting ************************************


***1. umap
#(combined.sct,"/biomedja01/disk1/lxu/viola/scRNAseq-TSC/scEPSC/wtko-e27t_SCT3000_intergrate_umap_RC5_xlcluster_res4.rds")
Idents(combined.sct)<- "xlcluster"
Idents(combined.sct)<- "condition"
wt<- subset(combined.sct,idents=c("WT"))
ko<- subset(combined.sct,idents=c("KO"))
Idents(wt)<- "xlcluster"
Idents(ko)<- "xlcluster"

scale_colour_gradient2(mid="#99cdce", low="lightgrey", high="#c63696") +theme_bw()+
colors<-c("#4583b3","#51c4c2","#0d8a8c","#4a9d47","#f5b375","#a38cbd","#f7afb9","#c63596")#used
p1 <- DimPlot(wt, reduction = "umap", label = FALSE,cols=colors,alpha=0.7)+
	tidydr::theme_dr(xlength=0.2,ylength=0.2,arrow=arrow(length=unit(0.1,"inches"),type="closed"))+
	theme(panel.grid=element_blank(),axis.title=element_text(face=1,hjust=0.03,size=7))+NoLegend()
ggsave("./manuscript_plot/wt_reduction_umap_color.pdf",plot=p1,width=5,height=5)

scale_colour_gradient2(mid="#99cdce", low="lightgrey", high="#c63696") +theme_bw()+
colors<-c("#4583b3","#51c4c2","#0d8a8c","#4a9d47","#f5b375","#a38cbd","#f7afb9","#c63596")#used
p1 <- DimPlot(ko, reduction = "umap", label = FALSE,cols=colors,alpha=0.7)+
	tidydr::theme_dr(xlength=0.2,ylength=0.2,arrow=arrow(length=unit(0.1,"inches"),type="closed"))+
	theme(panel.grid=element_blank(),axis.title=element_text(face=1,hjust=0.03,size=7))+NoLegend()
ggsave("./manuscript_plot/ko_reduction_umap_color.pdf",plot=p1,width=5,height=5)


Idents(combined.sct)<- "orig.ident"
colors<-c("#3477a9","#67a59b","#4a9d47","#f7afb9","#3477a9","#67a59b","#4a9d47","#f7afb9")#used
p1 <- DimPlot(combined.sct, reduction = "umap", label = FALSE,cols=colors,alpha=0.7)+
	tidydr::theme_dr(xlength=0.2,ylength=0.2,arrow=arrow(length=unit(0.1,"inches"),type="closed"))+
	theme(panel.grid=element_blank(),axis.title=element_text(face=1,hjust=0.03,size=7))+NoLegend()
ggsave("./manuscript_plot/timepoint_reduction_umap_color.pdf",plot=p1,width=5,height=5)


***2.feature_plot
markers4 <- c("POU5F1",	"NANOG", "HIST1H2AC",	"SESN3",	"FOXP1",	"CDX2",	"ISL1",	"BMP4",	"GABRP",	"IGFBP3",	"ERVMER34-1",	"HSD3B1",	"CCKBR",	"LRP2",	"RASA1",	"PEG10",	"FOS",	"GATA2",	"GATA3",	"ERVFRD-1",	"OVOL1",	"PGF")

p<-DotPlot(wt, features = markers4,col.min=0,col.max=1) + geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
scale_colour_gradient2(mid="#99cdce", low="lightgrey", high="#c63696") +theme_bw()+
 theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,colour = "black"),axis.text.y = element_text(colour = "black"))+
 guides(shape = guide_legend(override.aes = list(size = 0.5)))  +
        theme(legend.title = element_text(size = 6), 
              legend.text  = element_text(size = 6),legend.key.size = unit(0.3, 'cm'),legend.justification = "right",
    legend.direction = "vertical",panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_rect(colour = "black", size=1))

ggsave("./manuscript_plot/wt_all-clusters_related_markers4_genes-dotplot.pdf",p,width=8,height=3)

***3.cell_abandance_metadata

meta<- combined.sct@meta.data
xlcluster<-meta$xlcluster
names(xlcluster)<-meta$condition
table(xlcluster)

# Proportion / cell number composition per cluster
ggData = data.frame(prop.table(table(meta$condition, meta$xlcluster), margin = 2))
colnames(ggData) = c("library", "cluster", "value")
ggData=ggData[-(17:18),]###remove other

p1 <- ggplot(ggData, aes(cluster, value, fill = library)) +
  geom_col() + xlab("Cluster") + ylab("Cells fraction") +
  scale_fill_manual(values = c("red1","steelblue3"))  + coord_flip()

#ggsave("./manuscript_plot/all-clusters_cell_propotion.pdf",p1,width=3,height=3)

#ggData = data.frame(table(meta$condition, meta$xlcluster_res))
ggData = data.frame(prop.table(table(meta$condition, meta$xlcluster_res4), margin = 2))
colnames(ggData) = c("library", "cluster", "value")
ggData=ggData[-(19:20),]###remove other

p2 <- ggplot(ggData, aes(cluster, value, fill = library)) +
  geom_col() + xlab("Cluster") + ylab("Cells fraction") +
  scale_fill_manual(values = c("red1","steelblue3")) + coord_flip()
ggsave(p1 + p2 + plot_layout(guides = "collect"), 
       width = 6, height = 3, filename = "./manuscript_plot/all-clusters_cell_propotion.pdf")

***4.heatmap plot of clusters DEGs

dat<- read.csv("./manuscript_plot/cluster_wtko_DEGs.csv",row.names=1)
dat <- dat[match(dat[,6],rownames(dat)),1:5]
dat <- as.matrix(dat)
col_fun = colorRamp2(c(-3,0,3), c("#99cdce","white","#c63696"))
id<-c(rep("p53_signalling",6),rep("Trophoblast_markers",16),rep("Amnion-like",10))
row_ha = rowAnnotation(bar =id ,col = list(bar = c("p53_signalling" = "#4583b3","Trophoblast_markers"="#f2db96","Amnion-like"="tan3")))
pdf("all_cluster0-4_DEGs_noscale_adjusted_color.pdf", width = 3.5, height =6)
Heatmap(dat, na_col = "white",col = col_fun,left_annotation = row_ha,cluster_columns = FALSE,cluster_rows = FALSE,border=FALSE,row_names_gp=gpar(fontsize = 7))
dev.off()

dat<- read.csv("./manuscript_plot/cluster5-6_wtko_DEGs.csv",row.names=1)
dat <- dat[match(dat[,6],rownames(dat)),1:5]
dat <- as.matrix(dat)
col_fun = colorRamp2(c(-3,0,3), c("#99cdce","white","#c63696"))
pdf("all_cluster5-7_DEGs_noscale_adjusted_color.pdf", width = 2, height =6)
Heatmap(dat, na_col = "white",col = col_fun,cluster_columns = FALSE,cluster_rows = FALSE,border=FALSE,row_names_gp=gpar(fontsize = 7))
dev.off()

***5.dotplot
dat<- read.csv("./manuscript_plot/cluster_wtko_DEGs.csv",row.names=1)
gene=row.names(dat)
cIdents(combined.sct)<- "xlcluster"
combined<-subset(combined.sct,idents=c("C0-EPSC","C1-D2","C2-D2","C3-D7", "C4-D7","C5-TSC", "C6-TSC", "C7-TSC"))
Idents(combined)<- "condition"
wt<- subset(combined,idents=c("WT"))
ko<- subset(combined,idents=c("KO"))
Idents(wt)<- "xlcluster"
Idents(ko)<- "xlcluster"

p <-DotPlot(wt, features = gene,col.min=0,col.max=1) + geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
scale_colour_gradient2(mid="#99cdce", low="lightgrey", high="#c63696") +theme_bw()+
 theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,colour = "black"),axis.text.y = element_text(colour = "black"))+
 guides(shape = guide_legend(override.aes = list(size = 0.5)))  +
        theme(legend.title = element_text(size = 6), 
              legend.text  = element_text(size = 6),legend.key.size = unit(0.3, 'cm'),legend.justification = "right",
    legend.direction = "vertical",panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_rect(colour = "black", size=1))

ggsave("./manuscript_plot/revised-F4b-wt-dotplot.pdf",p,width=9,height=2.8)

p <-DotPlot(ko, features = gene,col.min=0,col.max=1) + geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
scale_colour_gradient2(mid="#99cdce", low="lightgrey", high="#c63696") +theme_bw()+
 theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,colour = "black"),axis.text.y = element_text(colour = "black"))+
 guides(shape = guide_legend(override.aes = list(size = 0.5)))  +
        theme(legend.title = element_text(size = 6), 
              legend.text  = element_text(size = 6),legend.key.size = unit(0.3, 'cm'),legend.justification = "right",
    legend.direction = "vertical",panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_rect(colour = "black", size=1))

ggsave("./manuscript_plot/revised-F4b-ko-dotplot.pdf",p,width=9,height=2.8)

