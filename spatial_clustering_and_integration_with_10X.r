library(Seurat)

library(qs)

n2_count = read.csv(file='/datg/wangmd_data/RE/ren_n_v3_raw_counts_matrix.csv')


rownames(n2_count) = n2_count$X


n2_count = n2_count[,-1]

n2_count = t(n2_count)


n2 = CreateSeuratObject(n2_count,min.cells = 0,min.features = 0,project = 'n2')

n2 = NormalizeData(n2)
n2 = FindVariableFeatures(n2)
n2 = ScaleData(n2)
n2 = RunPCA(n2)
n2 = RunUMAP(n2,dims = 1:15)

n2 = FindNeighbors(n2,dims = 1:15)

n2 = FindClusters(n2,resolution = 1)

DimPlot(n2)

FeaturePlot(n2,features = 'Chrna4')

FeaturePlot(n2,features = 'Trpc3')

FeaturePlot(n2,features = 'Ficd')



o2_count = read.csv(file='/datg/wangmd_data/RE/ren_o_v3_raw_counts_matrix.csv')

o2_count[1:5,1:5]

rownames(o2_count) = o2_count$X

o2_count = o2_count[,-1]

o2_count = t(o2_count)

o2_count[1:5,1:5]

o2 = CreateSeuratObject(o2_count,min.cells = 0,min.features = 0,project = 'o2')

o2 = NormalizeData(o2)
o2 = FindVariableFeatures(o2)
o2 = ScaleData(o2)
o2 = RunPCA(o2)
o2 = RunUMAP(o2,dims = 1:15)

o2 = FindNeighbors(o2,dims = 1:15)

o2 = FindClusters(o2,resolution = 1)

DimPlot(o2)


coord_n = read.csv(file='/datg/wangmd_data/RE/spatial_coord_n.csv')

rownames(coord_n) = coord_n$X

coord_n = coord_n[,-1]

head(coord_n)

all(rownames(coord_n)==colnames(n2))

colnames(coord_n) = c('spatial_1','spatial_2')

head(coord_n)

n2[['spatial']] <- CreateDimReducObject(embeddings = as.matrix(coord_n), key = 'spatial_', assay = 'RNA')

DimPlot(n2,reduction = 'spatial')

coord = as.data.frame(n2@reductions$spatial@cell.embeddings)

coord$type = n2$`RNA_snn_res.1`

head(coord)

write.csv(coord,file='/datg/wangmd_data/RE/n2_coord_for_tissumap.csv')



coord_o = read.csv(file='/datg/wangmd_data/RE/spatial_coord_o.csv')

rownames(coord_o) = coord_o$X

coord_o = coord_o[,-1]

head(coord_o)

all(rownames(coord_o)==colnames(o2))

colnames(coord_o) = c('spatial_1','spatial_2')

head(coord_o)

o2[['spatial']] <- CreateDimReducObject(embeddings = as.matrix(coord_o), key = 'spatial_', assay = 'RNA')

DimPlot(o2,reduction = 'spatial')

coord = as.data.frame(o2@reductions$spatial@cell.embeddings)

coord$type = o2$`RNA_snn_res.1`

head(coord)

write.csv(coord,file='/datg/wangmd_data/RE/o2_coord_for_tissumap.csv')



re_10x = readRDS(file='/datg/wangmd_data/RE/ReN_with_recallB_clustering.rds')

DimPlot(re_10x)

n2_tiss = read.csv(file='/datg/wangmd_data/RE/n2_tissueumap.csv')

head(n2_tiss)

library(ggplot2)

dim(n2_tiss)

n2_f = subset(n2,cells = n2_tiss$X)

n2_f

meta = read.csv('/datg/wangmd_data/RE/meta_n.csv')

head(meta)

pos = match(colnames(n2_f),meta$X)
grep('TRUE',is.na(pos))

n2_f$region = meta$region[pos]

table(n2_tiss$regionName)

pos = which(colnames(n2_f)%in%n2_tiss$X[n2_tiss$regionName=='region18'])

length(pos)

n2_f$region[pos]='region11'

DimPlot(n2_f,group.by = 'region',reduction = 'spatial',label=T)

n2_f$slice_seq = n2_f$region

n2_f$slice_seq = gsub('region1$','1',n2_f$slice_seq)
n2_f$slice_seq = gsub('region9','2',n2_f$slice_seq)
n2_f$slice_seq = gsub('region10','3',n2_f$slice_seq)
n2_f$slice_seq = gsub('region7','4',n2_f$slice_seq)
n2_f$slice_seq = gsub('region8','5',n2_f$slice_seq)
n2_f$slice_seq = gsub('region5','6',n2_f$slice_seq)
n2_f$slice_seq = gsub('region4','7',n2_f$slice_seq)
n2_f$slice_seq = gsub('region6','8',n2_f$slice_seq)
n2_f$slice_seq = gsub('region11','9',n2_f$slice_seq)

DimPlot(n2_f,group.by = 'slice_seq',reduction = 'spatial')

o2_tiss = read.csv(file='/datg/wangmd_data/RE/o2_tissueumap.csv')

dim(o2_tiss)

o2_f = subset(o2,cells = o2_tiss$X)

meta = read.csv(file='/datg/wangmd_data/RE/meta_o.csv')

pos = match(colnames(o2_f),meta$X)
grep('TRUE',is.na(pos))

o2_f$region = meta$region[pos]

DimPlot(o2_f,reduction = 'spatial',group.by = 'region',label=T)

o2_f$slice_seq = o2_f$region

o2_f$slice_seq = gsub('region13$','1',o2_f$slice_seq)
o2_f$slice_seq = gsub('region14','2',o2_f$slice_seq)
o2_f$slice_seq = gsub('region11','3',o2_f$slice_seq)
o2_f$slice_seq = gsub('region12','4',o2_f$slice_seq)
o2_f$slice_seq = gsub('region20','5',o2_f$slice_seq)
o2_f$slice_seq = gsub('region16','6',o2_f$slice_seq)
o2_f$slice_seq = gsub('region15','7',o2_f$slice_seq)
o2_f$slice_seq = gsub('region18','8',o2_f$slice_seq)
o2_f$slice_seq = gsub('region17','9',o2_f$slice_seq)

DimPlot(o2_f,reduction = 'spatial',group.by = 'slice_seq',label=T)

save(o2_f,n2_f,file='/datg/wangmd_data/RE/o2_n2_clean_RE_neuron.Robj')

load(file='/datg/wangmd_data/RE/o2_n2_clean_RE_neuron.Robj')

o2_f$batch = 'O'

n2_f$batch = 'N'

#spatial clustering
sp_merge = merge(o2_f,n2_f)

sp_merge = JoinLayers(sp_merge,assay = 'RNA')


sp_merge = NormalizeData(sp_merge)

table(sp_merge$batch)

sp_merge = FindVariableFeatures(sp_merge)
sp_merge = ScaleData(sp_merge)
sp_merge = RunPCA(sp_merge)
sp_merge = RunUMAP(sp_merge,dims = 1:15)

sp_merge[["RNA"]] <- split(sp_merge[["RNA"]], f = sp_merge$batch)

sp_merge <- FindVariableFeatures(sp_merge)

sp_merge <- ScaleData(sp_merge)
sp_merge <- RunPCA(sp_merge)


options(future.globals.maxSize=3e+12)

o2_f = FindVariableFeatures(o2_f)
n2_f = FindVariableFeatures(n2_f)


var_gene = intersect(VariableFeatures(n2_f),VariableFeatures(o2_f))
length(var_gene)


sp_merge <- IntegrateLayers(
  object = sp_merge, method = RPCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.rpca",
  verbose = FALSE,,k.anchor=80
)

sp_merge <- RunUMAP(sp_merge, reduction = "integrated.rpca", dims = 1:30, reduction.name = "umap.rpca")

sp_merge = FindNeighbors(sp_merge,dims = 1:30,reduction = 'integrated.rpca')

sp_merge = RunUMAP(sp_merge,dims = 1:30,reduction = 'integrated.rpca')

sp_merge <- RunUMAP(sp_merge, reduction = "integrated.rpca", dims = 1:30, reduction.name = "umap.rpca")


DimPlot(sp_merge,group.by = 'batch',reduction = 'umap.rpca')

FeaturePlot(sp_merge,features = 'Chrna4')

DimPlot(sp_merge,group.by  = 'proj',reduction = 'umap.rpca')

Idents(sp_merge)='proj'

sp_merge_dual = subset(sp_merge,idents = 'dual')

Idents(sp_merge_dual)='batch'

table(Idents(sp_merge_dual))

sp_merge_dual = JoinLayers(sp_merge_dual,assay = 'RNA')

dual_deg = FindAllMarkers(sp_merge_dual,only.pos = T)

table(dual_deg$cluster)

write.csv(dual_deg,file='/datg/wangmd_data/RE/N_O_dual_deg.csv')

saveRDS(sp_merge,file='/datg/wangmd_data/RE/N_O_RE_sp_integration.rds')

library(qs)

sp_merge = readRDS(file='/datg/wangmd_data/RE/N_O_RE_sp_integration.rds')


pos = grep('^Arc$',rownames(o2_f))
quantile(o2_f@assays$RNA@layers$data[pos,],0.9)
pos = grep('^Bdnf$',rownames(o2_f))
quantile(o2_f@assays$RNA@layers$data[pos,],0.9)
pos = grep('^Btg2$',rownames(o2_f))

quantile(o2_f@assays$RNA@layers$data[pos,],0.9)
pos = grep('^Fos$',rownames(o2_f))

quantile(o2_f@assays$RNA@layers$data[pos,],0.9)
pos = grep('^Fosl2$',rownames(o2_f))

quantile(o2_f@assays$RNA@layers$data[pos,],0.9)
pos = grep('^Homer1$',rownames(o2_f))

quantile(o2_f@assays$RNA@layers$data[pos,],0.9)
pos = grep('^Npas4$',rownames(o2_f))

quantile(o2_f@assays$RNA@layers$data[pos,],0.9)
pos = grep('^Nr4a1$',rownames(o2_f))

quantile(o2_f@assays$RNA@layers$data[pos,],0.9)

pos = grep('^Arc$',rownames(o2_f))
cell1 = pos = grep('TRUE',o2_f@assays$RNA@layers$data[pos,]>0)
length(cell1)
pos = grep('^Bdnf$',rownames(o2_f))
cell2 = pos = grep('TRUE',o2_f@assays$RNA@layers$data[pos,]>0)
length(cell2)
pos = grep('^Btg2$',rownames(o2_f))
cell3 = pos = grep('TRUE',o2_f@assays$RNA@layers$data[pos,]>0)
length(cell3)
pos = grep('^Fos$',rownames(o2_f))
cell4 = pos = grep('TRUE',o2_f@assays$RNA@layers$data[pos,]>1.33066991220973)
length(cell4)
pos = grep('^Fosl2$',rownames(o2_f))
cell5 = pos = grep('TRUE',o2_f@assays$RNA@layers$data[pos,]>0)
length(cell5)
pos = grep('^Homer1$',rownames(o2_f))
cell6 = pos = grep('TRUE',o2_f@assays$RNA@layers$data[pos,]>1.55409167281197)
length(cell6)
pos = grep('^Npas4$',rownames(o2_f))
cell7 = pos = grep('TRUE',o2_f@assays$RNA@layers$data[pos,]>0)
length(cell7)
pos = grep('^Nr4a1$',rownames(o2_f))
cell8 = pos = grep('TRUE',o2_f@assays$RNA@layers$data[pos,]> 1.30352676010436)
length(cell8)

cell_use = c(cell1,cell2,cell3,cell4,cell5,cell6,cell7,cell8)

cell_use = unique(cell_use)
length(cell_use)

o2_f$IEG_activate = 'no'
o2_f$IEG_activate[cell_use] ='yes'

table(o2_f$IEG_activate)

table(o2_f$proj[o2_f$IEG_activate=='yes'])/table(o2_f$proj)



pos = grep('^Arc$',rownames(n2_f))
quantile(n2_f@assays$RNA@layers$data[pos,],0.9)
pos = grep('^Bdnf$',rownames(n2_f))
quantile(n2_f@assays$RNA@layers$data[pos,],0.9)
pos = grep('^Btg2$',rownames(n2_f))

quantile(n2_f@assays$RNA@layers$data[pos,],0.9)
pos = grep('^Fos$',rownames(n2_f))

quantile(n2_f@assays$RNA@layers$data[pos,],0.9)
pos = grep('^Fosl2$',rownames(n2_f))

quantile(n2_f@assays$RNA@layers$data[pos,],0.9)
pos = grep('^Homer1$',rownames(n2_f))

quantile(n2_f@assays$RNA@layers$data[pos,],0.9)
pos = grep('^Npas4$',rownames(n2_f))

quantile(n2_f@assays$RNA@layers$data[pos,],0.9)
pos = grep('^Nr4a1$',rownames(n2_f))

quantile(n2_f@assays$RNA@layers$data[pos,],0.9)



pos = grep('^Arc$',rownames(n2_f))
cell1 = pos = grep('TRUE',n2_f@assays$RNA@layers$data[pos,]>0)
length(cell1)
pos = grep('^Bdnf$',rownames(n2_f))
cell2 = pos = grep('TRUE',n2_f@assays$RNA@layers$data[pos,]>0)
length(cell2)
pos = grep('^Btg2$',rownames(n2_f))
cell3 = pos = grep('TRUE',n2_f@assays$RNA@layers$data[pos,]>0)
length(cell3)
pos = grep('^Fos$',rownames(n2_f))
cell4 = pos = grep('TRUE',n2_f@assays$RNA@layers$data[pos,]>1.13099344823044)
length(cell4)
pos = grep('^Fosl2$',rownames(n2_f))
cell5 = pos = grep('TRUE',n2_f@assays$RNA@layers$data[pos,]>0)
length(cell5)
pos = grep('^Homer1$',rownames(n2_f))
cell6 = pos = grep('TRUE',n2_f@assays$RNA@layers$data[pos,]>1.48586916166111)
length(cell6)
pos = grep('^Npas4$',rownames(n2_f))
cell7 = pos = grep('TRUE',n2_f@assays$RNA@layers$data[pos,]>0)
length(cell7)
pos = grep('^Nr4a1$',rownames(n2_f))
cell8 = pos = grep('TRUE',n2_f@assays$RNA@layers$data[pos,]> 0)
length(cell8)

cell_use = c(cell1,cell2,cell3,cell4,cell5,cell6,cell7,cell8)

cell_use = unique(cell_use)
length(cell_use)

n2_f$IEG_activate = 'no'
n2_f$IEG_activate[cell_use] ='yes'

table(n2_f$IEG_activate)

table(n2_f$proj[n2_f$IEG_activate=='yes'])/table(n2_f$proj)

table(o2_f$proj[o2_f$IEG_activate=='yes'])/table(o2_f$proj)

save(o2_f,n2_f,file='/datg/wangmd_data/RE/o2_n2_clean_RE_neuron_with_IEG.Robj')

load(file='/datg/wangmd_data/RE/o2_n2_clean_RE_neuron_with_IEG.Robj')

table(ReN1$enrich_type)

on = merge(n2_f,o2_f)

on = JoinLayers(on,assay = 'RNA')

on = NormalizeData(on)

table(on$batch)

Idents(on)='batch'

on_deg = FindAllMarkers(on,only.pos = T,logfc.threshold = 0.01)

table(on_deg$cluster)s

Idents(on)='proj'

dual_on = subset(on,idents = 'dual')

Idents(dual_on)='batch'

on_deg = FindAllMarkers(dual_on,only.pos = T)

on_deg[on_deg$cluster=='O',][grep('Chr',on_deg[on_deg$cluster=='O',]$gene),]

re_10x$batch ='ReN'

re_10x

o2_f

n2_f

all_merge = merge(o2_f,n2_f)

all_merge = merge(all_merge,re_10x)

all_merge

table(all_merge$batch)

all_merge = JoinLayers(all_merge,assay = 'RNA')

all_merge = NormalizeData(all_merge)

pos = grep('ReN',all_merge$batch)
all_merge$batch[pos] = all_merge$orig.ident[pos]

table(re_10x$orig.ident)

table(all_merge$batch)

all_merge = FindVariableFeatures(all_merge)
all_merge = ScaleData(all_merge)
all_merge = RunPCA(all_merge)
all_merge = RunUMAP(all_merge,dims = 1:15)

all_merge[["RNA"]] <- split(all_merge[["RNA"]], f = all_merge$batch)

all_merge <- FindVariableFeatures(all_merge)

all_merge <- ScaleData(all_merge)
all_merge <- RunPCA(all_merge)

options(future.globals.maxSize=3e+12)

re_10x = FindVariableFeatures(re_10x)

o2_f = FindVariableFeatures(o2_f)
n2_f = FindVariableFeatures(n2_f)

var_gene = intersect(VariableFeatures(re_10x),VariableFeatures(o2_f))
var_gene = intersect(var_gene,VariableFeatures(n2_f))
length(var_gene)

inter_g = (intersect(VariableFeatures(re_10x),rownames(n2_f)))

all_merge <- IntegrateLayers(
  object = all_merge, method = RPCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.rpca",
  verbose = FALSE,features= inter_g,k.anchor=80
)

all_merge <- RunUMAP(all_merge, reduction = "integrated.rpca", dims = 1:30, reduction.name = "umap.rpca")

all_merge = FindNeighbors(all_merge,dims = 1:30,reduction = 'integrated.rpca')

all_merge = FindClusters(all_merge,resolution = 1)

all_merge = FindClusters(all_merge,resolution = 0.5)

all_merge = FindClusters(all_merge,resolution = 2)

DimPlot(all_merge,group.by='batch',reduction='umap.rpca',order='O')

DimPlot(all_merge,label=T,reduction='umap.rpca',group.by='RNA_snn_res.2')+NoLegend()

FeaturePlot(all_merge,features = 'Chrna4',order=T,pt.size = 1,reduction = 'umap.rpca')

DimPlot(all_merge,group.by='enrich_type',reduction='umap.rpca',order='O')

saveRDS(all_merge,file='/datg/wangmd_data/RE/sp_10x_integration.rds')



all_merge$proj = 'other'

pos = grep('^11$|^12$|^16$',all_merge$RNA_snn_res.0.5)
all_merge$proj[pos]='mPFC1'

pos = grep('^8$|^14$|^15$|^9$',all_merge$RNA_snn_res.0.5)
all_merge$proj[pos]='vHPC1'

pos = grep('^1$|^4$|^6$|^5$|^7$|^13$',all_merge$RNA_snn_res.0.5)
all_merge$proj[pos]='dual'

pos = grep('^32$',all_merge$RNA_snn_res.2)
all_merge$proj[pos]='vHPC1'

pos = grep('^33$|^10$|^18$',all_merge$RNA_snn_res.2)
all_merge$proj[pos]='dual'

pos = grep('^27$',all_merge$RNA_snn_res.2)
all_merge$proj[pos]='other'

DimPlot(all_merge,group.by = 'proj',reduction = 'umap.rpca')

table(all_merge$batch)

Idents(all_merge)='batch'

n_match = subset(all_merge,idents = 'N')

o_match = subset(all_merge,idents = 'O')

DimPlot(n_match,reduction = 'umap.rpca',group.by = 'proj')

DimPlot(o_match,reduction = 'umap.rpca',group.by = 'proj')



pos = grep('mPFC1',coord$proj)
coord = rbind(coord[-pos,],coord[pos,])

DimPlot(all_merge,group.by = 'RNA_snn_res.0.5',label=T,reduction = 'integrated.rpca')

pdf('/datg/wangmd_data/RE/revise_figure/3_proj_pca.pdf',height = 5,width = 8)
ggplot(coord,aes(x=x,y=y,color=proj))+geom_point(size=0.6)+theme_bw()+theme(panel.grid = element_blank())+
scale_color_manual(values=c('dual1'='#B6D8C8','mPFC1'='#40DECE','vPFC1'='#0065A1'))
dev.off()

save(all_merge,file='/datg/wangmd_data/RE/revise_figure/all_data_with_sp_IEG.rds')

count = all_merge@assays$RNA@layers$counts

count[1:5,1:5]

rownames(count) = rownames(all_merge)
colnames(count) = colnames(all_merge)
count[1:5,1:5]

coord = data.frame(x = all_merge@reductions$pca@cell.embeddings[,3],y = all_merge@reductions$pca@cell.embeddings[,4])
head(coord)

list_cluster = all_merge$state
table(list_cluster)

all_merge = FindNeighbors(all_merge,dims = c(2:10))
all_merge = FindClusters(all_merge,resolution = 2)

DimPlot(all_merge,group.by = 'RNA_snn_res.2',label=T,reduction = 'pca',dims = c(3,4))

New_matrix = count

length(deg1)

New_matrix <-New_matrix[deg1, ]

dim(New_matrix)

expression_matrix <- New_matrix

list_cluster <- all_merge@meta.data$state


names(list_cluster) <- colnames(all_merge)


colnames(coord) = c('pca_1','pca_2')


umap_emb = coord

gene_loading <- all_merge@reductions$pca@feature.loadings


head(gene_loading)

meta = all_merge@meta.data

gene_annotation <- as.data.frame(rownames(all_merge), row.names = rownames(all_merge))
head(gene_annotation)

cell_metadata = as.data.frame(colnames(all_merge))
colnames(cell_metadata) <- "barcode"

rownames(cell_metadata) = colnames(all_merge)

head(cell_metadata)

expression_matrix = as.matrix(expression_matrix)

dim(expression_matrix)

write.table(expression_matrix,file='/datg/wangmd_data/RE/revise_figure/all_data_expr.dat',sep = '\t')



save(meta,gene_loading,gene_annotation,cell_metadata,expression_matrix,list_cluster,umap_emb,file='/datg/wangmd_data/RE/revise_figure/all_data_IEG_activate_progress_sub_for_monocle3_dual.Robj')



DimPlot(all_merge,group.by='batch1',reduction='integrated.rpca;l')

Idents(all_merge)='batch'

table(all_merge$batch)

ReN1 

ReN1 = ScaleData(ReN1,features = g_use2$x)

ReN1 = RunPCA(ReN1,features = g_use2$x)
length(gene_use)

ReN1 = RunUMAP(ReN1,dims = 1:15)

DimPlot(ReN1,reduction = 'pca')

DimPlot(ReN1)





DimPlot(all_merge,reduction = 'integrated.rpca',group.by = 'batch1',dims=c(1,2))



pos = grep('FALSE',is.na(all_merge$enrich_type))
all_merge$proj[pos]=all_merge$enrich_type[pos]

all_merge <- RunUMAP(all_merge, reduction = "integrated.rpca", dims = 1:30,n.components = 3, reduction.name = "umap.rpca")

#all_merge = FindNeighbors(all_merge,dims = 1:30,reduction = 'integrated.rpca')

#all_merge = FindClusters(all_merge,resolution = 1)


Idents(all_merge)='batch'

all_merge$batch1 = all_merge$batch

all_merge$batch1 = gsub('N','recall_b',all_merge$batch1)

all_merge$batch1 = gsub('recall_b1|recall_b2','recall_b',all_merge$batch1)

?RunUMAP

all_merge <- RunUMAP(all_merge, reduction = "integrated.rpca", dims = 1:30, reduction.name = "umap.rpca",
                n.components=5)


DimPlot(all_merge,group.by = 'batch1',reduction = 'umap.rpca',dims=c(4,5))

Idents(all_merge)='batch1'

all_merge = JoinLayers(all_merge,assay = 'RNA')

batch_deg = FindAllMarkers(all_merge,only.pos = T)

library(dplyr)

library(reshape2)

batch_deg <-batch_deg%>%group_by(cluster)%>%arrange(desc(avg_log2FC),.by_group = T)

top100 <- batch_deg%>%group_by(cluster)%>%top_n(50,avg_log2FC)


table(top100$cluster)

length(unique(top100$gene))

?IntegrateLayers

g_use = (VariableFeatures(all_merge))

g_use = c(g_use,top100[top100$cluster=='O',]$gene)
length(g_use)

g_use = c(g_use,top100[top100$cluster=='recall_b',]$gene)
length(g_use)

g_use = unique(g_use)
length(g_use)

head(g_use)

all_merge



all_merge <- IntegrateLayers(
  object = all_merge, method = RPCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.rpca",
  verbose = FALSE,k.anchor=5,features=g_use2$x)

all_merge <- RunUMAP(all_merge, reduction = "integrated.rpca", dims = 1:30, reduction.name = "umap.rpca")

all_merge = FindNeighbors(all_merge,dims = 1:30,reduction = 'integrated.rpca')

all_merge = FindClusters(all_merge,resolution = 1)

DimPlot(all_merge,group.by = 'batch',reduction = 'umap.rpca')

all_merge <- SCTransform(all_merge, verbose = FALSE)

length(g_use)

table(top100$cluster)

length(unique(top100$gene))

g_use2 = read.csv(file='/datg/wangmd_data/RE/revise_figure/version1_g_use.csv')

pos = which(g_use2$x%in%rownames(all_merge))
length(pos)

all_merge <- RunPCA(all_merge,features=unique(g_use2$x))
all_merge <- RunUMAP(all_merge, dims = 1:10)

DimPlot(all_merge,group.by = 'batch',reduction = 'umap')



mean(ReN1$nFeature_RNA)

mean(n2_f$nFeature_RNA)

2596.0932783678/1558.10545905707

table(all_merge$batch)

all_merge$batch1 = all_merge$batch

all_merge$batch = gsub('recall_b1|recall_b2','recall_b',all_merge$batch)

all_merge$batch = gsub('N','recall_b',all_merge$batch)

all_merge$batch = gsub('N','recall_b',all_merge$batch)

pos = grep('FALSE',is.na(all_merge$enrich_type))
all_merge$proj[pos] = all_merge$enrich_type[pos]

all_merge$proj = gsub('^dual$','dual1',all_merge$proj)

library(ggplot2)



DimPlot(all_merge,group.by='proj',reduction='umap.rpca',dims = c(1,2))+scale_color_manual(values=c('dual1'='#B7D9C9','mPFC1'='#40DFCF','vHPC1'='#0065A2','vPFC1'='#0065A2'))

DimPlot(all_merge,group.by='batch',reduction='umap.rpca',dims = c(1,2))#+scale_color_manual(values=c('dual1'='#B7D9C9','mPFC1'='#40DFCF','vHPC1'='#0065A2','vPFC1'='#0065A2'))

Idents(all_merge)='batch'

all_merge = JoinLayers(all_merge,assay = 'RNA')

all_merge  = PrepSCTFindMarkers(all_merge)

deg_o = FindMarkers(all_merge,ident.1 = 'O',only.pos=T)

table(all_merge$batch)

all_deg = FindAllMarkers(all_merge,only.pos = T)

all_deg <-all_deg%>%group_by(cluster)%>%arrange(desc(avg_log2FC),.by_group = T)

top100 <- all_deg%>%group_by(cluster)%>%top_n(300,avg_log2FC)



table(top100$cluster)

table(top100$cluster)

dim(deg_o)

write.csv(deg_o,file='/datg/wangmd_data/RE/revise_figure/O_vs_all_other_deg.csv')

library(qs)

table(all_merge$IEG_activate)

qsave(all_merge,file='/datg/wangmd_data/RE/revise_figure/all_rna_sp_merged_with_and_whithout_ieg.qrds')

ls()

table(n2_f$IEG_activate)

table(n2_f$slice_seq)

table(o2_f$proj)

Idents(o2_f)='proj'

o2_f_dual = subset(o2_f,idents = 'dual')

Idents(n2_f)='proj1'

n2_f_dual = subset(n2_f,idents = 'dual')

n2_f_pfc = subset(n2_f,idents = 'mPFC1')

o2_f_pfc = subset(o2_f,idents = 'mPFC1')

n2_f_hpc = subset(n2_f,idents = 'vHPC1')

o2_f_hpc = subset(o2_f,idents = 'vHPC1')

Idents(o2_f)='proj1'

o2_f_dual = subset(o2_f,idents = 'dual')

table(n2_f_dual$slice_seq)

table(n2_f_dual$slice_seq[n2_f_dual$IEG_activate=='yes'])/table(n2_f_dual$slice_seq)

table(n2_f_dual$slice_seq[n2_f_dual$IEG_activate=='yes'])/table(n2_f_dual$slice_seq)

table(o2_f_dual$slice_seq[o2_f_dual$IEG_activate=='yes'])/table(o2_f_dual$slice_seq)

table(o2_f_pfc$slice_seq[o2_f_pfc$IEG_activate=='yes'])/table(o2_f_pfc$slice_seq)

table(n2_f_pfc$slice_seq[n2_f_pfc$IEG_activate=='yes'])/table(n2_f_pfc$slice_seq)

table(n2_f_hpc$slice_seq[n2_f_hpc$IEG_activate=='yes'])/table(n2_f_hpc$slice_seq)
table(o2_f_hpc$slice_seq[o2_f_hpc$IEG_activate=='yes'])/table(o2_f_hpc$slice_seq)

save(o2_f,n2_f,file='/datg/wangmd_data/RE/revise_figure/test_IEG_clustering.Robj')

table(o2_f_dual$slice_seq[o2_f_dual$IEG_activate=='yes'])
table(o2_f_dual$slice_seq)

Idents(n2_f)='slice_seq'

n1 = subset(n2_f,idents = '8')

Idents(o2_f)='slice_seq'

o3 = subset(o2_f,idents = '3')

table(o3$proj)

table(n1$slice_seq)

FeaturePlot(n1,features = 'Chrna4',order=T)

table(n1$`RNA_snn_res.0.5`)

n1 = FindClusters(n1,resolution = 2)

pos = grep('^1$|^2$',n1$RNA_snn_res.1)
n1$proj[pos]='other'
cell1 = colnames(n1)[pos]

pos = grep('^1$',n1$RNA_snn_res.1.5)
n1$proj[pos]='other'
cell3 = colnames(n1)[pos]

pos = grep('^1$',n1$RNA_snn_res.2)
n1$proj[pos]='other'
cell6 = colnames(n1)[pos]

pos = grep('^4$|^9$|^11$|^6$',n1$RNA_snn_res.2)
n1$proj[pos]='other'
cell8 = colnames(n1)[pos]

pos = grep('^0$',n1$RNA_snn_res.2)
n1$proj[pos]='other'
cell9 = colnames(n1)[pos]

pos = grep('^8$',o3$RNA_snn_res.2)
o3$proj[pos]='other'
cellon8 = colnames(o3)[pos]

pos = grep('^11$',o3$RNA_snn_res.2)
o3$proj[pos]='other'
cellon3 = colnames(o3)[pos]

pos = which(colnames(o2_f)%in%cello8)
table(o2_f$slice_seq[pos])

o3 = FindClusters(o3,resolution = 2)

FeaturePlot(o3,features = 'Chrna4',order=T)

DimPlot(o3,group.by = 'IEG_activate',label=T)

DimPlot(o3,group.by = 'proj',label=T)

DotPlot(o3,features = 'Chrna4',group.by = 'RNA_snn_res.2')

DimPlot(o3,group.by = 'RNA_snn_res.2',label=T)

DimPlot(n1,group.by = 'proj',label=T)

24/(24+54)

17/(17+29)

table(n1$IEG_activate[n1$proj=='dual'])

n2_f$proj1 = n2_f$proj

pos = which(colnames(n2_f)%in%c(cell1,cell3,cell6,cell8,cell9))

n2_f$proj1[pos]='other'

cello3

#o2_f$proj1 = o2_f$proj

pos = which(colnames(o2_f)%in%c(cellon3,cellon8))
length(pos)

o2_f$proj1[pos]='other'

table(o2_f$proj)
table(o2_f$proj1)

Idents(n2_f)='proj1'

n2_f_dual = subset(n2_f,idents = 'dual')

table(o3$IEG_activate[o3$proj=='dual'])

23/(23+37)

DimPlot(n1,group.by = 'proj')



FeaturePlot(o2_f,features = )

Idents(o2_f)='slice_seq'

table(Idents(o2_f))

s3 = subset(o2_f,idents = '3')

DimPlot(s3,group.by = 'IEG_activate')

DimPlot(s3,group.by = 'proj')

FeaturePlot(s3,features = 'Chrna4')





DimPlot(o2_f,group.by = 'IEG_activate')



Idents(n2_f)='slice_seq'

n1 = subset(n2_f,idents = '2')

FeaturePlot(n1,features = 'Trpc4',order=T)

n1 = FindClusters(n1,resolution = 2)

FeaturePlot(n1,features = 'Gpc3',order=T)

DimPlot(n1,group.by = 'proj1')
DimPlot(n1,group.by = 'IEG_activate')

DimPlot(n1,group.by = 'RNA_snn_res.2',label=T)

pos = grep('^0$|^1$|^11$|^6$',n1$RNA_snn_res.2)
pos1 = grep('^yes$',n1$IEG_activate[pos])
pos2 = grep('vHPC1',n1$proj1[pos[pos1]])
n1$proj1[pos[pos1[pos2]]]='other'
celln1= colnames(n1)[pos[pos1[pos2]]]

length(celln2)

pos = grep('^6$|^16$',n1$RNA_snn_res.2)
pos1 = grep('^yes$',n1$IEG_activate[pos])
pos2 = grep('vHPC1',n1$proj1[pos[pos1]])
n1$proj1[pos[pos1[pos2]]]='other'
celln2= colnames(n1)[pos[pos1[pos2]]]
length(celln2)

pos = grep('^4$|^7$|^5$|^11$',n1$RNA_snn_res.2)
pos1 = grep('^yes$',n1$IEG_activate[pos])
pos2 = grep('vHPC1',n1$proj1[pos[pos1]])
n1$proj1[pos[pos1[pos2]]]='other'
celln2= colnames(n1)[pos[pos1[pos2]]]

pos = grep('^1$|^3$|^14$|^15$|^11$',n1$RNA_snn_res.2)
pos1 = grep('^yes$',n1$IEG_activate[pos])
pos2 = grep('vHPC1',n1$proj1[pos[pos1]])
n1$proj1[pos[pos1[pos2]]]='other'
celln6= colnames(n1)[pos[pos1[pos2]]]

pos = grep('^13$|^15$|^11$|^9$',n1$RNA_snn_res.2)
pos1 = grep('^yes$',n1$IEG_activate[pos])
pos2 = grep('vHPC1',n1$proj1[pos[pos1]])
n1$proj1[pos[pos1[pos2]]]='other'
celln8= colnames(n1)[pos[pos1[pos2]]]

pos = grep('dual',n1$proj1)
n1$proj1[pos]='dual'

table(n1$IEG_activate[n1$proj1=='vHPC1'])

pos = which(colnames(n2_f)%in%c(celln1,celln2,celln6,celln8))
length(pos)

n2_f$proj1[pos]='other'

16/(35)

Idents(o2_f)='slice_seq'

o1 = subset(o2_f,idents = '1')

FeaturePlot(o1,features = 'Trpc4',order=T)

o1 = FindClusters(o1,resolution = 2)

DimPlot(o1,group.by = 'proj1')
DimPlot(o1,group.by = 'IEG_activate')

DimPlot(o1,group.by = 'RNA_snn_res.2',label=T)

pos = grep('^0$',o1$RNA_snn_res.2)
pos1 = grep('^no$',o1$IEG_activate[pos])
pos2 = grep('vHPC1',o1$proj1[pos[pos1]])
o1$proj1[pos[pos1[pos2]]]='other'
cello1= colnames(o1)[pos[pos1[pos2]]]

o1 = subset(o2_f,idents = '6')

FeaturePlot(o1,features = 'Trpc4',order=T)

o1 = FindClusters(o1,resolution = 2)

DimPlot(o1,group.by = 'proj1')
DimPlot(o1,group.by = 'IEG_activate')

DimPlot(o1,group.by = 'RNA_snn_res.2',label=T)

pos = grep('^3$',o1$RNA_snn_res.2)
pos1 = grep('^no$',o1$IEG_activate[pos])
pos2 = grep('vHPC1',o1$proj1[pos[pos1]])
o1$proj1[pos[pos1[pos2]]]='other'
cello6= colnames(o1)[pos[pos1[pos2]]]

FeaturePlot(o1,features = 'Gpc3',order=T)

o1 = subset(o2_f,idents = '9')

FeaturePlot(o1,features = 'Trpc4',order=T)

o1 = FindClusters(o1,resolution = 2)

DimPlot(o1,group.by = 'proj1')
DimPlot(o1,group.by = 'IEG_activate')



DimPlot(o1,group.by = 'RNA_snn_res.1',label=T)

DimPlot(o1,group.by = 'proj1')

pos = grep('^1$|^2$|^10$|^7$|^11$|^5$',o1$RNA_snn_res.2)
pos1 = grep('^no$',o1$IEG_activate[pos])
pos2 = grep('vHPC1',o1$proj1[pos[pos1]])
o1$proj1[pos[pos1[pos2]]]='other'
cello9= colnames(o1)[pos[pos1[pos2]]]

pos = grep('^6$',o1$RNA_snn_res.1)
pos1 = grep('^no$',o1$IEG_activate[pos])
pos2 = grep('vHPC1',o1$proj1[pos[pos1]])
o1$proj1[pos[pos1[pos2]]]='other'
cello9_1= colnames(o1)[pos[pos1[pos2]]]

o1 = subset(o2_f,idents = '7')

FeaturePlot(o1,features = 'Trpc4',order=T)

o1 = FindClusters(o1,resolution = 2)

DimPlot(o1,group.by = 'proj1')
DimPlot(o1,group.by = 'IEG_activate')

FeaturePlot(o1,features = 'Gpc3',order=T)

DimPlot(o1,group.by = 'RNA_snn_res.2',label=T)

pos = grep('^4$',o1$RNA_snn_res.2)
pos1 = grep('^no$',o1$IEG_activate[pos])
pos2 = grep('vHPC1',o1$proj1[pos[pos1]])
o1$proj1[pos[pos1[pos2]]]='other'
cello7= colnames(o1)[pos[pos1[pos2]]]

table(o1$IEG_activate[o1$proj1=='vHPC1'])

7/24

FeaturePlot(o1,features = 'Gpc3',order=T)

o1 = subset(o2_f,idents = '8')

FeaturePlot(o1,features = 'Trpc4',order=T)

o1 = FindClusters(o1,resolution = 2)

DimPlot(o1,group.by = 'proj1')
DimPlot(o1,group.by = 'IEG_activate')

DimPlot(o1,group.by = 'RNA_snn_res.1',label=T)

pos = grep('^6$',o1$RNA_snn_res.2)
pos1 = grep('^no$',o1$IEG_activate[pos])
pos2 = grep('vHPC1',o1$proj1[pos[pos1]])
o1$proj1[pos[pos1[pos2]]]='other'
cello8= colnames(o1)[pos[pos1[pos2]]]

#o2_f$proj1 = o2_f$proj
pos = which(colnames(o2_f)%in%c(cello1,cello6,cello7,cello8,cello9,cello9_1))
length(pos)

o2_f$proj1[pos]='other'

table(o1$IEG_activate[o1$proj1=='vHPC1'])

Idents(ReN1)='enrich_type'

table(ReN1$enrich_type)

ReN1$batch1 = ReN1$enrich_type

ReN1$batch1 = gsub('mPFC1|vPFC1','single',ReN1$batch1)
table(ReN1$batch1)

Idents(ReN1)='batch1'

deg_10x = FindMarkers(ReN1,ident.1 = 'dual1')

deg_10x = deg_10x[order(deg_10x$avg_log2FC,decreasing = T),]

pos = grep('TRUE',deg_10x$p_val<0.05)
length(pos)

deg_10x1 = deg_10x[pos,]

pos = which(rownames(deg_10x1)%in%g_name)
length(pos)

dim(deg_10x)

dim(deg_10x1)

write.csv(deg_10x1,file='/datg/wangmd_data/RE/spatial_dual_vs_other.csv')

pos

length(g_name)

pos = grep('TRUE',deg_10x$avg_log2FC>0)
pos1 = grep('TRUE',deg_10x$avg_log2FC<0)
length(pos)
length(pos1)

10442-3094

deg_10x1 = deg_10x[c(1:3094,7253:10442),]

dim(deg_10x)

deg_10x[pos,]

head(deg_10x)

g_name = c('Kcnip4','Vav3','Cadm2','Chrm3','Kcnc2','Cacna2d3','Kcnj3',
          'Kalrn','Gabra4','Grip1','Cntnap5c','Cacna1c','Cacnb4','Gabbr2',
          'Gabrb2','Oprm1','Grid1','Grin2b','Kcnd3','Kcnq5','Qrfpr','Npas2','Kcnq3','Serpini1',
          'Sgcz','Kcnt1','Tuba1a','Gria1','Lrrc4','Bc1','Calb2','Ly6h','Lrrc4c')

deg_10x1$logp = -log10(deg_10x1$p_val)

deg_10x1$group = 'single'
pos = grep('TRUE',deg_10x1$avg_log2FC>0)
deg_10x1$group[pos]='dual'

pos = which(rownames(deg_10x1)%in%g_name)
length(pos)

gene_label = deg_10x1[pos,]

gene_label$gene= rownames(gene_label)

library(ggrepel)

pdf('/datg/wangmd_data/RE/revise_figure/single_dual_IEG_deg_vol.pdf',height = 3.6,width = 4.2)
ggplot()+geom_point(data = deg_10x1,aes(x=avg_log2FC,y=logp,color=group))+
geom_text_repel(data =gene_label,aes(x = avg_log2FC,y = logp,label = gene), size = 3,
  box.padding = 0.5,
  point.padding = 0.3,
                max.overlaps=Inf,
  segment.color = "grey50",min.segment.length=0,
  segment.size =1)+scale_color_manual(values=c('single'='#006599','dual'='#AFEFCF'))+
theme_bw()+theme(panel.grid = element_blank())+geom_vline(xintercept = c(-0.25,0.25),lty='dashed')
dev.off()

ggplot(deg_10x1,aes(x=avg_log2FC,y= p_val))+geom_point()

DimPlot(all_merge,reduction = 'pca',group.by = 'RNA_snn_res.0.5',label=T,dims = c(3,4))

coord = data.frame(x= all_merge@reductions$pca@cell.embeddings[,3],
                  y = all_merge@reductions$pca@cell.embeddings[,4],
                  state = all_merge$state,
                  res = all_merge$RNA_snn_res.0.5)

pos = grep('^2$|^4$',coord$res)
coord1 = coord[pos,]
coord2 = coord[-pos,]

pos = grep('over_b',coord1$state)
length(pos)

pos_tmp = sample(length(pos),200)
length(pos_tmp)
coord1 = rbind(coord1[-pos[pos_tmp],],coord1[pos[pos_tmp],])

coord_all = rbind(coord2,coord1)
ggplot(coord_all,aes(x=x,y=y,color=state))+geom_point(size=0.5)+theme_bw()+theme(panel.grid = element_blank())

tiff('/datg/wangmd_data/RE/revise_figure/4_state_pca_dual.tif',height = 140,width = 160,units='mm',compression = 'lzw',res=600)
ggplot(coord_all,aes(x=x,y=y,color=state))+geom_point(size=1)+theme_bw()+theme(panel.grid = element_blank())+
scale_color_manual(values=c('training'='#61F2A2','recall_b'='#ADB2D8','recall_A'='#6534B2','over_b'='#AF5952'))
dev.off()

DimPlot(all_merge,reduction = 'pca',dims = c(3,4),group.by = 'state',order='over_b')

table(n2_f$IEG_activate)

table(n2_f$)

DimPlot(n2_f,group.by = 'proj1',reduction = 'spatial')

o2_f1 = o2_f
n2_f1 = n2_f

load(file='/datg/wangmd_data/RE/o2_n2_clean_RE_neuron.Robj')

ls()

DimPlot(o2_f,reduction = 'spatial')

table(o2_f1$proj1)

table(n2_f1$proj1)

table(all_merge$proj1)

table(all_merge$state)

Idents(all_merge)='state'

all_merge_b = subset(all_merge,idents = c('over_b','recall_b'))

no_deg = FindAllMarkers(all_merge_b,only.pos = T)

table(no_deg$cluster)

pos = grep('^Gm|^Rps|Rik$',no_deg$gene)
no_deg = no_deg[-pos,]

table(no_deg$cluster)

write.csv(no_deg,file='/datg/wangmd_data/RE/revise_figure/recall_over_b_deg_all_merged.csv')

DimPlot(ReN1)

ls()
