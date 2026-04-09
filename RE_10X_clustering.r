library(Seurat)

#read-in data
wt = Read10X('/gpfs2/wulab10/RE/wt_counts/')

wt= CreateSeuratObject(wt,min.cells = 3,min.features = 3)

rownames(wt)[grep('^mt',rownames(wt))]

#QC
wt[["percent.mt"]] <- PercentageFeatureSet(wt, pattern = "^mt")

wt <- subset(wt, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 1)

VlnPlot(wt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#preprocessing and normalization
wt = NormalizeData(wt)

wt = FindVariableFeatures(wt)

wt = ScaleData(wt)

#dimension reduction
wt = RunPCA(wt)
wt = RunUMAP(wt,dims=1:15)

#visualization
DimPlot(wt,label=T)

#clustering
wt = FindNeighbors(wt,dims = 1:15 )

wt = FindClusters(wt,resolution=1)

DimPlot(wt,label=T)


#read-in data
train = Read10X('/gpfs2/wulab10/RE/train_counts/')

dim(train)

train = CreateSeuratObject(train,min.cells = 3,min.features = 3,project = 'train')

#QC
train[["percent.mt"]] <- PercentageFeatureSet(train, pattern = "^mt")

VlnPlot(train, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


train <- subset(train, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 1)

#normalization and preprocessing
train = NormalizeData(train)

train = FindVariableFeatures(train)

train = ScaleData(train)

#dimension reduction
train = RunPCA(train)

train = RunUMAP(train,dims=1:15)

#visualization
DimPlot(train)

#clustering
train = FindNeighbors(train,dims = 1:15)

train = FindClusters(train,resolution=1)

DimPlot(train,label=T)


#read-in data
recall = Read10X('/gpfs2/wulab10/RE/recall_counts/')

dim(recall)

recall = CreateSeuratObject(recall,min.cells = 3,min.features = 3,project = 'recall')

#QC
recall[["percent.mt"]] <- PercentageFeatureSet(recall, pattern = "^mt")

#normalization and preprocessing
recall = NormalizeData(recall)

recall = FindVariableFeatures(recall)

recall = ScaleData(recall)

#dimension reduction
recall = RunPCA(recall)

recall = RunUMAP(recall,dims=1:15)

#visualization
DimPlot(recall)

#read-in data
recall_2 = Read10X('/gpfs2/wulab10/RE/recall_counts_2/')

dim(recall_2)

recall_2 = CreateSeuratObject(recall_2,min.cells = 3,min.features = 3,project = 'recall_2')

#QC
recall_2[["percent.mt"]] <- PercentageFeatureSet(recall_2, pattern = "^mt")

#normalization and preprocessing
recall_2 = NormalizeData(recall_2)

recall_2 = FindVariableFeatures(recall_2)

recall_2 = ScaleData(recall_2)

#dimension reduction
recall_2 = RunPCA(recall_2)

recall_2 = RunUMAP(recall_2,dims=1:15)

#clustering
recall_2 = FindNeighbors(recall_2,dims = 1:15)

recall_2 = FindClusters(recall_2,resolution = 1)


#merge all four datasets
wt$set = 'homecage'

train$set = 'training'

recall_2$set = 'recall'

Idents(wt) = 'DF.classifications_0.25_0.09_4279'

table(Idents(wt))

wt_f = subset(wt,idents='Singlet')

Idents(train) = 'DF.classifications_0.25_0.09_2914'

train_f = subset(train,idents='Singlet')

table(recall_2$DF.classifications_0.25_0.09_872)

Idents(recall_2) ='DF.classifications_0.25_0.09_872'

recall_2f = subset(recall_2,idents='Singlet')



all_merge = merge(wt_f,train_f)

all_merge = merge(all_merge,recall_2f)

all_merge = NormalizeData(all_merge)

all_merge = FindVariableFeatures(all_merge)

all_merge = ScaleData(all_merge)

all_merge = RunPCA(all_merge)

all_merge = RunUMAP(all_merge,dims=1:15)

DimPlot(all_merge)

all_merge = FindNeighbors(all_merge,dims=1:15)

all_merge = FindClusters(all_merge,resolution=1)


#celltype annotation
all_merge$major_type = as.character(all_merge$RNA_snn_res.1)

pos <- grep('^2$|^6$|^9$|^10$|^13$|^30$|^11$',all_merge$RNA_snn_res.1)
all_merge$major_type[pos] = 'RE_N'
pos <- grep('^1$|^4$|^21$|^29$',all_merge$RNA_snn_res.1)
all_merge$major_type[pos] = 'ZI_N'
pos <- grep('^12$|^5$|^18$|^34$|^20$|^19$|^8$|^26$',all_merge$RNA_snn_res.1)
all_merge$major_type[pos] = 'TH_N'

pos <- grep('^23$|^22$|^25$',all_merge$RNA_snn_res.1)
all_merge$major_type[pos] = 'Endothelial Cell'

pos <- grep('^35$|^37$|^16$',all_merge$RNA_snn_res.1)
all_merge$major_type[pos] = 'Microglia'

pos <- grep('^17$|^0$|^3$',all_merge$RNA_snn_res.1)
all_merge$major_type[pos] = 'Astrocyte'

pos <- grep('^15$|^36$|^28$|^31$',all_merge$RNA_snn_res.1)
all_merge$major_type[pos] = 'OPC'

pos <- grep('^28$',all_merge$RNA_snn_res.1)
all_merge$major_type[pos] = 'Oligodendrocyte'

pos <- grep('^27$|^7$',all_merge$RNA_snn_res.1)
all_merge$major_type[pos] = 'Oligodendrocyte'

pos <- grep('^24$',all_merge$RNA_snn_res.1)
all_merge$major_type[pos] = 'Immune cell'

pos <- grep('^14$',all_merge$RNA_snn_res.1)
all_merge$major_type[pos] = 'Macrophage'

pos <- grep('^33$',all_merge$RNA_snn_res.1)
all_merge$major_type[pos] = 'Ciliated cell'

pos <- grep('^32$',all_merge$RNA_snn_res.1)
all_merge$major_type[pos] = 'Epithelial cell'



tiff('/gpfs2/wulab10/RE/all_major_type_10x_umap.tif',height = 120,width = 150,units='mm',compression = 'lzw',res=600)
DimPlot(all_merge,group.by = 'major_type',label=F)+scale_color_manual(values=c('Astrocyte'='#EEDB6E','Ciliated cell'='#A9DE7D',
                                                                              'Endothelial Cell' = '#D02E29','Epithelial cell'='#5F6BC2',
                                                                              'Immune cell' = '#69BCE0','Macrophage'='#EFC1A8',
                                                                              'Microglia' = '#92E528','Oligodendrocyte'='#69BCE0',
                                                                              'OPC'='turquoise','RE_N'='#CA77DD','TH_N'='#F2ACD6','ZI_N'='brown'))
dev.off()

Idents(all_merge)<-'major_type'

DimPlot(all_merge)

major_deg = FindAllMarkers(all_merge,only.pos = T)

table(major_deg$cluster)

library(dplyr)

major_deg <- major_deg[order(major_deg$avg_log2FC,decreasing = T),]

head(major_deg)

major_deg <- arrange(major_deg,avg_log2FC,.by_group = T)

major_deg <- major_deg%>%group_by(cluster)%>%arrange(desc(avg_log2FC),.by_group = T)

top100 <-major_deg%>%group_by(cluster)%>%top_n(100,avg_log2FC)



tiff('/gpfs2/wulab10/RE/pdgfra.tif',height = 100,width = 110,units='mm',compression = 'lzw',res=600)
FeaturePlot(all_merge,features = 'Pdgfra')+scale_color_gradientn(colors = c('grey','red3'))
dev.off()

tiff('/gpfs2/wulab10/RE/nox4.tif',height = 100,width = 110,units='mm',compression = 'lzw',res=600)
FeaturePlot(all_merge,features = 'Nox4',order=T)+scale_color_gradientn(colors = c('grey','red3'))
dev.off()

tiff('/gpfs2/wulab10/RE/nova1.tif',height = 100,width = 110,units='mm',compression = 'lzw',res=600)
FeaturePlot(all_merge,features = 'Nova1',order=F)+scale_color_gradientn(colors = c('grey','red3'))
dev.off()

tiff('/gpfs2/wulab10/RE/prox1.tif',height = 100,width = 110,units='mm',compression = 'lzw',res=600)
FeaturePlot(all_merge,features = 'Prox1',order=F)+scale_color_gradientn(colors = c('grey','red3'))
dev.off()

tiff('/gpfs2/wulab10/RE/aqp4.tif',height = 100,width = 110,units='mm',compression = 'lzw',res=600)
FeaturePlot(all_merge,features = 'Aqp4',order=T,min.cutoff = 0.5)+scale_color_gradientn(colors = c('grey','red3'))
dev.off()

tiff('/gpfs2/wulab10/RE/cx3cr1.tif',height = 100,width = 110,units='mm',compression = 'lzw',res=600)
FeaturePlot(all_merge,features = 'Cx3cr1',order=T,min.cutoff = 0.5)+scale_color_gradientn(colors = c('grey','red3'))
dev.off()

tiff('/gpfs2/wulab10/RE/mog.tif',height = 100,width = 110,units='mm',compression = 'lzw',res=600)
FeaturePlot(all_merge,features = 'Mog',order=F,min.cutoff = 0.5)+scale_color_gradientn(colors = c('grey','red3'))
dev.off()

tiff('/gpfs2/wulab10/RE/fth1.tif',height = 100,width = 110,units='mm',compression = 'lzw',res=600)
FeaturePlot(all_merge,features = 'Fth1',order=T,min.cutoff = 2.5)+scale_color_gradientn(colors = c('grey','red3'))
dev.off()

tiff('/gpfs2/wulab10/RE/flt1.tif',height = 100,width = 110,units='mm',compression = 'lzw',res=600)
FeaturePlot(all_merge,features = 'Flt1',order=T,min.cutoff = 0)+scale_color_gradientn(colors = c('grey','red3'))
dev.off()

tiff('/gpfs2/wulab10/RE/ttr.tif',height = 100,width = 110,units='mm',compression = 'lzw',res=600)
FeaturePlot(all_merge,features = 'Ttr',order=F,min.cutoff = 1.5)+scale_color_gradientn(colors = c('grey','red3'))
dev.off()

tiff('/gpfs2/wulab10/RE/dnah12.tif',height = 100,width = 110,units='mm',compression = 'lzw',res=600)
FeaturePlot(all_merge,features = 'Dnah12',order=F,min.cutoff = 1.5)+scale_color_gradientn(colors = c('grey','red3'))
dev.off()

tiff('/gpfs2/wulab10/RE/rbfox3.tif',height = 100,width = 110,units='mm',compression = 'lzw',res=600)
FeaturePlot(all_merge,features = 'Rbfox3',order=T,min.cutoff = 1)+scale_color_gradientn(colors = c('grey','red3'))
dev.off()


save(all_merge,major_deg,file='/gpfs2/wulab10/RE/all_merged_major_type_deg.Robj')


tiff('/gpfs2/wulab10/RE/10x_nfeature_rna.tif',height = 160,width = 220,units='mm',compression = 'lzw',res=600)
VlnPlot(all_merge,features = 'nFeature_RNA',pt.size = 0)+geom_jitter(size=0.001,alpha=0.2)+scale_fill_manual(values=c('homecage'='#9ADCC5','training'='#DA3D39','recall'='#5680C0'))
dev.off()

Idents(all_merge)='set'

levels(all_merge) = as.factor(c('homecage','training','recall'))

tiff('/gpfs2/wulab10/RE/10x_ncount_rna.tif',height = 160,width = 220,units='mm',compression = 'lzw',res=600)
VlnPlot(all_merge,features = 'nCount_RNA',pt.size = 0)+geom_jitter(size=0.001,alpha=0.2)+scale_fill_manual(values=c('homecage'='#9ADCC5','training'='#DA3D39','recall'='#5680C0'))
dev.off()

table(major_deg$cluster)

#subset RE neurons
Idents(all_merge) <-'major_type'

ReN = subset(all_merge,idents='RE_N')

ReN = NormalizeData(ReN)

ReN = FindVariableFeatures(ReN)

ReN = ScaleData(ReN)

ReN = RunPCA(ReN)

ReN = RunUMAP(ReN,dims=1:15)

FeaturePlot(ReN,features = 'Nr4a1',order=T,pt.size = 1)

DimPlot(ReN,group.by = 'test',order='c1')

ReN = FindNeighbors(ReN,dims = 1:15)

ReN = FindClusters(ReN,resolution=1)

DimPlot(ReN,label=T)

FeaturePlot(ReN_wt,features = 'Gm15577',order=T)

#projection assignment
dual_gene1 = read.csv('/gpfs2/wulab10/RE/limma_dual_deg.csv')

vPFC_gene1 = read.csv('/gpfs2/wulab10/RE/limma_vPFC_deg.csv')
mPFC_gene1 = read.csv('/gpfs2/wulab10/RE/limma_mPFC_deg.csv')

dual_gene1 = dual_gene1[order(dual_gene1$logFC,decreasing = T),]

vPFC_gene1 = vPFC_gene1[order(vPFC_gene1$logFC,decreasing = T),]

mPFC_gene1 = mPFC_gene1[order(mPFC_gene1$logFC,decreasing = T),]

ReN_wt = AddModuleScore(ReN_wt,features = list(head(vPFC_gene1$X,100)),name='vPFC')

FeaturePlot(ReN_wt,features = 'dual1',min.cutoff = 0.17,order=T,pt.size = 1)+scale_color_gradientn(colors=c('gray','red2','brown'))

ReN = AddModuleScore(ReN,features = list(head(dual_gene1$X,100)),name='dual')

ReN = AddModuleScore(ReN,features = list(head(mPFC_gene1$X,100)),name='mPFC')

ReN = AddModuleScore(ReN,features = list(head(vPFC_gene1$X,100)),name='vPFC')

DimPlot(ReN,group.by = 'set')

ReN = FindClusters(ReN,resolution = 2.5)

DimPlot(ReN,group.by = 'RNA_snn_res.2.5',label=T)

p = DotPlot(ReN,features = c('vPFC1','mPFC1','dual1'),group.by = 'res_type')

enrich_list = p$data

enrich_list = enrich_list[,c(3,4,5)]

enrich_list = reshape(enrich_list, idvar = "id", timevar = "features.plot", direction = "wide")

for(i in 1:nrow(enrich_list)){
enrich_list$type[i] = names(which.max(enrich_list[i,c(2,3,4)]))
enrich_list$max_value[i] = enrich_list[i,c(2,3,4)][which.max(enrich_list[i,c(2,3,4)])][1,1]
    }

enrich_list$type = gsub('avg.exp.scaled.','',enrich_list$type)

min(enrich_list$max_value)

enrich_list[enrich_list$max_value<0.65,]$type = 'ns'

pos <- match(ReN$res_type,enrich_list$id)
grep('TRUE',is.na(pos))
length(pos)

ReN$enrich_type = 'other'


ReN$enrich_type = enrich_list$type[pos]

DimPlot(ReN,group.by = 'res_type1',label=T)

DimPlot(ReN,group.by = 'enrich_type',order='dual1')



pdf('/gpfs2/wulab10/RE/10x_ReN_subtype_deg_dotplot.pdf',height = 5,width = 12)
DotPlot(ReN,col.min = 0.5,scale.min = 5,features = c('Etv1','Lrtm1','Gpc5','Kcnab1','Lama1','Galntl6','Nkain3','Asic2',
                       'Hmcn1', 'Cntn6','Trpc3','Vwc2','Pdzrn3','Trpc4','Cemip',
                        'Lysmd4','P2ry14','Diaph2','Pex5l','Gnal',
                       'Pld5','Nfib','Erbb4', 'Angpt1','Cdh7','Serpini1',
                        'Nrbp2','Ly6h','Trhde','Pcp4',
                        'Fbn2','Vwc2l','Slit3','Gulp1','Grik1','Nwd2','Ldb2'))+theme(axis.text.x = element_text(angle=45,h=1))+
scale_color_gradientn(colors=c("gray80","blue"))
dev.off()

tiff('/gpfs2/wulab10/RE/10x_all_merge_ReN_subtype_umap.tif',height = 140,width = 150,units='mm',compression = 'lzw',res=600)
DimPlot(ReN,group.by = 'res_type',label=F,pt.size = 0.4)+scale_color_manual(values=c('0'='#5757ad','1'='#34661e','2'='#9af1ff','3'='#ff754c',
                                                               '4'='#ffe999',
                                                              '5'='#0e05f8','6'='#ff0172','7'='#cb885a','8'='#4bffb3','9'='#a82da0',
                                                              '10'='#ff73c9','11'='#b1acea','12'='#ffb473','13'='#5aad79','14'='#79b9ca',
                                                             '15'='#c0ff4c'))
dev.off()

Idents(ReN)<-'set'

table(Idents(ReN))

ReN_wt = subset(ReN,idents='homecage')

tiff('/gpfs2/wulab10/RE/10x_all_merge_ReN_homecage_subtype_umap.tif',height = 140,width = 150,units='mm',compression = 'lzw',res=600)
DimPlot(ReN_wt,group.by = 'res_type',label=F,pt.size = 0.4)+scale_color_manual(values=c('0'='#5757ad','1'='#34661e','2'='#9af1ff','3'='#ff754c',
                                                               '4'='#ffe999',
                                                              '5'='#0e05f8','6'='#ff0172','7'='#cb885a','8'='#4bffb3','9'='#a82da0',
                                                              '10'='#ff73c9','11'='#b1acea','12'='#ffb473','13'='#5aad79','14'='#79b9ca',
                                                             '15'='#c0ff4c'))
dev.off()

ReN_t = subset(ReN,idents='training')

tiff('/gpfs2/wulab10/RE/10x_all_merge_ReN_training_subtype_umap.tif',height = 140,width = 150,units='mm',compression = 'lzw',res=600)
DimPlot(ReN_t,group.by = 'res_type',label=F,pt.size = 0.4)+scale_color_manual(values=c('0'='#5757ad','1'='#34661e','2'='#9af1ff','3'='#ff754c',
                                                               '4'='#ffe999',
                                                              '5'='#0e05f8','6'='#ff0172','7'='#cb885a','8'='#4bffb3','9'='#a82da0',
                                                              '10'='#ff73c9','11'='#b1acea','12'='#ffb473','13'='#5aad79','14'='#79b9ca',
                                                             '15'='#c0ff4c'))
dev.off()

ReN_r = subset(ReN,idents='recall')

tiff('/gpfs2/wulab10/RE/10x_all_merge_ReN_recall_subtype_umap.tif',height = 140,width = 150,units='mm',compression = 'lzw',res=600)
DimPlot(ReN_r,group.by = 'res_type',label=F,pt.size = 0.4)+scale_color_manual(values=c('0'='#5757ad','1'='#34661e','2'='#9af1ff','3'='#ff754c',
                                                               '4'='#ffe999',
                                                              '5'='#0e05f8','6'='#ff0172','7'='#cb885a','8'='#4bffb3','9'='#a82da0',
                                                              '10'='#ff73c9','11'='#b1acea','12'='#ffb473','13'='#5aad79','14'='#79b9ca',
                                                             '15'='#c0ff4c'))
dev.off()

table(ReN$enrich_type)

tiff('/gpfs2/wulab10/RE/10x_all_merge_ReN_subtype_enrich_umap_240703.tif',height = 140,width = 150,units='mm',compression = 'lzw',res=600)
DimPlot(ReN,group.by = 'enrich_type',label=F,pt.size = 0.4)+scale_color_manual(values=c('dual1' = '#E88118','vPFC1' = '#92BFF4','mPFC1'='#0F01F7','ns'='gray'))
dev.off()

Idents(ReN) ='set'

table(ReN$set)

ReN_home = subset(ReN,idents = 'homecage')

ReN_recall = subset(ReN,idents = 'recall')

ReN_train = subset(ReN,idents = 'training')

tiff('/gpfs2/wulab10/RE/10x_all_merge_ReN_subtype_enrich_umap_home.tif',height = 140,width = 150,units='mm',compression = 'lzw',res=600)
DimPlot(ReN_home,group.by = 'enrich_type',label=F,pt.size = 0.4)+scale_color_manual(values=c('dual1'='#6F1D67','mPFC1'='#73B455','vPFC1'='#FCB96A','ns'='gray90'))
dev.off()

tiff('/gpfs2/wulab10/RE/10x_all_merge_ReN_subtype_enrich_umap_recall.tif',height = 140,width = 150,units='mm',compression = 'lzw',res=600)
DimPlot(ReN_recall,group.by = 'enrich_type',label=F,pt.size = 0.4)+scale_color_manual(values=c('dual1'='#6F1D67','mPFC1'='#73B455','vPFC1'='#FCB96A','ns'='gray90'))
dev.off()

tiff('/gpfs2/wulab10/RE/10x_all_merge_ReN_subtype_enrich_umap_train.tif',height = 140,width = 150,units='mm',compression = 'lzw',res=600)
DimPlot(ReN_train,group.by = 'enrich_type',label=F,pt.size = 0.4)+scale_color_manual(values=c('dual1'='#6F1D67','mPFC1'='#73B455','vPFC1'='#FCB96A','ns'='gray90'))
dev.off()



enrich_list_w = data.frame(id = enrich_list$id,mpfc = enrich_list$type,vhpc = enrich_list$type,dual = enrich_list$type)

enrich_list_w$mpfc = gsub('vPFC1|dual1','ns',enrich_list_w$mpfc)
enrich_list_w$vhpc = gsub('mPFC1|dual1','ns',enrich_list_w$vhpc)
enrich_list_w$dual = gsub('vPFC1|mPFC1','ns',enrich_list_w$dual)

table(enrich_list_w$mpfc)

head(enrich_list)

levels(enrich_list_w$id) = as.factor(c(seq(0,15)))

pdf('/gpfs2/wulab10/RE/ReN_subtype_enrich.pdf',height = 4,width = 6)
ggplot()+geom_point(data = enrich_list_w[,c(1,2)],aes(x=id,y=1,color=mpfc,shape=mpfc),size=3.5)+
geom_point(data = enrich_list_w[,c(1,3)],aes(x=id,y=1.05,color=vhpc,shape=vhpc),size=3.5)+
geom_point(data = enrich_list_w[,c(1,4)],aes(x=id,y=1.1,color=dual,shape=dual),size=3.5)+scale_y_continuous(breaks = c(1,1.05,1.1),limits = c(1,1.4))+
theme_bw()+scale_color_manual(values=c('ns'='gray80','dual1'='#4060b6','mPFC1'='#fce032',
                                                                                       'vPFC1'='#ef372f'))+scale_shape_manual(values=c('ns'=20,'dual1'=8,'mPFC1'=8,'vPFC1'=8))
dev.off()

#tiff('/gpfs2/wulab10/RE/10x_all_merge_ReN_subtype_trpc4.tif',height = 140,width = 150,units='mm',compression = 'lzw',res=600)
FeaturePlot(ReN,features = 'Gpc3',pt.size = 1,order=T,min.cutoff = 0)+scale_color_gradientn(colors=(c("gray80","#FB8861FF","#B63679FF","#51127CFF","#000004FF")))
#dev.off()

tiff('/gpfs2/wulab10/RE/10x_all_merge_ReN_subtype_trpc4.tif',height = 140,width = 150,units='mm',compression = 'lzw',res=600)
FeaturePlot(ReN,features = 'Trpc4',pt.size = 1,order=T,min.cutoff = 1)+scale_color_gradientn(colors=(c("gray80","#FB8861FF","#B63679FF","#51127CFF","#000004FF")))
dev.off()

#tiff('/gpfs2/wulab10/RE/10x_all_merge_ReN_subtype_chrna7.tif',height = 140,width = 150,units='mm',compression = 'lzw',res=600)
FeaturePlot(ReN,features = 'Nfix',pt.size = 1,order=T,min.cutoff = 0.2)+scale_color_gradientn(colors=(c("gray80","#FB8861FF","#B63679FF","#51127CFF","#000004FF")))
#dev.off()

coord = data.frame(type = ReN$enrich_type,value = ReN@assays$RNA@data['Gpc3',],x = ReN@reductions$umap@cell.embeddings[,1],
                  y = ReN@reductions$umap@cell.embeddings[,2])



tiff('/gpfs2/wulab10/RE/10x_all_merge_ReN_subtype_gpc3.tif',height = 140,width = 150,units='mm',compression = 'lzw',res=600)

ggplot(coord_all,aes(x=x,y=y,color=value))+geom_point(size=1)+scale_color_gradientn(colors=(c("gray80","#FB8861FF","#B63679FF","#51127CFF","#000004FF")))+
                                                                                theme_bw()+theme(panel.grid = element_blank())
dev.off()

coord = data.frame(type = ReN$enrich_type,value = ReN@assays$RNA@data['Gpc3',],x = ReN@reductions$umap@cell.embeddings[,1],
                  y = ReN@reductions$umap@cell.embeddings[,2])

tiff('/gpfs2/wulab10/RE/10x_all_merge_ReN_subtype_cdh13.tif',height = 140,width = 150,units='mm',compression = 'lzw',res=600)

ggplot(coord_all,aes(x=x,y=y,color=value))+geom_point(size=1)+scale_color_gradientn(colors=(c("gray80","#FB8861FF","#B63679FF","#51127CFF","#000004FF")))+
                                                                                theme_bw()+theme(panel.grid = element_blank())
dev.off()
coord = data.frame(type = ReN$enrich_type,value = ReN@assays$RNA@data['Chrna4',],x = ReN@reductions$umap@cell.embeddings[,1],
                  y = ReN@reductions$umap@cell.embeddings[,2])

tiff('/gpfs2/wulab10/RE/10x_all_merge_ReN_subtype_chrna4.tif',height = 140,width = 150,units='mm',compression = 'lzw',res=600)

ggplot(coord_all,aes(x=x,y=y,color=value))+geom_point(size=1)+scale_color_gradientn(colors=(c("gray80","#FB8861FF","#B63679FF","#51127CFF","#000004FF")))+
                                                                                theme_bw()+theme(panel.grid = element_blank())
dev.off()

head(enrich_list_w)




#projection enrichment output

Idents(ReN)='enrich_type'

ReN_hpc = subset(ReN,idents='vPFC1')

ReN_hpc

Idents(ReN_hpc)= 'set'

table(Idents(ReN_hpc))

ReN_hpc_h = subset(ReN_hpc,idents = 'homecage')

ReN_hpc_t = subset(ReN_hpc,idents = 'training')

ReN_hpc_r = subset(ReN_hpc,idents = 'recall')




pos = grep('recall',dual$set)
pos1 = grep('homecage',dual$set)

FeaturePlot(ReN,features = 'vPFC1',min.cutoff = 0)+scale_color_gradientn(colors = c('gray','red2'))

FeaturePlot(ReN,features = 'mPFC1',min.cutoff = 0,order=T)+scale_color_gradientn(colors = c('gray','red2'))

FeaturePlot(ReN,features = 'dual1',min.cutoff = 0.1)+scale_color_gradientn(colors = c('gray','red2'))

coord = data.frame(x=ReN@reductions$umap@cell.embeddings[,1],
                   y=ReN@reductions$umap@cell.embeddings[,2],
                  type = ReN$enrich_type,
                  dual = ReN$dual1,
                  mpfc = ReN$mPFC1,
                  vhpc = ReN$vPFC1)


library(scales)

tiff('/gpfs2/wulab10/RE/enrich_mpfc_umap.tif',height = 110,width = 120,units='mm',compression = 'lzw',res=600)
ggplot(coord_all,aes(x=x,y=y,color=mpfc))+geom_point(size=0.4)+scale_color_gradientn(colors=c('gray','red2'),limits=c(0.025,0.15),oob=scales::squish)+theme_void()
dev.off()

tiff('/gpfs2/wulab10/RE/enrich_vhpc_umap.tif',height = 110,width = 120,units='mm',compression = 'lzw',res=600)
ggplot(coord_all,aes(x=x,y=y,color=vhpc))+geom_point(size=0.4)+scale_color_gradientn(colors=c('gray','red2'),limits=c(0.025,0.25),oob=scales::squish)+theme_void()
dev.off()

tiff('/gpfs2/wulab10/RE/enrich_dual_umap.tif',height = 110,width = 120,units='mm',compression = 'lzw',res=600)
ggplot(coord_all,aes(x=x,y=y,color=dual))+geom_point(size=0.4)+scale_color_gradientn(colors=c('gray','red2'),limits=c(0.1,0.25),oob=scales::squish)+theme_void()
dev.off()

dim(condition_deg)

tiff('/gpfs2/wulab10/RE/gad1_expr.tif',height = 100,width = 130,units='mm',compression = 'lzw',res=600)
p =FeaturePlot(ReN,features = 'Gad1')
print(p)
dev.off()

ls()

save.image(file='/gpfs2/wulab10/RE/10x_all_image.rds')



