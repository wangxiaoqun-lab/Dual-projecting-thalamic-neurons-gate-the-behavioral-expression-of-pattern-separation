library(Seurat)

packageVersion('Matrix')

load('/gpfs2/wulab10/RE/clustering_all_data.Robj')


RE = ScaleData(RE,features =type_deg$gene)

RE = RunPCA(RE,features = pca_gene)

RE = RunUMAP(RE,dims=c(1:6),n.components = 2)

DimPlot(RE,pt.size = 3,dims = c(1,2),reduction = 'umap',group.by = 'proj_type')


library(limma)

DefaultAssay(RE)

RE1 = FindVariableFeatures(RE1,nfeatures = 22000)

exprs_mat = RE1@assays$RNA@data[VariableFeatures(RE1),]
dim(exprs_mat)

design <- model.matrix(~0+RE1$proj_type)

fit <- limma::lmFit(exprs_mat, design = design)


pos <- grep('TRUE',eBayes(fit)$coefficients[,1]>1)
length(pos)

min(eBayes(fit)$p.value[pos,1])

pos1 = (grep('TRUE',eBayes(fit)$p.value[pos,1]<=5e-2))

length(pos1)

dual_gene = rownames(eBayes(fit)$coefficients)[pos[pos1]]
length(dual_gene)

length(which(dual_gene%in%type_deg$gene))

pos <- grep('TRUE',eBayes(fit)$coefficients[,2]>1)
length(pos)

pos1 = (grep('TRUE',eBayes(fit)$p.value[pos,2]<=5e-2))

length(pos1)

mPFC_gene = rownames(eBayes(fit)$coefficients)[pos[pos1]]
length(mPFC_gene)

pos <- grep('TRUE',eBayes(fit)$coefficients[,3]>1)
length(pos)

pos1 = (grep('TRUE',eBayes(fit)$p.value[pos,3]<=5e-2))

vPFC_gene = rownames(eBayes(fit)$coefficients)[pos[pos1]]
length(vPFC_gene)

var_gene = unique(c(dual_gene,mPFC_gene,vPFC_gene))
length(var_gene)

p = which(rownames(eBayes(fit)$coefficients)%in%type_deg$gene)

grep('Chrna4|Ddhd1|Trpc3|Tanc1|Ptk2b|Kpna3|Tmem127|Cdh8|Ipo4|Lrp11',dual_gene)

dual_gene

VlnPlot(RE,features = 'Gm37180',group.by = 'proj_type')

coef= as.data.frame(eBayes(fit)$coefficients)

dim(coef)

coef_df = data.frame(gene='',type ='',value=0,idx=0)

for(i in 1:nrow(coef)){
    tmp = data.frame(gene='',type ='',value=0,idx=0)
   #if(sort((eBayes(fit)$coefficients)[i,],decreasing = T)[1]/(sort((eBayes(fit)$coefficients)[i,],decreasing = T))[2]>1.5){
   tmp$gene = rownames(eBayes(fit)$coefficients)[i]
    tmp$type = names(sort((eBayes(fit)$coefficients)[i,],decreasing = T)[1])
    tmp$value = sort((eBayes(fit)$coefficients)[i,],decreasing = T)[1]/(sort((eBayes(fit)$coefficients)[i,],decreasing = T))[2]
    tmp$idx = i
    #}
coef_df = rbind(coef_df,tmp)
}

for (i in 1:nrow(coef_df_f)){
    coef_df_f$pvalue[i] = eBayes(fit)$p.value[,coef_df_f$type[i]][coef_df_f$gene[i]]
    
}



pos <- grep('TRUE',coef_df_f$pvalue<5e-5)
length(pos)

coef_df_f1 = coef_df_f[pos,]



RE$res_type = as.character(RE$RNA_snn_res.2)

pos <- grep('^0$|^2$',RE$res_type)
RE$res_type[pos] = 'C1'

pos <- grep('^6$|^3$|^5$|^1$',RE$res_type)
RE$res_type[pos] = 'C2'

pos <- grep('^4$',RE$res_type)
RE$res_type[pos] = 'C3'

pos <- grep('^1$',RE$RNA_snn_res.3)
RE$res_type[pos] = 'C3'

pos <- grep('^22$|^18$|^7$|^19$',RE$RNA_snn_res.3)
RE$res_type[pos] = 'C3'

pos <- grep('^12$|^16$|^37$|^38$|^3$',RE$RNA_snn_res.4)
RE$res_type[pos] = 'C3'

RE1

DimPlot(RE,label=F,pt.size = 2,group.by = 'res_type')

table(RE$proj_type[RE$res_type=='C1'])

table(RE$proj_type[RE$res_type=='C2'])

table(RE$proj_type[RE$res_type=='C3'])



pdf('/gpfs2/wulab10/RE/RE_smart_seq_clustering_semi_supervised.pdf',height = 6,width = 6)
DimPlot(RE,label=F,pt.size = 3.5,group.by = 'res_type')
dev.off()

#pdf('/gpfs2/wulab10/RE/RE_smart_seq_clustering.pdf',height = 6,width = 6)
DimPlot(RE,label=F,pt.size = 3.5,group.by = 'res_type')
#dev.off()

#pdf('/gpfs2/wulab10/RE/RE_smart_seq_clustering.pdf',height = 6,width = 6)
DimPlot(RE,label=F,pt.size = 3.5,group.by = 'proj_type')
#dev.off()

DimPlot(RE,group.by = 'proj_type',pt.size = 2)

saveRDS(RE,file='/gpfs2/wulab10/RE/RE_smart_seq_clustering_v1.rds')


DimPlot(RE1,group.by = 'proj_type',pt.size = 2)

DimPlot(RE,group.by = 'res_type',pt.size = 2)

library(ggplot2)


coord = data.frame(x= RE1@reductions$umap@cell.embeddings[,1],y = RE1@reductions$umap@cell.embeddings[,2],type = RE1$proj_type)

pdf('/gpfs2/wulab10/RE/patch_proj_type_umap_250117.pdf',height = 5,width = 5.5)

ggplot(coord,aes(x=x,y=y))+geom_point(shape=21,color='black',aes(fill=type),size=7)+theme_void()+scale_fill_manual(values=c('dual' = '#8BC6BA','vPFC' = '#B2CCDE','mPFC'='royalblue'))
dev.off()

#pdf('/gpfs2/wulab10/RE/patch_proj_type_umap_240703.pdf',height = 5,width = 5.5)

ggplot(coord,aes(x=x,y=y))+geom_point(shape=21,color='black',aes(fill=type),size=7)+theme_void()+scale_fill_manual(values=c('dual' = '#E88118','vPFC' = '#92BFF4','mPFC'='#0F01F7'))
#dev.off()



#pdf('/gpfs2/wulab10/RE/patch_proj_type_umap.pdf',height = 5,width = 5.5)
DimPlot(RE,group.by = 'proj_type',pt.size = 3,order=c('dual','vPFC'))+scale_color_manual(values= c('dual'='#60308a','mPFC'='#63924d','vPFC'='#efd643'))
#dev.off()

pdf('/gpfs2/wulab10/RE/patch_res_type_umap.pdf',height = 5,width = 5.5)
DimPlot(RE,group.by = 'res_type',pt.size = 3)+scale_color_manual(values= c('C1'='#5b85ef','C3'='forestgreen','C2'='orange'))
dev.off()

density_plt = data.frame(cell = RE1@reductions$umap@cell.embeddings[,1],type = RE1$proj_type)

pdf('/gpfs2/wulab10/RE/density_plot_y_axis.pdf',height = 5,width = 5.5)
ggplot(density_plt,aes(cell,color=type,group=type,fill=type))+geom_density(alpha=0.5)+theme_bw()+theme(panel.grid = element_blank())+scale_color_manual(values= c('dual'='#C28018','mPFC'='#1000D2','vPFC'='#8099CF'))+scale_fill_manual(values= c('dual'='#C28018','mPFC'='#1000D2','vPFC'='#8099CF'))
dev.off()

pdf('/gpfs2/wulab10/RE/density_plot_x_axis.pdf',height = 5,width = 5.5)
ggplot(density_plt,aes(cell,color=type,group=type,fill=type))+geom_density(alpha=0.5)+theme_bw()+theme(panel.grid = element_blank())+scale_color_manual(values= c('dual'='#C28018','mPFC'='#1000D2','vPFC'='#8099CF'))+scale_fill_manual(values= c('dual'='#C28018','mPFC'='#1000D2','vPFC'='#8099CF'))
dev.off()

head(density_plt)

save(coef_df,coef_df_f,coef_df_f1,file='/gpfs2/wulab10/RE/coef_diff_gene.Robj')

RE = readRDS(file='/gpfs2/wulab10/RE/RE_smart_seq_clustering_v1.rds')

DimPlot(RE,group.by = 'proj_type',pt.size = 2)

FeaturePlot(RE,features = 'Prkcd',pt.size = 2,min.cutoff = 0.5)

Idents(RE)<-'proj_type'



library(DESeq2)

table(Idents(RE))

count = RE@assays$RNA@counts

dim(count)

meta = RE@meta.data

head(meta)

save(meta,count,file='/gpfs2/wulab10/RE/data_for_deseq2.Robj')

Idents(RE1)='proj_type'


library(edgeR)

count = as.matrix(RE1@assays$RNA@counts)


d0 = DGEList(counts = count)

d0 <- calcNormFactors(d0)


cutoff <- 1
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,] 
dim(d) # number of genes left


snames <- colnames(count) # Sample names

pos <- grep('^d|^m',snames)
length(pos)

snames[pos] = 'other'

pos <- grep('^v',snames)
length(pos)
group = rep('other',length(snames))
group[pos] = 'vHPC'


plotMDS(d, col = as.numeric(group))

mm <- model.matrix(~0 + group)

y <- voom(d, mm, plot = T)

tmp <- voom(d0, mm, plot = T)


fit <- lmFit(y, mm)
head(coef(fit))

contr <- makeContrasts(groupmPFC - groupother, levels = colnames(coef(fit)))


contr <- makeContrasts(groupvHPC - groupother, levels = colnames(coef(fit)))


contr <- makeContrasts(groupdual - groupother, levels = colnames(coef(fit)))


tmp <- contrasts.fit(fit, contr)

tmp <- eBayes(tmp)

dim(top.table[top.table$P.Value<0.05,])

fit = eBayes(fit)

top.table <- topTable(tmp, sort.by = "P", n = Inf)

top.table_mPFC = top.table[top.table$adj.P.Val<0.05,]

top.table_mPFC = top.table_mPFC[top.table_mPFC$logFC>0,]

top.table_mPFC


FeaturePlot(RE1,features = 'Trpc5',pt.size = 3,min.cutoff = 0)

dim(top.table_mPFC)

dim(top.table_mPFC)

pos <- grep('^Gm',rownames(top.table_mPFC))
length(pos)

mPFC_gene = top.table_mPFC[-pos,]

vHPC_gene = top.table_mPFC[-pos,]

dual_gene = top.table_mPFC[-pos,]

write.csv(dual_gene,file='/gpfs2/wulab10/RE/dual_limma_0713.csv')

write.csv(vHPC_gene,file='/gpfs2/wulab10/RE/vhpc_limma_0713.csv')

write.csv(mPFC_gene,file='/gpfs2/wulab10/RE/mPFC_limma_0713.csv')

dim(mPFC_gene)

rownames(vHPC_gene)

FeaturePlot(RE,pt.size = 3.5,features = 'Ficd',min.cutoff = 0.4,order=T)+scale_color_gradientn(colors=(c("gray","#FB8861FF","#B63679FF","#51127CFF","#000004FF")))


pdf('/gpfs2/wulab10/RE/patch_proj_type_umap_ficd.pdf',height = 5,width = 5.5)
FeaturePlot(RE,pt.size = 3.5,features = 'Ficd',min.cutoff = 0.5,order=T)+scale_color_gradientn(colors=(c("#8E8E8E","#FB8861FF","#B63679FF","#51127CFF","#000004FF")))

dev.off()



gene = c('Chrna4','Ddhd1','Trpc3','Cdh8','Tanc1','Ptk2b','Kpna3','Tmem127','Ipo4','Lrp11','Trpc4','Doc2b','Ficd','Snhg5','Gpc3','Mdm1','Zfp85','Hcn1','Sstr3','Hs3st4','Nfix','Ano3','Bcl11b','Etv3','Cnr1','Eml2','Tbrg4','Tamm41','Lmo3','Ppcs','Cdh13')

RE1 = ScaleData(RE1,features = gene)

RE = ScaleData(RE,features = c(dual_g,mpfc_g,vpfc_g))

avg_expr = AverageExpression(RE1,features =c(dual_g,mpfc_g,vpfc_g))

avg_expr = AverageExpression(RE1,features =gene)

Idents(RE1)='proj_type'

levels(RE1) = as.factor(c('dual','vPFC','mPFC'))

avg_expr = AverageExpression(RE1,features =gene)



#pdf('/gpfs2/wulab10/RE/patch_3_project_deg.pdf',height = 10,width = 5)
pheatmap(avg_expr$RNA[,],scale='row',cluster_rows = F,cluster_cols = F,color = colorRampPalette(colors=c('blue','white','white','red'))(100))
#dev.off()

pheatmap(avg_expr,scale='row',cluster_rows = F,cluster_cols = F,color = colorRampPalette(colors=c('white','steelblue4'))(100))


avg_expr = as.data.frame(avg_expr)



pdf('/gpfs2/wulab10/RE/patch_3_project_deg_1.pdf',height = 10,width = 5)
pheatmap(avg_expr[,],scale='row',cluster_rows = F,cluster_cols = F,color = colorRampPalette(colors=c('white','steelblue4'))(100))
dev.off()

avg_expr = avg_expr$RNA

saveRDS(avg_expr,file='/gpfs2/wulab10/RE/patch_3_project_deg.Rds')

gene = c('Chrna4','Ddhd1','Trpc3','Cdh8','Tanc1','Ptk2b','Kpna3','Tmem127',
        'Ipo4','Lrp11','Trpc4','Doc2b','Ficd','Snhg5','Gpc3','Mdm1','Zfp85',
        'Hcn1','Hs3st4','Nfix','Ano3','Bcl11b','Etv3','Cnr1','Eml2','Tbrg4',
        'Sstr3','Tamm41','Lmo3','Ppcs')

avg_expr = avg_expr$RNA

which(!gene%in%rownames(avg_expr))
gene[c(12,13)]

pos = which(rownames(avg_expr)%in%gene)
length(pos)

avg_expr = avg_expr[pos,]

library(pheatmap)

#pdf('/gpfs2/wulab10/RE/patch_3_project_deg_240703.pdf',height = 10,width = 5)
pheatmap(avg_expr,scale='row',cluster_rows = F,cluster_cols = F,border='black',color = colorRampPalette(colors=c('blue','white','white','#FF0000'))(100))
#dev.off()


library(ggplot2)

pdf('/gpfs2/wulab10/RE/patch_proj_type_umap_lmo3.pdf',height = 5,width = 5.5)
FeaturePlot(RE1,pt.size = 3.5,features = 'Lmo3',min.cutoff = 0.25,order=T)+scale_color_gradientn(colors=(c("#939393","#FB8861FF","#B63679FF","#51127CFF","#000004FF")))

dev.off()

pdf('/gpfs2/wulab10/RE/patch_proj_type_umap_ddiap1.pdf',height = 5,width = 5.5)
FeaturePlot(RE,pt.size = 3.5,features = 'Diap1',min.cutoff = 0.7,order=T)+scale_color_gradientn(colors=(c("gray","#FB8861FF","#B63679FF","#51127CFF","#000004FF")))

dev.off()



pdf('/gpfs2/wulab10/RE/patch_proj_type_umap_glut.pdf',height = 5,width = 5.5)
FeaturePlot(RE1,pt.size = 3.5,features = 'glut',order=T,min.cutoff = -1)+scale_color_gradientn(colors=(c("gray","#FB8861FF","#B63679FF","#51127CFF","#000004FF")))

dev.off()

pdf('/gpfs2/wulab10/RE/patch_proj_type_umap_gpc3.pdf',height = 5,width = 5.5)
FeaturePlot(RE1,pt.size = 3.5,features = 'Gpc3',min.cutoff = 0.35,order=T)+scale_color_gradientn(colors=(c("#939393","#FB8861FF","#B63679FF","#51127CFF","#000004FF")))

dev.off()

pdf('/gpfs2/wulab10/RE/patch_proj_type_umap_cdh13.pdf',height = 5,width = 5.5)
FeaturePlot(RE1,pt.size = 3.5,features = 'Cdh13',min.cutoff = 1,order=T)+scale_color_gradientn(colors=(c("#939393","#FB8861FF","#B63679FF","#51127CFF","#000004FF")))

dev.off()

pdf('/gpfs2/wulab10/RE/patch_proj_type_umap_trpc3.pdf',height = 5,width = 5.5)
FeaturePlot(RE1,pt.size = 3.5,features = 'Trpc3',min.cutoff = 0)+scale_color_gradientn(colors=(c("#939393","#FB8861FF","#B63679FF","#51127CFF","#000004FF")))

dev.off()

pdf('/gpfs2/wulab10/RE/patch_proj_type_umap_chrna4.pdf',height = 5,width = 5.5)
FeaturePlot(RE1,pt.size = 3.5,features = 'Chrna4',min.cutoff = 0.18)+scale_color_gradientn(colors=(c("#939393","#FB8861FF","#B63679FF","#51127CFF","#000004FF")))
dev.off()

pdf('/gpfs2/wulab10/RE/patch_proj_type_umap_ddhd1.pdf',height = 5,width = 5.5)
FeaturePlot(RE,pt.size = 3.5,features = 'Ddhd1',min.cutoff = 0.18)+scale_color_gradientn(colors=(c("gray","#FB8861FF","#B63679FF","#51127CFF","#000004FF")))
dev.off()

dual_g = c('Chrna4','Men1','Bap1','Peak1','Ddhd1','Trpc3','Tanc1','Pde4d','Ptk2b')

pdf('/gpfs2/wulab10/RE/patch_proj_type_umap_ptk2b.pdf',height = 5,width = 5.5)
FeaturePlot(RE,pt.size = 3.5,features = 'Ptk2b',min.cutoff = 0.8,order=F)+scale_color_gradientn(colors=(c("gray","#FB8861FF","#B63679FF","#51127CFF","#000004FF")))
dev.off()

pdf('/gpfs2/wulab10/RE/patch_proj_type_umap_nfix.pdf',height = 5,width = 5.5)
FeaturePlot(RE,pt.size = 3.5,features = 'Nfix',min.cutoff = 0.4,order=T)+scale_color_gradientn(colors=(c("gray","#FB8861FF","#B63679FF","#51127CFF","#000004FF")))
dev.off()

#pdf('/gpfs2/wulab10/RE/patch_proj_type_umap_hs3st4.pdf',height = 5,width = 5.5)
FeaturePlot(RE1,pt.size = 3.5,features = 'Hs3st4',min.cutoff = 0.5,order=T)+scale_color_gradientn(colors=(c("#939393","#FB8861FF","#B63679FF","#51127CFF","#000004FF")))
#dev.off()

pdf('/gpfs2/wulab10/RE/patch_proj_type_umap_tmem74b.pdf',height = 5,width = 5.5)
FeaturePlot(RE,pt.size = 3.5,features = 'Tmem74b',min.cutoff = 0.2,order=T)+scale_color_gradientn(colors=(c("gray","#FB8861FF","#B63679FF","#51127CFF","#000004FF")))
dev.off()

pdf('/gpfs2/wulab10/RE/patch_proj_type_umap_chrna7.pdf',height = 5,width = 5.5)
FeaturePlot(RE,pt.size = 3.5,features = 'Chrna7',min.cutoff = 0.1,order=T)+scale_color_gradientn(colors=(c("gray","#FB8861FF","#B63679FF","#51127CFF","#000004FF")))
dev.off()

pdf('/gpfs2/wulab10/RE/patch_proj_type_umap_chrm2.pdf',height = 5,width = 5.5)
FeaturePlot(RE,pt.size = 3.5,features = 'Chrm2',min.cutoff = 0.1,order=T)+scale_color_gradientn(colors=(c("gray","#FB8861FF","#B63679FF","#51127CFF","#000004FF")))
dev.off()

pdf('/gpfs2/wulab10/RE/patch_proj_type_umap_chrm3.pdf',height = 5,width = 5.5)
FeaturePlot(RE,pt.size = 3.5,features = 'Chrm3',min.cutoff = 0.1,order=T)+scale_color_gradientn(colors=(c("gray","#FB8861FF","#B63679FF","#51127CFF","#000004FF")))
dev.off()

pdf('/gpfs2/wulab10/RE/patch_proj_type_umap_chrm4.pdf',height = 5,width = 5.5)
FeaturePlot(RE,pt.size = 3.5,features = 'Chrm4',min.cutoff = 0.1,order=T)+scale_color_gradientn(colors=(c("gray","#FB8861FF","#B63679FF","#51127CFF","#000004FF")))
dev.off()

pdf('/gpfs2/wulab10/RE/patch_proj_type_umap_chrnb2.pdf',height = 5,width = 5.5)
FeaturePlot(RE,pt.size = 3.5,features = 'Chrnb2',min.cutoff = 0.1,order=T)+scale_color_gradientn(colors=(c("gray","#FB8861FF","#B63679FF","#51127CFF","#000004FF")))
dev.off()




tiff('/gpfs2/wulab10/RE/ncount_RNA.tif',height = 80,width = 120,units='mm',compression = 'lzw',res=600)
VlnPlot(RE1,features = 'nCount_RNA')+scale_fill_manual(values=c('dual' = '#E88118','vPFC' = '#92BFF4','mPFC'='#0F01F7'))
dev.off()

tiff('/gpfs2/wulab10/RE/ptc_ercc.tif',height = 80,width = 120,units='mm',compression = 'lzw',res=600)
VlnPlot(RE1,features = 'percent.ERCC')+scale_fill_manual(values=c('dual' = '#E88118','vPFC' = '#92BFF4','mPFC'='#0F01F7'))
dev.off()

tiff('/gpfs2/wulab10/RE/nFeature_rna.tif',height = 80,width = 120,units='mm',compression = 'lzw',res=600)
VlnPlot(RE1,features = 'nFeature_RNA')+scale_fill_manual(values=c('dual' = '#E88118','vPFC' = '#92BFF4','mPFC'='#0F01F7'))
dev.off()
