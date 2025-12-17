rm(list=ls()) #清空工作历史的所有变量
options(stringsAsFactors = F) #这两行命令可以清空环境
library(Seurat)
library(dplyr)
library(ggplot2)
library(magrittr)
library(gtools)
library(stringr)
library(Matrix)
library(tidyverse)
library(patchwork)
library(data.table)
library(RColorBrewer)
library(ggpubr)
library(ggsci)
dir_name=list.dirs('原始数据/',full.names = F,recursive = F)
dir_name
datalist=list()
#读取10x的数据创建CreateSeuratObject对象
for (i in 1:length(dir_name)){
  dir.10x = paste0("原始数据/",dir_name[i])
  my.data <- Read10X(data.dir = dir.10x) 
  #细胞增加标签
  colnames(my.data)=paste0(dir_name[i],colnames(my.data))
  datalist[[i]]=CreateSeuratObject(counts = my.data, project = dir_name[i], min.cells = 3, min.features = 250)
  datalist[[i]]$Samples=dir_name[i]
  datalist[[i]]$type=substr(dir_name[i],1,1)
}
names(datalist)=dir_name
#批量计算线粒体和rRNA的含量
for (i in 1:length(datalist)) {
  sce <- datalist[[i]]
  sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^MT-") 
  datalist[[i]] <- sce  
  rm(sce)  
}

#合并所有的数据
sce <- merge(datalist[[1]],y=datalist[2:length(datalist)])

#细胞数的统计
raw_cell=sce@meta.data
raw_count <- table(raw_cell$Samples)
raw_count
sum(raw_count)
pearplot_befor<-VlnPlot(sce,group.by ='Samples', 
                        features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.Ribo"),
                        pt.size = 0, ncol = 4)

pearplot_befor
ggsave('results/pearplot_befor.pdf',pearplot_befor,height = 5,width = 15)
ggsave('results/pearplot_befor.jpg',pearplot_befor,height = 5,width = 15,dpi = 300)

#样本的颜色
sample_color<-pal_nejm(alpha = 0.5)(8)[1:8]
#先把颜色调一下
sample_color <- c("#66c2a5","#fc8d62","#80b1d3","#ffdbb3","#f197a8","#b0e3ff","#bf847f","#b7d200",
                  
                  "#f68563","#af3e92")
sample_color
Feature_ber1<-FeatureScatter(sce,feature1 = 'nFeature_RNA',
                             feature2 = 'nCount_RNA',
                             group.by = 'Samples',
                             cols = sample_color)
Feature_ber2<-FeatureScatter(sce,feature1 = 'percent.mt',
                             feature2 = 'nCount_RNA',
                             group.by = 'Samples',
                             cols = sample_color)
Feature_ber3<-FeatureScatter(sce,feature1 = 'percent.mt',
                             feature2 = 'nFeature_RNA',
                             group.by = 'Samples',
                             cols = sample_color)
Feature_ber1=Feature_ber1+theme(legend.position = 'none')
Feature_ber2=Feature_ber2+theme(legend.position = 'none')

Feature_ber<-ggarrange(Feature_ber1,Feature_ber2,Feature_ber3,ncol = 3,nrow = 1,widths = c(1,1,1.2))
ggsave('results/Feature_cor.pdf',Feature_ber,height = 5,width = 17)
ggsave('results/Feature_cor.jpg',Feature_ber,height = 5,width = 17,dpi = 300)

#过滤
datalist <- lapply(X = datalist, FUN = function(x) {
  x<-subset(x,subset =nFeature_RNA < 5000 & percent.mt < 15)
})
#合并数据
sce <- merge(datalist[[1]],y=datalist[2:length(datalist)])
clean_cell=sce@meta.data

clean_count <- table(clean_cell$Samples)
clean_count
sum(clean_count)  

pearplot_after <- VlnPlot(sce,group.by ='Samples', 
                          features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.Ribo"),
                          pt.size = 0, 
                          ncol = 4)
pearplot_after
ggsave('results/pearplot_after.pdf',Feature_ber,height = 5,width = 15)
ggsave('results/pearplot_after.jpg',Feature_ber,height = 5,width = 15,dpi = 300)

#保存数据
save(datalist,file = 'datalist.RData')
#标准化
sce <- merge(datalist[[1]],y=datalist[2:length(datalist)])
sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 10000)
#筛选高变基因
sce <- FindVariableFeatures(sce, 
                            selection.method = "vst", 
                            nfeatures = 2000, #筛选前2000个高变，可修改的
                            mean.cutoff=c(0.0125,3),
                            dispersion.cutoff =c(1.5,Inf))
#对全部的数据进行scale
sce <- ScaleData(sce, features =  rownames(sce))
#PCA降维，将高维降到低维
sce <- RunPCA(sce, features = VariableFeatures(sce)) 

elbowplot <- ElbowPlot(sce, ndims=50, reduction="pca") 
elbowplot
ggsave('results/elbowplot.pdf',height = 5,width = 5)

#可修改，选择合适的PC进行后续的聚类和降维
Dims <- 30
sce <- RunUMAP(sce, dims=1:Dims, reduction="pca")
sce <- RunTSNE(sce,dims=1:Dims, reduction="pca")
raw.umap<-DimPlot(sce,group.by='Samples',
                  reduction="umap",
                  label = "T", 
                  pt.size = 0.2,
                  label.size = 0)+
  ggtitle('')
raw.umap
ggsave('results/raw.umap.pdf',raw.umap,height = 7,width = 7)

raw.tsne<-DimPlot(sce,group.by='Samples',
                  reduction="tsne",
                  label = "T", 
                  pt.size = 0.2,
                  label.size = 0)+
  ggtitle('')
raw.tsne
ggsave('results/raw.tsne.pdf',raw.tsne,height = 7,width = 7)
#聚类
library(clustree)
sce <- FindNeighbors(sce, dims = 1:Dims)
sce <- FindClusters(
  object = sce,
  resolution = c(seq(.1,1,.1))
)
colnames(sce@meta.data)
clustree(sce@meta.data, prefix = "RNA_snn_res.")

pdf('results/clust.snn_res.pdf',he=15,wi=15)

dev.off()

#注释NK细胞
#聚类分析
Resolution <- 1
sce <- FindNeighbors(object = sce, dims = 1:Dims)
sce <- FindClusters(object = sce, resolution = Resolution)
#亚群注释
DefaultAssay(sce) <- "RNA"
#NK:5,8,10,16,17,21,20,25 

VlnPlot(sce,features = c('KLRK1','NCAM1','FCGR3A','KLRG1'),pt.size = 0,group.by = 'seurat_clusters',ncol = 2)
library(randomcoloR)
allcolour <- c(pal_npg(alpha = 0.8)(9),
               pal_igv(alpha = 0.8)(9),
               pal_jama(alpha = 0.8)(7),
               pal_jco(alpha = 0.8)(9),
               pal_nejm(alpha = 0.8)(8))
length(table(sce@active.ident))
#32
mycolor1 = allcolour[1:length(table(sce$seurat_clusters))]

figs2b<-FeaturePlot(sce,
                    features =  c('KLRK1','NCAM1','FCGR3A','KLRG1'),
                    pt.size = 0.3,reduction = 'tsne',ncol = 2)
figs2a<-DimPlot(sce,cols =mycolor1 ,group.by = 'seurat_clusters',
                reduction="tsne",
                label = "T", 
                pt.size = 0.3,
                label.size = 3) +
  theme(axis.line = element_line(size=0.1, colour = "black"), 
        #axis.text = element_blank(), 
        #axis.title = element_blank(),
        axis.ticks = element_blank()
  ) +ggtitle('')

figs2ab<-ggarrange(figs2a,figs2b,nrow = 1,ncol = 2,widths = c(1,1),labels = c('A','B'))
figs2ab
ggsave('results/figs2ab.pdf',figs2ab,height = 7,width = 17)
table(sce$seurat_clusters)
load('sce1.RData')
save(sce,file = 'sce.RData')
DimPlot(sce,cols =mycolor1 ,group.by = 'seurat_clusters',
        reduction="tsne",
        label = "T", 
        pt.size = 0.4,
        label.size = 3) +
  theme(axis.line = element_line(size=0.1, colour = "black"), 
        #axis.text = element_blank(), 
        #axis.title = element_blank(),
        axis.ticks = element_blank()
  ) +ggtitle('')

#提取9,11,21,重新聚类
load('sce1.RData')
Idents(sce)='seurat_clusters'
sce<-subset(sce,idents =c(6,8,11,12,14,16,31,32))
#二次聚类
Resolution <- 0.1
DefaultAssay(sce) <- "RNA"

sce <- FindNeighbors(object = sce, dims = 1:10)
sce <- FindClusters(object = sce, resolution = Resolution)

DefaultAssay(sce) <- "RNA"

VlnPlot(sce,features = c('KLRK1','NCAM1','FCGR3A','KLRG1'),pt.size = 0,group.by = 'seurat_clusters',ncol = 2)

#重新降维
sce <- RunUMAP(sce, 
               dims=1:10, 
               reduction="pca",
               perplexity=30,
               max_iter=1000)
sce <- RunTSNE(sce, 
               dims=1:Dims, 
               reduction="pca",
               perplexity=30,
               max_iter=1000)
figs2c<-DimPlot(sce,cols =mycolor1 ,group.by = 'seurat_clusters',
                reduction="umap",
                label = "T", 
                pt.size = 0.2,
                label.size = 3) +
  theme(axis.line = element_line(size=0.1, colour = "black"), 
        #axis.text = element_blank(), 
        #axis.title = element_blank(),
        axis.ticks = element_blank()
  ) +ggtitle('')

figs2c <- DimPlot(sce,cols =mycolor1 ,group.by = 'seurat_clusters',
                  reduction="tsne",
                  label = "T", 
                  pt.size = 0.4,
                  label.size = 5) +
  theme(axis.line = element_line(size=0.1, colour = "black"), 
        #axis.text = element_blank(), 
        #axis.title = element_blank(),
        axis.ticks = element_blank()
  ) +ggtitle('')

figs2d<-FeaturePlot(sce,
                    features =  c('KLRK1','NCAM1','FCGR3A','KLRG1'),
                    pt.size = 0.1,reduction = 'tsne',ncol = 2)

fig2cd <- ggarrange(figs2c,figs2d,nrow = 1,ncol = 2,widths = c(1,1),labels = c('C','D'))
fig2cd

figs2<-ggarrange(figs2ab,fig2cd,nrow = 2,ncol = 1)
ggsave('results/FigS3.pdf',figs2,height =15,width = 15)
FeaturePlot(sce,
            features =  c('KLRK1','NCAM1','FCGR3A','KLRG1'),
            pt.size = 0.1,reduction = 'tsne',ncol = 2)

#寻找差异基因时的差异倍数
Logfc = 0.5
#差异基因时最小的表达比例
Minpct = 0.35
DefaultAssay(sce) <- 'RNA'

Idents(sce)<-'seurat_clusters'
sce.markers <- FindAllMarkers(object = sce,logfc.threshold = Logfc, min.pct = Minpct,only.pos = T)
sce.markers$pct.diff=sce.markers$pct.1-sce.markers$pct.2
sce.markers <- sce.markers[sce.markers$p_val_adj<0.05,]
length(unique(sce.markers$gene))
head(sce.markers)
write.table(sce.markers,'results/scRNA_marker_gene.txt',quote = F,row.names = F,sep='\t')

### 选择前5个marker基因
Top5 <- sce.markers %>% group_by(cluster) %>% slice_max(n =5, order_by = avg_log2FC)  

Top5 <- intersect(unique(Top5$gene),rownames(sce@assays$RNA@meta.features))

sc_marker_dotplot <- DotPlot(object = sce, features = Top5,cols=c("blue", "red"), scale = T)+ 
  RotatedAxis()+ ggtitle("Top 5 Marker Genes")+ 
  theme(plot.title = element_text(hjust = 0.5)) +xlab('')

sc_marker_dotplot
ggsave('results/sc_marker_dotplot.pdf',sc_marker_dotplot,height = 7,width = 9)
#第二种
bubble.df=as.matrix(sce[["RNA"]]@data[Top5,])
bubble.df=t(bubble.df)
bubble.df=as.data.frame(scale(bubble.df))
bubble.df$CB=rownames(bubble.df)
bubble.df=merge(bubble.df,
                data.frame(CB=rownames(sce@meta.data),
                           celltype=sce@meta.data$seurat_clusters),
                by = "CB")
bubble.df$CB=NULL

celltype_v=c()
gene_v=c()
mean_v=c()
ratio_v=c()
for (i in unique(bubble.df$celltype)) {
  bubble.df_small=bubble.df%>%filter(celltype==i)
  for (j in Top5) {
    exp_mean=mean(bubble.df_small[,j])
    exp_ratio=sum(bubble.df_small[,j] > min(bubble.df_small[,j])) / length(bubble.df_small[,j])
    celltype_v=append(celltype_v,i)
    gene_v=append(gene_v,j)
    mean_v=append(mean_v,exp_mean)
    ratio_v=append(ratio_v,exp_ratio)
  }
}
plotdf=data.frame(
  celltype=celltype_v,
  gene=gene_v,
  exp=mean_v,
  ratio=ratio_v)
plotdf$celltype=factor(plotdf$celltype,levels = unique(as.character(sce.markers$cluster)))
plotdf$gene=factor(plotdf$gene,levels = rev(as.character(Top5)))
plotdf$exp=ifelse(plotdf$exp>3,3,plotdf$exp)
sc_marker_dotplot1<-plotdf%>%ggplot(aes(x=celltype,y=gene,size=ratio,color=exp))+geom_point()+
  scale_x_discrete("")+scale_y_discrete("")+
  scale_color_gradientn(colours = rev(c("#FFD92F","#FEE391",brewer.pal(11, "Spectral")[7:11])))+
  scale_size_continuous(limits = c(0,1))+theme_bw()+
  theme(
    axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45)
  )
sc_marker_dotplot1
ggsave('results/sc_marker_dotplot1.pdf',sc_marker_dotplot1,height = 7,width = 9)

#绘图
save(sce,file = 'sce1.rdata')
load('sce1.rdata')
mycolor =pal_npg('nrc')(9)

fig1a = DimPlot(sce,group.by = 'Samples',
                reduction="tsne",
                label = "F", 
                pt.size = 0.5,
                label.size = 5) +
  theme(axis.line = element_line(size=0.1, colour = "black"), 
        #axis.text = element_blank(), 
        #axis.title = element_blank(),
        axis.ticks = element_blank()
  ) +ggtitle('')+guides(colour = guide_legend(ncol = 1))

fig1a
ggsave('results/fig1a.pdf',fig1a,height = 7,width = 7)

#mycolor <- rainbow(41)  # Generates 41 distinct colors
#library(RColorBrewer)
#mycolor <- brewer.pal(n = 8, name = "Set1")  # Change `n` to a larger number if needed
#mycolor <- rep(mycolor, length.out = 41)  # Repeat colors if fewer than 41
table(sce@meta.data$orig.ident)
sce@meta.data$Group <- ifelse(
  grepl(pattern = "NT$", x = sce@meta.data$orig.ident),  # 匹配以NT结尾的样本
  yes = "Normal",                                      # NT → Normal组
  no = "Tumor"                                         # T → Tumor组
)
fig1b<-DimPlot(sce,cols=mycolor,group.by = 'seurat_clusters',
               reduction="tsne",split.by = 'Group',
               label = "F", 
               pt.size = 0.5,
               label.size = 5) +
  theme(axis.line = element_line(size=0.1, colour = "black"), 
        #axis.text = element_blank(), 
        #axis.title = element_blank(),
        axis.ticks = element_blank()
  ) +ggtitle('')

fig1b
ggsave('results/fig1b.pdf',fig1b,height = 7,width = 7)
#
Idents(sce)='seurat_clusters'
library("ggplot2")
sample_clust<-as.matrix(table(sce$Samples,sce$seurat_clusters))
sample_clust=apply(sample_clust,1,function(x){return(x/sum(x))})
sample_clust=reshape2::melt(sample_clust)
colnames(sample_clust)<-c("cluster","Samples","proportion")
sample_clust$cluster=paste0('NK_',sample_clust$cluster)
write.table(sample_clust,'results/sample_clust1.txt',quote = F,row.names = T,sep='\t')

clust_freq<-as.data.frame(table(sce$Samples))
colnames(clust_freq)=c('Samples','cell_num')
clust_freq=clust_freq[order(clust_freq$cell_num,decreasing = T),]
clust_freq$Samples=factor(clust_freq$Samples,levels = clust_freq$Samples)
sample_clust$Samples=factor(sample_clust$Samples,levels =clust_freq$Samples)

fig1e1<-ggplot(sample_clust,aes(x = Samples,y = proportion,fill=cluster))+
  geom_bar(stat = "identity", position = "fill") +
  ggtitle("") +scale_fill_manual(values = mycolor[1:9])+
  theme_bw() + 
  theme(axis.ticks.length = unit(0.1, 'cm'),
        legend.position = "left") +xlab('')+
  coord_flip()+scale_y_continuous(expand = expand_scale(mult = c(0, 0)))
fig1e1
fig1e1 <- ggplot(sample_clust, aes(x = Samples, y = proportion, fill = cluster)) +
  geom_bar(stat = "identity", position = "fill") +
  ggtitle("") +
  scale_fill_manual(values = mycolor[1:9]) +
  theme_bw() +
  # 核心：统一调整字体大小 + 保留原主题设置
  theme(
    # 全局字体大小（基础参考，可按需调整）
    text = element_text(size = 14),  # 全局默认字体大小（建议10-14）
    
    # 坐标轴相关字体
    axis.text = element_text(size = 10),    # 坐标轴刻度文本大小
    axis.title = element_text(size = 12),   # 坐标轴标题大小（当前xlab为空，可保留）
    axis.ticks.length = unit(0.1, 'cm'),    # 保留原刻度长度
    
    # 图例相关字体
    legend.position = "left",
    legend.title = element_text(size = 11), # 图例标题大小
    legend.text = element_text(size = 9),   # 图例项文本大小
    
    # 面板样式（保留原theme_bw基础）
    panel.grid = element_blank(),           # 可选：移除网格线，更简洁
    panel.border = element_rect(linewidth = 0.5) # 面板边框粗细
  ) +
  xlab('') +
  coord_flip() +
  scale_y_continuous(expand = expand_scale(mult = c(0, 0)))

# 查看调整后的图表
print(fig1e1)

fig1e2<-ggplot(clust_freq,aes(x = Samples,y = cell_num,fill=Samples))+
  geom_bar(stat="identity")+ggtitle("") +
  theme_bw() + scale_fill_manual(values = sample_color)+
  theme(axis.ticks.length = unit(0, 'cm'),
        axis.text.y = element_blank(),
        axis.title.y = element_blank()) +coord_flip()+
  scale_y_continuous(expand = expand_scale(mult = c(0, 0)))+ylim(0,max(clust_freq$cell_num)+10)
fig1e2

fig1e3<-ggpubr::ggarrange(fig1e1,fig1e2,nrow = 1,ncol = 2,widths = c(2,1))
fig1e3
ggsave('results/fig1e3.pdf',fig1e3,height = 7,width = 7)

#marker基因进行注释
library(clusterProfiler)
library(org.Hs.eg.db)
ids=bitr(sce.markers$gene,'SYMBOL','ENTREZID','org.Hs.eg.db') ## 将SYMBOL转成ENTREZID
sce.markers2=merge(sce.markers,ids,by.x='gene',by.y='SYMBOL')

## 函数split()可以按照分组因子，把向量，矩阵和数据框进行适当的分组。
## 它的返回值是一个列表，代表分组变量每个水平的观测。
gcSample=split(sce.markers2$ENTREZID, sce.markers2$cluster) 
## KEGG
sce.markers2.enrich.res <- compareCluster(gcSample,
                                          fun = "enrichKEGG",
                                          organism = "hsa", pvalueCutoff = 0.05)
fig1f<-dotplot(sce.markers2.enrich.res)+ 
  theme(axis.text.x = element_text(angle = 45,hjust = 1,size=10),
        axis.text.y=element_text(size=10))
fig1f
ggsave('results/fig1f.pdf',fig1f,height = 7,width = 7)
save(sce,file = 'sce2.RData')

#恶性和非恶性的区别
load('sce2.RData')
library(copykat)
library(Seurat)
library(gplots)
library(ggplot2)
library(qs)
library(SingleCellExperiment)
library(Matrix)

exp.rawdata <- as.matrix(sce@assays$RNA@counts)

copykat.test <- copykat(rawmat=exp.rawdata,
                        id.type='S',
                        cell.line='no',
                        ngene.chr=5, 
                        #每个染色体中至少有 5 个基因来计算 DNA 拷贝数
                        win.size=25, 
                        #每个片段至少取 25 个基因
                        KS.cut=0.15, 
                        #0-1,值越大灵敏度越低
                        sam.name = 'HCC',  
                        #随意固定一个名
                        distance='euclidean', 
                        #并行计算
                        n.cores=8)

save(copykat.test,file = 'copykat.test.RData')

#读取CNV的结果
copykat.test<-read.delim('HCC_copykat_prediction.txt',sep='\t',header = T)
head(copykat.test)
table(copykat.test$copykat.pred)
rownames(copykat.test)=copykat.test$cell.names
copykat.test=copykat.test[rownames(sce@meta.data),]
#添加分组
sce <- AddMetaData(sce, copykat.test$copykat.pred,col.name = 'copykat.pred')
sce$copykat.pred[is.na(sce$copykat.pred)]<-'Unknown'
table(sce$copykat.pred)

sce$copykat.pred=ifelse(sce$copykat.pred=='aneuploid','malignant','no_malignant')

save(sce,file = 'sce2.RData')

fig1h<-DimPlot(sce,cols=c('red','blue'),group.by = 'copykat.pred',
               reduction='tsne',
               label = 'F', 
               pt.size = 0.5,
               label.size = 5) +
  theme(axis.line = element_line(size=0.1, colour = 'black'),
        #axis.text = element_blank(), 
        #axis.title = element_blank(),
        axis.ticks = element_blank()
  ) + ggtitle('')

fig1h
ggsave(filename = 'results/fig1h.pdf',plot = fig1,he=15,wi=18)
fig1ef<-ggarrange(fig1e3,fig1f,fig1h,labels = c('D','E','F'),nrow = 1,ncol = 3,widths = c(1.2,1.3,1))

fig1ab<-ggarrange(fig1a,fig1b,nrow = 1,ncol=2,labels = c('A','B'),widths = c(1,1.5))
fig1=ggarrange(fig1ab,sc_marker_dotplot,fig1ef,labels = c('','C',''),nrow = 3,ncol = 1,heights = c(2,1,1))

ggsave(filename = 'results/Fig1.pdf',plot = fig1,he=15,wi=18)
ggsave(filename = 'results/Fig1.jpg',plot = fig1,he=15,wi=18)

#### 计算细胞恶性积分####
rm(list=ls()) #清空工作历史的所有变量
options(stringsAsFactors = F,check.bounds = F) #这两行命令可以清空环境
setwd("~/LIHC/GSE242889/")
dir.create('results2')
setwd("results2")

library(Seurat)
library(dplyr)
library(data.table) #data.table包是自带包data.frame的升级版，用于数据框格式数据的处理，最大的特点是快
library(R.utils)
library(tidyverse)
library(patchwork)
library(ggplot2)
library(cowplot)
library(hdf5r)
library(devtools)
library(magrittr)
library(gtools)
library(stringr)
library(Matrix)
library(patchwork)
library(data.table)
library(RColorBrewer)
library(ggpubr)
library(ggsci)
#肿瘤发生的10个关键通路
load('sce2.RData')
pathway.score<-function(exp,gene){
  ssGSEAScore_by_genes<-function(gene.exp,genes){
    gs=GSEABase::GeneSet(setName='GeneSet', setIdentifier=paste0('101'), geneIds=unique(genes), GSEABase::SymbolIdentifier()) 
    
    gsc <- GSEABase::GeneSetCollection(list(gs))
    fl <- tempfile()
    GSEABase::toGmt(gsc, fl)
    cgeneset=GSEABase::getGmt(fl)
    ssGSEA.geneset <- GSVA::gsva(as.matrix(gene.exp), cgeneset,method='ssgsea',
                                 min.sz=1, max.sz=5000, verbose=TRUE)
    return(ssGSEA.geneset)
  }
  
  pathway_score<-data.frame()
  for (i in unique(gene[,2])){
    gene_set=gene[gene[,2]==i,1]
    score=ssGSEAScore_by_genes(exp,gene_set)
    rownames(score)=i
    pathway_score=rbind.data.frame(pathway_score,score)
  }
  return(t(pathway_score))
}
#_pmid_29625050
tumor.pathway=read.delim('~/LIHC/GSE242889/pmid_29625050_pathway.txt',sep='\t',header = T)
head(tumor.pathway)
tumor.pathway=tumor.pathway[,c('Gene','OG.TSG')]

#每一个细胞计算得分
tumor.pathway.score<-pathway.score(exp = as.matrix(sce@assays$RNA@counts),gene = tumor.pathway)
head(tumor.pathway.score)

tumor.pathway.score.group<-merge(data.frame(cell.names=rownames(sce@meta.data),
                                            sce@meta.data),
                                 data.frame(cell.names=rownames(tumor.pathway.score),
                                            tumor.pathway.score),
                                 by='cell.names')
rownames(tumor.pathway.score.group)=tumor.pathway.score.group$cell.names
head(tumor.pathway.score.group)
tumor.pathway.score.group=tumor.pathway.score.group[,-1]

#读取CNV的结果
copykat.test=data.frame(cell.names=rownames(sce@meta.data),
                        copykat.pred=sce@meta.data$copykat.pred)

head(copykat.test)

table(copykat.test$copykat.pred)
tumor.score.copy<-cbind.data.frame(tumor.pathway.score.group[copykat.test$cell.names,],
                                   copykat.pred=copykat.test$copykat.pred)
head(tumor.score.copy)
table(tumor.score.copy$copykat.pred)
tumor.score.copy$seurat_clusters=paste0('NK_',tumor.score.copy$seurat_clusters)
library(pheatmap)
head(tumor.score.copy)
table(tumor.score.copy$copykat.pred)

tumor.score.copy=tumor.score.copy[,c('Samples','seurat_clusters',
                                     'CellCyle','HIPPO','MYC','NOTCH','NRF1','PI3K','TGF.Beta','RAS','TP53','WNT',
                                     'copykat.pred')]
colnames(tumor.score.copy)[9]
colnames(tumor.score.copy)[9]='TGF-Beta'
mat=tumor.score.copy[,as.character(unique(tumor.pathway$OG.TSG))]
anno_col<-tumor.score.copy[,c('seurat_clusters','copykat.pred')]
anno_col=anno_col[order(anno_col$copykat.pred,anno_col$seurat_clusters),]
pdf('results/Fig2a.pdf',height = 9,width = 12)

pheatmap(
  t(mat),  # Transpose if necessary
  scale = 'row',
  show_colnames = FALSE,
  annotation_col = anno_col,
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  color = colorRampPalette(c("#FF7F00", "white", "#E41A1C"))(100),
  annotation_names_row = FALSE,  # Or TRUE based on your requirement
  annotation_colors = list(copykat.pred = c('malignant' = "red", 'no_malignant' = "blue")),
  breaks = unique(c(seq(-2, 2, length = 100)))
)


dev.off()  # 关闭PDF设备

write.table(tumor.score.copy,'results/tumor.score.copy.txt',quote = F,row.names = T,sep='\t')
#恶性细胞和非恶性细胞中CAF的比例
clust.malig<-table(tumor.score.copy$copykat.pred,tumor.score.copy$seurat_clusters)
clust.malig
write.table(clust.malig,'results/clust.malig.txt',quote = F,sep='\t',row.names = T)

plotiBarplot<-function(dat,palette,ist=F,margin=T,lineCol='black',legTitle='Group',showValue=F,showLine=T){
  library(ggplot2)
  xlb='';ylb='';lineW=0.5;xangle=0;isAuto=T;
  library(tidyverse)
  library(reshape2)
  library(optparse)
  if(ist){
    dat=t(dat)
  }
  lbc=colnames(dat)
  lbr=row.names(dat)
  bk_dat=dat
  if(margin){
    dat=dat%*%diag(1/c(apply(t(dat), 1, sum)))
  }
  row.names(dat)=paste0('R',1:(nrow(dat)))
  colnames(dat)=paste0('C',1:(ncol(dat)))
  row.names(bk_dat)=paste0('R',1:(nrow(bk_dat)))
  colnames(bk_dat)=paste0('C',1:(ncol(bk_dat)))
  #df=cbind(bg=paste0('R',1:nrow(dat)),dat)
  #colnames(df)=c('bg',paste0('C',1:(ncol(dat))))
  tp.dat=as.data.frame(cbind(bg=row.names(dat),dat))
  tp.dat[,1]=as.character(tp.dat[,1])
  for(i in 2:ncol(tp.dat)){
    tp.dat[,i]=as.numeric(as.character(tp.dat[,i]))
  }
  mt.df=reshape2::melt(tp.dat)
  colnames(mt.df)=c('bg','variable','value')
  
  pg=ggplot(mt.df, aes(x=variable, y=value, fill=bg))+
    geom_bar(stat = "identity", width=lineW, col=lineCol)
  if(showLine){
    for (i in 2:(ncol(tp.dat)-1)) {
      tmp=tp.dat[order(tp.dat[,1],decreasing = T),]
      tmp[,i]=base::cumsum(tmp[,i])
      tmp[,i+1]=base::cumsum(tmp[,i+1])
      colnames(tmp)[c(i,i+1)]=c('STY','ED')
      tmp1=cbind(tmp,STX=rep(i-1+lineW/2,nrow(tmp))
                 ,EDX=rep(i-lineW/2,nrow(tmp)))
      pg=pg+geom_segment(data=tmp1,aes(x=STX, xend=EDX, y=STY, yend=ED))
    }
  }
  
  if(showValue){
    pg=pg+geom_text(data=mt.df,aes(label=sprintf("%0.2f", round(value, digits = 2))),position=position_stack(vjust=0.5))
  }
  pg=pg+scale_x_discrete(breaks = paste0('C',1:(ncol(dat))),label = lbc)
  pg=pg+labs(x=xlb, y=ylb)+theme(legend.position = "bottom")
  #pg=pg+scale_fill_discrete(breaks = paste0('R',1:nrow(dat)),label = lbr,name=legTitle)
  pg=pg+scale_fill_manual(breaks = paste0('R',1:nrow(dat)),label = lbr,name=legTitle,values=palette)
  if(xangle>0){
    pg=pg+theme(axis.text.x = element_text(angle = xangle, hjust = 1),legend.position = "bottom")
  }
  
  g.tb=matrix(0,nrow=ncol(dat),ncol=ncol(dat))
  for(i in 1:(ncol(dat))){
    for(j in 1:ncol(dat)){
      if(i!=j){
        g.tb[i,j]=round(-log10((chisq.test(bk_dat[,c(i,j)])$p.value)),2)
      }
    }
  }
  colnames(g.tb)=lbc
  row.names(g.tb)=lbc
  g.tb=reshape2::melt(g.tb) 
  colnames(g.tb)=c('A1','A2','A3')
  g.tb$A4=paste0(g.tb[,3],ifelse(g.tb[,3]>-log10(0.05),'(*)',''))
  stable.p=ggplot(g.tb, aes(A1, A2)) + geom_tile(aes(fill = A3),colour = "white") +xlab('')+ylab('')+ scale_fill_gradient(low = "white",high = "steelblue")+geom_text(aes(x=A1,y=A2,label=A4))+theme(legend.position="none",axis.title.x=element_blank(),axis.text.x=element_blank())
  stable.p=stable.p+ggtitle('-log10(anova p value)')
  if(isAuto){
    g1=ggpubr::ggarrange(stable.p,pg, ncol = 1, nrow = 2,heights = c(0.5,1),align = "hv")
    return(g1)
  }else{
    return(list(Bar=pg,Table=stable.p))
  }
}

#cellcolor =ggsci::pal_jama()(9) #最多7个颜色
cellcolor = ggsci::pal_npg()(9)

fig2b<-plotiBarplot(dat = clust.malig,palette=c('#E41A1C','#FF7F00'),ist = F,margin=T,lineCol='black',
                    legTitle = 'Predict',showValue=T,showLine=T)
fig2b
ggsave('results/Fig2b.pdf',fig2b,height = 7,width = 10)
head(tumor.pathway.score.group)
table(tumor.score.copy$seurat_clusters)

#install.packages("ggplot2")
#install.packages("ggpubr")
#install.packages("reshape2")
##install.packages("ggdist") 

library(reshape2)
library(ggplot2)  # 确保 ggplot2 已加载
library(ggdist)
library(ggpubr)
#C0
tumor.pathway.score.group1=tumor.score.copy[tumor.score.copy$seurat_clusters=='NK_0',]
#箱线图
Muti_Boxplot<-function(dat,group,group_cols,leg,
                       test_method = 'wilcox.test',ylabs){
  library(ggpubr)
  library(ggplot2)
  library(reshape2)
  library(ggdist) 
  library(ggsci)
  dat1=reshape2::melt(cbind.data.frame(dat,group))
  p=ggplot(dat1,aes(x = variable,y = value,fill = group)) +
    # split violin
    geom_violin(trim = TRUE, scale = "area", na.rm = TRUE) +
    # mean point
    stat_summary(fun = 'mean', geom = 'point',position = position_dodge(0.2)) +
    # errorbar
    stat_summary(fun.data = 'mean_sd', geom = 'errorbar', width = .15,
                 size = 0.3,
                 position = position_dodge(0.2)) +
    theme_bw(base_size = 16) +
    theme(axis.text.x = element_text(angle = 90,color = 'black',hjust = 1),
          legend.position = 'top') +
    # scale_fill_brewer(palette = 'Set1') +
    scale_fill_jco(name = '') +
    ylim(0,4) + #修正Y轴限制
    stat_compare_means(aes(group=group),method = test_method,
                       symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                                        symbols = c('***', '**', '*', 'ns')),
                       label = 'p.signif')+
    theme(axis.text.x = element_text(angle = 30,hjust = 1))+labs(fill=leg)
  return(p)
}

fig2d<-Muti_Boxplot(dat =tumor.pathway.score.group1[,as.character(unique(tumor.pathway$OG.TSG))],
                    group = tumor.pathway.score.group1$copykat.pred,
                    group_cols = ggsci::pal_lancet()(9)[c(2,1)],
                    test_method = 'wilcox.test',
                    leg = 'NK_0',ylab = 'GSVA Score')
fig2d
ggsave('results2/NK_0.pdf',fig2d,height = 5,width = 9)
1
#C1
tumor.pathway.score.group2=tumor.score.copy[tumor.score.copy$seurat_clusters=='NK_1',]
fig2e<-Muti_Boxplot(dat =tumor.pathway.score.group2[,as.character(unique(tumor.pathway$OG.TSG))],
                    group = tumor.pathway.score.group2$copykat.pred,
                    group_cols = ggsci::pal_lancet()(9)[c(2,1)],
                    test_method = 'wilcox.test',
                    leg = 'NK_1',ylab = 'GSVA Score')
fig2e
ggsave('results2/NK_1.pdf',fig2e,height = 5,width = 9)

#C2
tumor.pathway.score.group3=tumor.score.copy[tumor.score.copy$seurat_clusters=='NK_2',]
fig2f<-Muti_Boxplot(dat =tumor.pathway.score.group3[,as.character(unique(tumor.pathway$OG.TSG))],
                    group = tumor.pathway.score.group3$copykat.pred,
                    group_cols = ggsci::pal_lancet()(9)[c(2,1)],
                    test_method = 'wilcox.test',
                    leg = 'NK_2',ylab = 'GSVA Score')
fig2f
ggsave('results2/NK_2.pdf',fig2f,height = 5,width = 9)
0
#C3
tumor.pathway.score.group4=tumor.score.copy[tumor.score.copy$seurat_clusters=='NK_3',]
fig2g<-Muti_Boxplot(dat =tumor.pathway.score.group4[,as.character(unique(tumor.pathway$OG.TSG))],
                    group = tumor.pathway.score.group4$copykat.pred,
                    group_cols = ggsci::pal_lancet()(9)[c(2,1)],
                    test_method = 'wilcox.test',
                    leg = 'NK_3',ylab = 'GSVA Score')
fig2g
ggsave('results2/NK_3.pdf',fig2g,height = 5,width = 9)

#C4
tumor.pathway.score.group5=tumor.score.copy[tumor.score.copy$seurat_clusters=='NK_4',]
fig2h<-Muti_Boxplot(dat =tumor.pathway.score.group5[,as.character(unique(tumor.pathway$OG.TSG))],
                    group = tumor.pathway.score.group5$copykat.pred,
                    group_cols = ggsci::pal_lancet()(9)[c(2,1)],
                    test_method = 'wilcox.test',
                    leg = 'CAF_4',ylab = 'GSVA Score')
fig2h
ggsave('results2/NK_4.pdf',fig2h,height = 5,width = 9)


tumor.pathway.score.group6=tumor.score.copy[tumor.score.copy$seurat_clusters=='NK_5',]
fig2h<-Muti_Boxplot(dat =tumor.pathway.score.group6[,as.character(unique(tumor.pathway$OG.TSG))],
                    group = tumor.pathway.score.group6$copykat.pred,
                    group_cols = ggsci::pal_lancet()(9)[c(2,1)],
                    test_method = 'wilcox.test',
                    leg = 'CAF_5',ylab = 'GSVA Score')
fig2h
ggsave('results/CAF_5.pdf',fig2h,height = 5,width = 9)

tumor.pathway.score.group7=tumor.score.copy[tumor.score.copy$seurat_clusters=='CAF_6',]
fig2h<-Muti_Boxplot(dat =tumor.pathway.score.group7[,as.character(unique(tumor.pathway$OG.TSG))],
                    group = tumor.pathway.score.group7$copykat.pred,
                    group_cols = ggsci::pal_lancet()(9)[c(2,1)],
                    test_method = 'wilcox.test',
                    leg = 'CAF_6',ylab = 'GSVA Score')
fig2h
ggsave('results/CAF_6.pdf',fig2h,height = 5,width = 9)

tumor.pathway.score.group8=tumor.score.copy[tumor.score.copy$seurat_clusters=='CAF_7',]
fig2h<-Muti_Boxplot(dat =tumor.pathway.score.group8[,as.character(unique(tumor.pathway$OG.TSG))],
                    group = tumor.pathway.score.group8$copykat.pred,
                    group_cols = ggsci::pal_lancet()(9)[c(2,1)],
                    test_method = 'wilcox.test',
                    leg = 'CAF_7',ylab = 'GSVA Score')
fig2h

tumor.pathway.score.group9=tumor.score.copy[tumor.score.copy$seurat_clusters=='CAF_8',]
fig2h<-Muti_Boxplot(dat =tumor.pathway.score.group9[,as.character(unique(tumor.pathway$OG.TSG))],
                    group = tumor.pathway.score.group9$copykat.pred,
                    group_cols = ggsci::pal_lancet()(9)[c(2,1)],
                    test_method = 'wilcox.test',
                    leg = 'CAF_8',ylab = 'GSVA Score')
fig2h
ggsave('results/CAF_8.pdf',fig2h,height = 5,width = 9)

save.image('all.RData')



####04.NK.select####
setwd("/home/liangyu-kaikai/LIHC/GEO/04.NK.select/")
dir.create('results')
bioSurvival=function(OS,OS.time,riskscore,labs,palette,leg){
  dat=data.frame(OS=OS,OS.time=OS.time,riskscore=riskscore)
  library(survival)
  library(survminer)
  res.cut <- surv_cutpoint(dat,
                           time = "OS.time", 
                           event = "OS", 
                           variables = c("riskscore"))
  cut_va=as.numeric(res.cut$cutpoint[1])
  print(cut_va)
  dat$risk=ifelse(dat$riskscore>=cut_va,'High','Low')
  ggplotKM<-function(time,status,group,labs,palette,leg){
    library(ggplot2)
    library(survival)
    dat1=data.frame(time=time,status=status,group=group)
    colnames(dat1)=c('time','status','groups')
    sf<-survival::survfit(Surv(time,status) ~ groups,data=dat1)
    surv=survminer::ggsurvplot(sf, data = dat1, 
                               palette = palette, 
                               pval = TRUE,
                               surv.median.line='hv'
                               #,conf.int = T
                               ,conf.int.style ='step'
                               , pval.coord=c(0, 0.2), #Add p-value 
                               risk.table = TRUE, 
                               legend.title = leg
                               ,legend.labs =labs
                               ,conf.int=T
    )
    p1=surv$plot+theme_pubr()+
      theme(axis.text.y=element_text(family="Times",face="plain"),
            #axis.text.x=element_blank(),
            axis.title.x=element_blank(),
            plot.margin=unit(c(0.2, 0.2, 0, 0.1), "inches"),
            legend.position=c(1,1),
            legend.justification=c(1,1),
            legend.background = element_rect(fill = NA, colour = NA),
            legend.title = element_text(family="Times",face="plain"),
            legend.text = element_text(family="Times",face="plain"))
    p2=surv$table+theme_pubr()+
      theme(axis.text.y=element_text(family="Times",face="plain"),
            plot.margin=unit(c(0, 0.2, 0.2, 0.1), "inches"),
            plot.title=element_blank(),
            legend.position=c(1,1), 
            legend.justification=c(1,1),
            legend.title = element_text(family="Times",face="plain"),
            legend.text = element_text(family="Times",face="plain"))+
      xlab('Year')
    
    g2=ggpubr::ggarrange(p1,p2, ncol = 1, nrow = 2,heights = c(1,0.3),align = "v")
    return(p1)
  }
  p=ggplotKM(time = dat$OS.time,status=dat$OS,group=dat$risk,labs=labs,palette=palette,leg=leg)
  return(p)
}
#
tcga.dat<-read.delim('../03.data.pre/GSE53624/results/geo_dat.txt',sep='\t',header = T,row.names = 1,check.names = F)
tcga.group<-read.delim('/home/liangyu-kaikai/LIHC/GEO/geogroup.txt',sep='\t',header = T)
#CAF_0的marke基因
cell.module.gene=read.delim('/home/liangyu-kaikai/LIHC/GSE242889/results/scRNA_marker_gene.txt',sep='\t',header = T)
cell.module.gene=cell.module.gene[,c("gene","cluster")]
head(cell.module.gene)
cell.module.gene$cluster=paste0('NK_',cell.module.gene$cluster)
#ssGSEA计算TCGA样本的NK评分
pathway.score<-function(exp,gene){
  ssGSEAScore_by_genes<-function(gene.exp,genes){
    #library('GSVA')
    #library(GSEABase)
    #all.list=list()
    gs=GSEABase::GeneSet(setName='GeneSet', 
                         setIdentifier=paste0("101"),
                         geneIds=unique(genes),
                         GSEABase::SymbolIdentifier()) 
    gsc <- GSEABase::GeneSetCollection(list(gs))
    fl <- tempfile()
    GSEABase::toGmt(gsc, fl)
    cgeneset=GSEABase::getGmt(fl)
    ssGSEA.geneset <- GSVA::gsva(as.matrix(gene.exp),
                                 cgeneset,method='ssgsea',
                                 min.sz=1, max.sz=5000, 
                                 verbose=TRUE)
    #detach('package:GSVA')
    #detach('package:GSEABase')
    #row.names(ssGSEA.geneset)
    return(ssGSEA.geneset)
  }
  pathway_score<-data.frame()
  for (i in unique(gene[,2])){
    gene_set=gene[gene[,2]==i,1]
    score=ssGSEAScore_by_genes(exp,gene_set)
    rownames(score)=i
    pathway_score=rbind.data.frame(pathway_score,score)
  }
  return(t(pathway_score))
}
tcga.nk.score<-pathway.score(exp = tcga.dat,gene = cell.module.gene)
head(tcga.nk.score)
library(ggplot2)
library(ggsignif)
library(ggdist) 
library(ggpubr)
sig_violin<-function(dat,leg,ylab,palette=ggsci::pal_lancet()(10)[3:4]){
  library(ggpubr)
  dat=na.omit(dat)
  colnames(dat)=c('group','gene')
  dat=dat[order(dat$group),]
  all.combn=combn(as.character(unique(dat$group)),2)
  my_comparisons=lapply(seq_len(ncol(all.combn)), function(i) all.combn[,i])
  p<-ggplot(dat, aes(group, gene, fill = group, color = group)) +
    geom_rain(alpha = .5, rain.side = 'l',
              boxplot.args = list(color = "black", outlier.shape = NA),
              boxplot.args.pos = list(
                position = ggpp::position_dodgenudge(x = .1), width = 0.1
              )) +
    theme_classic() +
    scale_fill_brewer(palette = 'Dark2') +
    scale_color_brewer(palette = 'Dark2') +
    guides(fill = 'none', color = 'none') +
    stat_compare_means(comparisons=my_comparisons,method="wilcox.test",label = "p.signif")+
    ylab(ylab)+xlab('')+labs(color=leg,fill=leg)+scale_color_manual(values = palette)+theme_classic()+
    scale_fill_manual(values = palette)
  return(p)
}
tcga.nk.score.group<-merge(geogroup,
                           data.frame(Samples=rownames(tcga.nk.score),tcga.nk.score),
                           by='Samples')
write.table(tcga.nk.score.group,'results/tcga.nk.score.txt',quote = F,row.names = T,sep='\t')

fig1a<-sig_violin(dat = tcga.nk.score.group[,c("Type","NK_0")],
                  leg = 'Groups',ylab = 'NK_0 Score',palette = c('#377EB8','indianred'))

fig1a <- ggplot(tcga.nk.score.group, aes(x = Type, y = NK_0, fill = Type)) +
  geom_violin(trim = FALSE, alpha = 0.5) +  # 绘制小提琴图
  geom_boxplot(width = 0.1, color = "black", alpha = 0.5) +  # 添加箱线图
  labs(title = "NK_0 Score by Groups", x = "Groups", y = "NK_0 Score") +
  scale_fill_manual(values = c('#377EB8', 'indianred')) +
  theme_minimal() +
  stat_compare_means(comparisons = list(c("N", "T")),
                     label = "p.signif",map_signif_level = TRUE, 
                     textsize = 8,  # 设置显著性文本的字体大小
                     size = 8)   # 设置显著性线条的大小) +  # 使用显著性标记
theme(
  plot.title = element_text(size = 16, face = "bold"),  # 标题字体大小
  axis.title.x = element_text(size = 14),                # x轴标题字体大小
  axis.title.y = element_text(size = 14),                # y轴标题字体大小
  axis.text.x = element_text(size = 12),                 # x轴刻度字体大小
  axis.text.y = element_text(size = 12),                 # y轴刻度字体大小
  legend.title = element_text(size = 14),                 # 图例标题字体大小
  legend.text = element_text(size = 12)                   # 图例文本字体大小
)
fig1a

fig1b<-sig_violin(dat = tcga.nk.score.group[,c("Type","NK_1")],
                  leg = 'Groups',ylab = 'NK_1 Score',palette = c('#377EB8','indianred'))

fig1b <- ggplot(tcga.nk.score.group, aes(x = Type, y = NK_1, fill = Type)) +
  geom_violin(trim = FALSE, alpha = 0.5) +  # 绘制小提琴图
  geom_boxplot(width = 0.1, color = "black", alpha = 0.5) +  # 添加箱线图
  labs(title = "NK_1 Score by Groups", x = "Groups", y = "NK_1 Score") +
  scale_fill_manual(values = c('#377EB8', 'indianred')) +
  theme_minimal() +
  stat_compare_means(comparisons = list(c("N", "T")),
                     label = "p.signif",map_signif_level = TRUE, 
                     textsize = 8,  # 设置显著性文本的字体大小
                     size = 8)   # 设置显著性线条的大小) +  # 使用显著性标记
theme(
  plot.title = element_text(size = 16, face = "bold"),  # 标题字体大小
  axis.title.x = element_text(size = 14),                # x轴标题字体大小
  axis.title.y = element_text(size = 14),                # y轴标题字体大小
  axis.text.x = element_text(size = 12),                 # x轴刻度字体大小
  axis.text.y = element_text(size = 12),                 # y轴刻度字体大小
  legend.title = element_text(size = 14),                 # 图例标题字体大小
  legend.text = element_text(size = 12)                   # 图例文本字体大小
)


fig1b




fig1c<-sig_violin(dat = tcga.caf.score.group[,c("Type","CAF_2")],
                  leg = 'Groups',ylab = 'CAF_2 Score',palette = c('#377EB8','indianred'))
fig1c <- ggplot(tcga.nk.score.group, aes(x = Type, y = NK_2, fill = Type)) +
  geom_violin(trim = FALSE, alpha = 0.5) +  # 绘制小提琴图
  geom_boxplot(width = 0.1, color = "black", alpha = 0.5) +  # 添加箱线图
  labs(title = "NK_2 Score by Groups", x = "Groups", y = "NK_2 Score") +
  scale_fill_manual(values = c('#377EB8', 'indianred')) +
  theme_minimal() +
  stat_compare_means(comparisons = list(c("N", "T")),
                     label = "p.signif",map_signif_level = TRUE, 
                     textsize = 8,  # 设置显著性文本的字体大小
                     size = 8)   # 设置显著性线条的大小) +  # 使用显著性标记
theme(
  plot.title = element_text(size = 16, face = "bold"),  # 标题字体大小
  axis.title.x = element_text(size = 14),                # x轴标题字体大小
  axis.title.y = element_text(size = 14),                # y轴标题字体大小
  axis.text.x = element_text(size = 12),                 # x轴刻度字体大小
  axis.text.y = element_text(size = 12),                 # y轴刻度字体大小
  legend.title = element_text(size = 14),                 # 图例标题字体大小
  legend.text = element_text(size = 12)                   # 图例文本字体大小
)

fig1c


fig1d<-sig_violin(dat = tcga.caf.score.group[,c("Type","CAF_3")],
                  leg = 'Groups',ylab = 'CAF_3 Score',palette = c('#377EB8','indianred'))

fig1d <- ggplot(tcga.nk.score.group, aes(x = Type, y = NK_3, fill = Type)) +
  geom_violin(trim = FALSE, alpha = 0.5) +  # 绘制小提琴图
  geom_boxplot(width = 0.1, color = "black", alpha = 0.5) +  # 添加箱线图
  labs(title = "NK_3 Score by Groups", x = "Groups", y = "NK_3 Score") +
  scale_fill_manual(values = c('#377EB8', 'indianred')) +
  theme_minimal() +
  stat_compare_means(comparisons = list(c("N", "T")),
                     label = "p.signif",map_signif_level = TRUE, 
                     textsize = 8,  # 设置显著性文本的字体大小
                     size = 8) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),  # 标题字体大小
    axis.title.x = element_text(size = 14),                # x轴标题字体大小
    axis.title.y = element_text(size = 14),                # y轴标题字体大小
    axis.text.x = element_text(size = 12),                 # x轴刻度字体大小
    axis.text.y = element_text(size = 12),                 # y轴刻度字体大小
    legend.title = element_text(size = 14),                 # 图例标题字体大小
    legend.text = element_text(size = 12)                   # 图例文本字体大小
  )

fig1d

fig1e<-sig_violin(dat = tcga.caf.score.group[,c("Type","CAF_4")],
                  leg = 'Groups',ylab = 'CAF_4 Score',palette = c('#377EB8','indianred'))
fig1e <- ggplot(tcga.nk.score.group, aes(x = Type, y = NK_4, fill = Type)) +
  geom_violin(trim = FALSE, alpha = 0.5) +  # 绘制小提琴图
  geom_boxplot(width = 0.1, color = "black", alpha = 0.5) +  # 添加箱线图
  labs(title = "NK_4 Score by Groups", x = "Groups", y = "NK_4 Score") +
  scale_fill_manual(values = c('#377EB8', 'indianred')) +
  theme_minimal() +
  stat_compare_means(comparisons = list(c("N", "T")),
                     label = "p.signif",map_signif_level = TRUE, 
                     textsize = 8,  # 设置显著性文本的字体大小
                     size = 8)   # 设置显著性线条的大小) +  # 使用显著性标记
theme(
  plot.title = element_text(size = 16, face = "bold"),  # 标题字体大小
  axis.title.x = element_text(size = 14),                # x轴标题字体大小
  axis.title.y = element_text(size = 14),                # y轴标题字体大小
  axis.text.x = element_text(size = 12),                 # x轴刻度字体大小
  axis.text.y = element_text(size = 12),                 # y轴刻度字体大小
  legend.title = element_text(size = 14),                 # 图例标题字体大小
  legend.text = element_text(size = 12)                   # 图例文本字体大小
)

fig1e

fig1<-ggarrange(fig1a,fig1b,fig1c,fig1d,fig1e,nrow = 2,ncol = 3,common.legend = T)
fig1
ggsave('results/Fig1.pdf',fig1,height = 7,width = 12)
library(ggrain)
fig1a+fig1b+fig1c+fig1d+fig1e+ plot_layout(guides = 'collect') &
  theme_bw()

#KM曲线
tcga.cli<-read.delim('/home/liangyu-kaikai/LIHC/GEO/geo_cli.txt',sep='\t',header = T)
tcga.nk.cli<-merge(tcga.cli,
                   data.frame(Sample=rownames(tcga.nk.score),tcga.nk.score),
                   by='Sample')

write.table(tcga.nk.cli,'results/tcga.nk.cli.txt',quote = F,row.names = F,sep='\t')
fig2a<-bioSurvival(OS = tcga.nk.cli$OS,
                   OS.time = tcga.nk.cli$OS.time/365,
                   riskscore = tcga.nk.cli$NK_0,labs = c('High','Low'),
                   palette = ggsci::pal_jco(alpha = 0.5)(9),leg = 'TCGA NK_0')
fig2a
fig2b<-bioSurvival(OS = tcga.nk.cli$OS,
                   OS.time = tcga.nk.cli$OS.time/365,
                   riskscore = tcga.nk.cli$NK_1,labs = c('High','Low'),
                   palette = ggsci::pal_jco(alpha = 0.5)(9),leg = 'TCGA NK_1')
fig2b
fig2c<-bioSurvival(OS = tcga.nk.cli$OS,
                   OS.time = tcga.nk.cli$OS.time/365,
                   riskscore = tcga.nk.cli$NK_2,labs = c('High','Low'),
                   palette = ggsci::pal_jco(alpha = 0.5)(9),leg = 'TCGA NK_2')
fig2c
fig2d<-bioSurvival(OS = tcga.nk.cli$OS,
                   OS.time = tcga.nk.cli$OS.time/365,
                   riskscore = tcga.nk.cli$NK_3,labs = c('High','Low'),
                   palette = ggsci::pal_jco(alpha = 0.5)(9),leg = 'TCGA NK_3')
fig2d
fig2e<-bioSurvival(OS = tcga.nk.cli$OS,
                   OS.time = tcga.nk.cli$OS.time/365,
                   riskscore = tcga.nk.cli$NK_4,labs = c('High','Low'),
                   palette = ggsci::pal_jco(alpha = 0.5)(9),leg = 'TCGA NK_4')
fig2e


fig2<-ggarrange(fig2a,fig2b,fig2c,fig2d,fig2e,nrow = 2,ncol = 3)
fig2
ggsave('results/tcga.NK.KM.pdf',fig2,height = 7,width = 12)

#临床特征表达的比较

library(ggplot2)
library(ggsci)
fig3a=list()


fig3a[[1]] <- ggplot(tcga.nk.cli, aes(x = TNM.staging, y = NK_0)) +
  geom_violin(fill = ggsci::pal_nejm(alpha = 0.8)(8)[1]) +
  geom_jitter(alpha = 0.5, width = 0.1) +
  labs(x = 'TNM.Stage', y = 'NK_0 Score') +
  stat_compare_means(comparisons = list(c("I", "II"), 
                                        c("II", "III"), 
                                        c("I", "III")), 
                     label = "p.signif",map_signif_level = TRUE,
                     textsize = 8,  # 设置显著性文本的字体大小
                     size = 8) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),  # 标题字体大小
    axis.title.x = element_text(size = 14),                # x轴标题字体大小
    axis.title.y = element_text(size = 14),                # y轴标题字体大小
    axis.text.x = element_text(size = 12),                 # x轴刻度字体大小
    axis.text.y = element_text(size = 12),                 # y轴刻度字体大小
    legend.title = element_text(size = 14),                 # 图例标题字体大小
    legend.text = element_text(size = 12)                   # 图例文本字体大小
  )
fig3a[[1]]

fig3a[[2]] <- ggplot(tcga.nk.cli, aes(x = TNM.staging, y = NK_1)) +
  geom_violin(fill = ggsci::pal_nejm(alpha = 0.8)(8)[1]) +
  geom_jitter(alpha = 0.5, width = 0.1) +
  labs(x = 'TNM.Stage', y = 'NK_1 Score') +
  stat_compare_means(comparisons = list(c("I", "II"), 
                                        c("II", "III"), 
                                        c("I", "III")), 
                     label = "p.signif",map_signif_level = TRUE,
                     textsize = 8,  # 设置显著性文本的字体大小
                     size = 8) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),  # 标题字体大小
    axis.title.x = element_text(size = 14),                # x轴标题字体大小
    axis.title.y = element_text(size = 14),                # y轴标题字体大小
    axis.text.x = element_text(size = 12),                 # x轴刻度字体大小
    axis.text.y = element_text(size = 12),                 # y轴刻度字体大小
    legend.title = element_text(size = 14),                 # 图例标题字体大小
    legend.text = element_text(size = 12)                   # 图例文本字体大小
  )
fig3a[[2]]

fig3a[[3]] <- ggplot(tcga.nk.cli, aes(x = TNM.staging, y = NK_2)) +
  geom_violin(fill = ggsci::pal_nejm(alpha = 0.8)(8)[1]) +
  geom_jitter(alpha = 0.5, width = 0.1) +
  labs(x = 'TNM.Stage', y = 'NK_2 Score') +
  stat_compare_means(comparisons = list(c("I", "II"), 
                                        c("II", "III"), 
                                        c("I", "III")), 
                     label = "p.signif",map_signif_level = TRUE,
                     textsize = 8,  # 设置显著性文本的字体大小
                     size = 8) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),  # 标题字体大小
    axis.title.x = element_text(size = 14),                # x轴标题字体大小
    axis.title.y = element_text(size = 14),                # y轴标题字体大小
    axis.text.x = element_text(size = 12),                 # x轴刻度字体大小
    axis.text.y = element_text(size = 12),                 # y轴刻度字体大小
    legend.title = element_text(size = 14),                 # 图例标题字体大小
    legend.text = element_text(size = 12)                   # 图例文本字体大小
  )
fig3a[[3]]

fig3a[[4]] <- ggplot(tcga.nk.cli, aes(x = TNM.staging, y = NK_3)) +
  geom_violin(fill = ggsci::pal_nejm(alpha = 0.8)(8)[1]) +
  geom_jitter(alpha = 0.5, width = 0.1) +
  labs(x = 'TNM.Stage', y = 'NK_3 Score') +
  stat_compare_means(comparisons = list(c("I", "II"), 
                                        c("II", "III"), 
                                        c("I", "III")), 
                     label = "p.signif",map_signif_level = TRUE,
                     textsize = 8,  # 设置显著性文本的字体大小
                     size = 8) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),  # 标题字体大小
    axis.title.x = element_text(size = 14),                # x轴标题字体大小
    axis.title.y = element_text(size = 14),                # y轴标题字体大小
    axis.text.x = element_text(size = 12),                 # x轴刻度字体大小
    axis.text.y = element_text(size = 12),                 # y轴刻度字体大小
    legend.title = element_text(size = 14),                 # 图例标题字体大小
    legend.text = element_text(size = 12)                   # 图例文本字体大小
  )
fig3a[[4]]

fig3a[[5]] <- ggplot(tcga.nk.cli, aes(x = TNM.staging, y = NK_4)) +
  geom_violin(fill = ggsci::pal_nejm(alpha = 0.8)(8)[1]) +
  geom_jitter(alpha = 0.5, width = 0.1) +
  labs(x = 'TNM.Stage', y = 'NK_4 Score') +
  stat_compare_means(comparisons = list(c("I", "II"), 
                                        c("II", "III"), 
                                        c("I", "III")), 
                     label = "p.signif",map_signif_level = TRUE,
                     textsize = 8,  # 设置显著性文本的字体大小
                     size = 8) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),  # 标题字体大小
    axis.title.x = element_text(size = 14),                # x轴标题字体大小
    axis.title.y = element_text(size = 14),                # y轴标题字体大小
    axis.text.x = element_text(size = 12),                 # x轴刻度字体大小
    axis.text.y = element_text(size = 12),                 # y轴刻度字体大小
    legend.title = element_text(size = 14),                 # 图例标题字体大小
    legend.text = element_text(size = 12)                   # 图例文本字体大小
  )
fig3a[[5]]

if (length(fig3a) >= 4) {
  fig3 <- ggarrange(fig3a[[1]], fig3a[[2]], fig3a[[3]], fig3a[[4]], fig3a[[5]], nrow = 2, ncol = 2)
} else {
  fig3 <- ggarrange(fig3a[[1]], fig3a[[2]], nrow = 1, ncol = 2)  # Example for 2 plots
}
fig3<-ggarrange(fig3a[[1]],fig3a[[2]],fig3a[[3]],fig3a[[4]],fig3a[[5]],nrow = 3,ncol = 2)
fig3

ggsave('results/tcga.NK.KM.pdf',fig3,height = 7,width = 12)


#
fig3b=list()
fig3b[[1]]<-sig_violin(dat = tcga.nk.cli[,c("T.Stage","NK_1")],
                       leg = 'T.Stage',ylab = 'NK_1 Score',
                       palette = ggsci::pal_nejm(alpha = 0.8)(8))
fig3b[[1]]
fig3b[[2]]<-sig_violin(dat = tcga.nk.cli[,c("N.Stage","NK_1")],
                       leg = 'N.Stage',ylab = 'NK_1 Score',
                       palette = ggsci::pal_nejm(alpha = 0.8)(8))
fig3b[[2]]

fig3b[[3]]<-sig_violin(dat = tcga.nk.cli[,c("Stage","NK_1")],
                       leg = 'Stage',ylab = 'NK_1 Score',
                       palette = ggsci::pal_nejm(alpha = 0.8)(8))
fig3b[[3]]
fig3b=ggarrange(plotlist = fig3b,ncol = 3,nrow = 1,widths = c(1,1,1))
fig3b
#
fig3c=list()
fig3c[[1]]<-sig_violin(dat = tcga.nk.cli[,c("T.Stage","NK_2")],
                       leg = 'T.Stage',ylab = 'NK_2 Score',
                       palette = ggsci::pal_nejm(alpha = 0.8)(8))
fig3c[[1]]
fig3c[[2]]<-sig_violin(dat = tcga.nk.cli[,c("N.Stage","NK_2")],
                       leg = 'N.Stage',ylab = 'CAF_2 Score',
                       palette = ggsci::pal_nejm(alpha = 0.8)(8))
fig3c[[2]]

fig3c[[3]]<-sig_violin(dat = tcga.nk.cli[,c("Stage","NK_2")],
                       leg = 'Stage',ylab = 'NK_2 Score',
                       palette = ggsci::pal_nejm(alpha = 0.8)(8))
fig3c[[3]]
fig3c=ggarrange(plotlist = fig3c,ncol = 3,nrow = 1,widths = 1)
fig3c
#
fig3d=list()
fig3d[[1]]<-sig_violin(dat = tcga.nk.cli[,c("T.Stage","NK_3")],
                       leg = 'T.Stage',ylab = 'NK_3 Score',
                       palette = ggsci::pal_nejm(alpha = 0.8)(8))
fig3d[[1]]
fig3d[[2]]<-sig_violin(dat = tcga.nk.cli[,c("N.Stage","NK_3")],
                       leg = 'N.Stage',ylab = 'NK_3 Score',
                       palette = ggsci::pal_nejm(alpha = 0.8)(8))
fig3d[[2]]

fig3d[[3]]<-sig_violin(dat = tcga.nk.cli[,c("Stage","NK_3")],
                       leg = 'Stage',ylab = 'NK_3 Score',
                       palette = ggsci::pal_nejm(alpha = 0.8)(8))
fig3d[[3]]
fig3d=ggarrange(plotlist = fig3d,ncol = 3,nrow = 1,widths = 1)
fig3d
#
fig3e=list()
fig3e[[1]]<-sig_violin(dat = tcga.nk.cli[,c("T.Stage","NK_4")],
                       leg = 'T.Stage',ylab = 'NK_4 Score',
                       palette = ggsci::pal_nejm(alpha = 0.8)(8))
fig3e[[1]]
fig3e[[2]]<-sig_violin(dat = tcga.nk.cli[,c("N.Stage","NK_4")],
                       leg = 'N.Stage',ylab = 'NK_4 Score',
                       palette = ggsci::pal_nejm(alpha = 0.8)(8))
fig3e[[2]]

fig3e[[3]]<-sig_violin(dat = tcga.nk.cli[,c("Stage","NK_4")],
                       leg = 'Stage',ylab = 'NK_4 Score',
                       palette = ggsci::pal_nejm(alpha = 0.8)(8))
fig3e[[3]]
fig3e=ggarrange(plotlist = fig3e,ncol = 3,nrow = 1,widths = 1)
fig3e
fig3f=list()
fig3f[[1]]<-sig_violin(dat = tcga.nk.cli[,c("T.Stage","NK_5")],
                       leg = 'T.Stage',ylab = 'NK_5 Score',
                       palette = ggsci::pal_nejm(alpha = 0.8)(8))
fig3f[[1]]
fig3f[[2]]<-sig_violin(dat = tcga.nk.cli[,c("N.Stage","NK_5")],
                       leg = 'N.Stage',ylab = 'NK_5 Score',
                       palette = ggsci::pal_nejm(alpha = 0.8)(8))
fig3f[[2]]

fig3f[[3]]<-sig_violin(dat = tcga.nk.cli[,c("Stage","NK_5")],
                       leg = 'Stage',ylab = 'NK_5 Score',
                       palette = ggsci::pal_nejm(alpha = 0.8)(8))
fig3f[[3]]
fig3f=ggarrange(plotlist = fig3f,ncol = 3,nrow = 1,widths = 1)
fig3f
fig3<-ggarrange(fig3a,fig3b,fig3c,fig3d,fig3e,fig3f,nrow = 3,ncol = 2,labels = '')
ggsave('results/tcga.caf.cli.pdf',fig3,height = 20,width = 25)



####05.CAF.gene.set.select####


setwd("/home/liangyu-kaikai/LIHC/GEO/05.CAF.gene.set.select/")

dir.create('results')
tcga.dat<-read.delim("/home/liangyu-kaikai/LIHC/GEO/GSE14520-1.txt",sep='\t',header = T,row.names = 1,check.names = F)
tcga.nk.score<-read.delim('/home/liangyu-kaikai/LIHC/GEO/04.CAF.select/results/tcga.nk.score.txt',sep='\t',header = T)
tcga.cli<-read.delim('/home/liangyu-kaikai/LIHC/GEO/geo_cli.txt',sep='\t',header = T)
#差异分析
limma_DEG=function(exp,group,ulab,dlab){
  library(limma)
  ind1=which(group==ulab)
  ind2=which(group==dlab)
  sml <- c(rep('A',length(ind1)),rep('B',length(ind2)))    # set group names
  eset=exp[,c(ind1,ind2)]
  fl <- as.factor(sml)
  
  design <- model.matrix(~fl+0)
  colnames(design) <- levels(fl)
  cont.matrix<-makeContrasts(contrasts='A-B',levels=design)
  #print(head(eset))
  fit<-lmFit (eset,design)
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit2)
  #print(sml)
  tT <- topTable(fit2, adjust="fdr", sort.by="B", number=nrow(eset))
  return(tT)
}
#火山图
volcano_plot=function(logfc,pvalue,cutFC=1,cutPvalue=0.05
                      ,colors=c('red','grey','blue')
                      ,ylab='-log10(FDR)',
                      leg='Group',
                      xlab='log2(FoldChange)'){
  library(ggplot2)
  diff.name=rep('None',length(logfc))
  diff.name[which(logfc>cutFC&pvalue<cutPvalue)]='Up'
  diff.name[which(logfc< -cutFC&pvalue<cutPvalue)]='Down'
  dat.input=data.frame(logFC=logfc,FDR=pvalue,change=diff.name)
  p1 <- ggplot(data = dat.input, 
               aes(x = logFC, 
                   y = -log10(FDR)))+
    geom_point(alpha=0.4, size=3.5, aes(color=change))+
    scale_color_manual(values=colors,limits = c("Down",'None', "Up"),name=leg)+
    geom_vline(xintercept=c(-cutFC,cutFC),lty=2,col="black",lwd=0.8)+
    geom_hline(yintercept = -log10(cutPvalue),lty=2,col="black",lwd=0.8)+
    ylab(ylab)+xlab(xlab)+
    theme_get()+
    theme(legend.background = element_rect(fill = NA, colour = NA))
  return(p1)
}
#差异分析
library(dplyr)
tcga.diff<-limma_DEG(exp = tcga.dat[,tcga.nk.score$Samples],
                     group = tcga.nk.score$Type,
                     ulab = 'T',dlab = 'N')
fc.fit=1;p.fit=0.05
# Set the threshold for significance
p_threshold <- 0.05  # Adjust if needed
fc_threshold <- 1.5   # Adjust this for fold change threshold

# Create a data frame with logFC and adj.P.Val
gene_data <- tcga.diff %>%
  mutate(Significant = adj.P.Val < p_threshold) %>%
  arrange(desc(logFC))

# Get top 5 upregulated genes
top_upregulated <- gene_data %>%
  filter(Significant & logFC > fc_threshold) %>%
  top_n(5, logFC)

# Get top 5 downregulated genes
top_downregulated <- gene_data %>%
  filter(Significant & logFC < -fc_threshold) %>%
  top_n(5, -logFC)  # Using negative logFC to get the top downregulated
# Combine the top genes
top_genes <- rbind(top_upregulated, top_downregulated)
top_genes$Gene <- rownames(top_genes)
tcga.diff$Gene <- rownames(tcga.diff)
fig1a <- volcano_plot(logfc = tcga.diff$logFC,
                      pvalue = tcga.diff$adj.P.Val,
                      cutFC = fc.fit,
                      cutPvalue = p.fit,
                      colors = c('blue','grey','red'),
                      ylab = '-log10(FDR)',
                      leg='GEO TvsN',
                      xlab='log2(FoldChange)') +
  geom_text(data = top_genes, aes(x = logFC, y = -log10(adj.P.Val),label = Gene), 
            vjust = -0.5, color = "black", size = 4)
fig1a


ggsave('results/volcano.plot.pdf',fig1a,height = 5,width = 6)
#
write.table(tcga.diff,'results/tcga.diff.txt',quote = F,row.names = T,sep='\t')
tcga.diff.fit<-tcga.diff[which(abs(tcga.diff$logFC)>fc.fit & tcga.diff$adj.P.Val<p.fit),]
table(tcga.diff.fit$logFC>0)
#FALSE  TRUE 
#1006   725
unique(tcga.nk.score$Type)
tcga.nk.score <- na.omit(tcga.nk.score)
sam_T<-tcga.nk.score[which(tcga.nk.score$Type=='T'),"Samples"]
tcga.gene.exp<-t(tcga.dat[rownames(tcga.diff.fit),sam_T])
tcga.nk.exp<-merge(data.frame(Samples=rownames(tcga.gene.exp),tcga.gene.exp),
                   tcga.nk.score[,c("Samples","NK_0","NK_1","NK_2","NK_3","NK_4")],
                   by='Samples')
write.table(tcga.nk.exp,'results/tcga.nk.exp.txt',quote = F,row.names = F,sep='\t')
#相关性分析
tcga.nk.exp=tcga.nk.exp[,-1]
tcga.nk.gene.cor<-Hmisc::rcorr(as.matrix(tcga.nk.exp))
tcga.nk.gene.cor.r<-reshape2::melt(tcga.nk.gene.cor$r)
tcga.nk.gene.cor.p<-reshape2::melt(tcga.nk.gene.cor$P)
colnames(tcga.nk.gene.cor.r)=c('nk','gene','cor')
colnames(tcga.nk.gene.cor.p)=c('nk','gene','p')
tcga.nk.gene.cor.r=tcga.nk.gene.cor.r[tcga.nk.gene.cor.r$nk %in% c("NK_0","NK_1","NK_2","NK_3","NK_4") & tcga.nk.gene.cor.r$gene %in% rownames(tcga.diff.fit),]
tcga.nk.gene.cor.p=tcga.nk.gene.cor.p[tcga.nk.gene.cor.p$nk %in% c("NK_0","NK_1","NK_2","NK_3","NK_4") & tcga.nk.gene.cor.p$gene %in% rownames(tcga.diff.fit),]
tcga.nk.gene.cor.res<-merge(tcga.nk.gene.cor.r,tcga.nk.gene.cor.p,by=c('nk','gene'))
head(tcga.nk.gene.cor.res)
tcga.nk.gene.cor.res.fit<-tcga.nk.gene.cor.res[abs(tcga.nk.gene.cor.res$cor)>0.4 & tcga.nk.gene.cor.res$p<0.05,]
dim(tcga.nk.gene.cor.res.fit)
#867
tcga.nk.gene.cor.res.fit$nk=as.character(tcga.nk.gene.cor.res.fit$nk)
tcga.nk.gene.cor.res.fit$gene=as.character(tcga.nk.gene.cor.res.fit$gene)
table(tcga.nk.gene.cor.res.fit$nk,tcga.nk.gene.cor.res.fit$cor>0)
#         neg pos
# CAF_1     4  196
# CAF_2   167  500
write.table(tcga.nk.gene.cor.res.fit,'results/tcga.nk.gene.cor.res.fit.txt',quote = F,row.names = F,sep='\t')
#富集分析
library("clusterProfiler")
library(AnnotationHub)
library(AnnotationDbi)
library(GOplot)
library(ggplot2)

module_genes=as.character(unique(tcga.nk.gene.cor.res.fit$gene))
gene = bitr(module_genes, #数据集
            fromType="SYMBOL", #输入为SYMBOL格式
            toType="ENTREZID",  # 转为ENTERZID格式
            OrgDb="org.Hs.eg.db") #人类 数据库
head(gene,2)

#GO富集分析
library(org.Hs.eg.db)
egoBP <- enrichGO(gene = gene$ENTREZID,
                  OrgDb = org.Hs.eg.db, 
                  pvalueCutoff =0.05, 
                  qvalueCutoff = 0.05,
                  ont="BP",
                  readable =T)
dotplot(egoBP)
egoMF <- enrichGO(gene = gene$ENTREZID,
                  OrgDb = org.Hs.eg.db, 
                  pvalueCutoff =0.05, 
                  qvalueCutoff = 0.05,
                  ont="MF",
                  readable =T)
dotplot(egoMF)

egoCC <- enrichGO(gene = gene$ENTREZID,
                  OrgDb = org.Hs.eg.db, 
                  pvalueCutoff =0.05, 
                  qvalueCutoff = 0.05,
                  ont="CC",
                  readable =T)
dotplot(egoCC)

egoALL <- enrichGO(gene = gene$ENTREZID,
                   OrgDb = org.Hs.eg.db, 
                   pvalueCutoff =0.05, 
                   qvalueCutoff = 0.05,
                   ont="ALL",
                   readable =T)
dotplot(egoALL)


#KEGG富集
R.utils::setOption("clusterProfiler.download.method",'auto') 
kk <- enrichKEGG(gene = gene$ENTREZID,
                 organism = 'hsa', #KEGG可以用organism = 'hsa'
                 pvalueCutoff = 0.05)
write.table(kk@result,file="results/kegg.txt",sep="\t",quote=F,row.names = F)

dotplot(kk)


####06. CAF moldue####
setwd('~/LIHC/GEO/06.nk.module/')
dir.create('results')

tcga.cli<-read.delim('C:\\Users\\ZZD\\Desktop\\套路七单细胞+成纤维+免疫治疗\\03.data.pre\\GSE53624\\results\\geo_cli.txt',sep='\t',header = T)
tcga.dat<-read.delim('~/LIHC/GEO/GSE14520-1.txt',sep='\t',header = T,row.names = 1,check.names = F)
tcga.nk.gene.cor.res.fit<-read.delim('C:\\Users\\ZZD\\Desktop\\套路七单细胞+成纤维+免疫治疗\\05.CAF.gene.set.select\\results\\tcga.nk.gene.cor.res.fit.txt',sep='\t',header = T)
library(ggsci)
mycols <- pal_npg('nrc',alpha = 0.6)(9)
mycols <- rep(mycols,14)
cox_batch<-function(dat,time,event){
  coxRun<-function(dat){
    library(survival)
    colnames(dat)=c('time','status','AS')  
    dat=dat[which(!is.na(dat[,1])&!is.na(dat[,3])&!is.na(dat[,2])),]
    #print(nrow(dat))
    if(nrow(dat)<10){
      print(paste0('Sample Num is small:',nrow(dat)))
      return(c(NA,NA,NA,NA))
    }
    #if(quantile(dat[,3])['25%']==quantile(dat[,3])['50%']) return(c(NA,NA,NA,NA))
    fmla <- as.formula("Surv(time, status) ~AS")
    if(table(dat[,2])[1]>1&table(dat[,2])[2]>1){
      cox <- survival::coxph(fmla, data = dat)
      re=c(summary(cox)[[7]][5],summary(cox)[[7]][2],summary(cox)[[8]][3],summary(cox)[[8]][4])
      return(re)
    }else{
      return(c(NA,NA,NA,NA))
    }
  }
  t.inds=which(!is.na(time)&!is.na(event))
  dat1=dat[,t.inds]
  os=time[t.inds]
  ev=event[t.inds]
  
  ct=sum(ev%in%c(0,1))
  if(ct!=length(ev)){
    print('event must be 0(alive) or 1(dead)')
    return(NULL)
  }
  
  res=t(apply(dat1, 1, function(x){
    ep=as.numeric(as.character(x))
    ind2=which(!is.na(ep))
    print(length(ind2))
    if(length(ind2)>1){
      os1=os[ind2]
      ev1=ev[ind2]
      ep1=ep[ind2]
      return(coxRun(data.frame(os1,ev1,ep1)))
    }else{
      return(c(NA,NA,NA,NA))
    }
  }))
  colnames(res)=c('p.value','HR','Low 95%CI','High 95%CI')
  row.names(res)=row.names(dat1)
  return(as.data.frame(res))
}

tcga.cox<-cox_batch(dat = tcga.dat[as.character(unique(tcga.nk.gene.cor.res.fit$gene)),tcga.cli$Sample],
                    time = tcga.cli$OS.time,event = tcga.cli$OS)
head(tcga.dat)
colnames(tcga.dat)
rownames(tcga.dat)
# 获取唯一基因名称
unique_genes <- as.character(unique(tcga.nk.gene.cor.res.fit$gene))

# 找出存在的基因名称
existing_genes <- intersect(unique_genes, rownames(tcga.dat))

# 列出不存在的基因名称
missing_genes <- setdiff(unique_genes, rownames(tcga.dat))
# 查看 tcga.cli$Sample 的内容
unique_samples <- unique(tcga.cli$Sample)
print(unique_samples)

# 确认这些样本是否在 tcga.dat 的列名中存在
missing_samples <- setdiff(unique_samples, colnames(tcga.dat))

# 确保选择的基因和样本都是存在的
valid_genes <- intersect(unique_genes, rownames(tcga.dat))
valid_samples <- intersect(tcga.cli$Sample, colnames(tcga.dat))

# 进行 Cox 分析
tcga.cox <- cox_batch(dat = tcga.dat[valid_genes, valid_samples],
                      time = tcga.cli$OS.time[tcga.cli$Sample %in% valid_samples],
                      event = tcga.cli$OS[tcga.cli$Sample %in% valid_samples])

cox.p<-0.01#可修改
table(tcga.cox$p.value<cox.p)
# FALSE  TRUE 
# 639   26 
tcga.cox.fit<-tcga.cox[which(tcga.cox$p.value<cox.p),]
dim(tcga.cox.fit)#219
dim(tcga.cox)
tcga.cox$coef=log(tcga.cox$HR)
tcga.cox$Gene=rownames(tcga.cox)
tcga.cox$type=rep('None',nrow(tcga.cox))
tcga.cox$type[which(tcga.cox$p.value<cox.p & tcga.cox$coef>0)]='Risk'
tcga.cox$type[which(tcga.cox$p.value<cox.p & tcga.cox$coef<0)]='Protective'
table(tcga.cox$type)
# None Protective       Risk 
# 639        8        18 

write.table(tcga.cox,'results/tcga.cox.txt',quote = F,row.names = T,sep='\t')

num_na_needed <- nrow(tcga.cox) - length(names(tcga_risk$module.gene))
tcga.cox$label <- c(names(tcga_risk$module.gene), rep(NA, num_na_needed))
#tcga.cox$label <- c(names(tcga_risk$module.gene),rep(NA,nrow(tcga.cox)-9))

#devtools::install_github("BioSenior/ggVolcano")
library(ggVolcano)
library(latex2exp)
library(ggrepel)
library(ggplot2)
fig1a <- ggplot(data = tcga.cox,
                aes(x = coef,
                    y = -log10(p.value)))+
  geom_point(alpha=0.4, size=3.5, aes(color=type))+
  scale_color_manual(values=c('#4DAF4A','grey','#E41A1C'),
                     limits = c("Protective",'None', "Risk"),name='State')+
  geom_hline(yintercept = -log10(cox.p),lty=4,col="black",lwd=0.8)+
  ylab('-log10(pvalue)')+xlab('Cox coefficient')+
  theme_bw()+
  theme(axis.text.y=element_text(family="Times",face="plain"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.y=element_text(family="Times",face="plain"), #设置y轴标题的字体属性
        legend.text=element_text(face="plain", family="Times", colour="black"),  #设置图例的子标题的字体属性
        legend.title=element_text(face="plain", family="Times", colour="black" ),#设置图例的总标题的字体属性
        legend.justification=c(1,1), legend.position=c(1,1),
        legend.background = element_rect(fill = NA, colour = NA)
  )
data_selected <- tcga.cox[names(tcga_risk$module.gene),]
fig <- fig1a + geom_label_repel(data=data_selected,
                                aes(label=rownames(data_selected)))

fig
ggsave('results/tcga.sig.cox.pdf',fig,height = 5,width = 5)


#
#lasso进一步压缩
get_riskscore.lasso<-function(dat,os,os.time,labels=c('A','B')){
  crbind2DataFrame=function(dat){
    print(class(dat))
    if(class(dat)=='table'){
      if(!is.na(ncol(dat))){
        dat=apply(dat,2,function(x){
          return(x)
        })
      }
    }
    if(class(dat)!='data.frame'){
      dat1=as.data.frame(as.matrix(dat))
    }else{
      dat1=dat
    }
    #print(head(dat1))
    for(i in 1:ncol(dat1)){
      dat1[,i]=as.character(dat1[,i])
      dt=dat1[which(gsub(' ','',dat1[,i])!=''&!is.na(dat1[,i])),i]
      dt=dt[which(dt!='Inf'&dt!='NaN'&dt!='NA')]
      if(sum(is.na(as.numeric(dt)))==0){
        #print(dat1[,i])
        dat1[,i]=as.numeric(dat1[,i])
      }
    }
    return(dat1)  
  }
  library(glmnet)
  set.seed(2021)
  options(ggrepel.max.overlaps=Inf)
  fit1=glmnet(as.matrix(dat)
              ,cbind(time=os.time,
                     status=os)
              ,family="cox"
              ,nlambda=100
              , alpha=1) 
  
  cv.fit<-cv.glmnet(as.matrix(dat)
                    ,cbind(time=os.time,
                           status=os)
                    ,family="cox"
                    ,nfolds = 10,
                    ,nlambda=100
                    , alpha=1)
  sig.coef <- coefficients(cv.fit,s=cv.fit$lambda.min)[which(coefficients(cv.fit,s=cv.fit$lambda.min)[,1]!=0),1]
  #print(cv.fit$lambda.min)
  #length(names(sig.coef))
  #10
  mg_plot_lasso <- function(fit,cv_fit,lambda=NULL,show_text=T,figLabels=c('A','B')){
    if(is.null(lambda)){
      lmda=cv_fit$lambda.min
    }else{
      lmda=lambda
    }
    fit.coef=fit$beta[(apply(fit$beta,1,function(x){
      return(sum(x!=0))
    })>0),]
    
    fit.coef=as.matrix(fit.coef)
    colnames(fit.coef)=fit$lambda
    #fit$lambda==cv_fit$lambda
    library(ggplot2)
    dat=data.table::melt(t(as.matrix(fit.coef)))
    dat_z=dat[which(dat$value==0),]
    dat=dat[which(dat$value!=0),]
    dat.sv=rbind()
    for (u in unique(dat_z[,2])) {
      t.z=dat_z[which(dat_z[,2]==u),1]
      t.zx=max(t.z)
      dat.sv=rbind(dat.sv,c(t.zx,u,0))
      t.zn=min(t.z)
      if(t.zx!=t.zn){
        dat.sv=rbind(dat.sv,c(t.zn,u,0))
      }
    }
    colnames(dat.sv)=colnames(dat_z)
    #dat_z=dat_z[dat_z[,2]%in%names(which(fit.coef[,which(fit$lambda==lmda)]!=0)),]
    dat=crbind2DataFrame(rbind(dat,dat.sv))
    mn=min(-log(dat$Var1))
    mx=max(-log(dat$Var1))
    if(show_text){
      mx=(mx-mn)*0.1+mx
    }
    p=ggplot(dat, aes(x=-log(Var1), y=value,colour=Var2))+geom_line()+theme_bw()+theme(legend.position = "none")
    p=p+coord_cartesian(xlim=c(mn, mx))+xlab('-ln(lambda)')+ylab('Coefficients')
    if(show_text){
      fl=fit.coef[which(fit.coef[,which(fit$lambda==lmda)]!=0),ncol(fit.coef)]
      for_label=data.frame(Var1=rep(min(dat$Var1),length(fl)),Var2=names(fl),value=fl)
      p=p+ggrepel::geom_label_repel(
        aes(label = Var2,color=Var2),
        data = for_label,hjust = 0
      )
    }
    p=p+geom_vline(aes(xintercept=-log(lmda)), colour="#BB0000", linetype="dashed")
    p=p+annotate('text',x=-log(lmda),y=min(dat[,3]),label=paste0('lambda=',round(lmda,4)))
    tgc=data.frame(lambda=cv_fit$lambda,cvm=cv_fit$cvm,cvup=cv_fit$cvup,cvlo=cv_fit$cvlo,cvsd=cv_fit$cvsd
                   ,col=ifelse(cv_fit$lambda>=cv_fit$lambda.min&cv_fit$lambda<=cv_fit$lambda.1se,ifelse(cv_fit$lambda==lmda,'A','C'),'B'))
    p1=ggplot(tgc, aes(x=log(lambda), y=cvm)) + xlab('ln(lambda)')+ ylab('Parial Likelihood Deviance')+
      geom_errorbar(aes(ymin=cvm-cvsd, ymax=cvm+cvsd)) +
      geom_point(aes(colour=col))
    p1=p1+theme_bw()+theme(legend.position = "none")
    gal=ggpubr::ggarrange(p,p1, ncol = 2, nrow = 1
                          #,align = "hv"
                          ,labels = figLabels)
    return(gal)
  }
  lasso.pdf <- mg_plot_lasso(fit1,
                             cv.fit,
                             show_text=T,
                             figLabels=labels)
  return(list(lasso.gene=sig.coef,lambda.min=cv.fit$lambda.min,plot=lasso.pdf))
}
get_riskscore<-function(dat,os,os.time,step=T,direction=c("both", "backward", "forward")[1]){
  crbind2DataFrame=function(dat){
    print(class(dat))
    if(class(dat)=='table'){
      if(!is.na(ncol(dat))){
        dat=apply(dat,2,function(x){
          return(x)
        })
      }
    }
    if(class(dat)!='data.frame'){
      dat1=as.data.frame(as.matrix(dat))
    }else{
      dat1=dat
    }
    #print(head(dat1))
    for(i in 1:ncol(dat1)){
      dat1[,i]=as.character(dat1[,i])
      dt=dat1[which(gsub(' ','',dat1[,i])!=''&!is.na(dat1[,i])),i]
      dt=dt[which(dt!='Inf'&dt!='NaN'&dt!='NA')]
      if(sum(is.na(as.numeric(dt)))==0){
        #print(dat1[,i])
        dat1[,i]=as.numeric(dat1[,i])
      }
    }
    return(dat1)  
  }
  
  tcga_dat1 <- cbind(time=os.time,
                     status=os,
                     dat)
  tcga_dat1=crbind2DataFrame(tcga_dat1)
  colnames(tcga_dat1)=gsub('-','__',colnames(tcga_dat1))
  gene111=gsub('-','__',colnames(dat))
  fmla <- as.formula(paste0("Surv(time, status) ~"
                            ,paste0(gene111,collapse = '+')))
  cox <- coxph(fmla, data =as.data.frame(tcga_dat1))
  if(step==T){
    cox1 <- step(cox,direction =direction)
  }else{
    cox1=cox
  }
  lan <- coef(cox1)
  #round(lan, 3)
  genes <- names(cox1$coefficients)
  mult_results=paste0(round(lan, 3), '*', names(lan),collapse = '+')
  risk.tcga <- as.numeric(lan%*%as.matrix(t(tcga_dat1[,genes])))
  
  data_gene_score_final<-tcga_dat1
  data_gene_score_final$Samples<-rownames(data_gene_score_final)
  data_gene_score_final$riskscore=risk.tcga
  data_gene_score_final$riskscorez=mosaic::zscore(risk.tcga)
  optimalCutoff <- survminer::surv_cutpoint(data.frame(time=data_gene_score_final$time/365,
                                                       event=data_gene_score_final$status,
                                                       risk=data_gene_score_final$riskscore), 
                                            time = "time", event = "event",variables = c("risk"))
  optimalCutoff=optimalCutoff$cutpoint$cutpoint[1]
  #print(optimalCutoff)
  #optimalCutoff=median(data_gene_score_final$riskscore)
  #optimalCutoff=0
  data_gene_score_final$Risk=ifelse(data_gene_score_final$riskscore>optimalCutoff,'High','Low')
  table(data_gene_score_final$Risk)
  data_gene_score_final$cutoff=optimalCutoff
  return(list(result=data_gene_score_final,module.gene=cox1$coefficients,model=mult_results))
}

colnames(tcga.dat) 
unique(tcga.cli$Sample)
# 检查哪个样本名不存在
missing_samples <- tcga.cli$Sample[!tcga.cli$Sample %in% colnames(tcga.dat)]
if (length(missing_samples) > 0) {
  print(paste("以下样本名未在 tcga.dat 中找到:", paste(missing_samples, collapse = ", ")))
}
# 检查列名是否存在
valid_samples <- tcga.cli$Sample[tcga.cli$Sample %in% colnames(tcga.dat)]
if (length(valid_samples) == 0) {
  stop("没有有效的样本名可以用于选择 tcga.dat 中的列。")
}

# 进行数据子集选择
tcga.module.dat <- tcga.dat[rownames(tcga.cox.fit), valid_samples, drop = FALSE]

rownames(tcga.module.dat)=gsub('-','__',rownames(tcga.module.dat))

# 检查哪些样本名在 tcga.module.dat 中不存在
missing_samples <- tcga.cli$Sample[!tcga.cli$Sample %in% colnames(tcga.module.dat)]
if (length(missing_samples) > 0) {
  print(paste("以下样本名在 tcga.module.dat 中未找到:", paste(missing_samples, collapse = ", ")))
}
# 选择存在的样本名
valid_samples <- tcga.cli$Sample[tcga.cli$Sample %in% colnames(tcga.module.dat)]

# 确保有效样本的存在
if (length(valid_samples) == 0) {
  stop("没有有效的样本名可以用于选择 tcga.module.dat 中的列。")
}

# 进行数据子集选择
tcga.module.dat.valid <- tcga.module.dat[, valid_samples, drop = FALSE]
# 转置数据并调用函数
tcga.lasso <- get_riskscore.lasso(dat = t(tcga.module.dat.valid),
                                  os = tcga.cli$OS,
                                  os.time = tcga.cli$OS.time,
                                  labels = '')

tcga.lasso<-get_riskscore.lasso(dat = t(tcga.module.dat[,geo_cli$Sample]),
                                os=geo_cli$OS,
                                os.time = geo_cli$OS.time,
                                labels='')

length(tcga.lasso$lasso.gene) 
tcga.lasso$plot
ggsave('results/lasso.pdf',tcga.lasso$plot,height = 5,width = 10)
tcga.lasso$lambda.min
#0.06533129
library(survival)
library(mosaic)
tcga_risk<-get_riskscore(dat=as.data.frame(t(tcga.module.dat[names(tcga.lasso$lasso.gene),geo_cli$Sample])),
                         os=geo_cli$OS,
                         os.time=geo_cli$OS.time,
                         step=F,#不进行逐步回归
                         direction=c("both", "backward", "forward")[1])

length(tcga_risk$module.gene) #9
names(tcga_risk$module.gene)

#"ANGPTL7" "C6"      "CSRP1"   "EXPH5"   "F2RL2"   "KCNMA1"  "MAGEC3"  "MAMDC2"  "SLC4A9" 
tcga_risk$model
#"0.093*ANGPTL7+0.15*C6+0.121*CSRP1+-0.08*EXPH5+0.12*F2RL2+0.014*KCNMA1+-0.373*MAGEC3+0.143*MAMDC2+-0.188*SLC4A9"

#柱状图
gene.coef=as.data.frame(tcga_risk$module.gene)
gene.coef=data.frame(Gene=rownames(gene.coef),Coef=gene.coef[,1])
gene.coef$Type=ifelse(gene.coef$Coef>0,'Risk','Protective')
gene.coef$Type=factor(gene.coef$Type,levels=c('Risk','Protective'))
library(dplyr)
fig1d=gene.coef %>% 
  ggplot(aes(reorder(Gene, Coef), Coef)) +
  geom_col(aes(fill = Type)) +
  coord_flip() +
  scale_fill_manual(values=ggsci::pal_lancet('lanonc',alpha =0.6)(9)[c(7,1)]) +
  coord_flip() +
  labs(x = "") +
  labs(y = "Cox coefficient") +
  theme_bw()+
  theme(axis.text.y = element_text(angle = 0, hjust = 1),legend.position="top")

fig1d
ggsave('results/fig1d.pdf',fig1d,height = 5,width = 5)
#AUC和ROC
risk_plot<-function(time,event,riskscore,
                    group,mk,labs=c('High','Low'),
                    palette=brewer.pal(3,'Set2')){
  dat=data.frame(time=time,status=event,riskscore=riskscore,group=group)
  #ROC曲线
  ROC_rt=timeROC::timeROC(T=dat$time, 
                          delta=dat$status,
                          marker=dat$riskscore, cause=1,
                          weighting='marginal',
                          times=mk, 
                          ROC=TRUE,iid = T)
  p.dat=data.frame()
  for(i in which(ROC_rt$times>0)){
    lbs=paste0('AUC at ',mk[i],' years: ',round(ROC_rt$AUC[i],2))
    p.dat=rbind.data.frame(p.dat,data.frame(V1=ROC_rt$FP[,i],V2=ROC_rt$TP[,i],Type=lbs))
  }
  #colnames(p.dat)=c('V1','V2','Type')
  p.dat=as.data.frame(p.dat)
  p.dat$Type=as.factor(p.dat$Type)
  roc_plot=ggplot(p.dat, aes(x=V1,y=V2, fill=Type))+
    stat_smooth(aes(colour=Type),se = FALSE, size = 1)+
    theme_bw()+xlab('False positive fraction')+
    ylab('True positive fraction') +
    theme(axis.text.y=element_text(family="Times",face="plain"),
          axis.text.x=element_text(family="Times",face="plain"),
          axis.title.x=element_text(family="Times",face="plain"),
          axis.title.y=element_text(family="Times",face="plain"),
          plot.title=element_blank(),
          plot.margin=unit(c(0.2, 0.2, 0.2, 0.2), "inches"),
          legend.position=c(1,0),
          legend.justification=c(1,0),
          legend.background = element_rect(fill = NA, colour = NA),
          legend.title = element_text(family="Times",face="plain"),
          legend.text = element_text(family="Times",face="plain"))
  #KM曲线
  library(ggplot2)
  library(survival)
  dat1=dat[,c("time","status","group")]
  colnames(dat1)=c('time','status','groups')
  sf<-survival::survfit(Surv(time,status) ~ groups,data=dat1)
  surv=survminer::ggsurvplot(sf, data = dat1, 
                             palette = c("indianred1","#377eb8" ), 
                             pval = TRUE,
                             surv.median.line='hv'
                             #,conf.int = T
                             ,conf.int.style ='step'
                             , pval.coord=c(0, 0.2), #Add p-value 
                             risk.table = TRUE, 
                             linetype = 'solid',
                             legend.title = 'Group'
                             ,legend.labs =labs
                             ,conf.int=T
  )
  p1=surv$plot+theme_bw()+
    theme(axis.text.y=element_text(family="Times",face="plain"),
          axis.text.x=element_blank(),
          axis.title.x=element_blank(),
          plot.margin=unit(c(0.2, 0.2, 0, 0.1), "inches"),
          legend.position=c(1,1),
          legend.justification=c(1,1),
          legend.background = element_rect(fill = NA, colour = NA),
          legend.title = element_text(family="Times",face="plain"),
          legend.text = element_text(family="Times",face="plain"))
  p2=surv$table+theme_bw()+
    theme(axis.text.y=element_text(family="Times",face="plain"),
          plot.margin=unit(c(0, 0.2, 0.2, 0.1), "inches"),
          plot.title=element_blank(),
          legend.position=c(1,1), 
          legend.justification=c(1,1),
          legend.title = element_text(family="Times",face="plain"),
          legend.text = element_text(family="Times",face="plain"))
  
  g2=ggpubr::ggarrange(p1,p2, ncol = 1, nrow = 2,heights = c(1,0.3),align = "v")
  gg<-ggpubr::ggarrange(g2,roc_plot,ncol = 2,nrow = 1,labels = '')
  return(gg)
}

#risk_plot <- function(time, event, riskscore, group, mk, labs = c('High', 'Low'), 
palette = brewer.pal(3, 'Set2')) 
{
  # Create data frame
  dat <- data.frame(time = time, status = event, riskscore = riskscore, group = group)
  
  # Generate ROC curve data
  ROC_rt <- timeROC::timeROC(T = dat$time, 
                             delta = dat$status,
                             marker = dat$riskscore, 
                             cause = 1,
                             weighting = 'marginal',
                             times = mk, 
                             ROC = TRUE, 
                             iid = TRUE)
  
  # Create data frame for ROC plot
  p.dat <- data.frame()
  for (i in which(ROC_rt$times > 0)) {
    lbs <- paste0('AUC at ', mk[i], ' years: ', round(ROC_rt$AUC[i], 2))
    p.dat <- rbind.data.frame(p.dat, data.frame(V1 = ROC_rt$FP[, i], V2 = ROC_rt$TP[, i], Type = lbs))
  }
  
  p.dat$Type <- as.factor(p.dat$Type)
  
  # ROC plot
  roc_plot <- ggplot(p.dat, aes(x = V1, y = V2, fill = Type)) +
    stat_smooth(aes(colour = Type), se = FALSE, size = 1) +
    theme_bw() +
    xlab('False Positive Fraction') +
    ylab('True Positive Fraction') +
    theme(axis.text.y = element_text(family = "Times", face = "plain"),
          axis.text.x = element_text(family = "Times", face = "plain"),
          axis.title.x = element_text(family = "Times", face = "plain"),
          axis.title.y = element_text(family = "Times", face = "plain"),
          plot.title = element_blank(),
          plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "inches"),
          legend.position = c(1, 0),
          legend.justification = c(1, 0),
          legend.background = element_rect(fill = NA, colour = NA),
          legend.title = element_text(family = "Times", face = "plain"),
          legend.text = element_text(family = "Times", face = "plain"))
  
  # Survival analysis
  dat1 <- dat[, c("time", "status", "group")]
  colnames(dat1) <- c('time', 'status', 'groups')
  sf <- survival::survfit(Surv(time, status) ~ groups, data = dat1)
  
  surv <- survminer::ggsurvplot(sf, data = dat1, 
                                palette = c("indianred1", "#377eb8"), 
                                pval = TRUE,
                                surv.median.line = 'hv',
                                pval.coord = c(0, 0.2), 
                                risk.table = TRUE, 
                                linetype = 'solid',
                                legend.title = 'Group',
                                legend.labs = labs,
                                conf.int = TRUE)
  
  # Customize survival plot
  p1 <- surv$plot + theme_bw() +
    theme(axis.text.y = element_text(family = "Times", face = "plain"),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          plot.margin = unit(c(0.2, 0.2, 0, 0.1), "inches"),
          legend.position = c(1, 1),
          legend.justification = c(1, 1),
          legend.background = element_rect(fill = NA, colour = NA),
          legend.title = element_text(family = "Times", face = "plain"),
          legend.text = element_text(family = "Times", face = "plain"))
  
  p2 <- surv$table + theme_bw() +
    theme(axis.text.y = element_text(family = "Times", face = "plain"),
          plot.margin = unit(c(0, 0.2, 0.2, 0.1), "inches"),
          plot.title = element_blank(),
          legend.position = c(1, 1), 
          legend.justification = c(1, 1),
          legend.title = element_text(family = "Times", face = "plain"),
          legend.text = element_text(family = "Times", face = "plain"))
  
  # Combine survival plot and risk table
  g2 <- ggpubr::ggarrange(p1, p2, ncol = 1, nrow = 2, heights = c(1, 0.3), align = "v")
  
  # Combine ROC plot and survival plots into one final figure
  final_plot <- ggpubr::ggarrange(g2, roc_plot, ncol = 2, nrow = 1, labels = '')
  
  return(final_plot)
}



#开始分析
library(ggplot2)
library(pROC)
library(ggplot2)

tcga_risk$result$Risk=ifelse(tcga_risk$result$riskscorez>0,'High','Low')
tcga_roc_km<-risk_plot(time=tcga_risk$result$time/365,
                       event=tcga_risk$result$status,
                       riskscore=tcga_risk$result$riskscorez,
                       group=tcga_risk$result$Risk,
                       mk=c(1,3,5),labs=c('High','Low'),
                       palette=brewer.pal(3,'Set2'))
tcga_roc_km
ggsave('results/gse53624.auc.km.pdf',tcga_roc_km,height = 5,width = 10)


#gse31210
gse31210.dat<-read.delim('/home/liangyu-kaikai/LIHC/GEO/06.nk.module/TCGA.txt',sep='\t',header = T,row.names = 1)
gse31210.cli<-read.delim('/home/liangyu-kaikai/LIHC/GEO/06.nk.module/tcga_cli.txt',sep='\t',header = T,row.names = 1)
gse31210.dat <- as.data.frame(t(gse31210.dat))
colnames(gse31210.dat) <- rownames(gse31210.cli)
gse31210_dat_m<-cbind.data.frame(OS.time=gse31210.cli$OS.time,
                                 OS=gse31210.cli$OS,
                                 t(gse31210.dat[intersect(names(tcga_risk$module.gene),rownames(gse31210.dat)),]))
gse31210_dat_m[1:4,1:4]
gse31210_risk<-get_riskscore(dat=gse31210_dat_m[,-c(1,2)],
                             os=gse31210_dat_m$OS,
                             os.time=gse31210_dat_m$OS.time,
                             step=F,#不进行逐步回归
                             direction=c("both", "backward", "forward")[1])
gse31210_risk$result$Risk=ifelse(gse31210_risk$result$riskscorez>0,'High','Low')
gse31210_roc_km<-risk_plot(time=gse31210_risk$result$time/365,
                           event=gse31210_risk$result$status,
                           riskscore=gse31210_risk$result$riskscorez,
                           group=gse31210_risk$result$Risk,
                           mk=c(1,2,3),labs=c('High','Low'),
                           palette=ggsci::pal_lancet()(10)[c(2,1)])
gse31210_roc_km
ggsave('results/tcga.auc.km.pdf',gse31210_roc_km,height = 5,width = 10)
#
gse3141.dat<-read.delim('../03.data.pre/GSE3141/results/gse3141.dat.txt',sep='\t',header = T,row.names = 1)
gse3141.cli<-read.delim('../03.data.pre/GSE3141/results/gse3141.cli.txt',sep='\t',header = T)
gse3141_dat_m<-cbind.data.frame(OS.time=gse3141.cli$OS.time,
                                OS=gse3141.cli$OS,
                                t(gse3141.dat[intersect(names(tcga_risk$module.gene),rownames(gse3141.dat)),gse3141.cli$Accession]))
gse3141_dat_m[1:4,1:4]

gse3141_risk<-get_riskscore(dat=gse3141_dat_m[,-c(1,2)],
                            os=gse3141_dat_m$OS,
                            os.time=gse3141_dat_m$OS.time,
                            step=F,#不进行逐步回归
                            direction=c("both", "backward", "forward")[1])
gse3141_risk$result$Risk=ifelse(gse3141_risk$result$riskscorez>0,'High','Low')
gse3141_roc_km<-risk_plot(time=gse3141_risk$result$time/365,
                          event=gse3141_risk$result$status,
                          riskscore=gse3141_risk$result$riskscorez,
                          group=gse3141_risk$result$Risk,
                          mk=c(1,3,5),labs=c('High','Low'),
                          palette=ggsci::pal_lancet()(10)[c(2,1)])
gse3141_roc_km
ggsave('results/gse3141.auc.km.pdf',gse3141_roc_km,height = 5,width = 10)
#GSE37745
gse37745.dat<-read.delim('../03.data.pre/GSE37745/results/GSE37745.exp.txt',sep='\t',header = T,row.names = 1)
gse37745.cli<-read.delim('../03.data.pre/GSE37745/results/GSE37745.cli.txt',sep='\t',header = T)
gse37745_dat_m<-cbind.data.frame(OS.time=gse37745.cli$OS.time,
                                 OS=gse37745.cli$OS,
                                 t(gse37745.dat[intersect(names(tcga_risk$module.gene),rownames(gse37745.dat)),gse37745.cli$Samples]))
gse37745_dat_m[1:4,1:4]

gse37745_risk<-get_riskscore(dat=gse37745_dat_m[,-c(1,2)],
                             os=gse37745_dat_m$OS,
                             os.time=gse37745_dat_m$OS.time,
                             step=F,#不进行逐步回归
                             direction=c("both", "backward", "forward")[1])
gse37745_risk$result$Risk=ifelse(gse37745_risk$result$riskscorez>0,'High','Low')
gse37745_roc_km<-risk_plot(time=gse37745_risk$result$time/365,
                           event=gse37745_risk$result$status,
                           riskscore=gse37745_risk$result$riskscorez,
                           group=gse37745_risk$result$Risk,
                           mk=c(1,3,5),labs=c('High','Low'),
                           palette=ggsci::pal_lancet()(10)[c(2,1)])
gse37745_roc_km
ggsave('results/gse37745.auc.km.pdf',gse37745_roc_km,height = 5,width = 10)
#GSE50081
gse50081.dat<-read.delim('../03.data.pre/GSE50081//results/GSE50081.exp.txt',sep='\t',header = T,row.names = 1)
gse50081.cli<-read.delim('../03.data.pre/GSE50081/results/GSE50081.cli.txt',sep='\t',header = T)
gse50081_dat_m<-cbind.data.frame(OS.time=gse50081.cli$OS.time,
                                 OS=gse50081.cli$OS,
                                 t(gse50081.dat[intersect(names(tcga_risk$module.gene),rownames(gse50081.dat)),gse50081.cli$Samples]))
gse50081_dat_m[1:4,1:4]

gse50081_risk<-get_riskscore(dat=gse50081_dat_m[,-c(1,2)],
                             os=gse50081_dat_m$OS,
                             os.time=gse50081_dat_m$OS.time,
                             step=F,#不进行逐步回归
                             direction=c("both", "backward", "forward")[1])
gse50081_risk$result$Risk=ifelse(gse50081_risk$result$riskscorez>0,'High','Low')
gse50081_roc_km<-risk_plot(time=gse50081_risk$result$time,
                           event=gse50081_risk$result$status,
                           riskscore=gse50081_risk$result$riskscorez,
                           group=gse50081_risk$result$Risk,
                           mk=c(1,3,5),labs=c('High','Low'),
                           palette=ggsci::pal_lancet()(10)[c(2,1)])
gse50081_roc_km
ggsave('results/gse50081.auc.km.pdf',gse50081_roc_km,height = 5,width = 10)
#GSE68465
gse68465.dat<-read.delim('../03.data.pre/GSE68465/results/GSE68465.exp.txt',sep='\t',header = T,row.names = 1)
gse68465.cli<-read.delim('../03.data.pre/GSE68465/results/GSE68465.cli.txt',sep='\t',header = T)
gse68465_dat_m<-cbind.data.frame(OS.time=gse68465.cli$OS.time,
                                 OS=gse68465.cli$OS,
                                 t(gse68465.dat[intersect(names(tcga_risk$module.gene),rownames(gse68465.dat)),gse68465.cli$Samples]))
gse68465_dat_m[1:4,1:4]

gse68465_risk<-get_riskscore(dat=gse68465_dat_m[,-c(1,2)],
                             os=gse68465_dat_m$OS,
                             os.time=gse68465_dat_m$OS.time,
                             step=F,#不进行逐步回归
                             direction=c("both", "backward", "forward")[1])
gse68465_risk$result$Risk=ifelse(gse68465_risk$result$riskscorez>0,'High','Low')
gse68465_roc_km<-risk_plot(time=gse68465_risk$result$time/365,
                           event=gse68465_risk$result$status,
                           riskscore=gse68465_risk$result$riskscorez,
                           group=gse68465_risk$result$Risk,
                           mk=c(1,3,5),labs=c('High','Low'),
                           palette=ggsci::pal_lancet()(10)[c(2,1)])
gse68465_roc_km
ggsave('results/gse68465.auc.km.pdf',gse68465_roc_km,height = 5,width = 10)
tcga_risk1=tcga_risk$result
gse31210_risk1=gse31210_risk$result
gse3141_risk1=gse3141_risk$result
gse37745_risk1=gse37745_risk$result
gse50081_risk1=gse50081_risk$result
gse68465_risk1=gse68465_risk$result

write.table(tcga_risk,'results/tcga.risk.txt',quote = F,row.names = F,sep='\t')
write.table(gse31210_risk1,'results/gse31210_risk.txt',quote = F,row.names = F,sep='\t')
write.table(gse3141_risk1,'results/gse3141_risk.txt',quote = F,row.names = F,sep='\t')
write.table(gse37745_risk1,'results/gse37745_risk.txt',quote = F,row.names = F,sep='\t')
write.table(gse50081_risk1,'results/gse50081_risk.txt',quote = F,row.names = F,sep='\t')
write.table(gse68465_risk1,'results/gse68465_risk.txt',quote = F,row.names = F,sep='\t')

#07.CAF.module.cli/
####07.CAF.module.cli####
setwd('~/LIHC/GEO/07.CAF.module.cli/')
dir.create('results')

tcga.risk<-read.delim('~/LIHC/GEO/07.CAF.module.cli/tcga_risk.csv',sep='\t',header = T)

head(tcga.risk)
tcga.cli<-read.delim('~/LIHC/GEO/06.nk.module/geo_cli.txt',sep='\t',header = T)
tcga.risk.cli<-merge(tcga_risk,geo_cli,by='Sample')
write.table(tcga.risk.cli,'results/tcga.risk.cli.txt',quote = F,row.names = F,sep='\t')
write.table(tcga_risk$result,'results/tcga_risk.txt',quote = F,row.names = F,sep='\t')

install.packages('forestplot')
library(forestplot)
library(survcomp)
tcga_cox_datas=tcga.risk2.cli
table(tcga_cox_datas$ALT.....50U.L.)
tcga_cox_datas$T.Stage[tcga_cox_datas$T.Stage == 'T1' | tcga_cox_datas$T.Stage == 'T2'] <- 'T1+T2'
tcga_cox_datas$T.Stage[tcga_cox_datas$T.Stage == 'T3' | tcga_cox_datas$T.Stage == 'T4'] <- 'T3+T4'

table(tcga_cox_datas$N.Stage)
tcga_cox_datas$N.Stage[tcga_cox_datas$N.Stage == 'N1' | tcga_cox_datas$N.Stage == 'N2' | tcga_cox_datas$N.Stage == 'N3'] <- 'N1+N2+N3'

table(tcga_cox_datas$M.Stage)

table(tcga_cox_datas$Stage)
tcga_cox_datas$Stage[tcga_cox_datas$Stage == 'I' | tcga_cox_datas$Stage == 'II'] <- 'I+II'
tcga_cox_datas$Stage[tcga_cox_datas$Stage == 'III' | tcga_cox_datas$Stage == 'IV'] <- 'III+IV'

table(tcga_cox_datas$Gender)

#
# 单因素分析结果 #######
Age_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~Age,
                             data=tcga_cox_datas))
Age_sig_cox_dat <- data.frame(Names=rownames(Age_sig_cox[[8]]),
                              HR = round(Age_sig_cox[[7]][,2],3),
                              lower.95 = round(Age_sig_cox[[8]][,3],3),
                              upper.95 = round(Age_sig_cox[[8]][,4],3),
                              p.value=round(Age_sig_cox[[7]][,5],3))
Age_sig_cox_dat



Gender_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~Gender,
                                data=tcga_cox_datas))
Gender_sig_cox_dat <- data.frame(Names=rownames(Gender_cox[[8]]),
                                 HR = round(Gender_sig_cox[[7]][,2],3),
                                 lower.95 = round(Gender_sig_cox[[8]][,3],3),
                                 upper.95 = round(Gender_sig_cox[[8]][,4],3),
                                 p.value=round(Gender_sig_cox[[7]][,5],3))
Gender_sig_cox_dat



Main.Tumor.Size_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~Main.Tumor.Size,
                                         data=tcga_cox_datas))
Main.Tumor.Size_sig_cox_dat <- data.frame(Names=rownames(Main.Tumor.Size_sig_cox[[8]]),
                                          HR = round(Main.Tumor.Size_sig_cox[[7]][,2],3),
                                          lower.95 = round(Main.Tumor.Size_sig_cox[[8]][,3],3),
                                          upper.95 = round(Main.Tumor.Size_sig_cox[[8]][,4],3),
                                          p.value=round(Main.Tumor.Size_sig_cox[[7]][,5],3))
Main.Tumor.Size_sig_cox_dat


ALT_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~ALT,
                             data=tcga_cox_datas))
ALT_sig_cox_dat <- data.frame(Names=rownames(ALT_sig_cox[[8]]),
                              HR = round(ALT_sig_cox[[7]][,2],3),
                              lower.95 = round(ALT_sig_cox[[8]][,3],3),
                              upper.95 = round(ALT_sig_cox[[8]][,4],3),
                              p.value=round(ALT_sig_cox[[7]][,5],3))
ALT_sig_cox_dat


TNM.staging_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~TNM.staging,
                                     data=tcga_cox_datas))
TNM.staging_sig_cox_dat <- data.frame(Names=rownames(TNM.staging_sig_cox[[8]]),
                                      HR = round(TNM.staging_sig_cox[[7]][,2],3),
                                      lower.95 = round(TNM.staging_sig_cox[[8]][,3],3),
                                      upper.95 = round(TNM.staging_sig_cox[[8]][,4],3),
                                      p.value=round(TNM.staging_sig_cox[[7]][,5],3))
TNM.staging_sig_cox_dat


Risk_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~Risk,
                              data=tcga_cox_datas))
Risk_sig_cox_dat <- data.frame(Names=rownames(Risk_sig_cox[[8]]),
                               HR = round(Risk_sig_cox[[7]][,2],3),
                               lower.95 = round(Risk_sig_cox[[8]][,3],3),
                               upper.95 = round(Risk_sig_cox[[8]][,4],3),
                               p.value=round(Risk_sig_cox[[7]][,5],3))
Risk_sig_cox_dat

Multinodular_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~Multinodular,
                                      data=tcga_cox_datas))
Multinodular_sig_cox_dat <- data.frame(Names=rownames(Multinodular_sig_cox[[8]]),
                                       HR = round(Multinodular_sig_cox[[7]][,2],3),
                                       lower.95 = round(Multinodular_sig_cox[[8]][,3],3),
                                       upper.95 = round(Multinodular_sig_cox[[8]][,4],3),
                                       p.value=round(Multinodular_sig_cox[[7]][,5],3))
Multinodular_sig_cox_dat

Cirrhosis_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~Cirrhosis,
                                   data=tcga_cox_datas))
Cirrhosis_sig_cox_dat <- data.frame(Names=rownames(Cirrhosis_sig_cox[[8]]),
                                    HR = round(Cirrhosis_sig_cox[[7]][,2],3),
                                    lower.95 = round(Cirrhosis_sig_cox[[8]][,3],3),
                                    upper.95 = round(Cirrhosis_sig_cox[[8]][,4],3),
                                    p.value=round(Cirrhosis_sig_cox[[7]][,5],3))
Cirrhosis_sig_cox_dat

BCLC.staging_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~BCLC.staging,
                                      data=tcga_cox_datas))
BCLC.staging_sig_cox_dat <- data.frame(Names=rownames(BCLC.staging_sig_cox[[8]]),
                                       HR = round(BCLC.staging_sig_cox[[7]][,2],3),
                                       lower.95 = round(BCLC.staging_sig_cox[[8]][,3],3),
                                       upper.95 = round(BCLC.staging_sig_cox[[8]][,4],3),
                                       p.value=round(BCLC.staging_sig_cox[[7]][,5],3))
BCLC.staging_sig_cox_dat

CLIP.staging_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~CLIP.staging,
                                      data=tcga_cox_datas))
CLIP.staging_sig_cox_dat <- data.frame(Names=rownames(CLIP.staging_sig_cox[[8]]),
                                       HR = round(CLIP.staging_sig_cox[[7]][,2],3),
                                       lower.95 = round(CLIP.staging_sig_cox[[8]][,3],3),
                                       upper.95 = round(CLIP.staging_sig_cox[[8]][,4],3),
                                       p.value=round(CLIP.staging_sig_cox[[7]][,5],3))
CLIP.staging_sig_cox_dat


AFP_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~AFP,
                             data=tcga_cox_datas))
AFP_sig_cox_dat <- data.frame(Names=rownames(AFP_sig_cox[[8]]),
                              HR = round(AFP_sig_cox[[7]][,2],3),
                              lower.95 = round(AFP_sig_cox[[8]][,3],3),
                              upper.95 = round(AFP_sig_cox[[8]][,4],3),
                              p.value=round(AFP_sig_cox[[7]][,5],3))
AFP_sig_cox_dat



sig_cox_dat <- rbind(Age_sig_cox_dat,
                     Gender_sig_cox_dat,
                     Main.Tumor.Size_sig_cox_dat,
                     ALT_sig_cox_dat,
                     TNM.staging_sig_cox_dat,
                     Multinodular_sig_cox_dat,
                     Cirrhosis_sig_cox_dat,
                     BCLC.staging_sig_cox_dat,
                     CLIP.staging_sig_cox_dat,
                     AFP_sig_cox_dat,
                     Risk_sig_cox_dat)

data.sig <- data.frame(Names=sig_cox_dat$Names,
                       p.value=sig_cox_dat$p.value,
                       sig_cox_dat$HR,
                       sig_cox_dat$lower.95,
                       sig_cox_dat$upper.95)

rownames(data.sig) <- c("Age",
                        "Gender",
                        "Main.Tumor.Size",
                        "ALT",
                        "TNM.staging",
                        "Multinodular",
                        "Cirrhosis",
                        "BCLC.staging",
                        "CLIP.staging",
                        "AFP",
                        "Riskscore")
data.sig$Names <- rownames(data.sig)
data.sig$p.value=ifelse(data.sig$p.value==0.000,'<0.001',data.sig$p.value)
write.table(data.sig,'results/data.sig.txt',sep='\t',row.names = T,quote = F)
vio_fores<-function(dat){
  library(forestplot)
  col=fpColors(box='red',summary='#8B008B',lines = 'grey',zero = 'grey')
  smary=rep(F,length(as.numeric(dat[,ncol(dat)-2])))
  lower=as.numeric(dat[,ncol(dat)-1])
  upper=as.numeric(dat[,ncol(dat)])
  nind=which(is.na(lower)|is.na(upper)|is.na(as.numeric(dat[,ncol(dat)-2])))
  smary[nind]=T
  labeltext=as.matrix(dat[,1:(ncol(dat)-3)])
  adt=paste0(dat[,3],'(',dat[,4],',',dat[,5],')')
  
  labeltext=cbind(labeltext,adt)
  colnames(labeltext)[3]=c('Hazard Ratio(95% CI)')
  
  hz_list=list('2'=gpar(lty=1,col='#8B008B'),
               '3'=gpar(lty=1,col='#8B008B'))
  names(hz_list)=c(2,nrow(labeltext)+2)
  forestplot(labeltext = rbind(colnames(labeltext),labeltext),
             hrzl_lines = hz_list,
             mean = c(NA,as.numeric(dat[,ncol(dat)-2])),
             lower =c(NA,lower),
             upper = c(NA,upper),
             is.summary=c(T,smary),
             zero = 1,
             boxsize = 0.4, 
             lineheight = unit(10,'mm'),
             colgap = unit(8,'mm'),
             lwd.zero = 2,
             lwd.ci = 2,
             col=col,
             xlab='HR',
             lwd.xaxis=2,
             lty.ci = 'solid',
             clip = c(min(lower,na.rm = T),max(upper,na.rm = T)),
             xlog=T,
             graph.pos = 2)
  
  
}
dat=data.sig

pdf('results/sig_cox.pdf', width = 8, height = 5)
vio_fores(data.sig)
dev.off()
#多因素分析结果
muti_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~Gender+Main.Tumor.Size+TNM.staging+Multinodular+Cirrhosis+BCLC.staging+CLIP.staging+AFP+Risk, 
                              data=tcga_cox_datas))
muti_cox_dat <- data.frame(Names=rownames(muti_sig_cox[[8]]),
                           HR = round(muti_sig_cox[[7]][,2],3),
                           lower.95 = round(muti_sig_cox[[8]][,3],3),
                           upper.95 = round(muti_sig_cox[[8]][,4],3),
                           p.value=round(muti_sig_cox[[7]][,5],3))
data.muti <- data.frame(Names=muti_cox_dat$Names,
                        p.value=muti_cox_dat$p.value,
                        muti_cox_dat$HR,
                        muti_cox_dat$lower.95,
                        muti_cox_dat$upper.95)
rownames(data.muti) <- c("Gender","Main.Tumor.Size","TNM.staging","Multinodular","Cirrhosis","BCLC.staging","CLIP.staging","AFP","RiskScore")
data.muti$Names <- rownames(data.muti)
data.muti$p.value=ifelse(data.muti$p.value==0.000,'<0.001',data.muti$p.value)

write.table(data.muti,'results/data.muti.txt',sep='\t',row.names = T,quote = F)
pdf('results/muti_cox.pdf', width = 8, height = 5)
vio_fores(data.muti)
dev.off()


# 列线图以及 DCA 图 #############
#install_github("zzawadz/customLayout")
library(devtools)
library(customLayout)

getC_index=function(riskscore,os,status){
  inds=which(!is.na(riskscore)&!is.na(os)&!is.na(status))
  riskscore=riskscore[inds]
  os=os[inds]
  status=status[inds]
  c1=survcomp::concordance.index(x=riskscore, surv.time=os, surv.event=status,
                                 method="noether")
  c2=concordance.index(x=riskscore[order(rnorm(length(riskscore)))], surv.time=os, surv.event=status,
                       method="noether")
  p=min(cindex.comp(c1, c2)$p.value,cindex.comp(c2, c1)$p.value)  
  return(c1)
}
crbind2DataFrame=function(dat,full=F){
  print(class(dat))
  if(class(dat)=='table'){
    if(!is.na(ncol(dat))){
      dat=apply(dat,2,function(x){
        return(x)
      })
    }
  }
  if(class(dat)!='data.frame'){
    dat1=as.data.frame(as.matrix(dat))
  }else{
    dat1=dat
  }
  dat1.class=apply(dat1, 2, class)
  #which(dat1.class!='numeric')
  #print(head(dat1))
  for(i in which(dat1.class!='numeric')){
    dat1[,i]=as.character(dat1[,i])
    if(full){
      dat1[,i]=as.numeric(dat1[,i])
    }else{
      dt=dat1[which(gsub(' ','',dat1[,i])!=''&!is.na(dat1[,i])),i]
      dt=dt[which(dt!='Inf'&dt!='NaN'&dt!='NA')]
      #dt[which(is.na(as.numeric(dt)))]
      if(sum(is.na(as.numeric(dt)))<length(dt)*0.1){
        #print(dat1[,i])
        dat1[,i]=as.numeric(dat1[,i])
      }
    }
  }
  return(dat1)  
}
mg_plotDCA=function(status,fmlas,modelNames,data){
  set.seed(123)
  all.mod=list()
  for(i in 1:length(fmlas)){
    fmla <- as.formula(paste0("status~",fmlas[i]))
    model<-rmda::decision_curve(fmla,
                                data=data,
                                bootstraps=500)
    all.mod=c(all.mod,list(model))
  }
  rmda::plot_decision_curve(all.mod,
                            curve.names=modelNames,
                            xlim=c(0,1),legend.position="topright",
                            confidence.intervals=FALSE)
}

nomogram_plot=function(clinical_riskscore,
                       os,
                       status,
                       title='Nomogram',
                       quick=T,
                       mks = c(1,3,5)){
  mg_colors= ggsci::pal_npg("nrc")(10)
  norm.stat.al=data.frame(clinical_riskscore,time=os,status=status)
  norm.stat.al=as.data.frame(norm.stat.al)
  library(rms)
  env <- globalenv()
  env$MG_Grobal_DDSet <- rms::datadist(norm.stat.al) 
  options(datadist='MG_Grobal_DDSet')
  fmla <- as.formula(paste0("Surv(time, status) ~",paste0(colnames(clinical_riskscore),collapse = '+')))
  cox2 <- cph(fmla, data = norm.stat.al,surv = T,x = T,y = T)
  
  fp <- predict(cox2)
  cindex=getC_index(fp,norm.stat.al$time,norm.stat.al$status)
  cut.time=c()
  if(quantile(os[!is.na(os)])['75%']<12){
    cut.time=mks
  }else if(quantile(os[!is.na(os)])['75%']<365){
    cut.time=c(12*mks[1],12*mks[2],12*mks[3])
  }else{
    cut.time=c(365*mks[1],365*mks[2],365*mks[3])
  }
  cut.time=cut.time[which(cut.time<quantile(os,seq(0,1,0.01))['95%'])]
  print(cut.time)
  surv=Survival(cox2)
  survs=list()
  cal_all=list()
  for(i in 1:length(cut.time)){
    f1<-cph(formula = fmla,data=norm.stat.al,x=T,y=T,surv = T,na.action=na.delete,time.inc = cut.time[i]) 
    cal1<-calibrate(f1, cmethod="KM", method="boot",u=cut.time[i],m=floor(sum(f1$n)/3)) 
    cal_all=c(cal_all,list(cal1))
  }
  if(length(cut.time)==1){
    surv1 <- function(x)surv(cut.time[1],lp=x) 
    survs=list(surv1)
  }else if(length(cut.time)==2){
    surv1 <- function(x)surv(cut.time[1],lp=x) 
    surv2 <- function(x)surv(cut.time[2],lp=x) 
    survs=list(surv1,surv2)
  }else if(length(cut.time)==3){
    surv1 <- function(x)surv(cut.time[1],lp=x) 
    surv2 <- function(x)surv(cut.time[2],lp=x) 
    surv3 <- function(x)surv(cut.time[3],lp=x) 
    survs=list(surv1,surv2,surv3)
  }
  nom=nomogram(cox2,fun=survs,lp= F
               ,funlabel=c(paste0(mks[1], '-Year Survival'),
                           paste0(mks[2], '-Year Survival'),
                           paste0(mks[3], '-Year Survival'))[1:length(cut.time)]
               ,maxscale=100
               ,fun.at=seq(0,1,0.2)
  )
  
  if(!quick){
    cal_all=list()
    for(i in 1:length(cut.time)){
      cal1=get_best_calibrate(cox2,cut.time[i])
      cal_all=c(cal_all,list(cal1))
    }
  }
  lay2 <- customLayout::lay_new(matrix(1:2))
  lay1 <- customLayout::lay_new(matrix(1:1))
  cl <- customLayout::lay_bind_col(lay2, lay1, widths = c(1, 1.5)) 
  customLayout::lay_set(cl) 
  
  plot(cal_all[[1]],lwd = 2,lty = 1,errbar.col = mg_colors[1], bty = "l",xlim = c(0,1),ylim= c(0,1),xlab = "Nomogram-prediced OS (%)"
       ,ylab = "Observed OS (%)",col = mg_colors[1],cex.lab=1.2,cex.axis=1, cex.main=1.2, cex.sub=0.6,subtitles = F)
  mtext("")
  if(length(cal_all)>1){
    for(i in 2:length(cal_all)){
      plot(cal_all[[i]],lwd = 2,lty = 1,errbar.col = mg_colors[i],xlim = c(0,1),ylim= c(0,1),col = mg_colors[i],add = T)
    }
  }
  abline(0,1, lwd = 2, lty = 3, col = 'black')
  legend("topleft", legend = c(paste0(mks[1], '-Year'),
                               paste0(mks[2], '-Year'),
                               paste0(mks[3], '-Year'))[1:length(cut.time)],col =mg_colors[1:length(cut.time)],lwd = 2,cex = 1.2,bty = "n")
  
  fp <- predict(cox2)
  cindex=getC_index(fp,norm.stat.al$time,norm.stat.al$status)
  dca_dat=data.frame(Nomogram=fp,status=norm.stat.al$status)
  fp.al=cbind(fp)
  for(i in 1:ncol(clinical_riskscore)){
    fmla1 <- as.formula(paste0("Surv(time, status) ~",colnames(clinical_riskscore)[i]))
    cox1 <- cph(fmla1, data = norm.stat.al,surv = T,x = T,y = T)
    fp1 <- predict(cox1)
    fp.al=cbind(fp.al,fp1)
  }
  colnames(fp.al)=c('Nomogram',colnames(clinical_riskscore))
  fp.al=crbind2DataFrame(fp.al)
  fp.al$status=norm.stat.al$status
  mg_plotDCA(fp.al$status
             ,c('Nomogram',colnames(clinical_riskscore))
             ,c('Nomogram',colnames(clinical_riskscore)),fp.al)
  plot(nom)
  options(datadist=NULL)
  return(list(Mod=cox2,Cindex=cindex,CutTime=cut.time))
}


nomogram_plot <- function(clinical_riskscore,
                          os,
                          status,
                          title='Nomogram',
                          quick=TRUE,
                          mks = c(1, 3, 5)) {
  # Load necessary libraries
  library(rms)
  library(ggsci)
  
  mg_colors = ggsci::pal_npg("nrc")(10)
  
  # Prepare data frame
  norm.stat.al = data.frame(clinical_riskscore, time = os, status = status)
  
  # Check for NA values and remove if necessary
  norm.stat.al = na.omit(norm.stat.al)
  
  # Set up the datadist object
  env <- globalenv()
  env$MG_Grobal_DDSet <- rms::datadist(norm.stat.al)
  options(datadist = 'MG_Grobal_DDSet')
  
  # Create formula for Cox model
  fmla <- as.formula(paste0("Surv(time, status) ~ ", paste(colnames(clinical_riskscore), collapse = '+')))
  
  # Fit Cox model
  cox2 <- cph(fmla, data = norm.stat.al, surv = TRUE, x = TRUE, y = TRUE)
  
  # Predictions and C-index
  fp <- predict(cox2)
  cindex <- getC_index(fp, norm.stat.al$time, norm.stat.al$status)
  
  # Determine cut-off times
  cut.time = numeric()
  if (quantile(os[!is.na(os)], 0.75) < 12) {
    cut.time = mks
  } else if (quantile(os[!is.na(os)], 0.75) < 365) {
    cut.time = c(12 * mks[1], 12 * mks[2], 12 * mks[3])
  } else {
    cut.time = c(365 * mks[1], 365 * mks[2], 365 * mks[3])
  }
  
  cut.time = cut.time[cut.time < quantile(os, 0.95)]
  print(cut.time)
  
  surv=Survival(cox2)
  survs=list()
  cal_all=list()
  for(i in 1:length(cut.time)){
    f1<-cph(formula = fmla,data=norm.stat.al,x=T,y=T,surv = T,na.action=na.delete,time.inc = cut.time[i]) 
    cal1<-calibrate(f1, cmethod="KM", method="boot",u=cut.time[i],m=floor(sum(f1$n)/3)) 
    cal_all=c(cal_all,list(cal1))
  }
  if(length(cut.time)==1){
    surv1 <- function(x)surv(cut.time[1],lp=x) 
    survs=list(surv1)
  }else if(length(cut.time)==2){
    surv1 <- function(x)surv(cut.time[1],lp=x) 
    surv2 <- function(x)surv(cut.time[2],lp=x) 
    survs=list(surv1,surv2)
  }else if(length(cut.time)==3){
    surv1 <- function(x)surv(cut.time[1],lp=x) 
    surv2 <- function(x)surv(cut.time[2],lp=x) 
    surv3 <- function(x)surv(cut.time[3],lp=x) 
    survs=list(surv1,surv2,surv3)
  }
  nom=nomogram(cox2,fun=survs,lp= F
               ,funlabel=c(paste0(mks[1], '-Year Survival'),
                           paste0(mks[2], '-Year Survival'),
                           paste0(mks[3], '-Year Survival'))[1:length(cut.time)]
               ,maxscale=100
               ,fun.at=seq(0,1,0.2)
  )
  
  if(!quick){
    cal_all=list()
    for(i in 1:length(cut.time)){
      cal1=get_best_calibrate(cox2,cut.time[i])
      cal_all=c(cal_all,list(cal1))
    }
  }
  lay2 <- customLayout::lay_new(matrix(1:2))
  lay1 <- customLayout::lay_new(matrix(1:1))
  cl <- customLayout::lay_bind_col(lay2, lay1, widths = c(1, 1.5)) 
  customLayout::lay_set(cl) 
  
  plot(cal_all[[1]],lwd = 2,lty = 1,errbar.col = mg_colors[1], bty = "l",xlim = c(0,1),ylim= c(0,1),xlab = "Nomogram-prediced OS (%)"
       ,ylab = "Observed OS (%)",col = mg_colors[1],cex.lab=1.2,cex.axis=1, cex.main=1.2, cex.sub=0.6,subtitles = F)
  mtext("")
  if(length(cal_all)>1){
    for(i in 2:length(cal_all)){
      plot(cal_all[[i]],lwd = 2,lty = 1,errbar.col = mg_colors[i],xlim = c(0,1),ylim= c(0,1),col = mg_colors[i],add = T)
    }
  }
  abline(0,1, lwd = 2, lty = 3, col = 'black')
  legend("topleft", legend = c(paste0(mks[1], '-Year'),
                               paste0(mks[2], '-Year'),
                               paste0(mks[3], '-Year'))[1:length(cut.time)],col =mg_colors[1:length(cut.time)],lwd = 2,cex = 1.2,bty = "n")
  
  fp <- predict(cox2)
  cindex=getC_index(fp,norm.stat.al$time,norm.stat.al$status)
  dca_dat=data.frame(Nomogram=fp,status=norm.stat.al$status)
  fp.al=cbind(fp)
  for(i in 1:ncol(clinical_riskscore)){
    fmla1 <- as.formula(paste0("Surv(time, status) ~",colnames(clinical_riskscore)[i]))
    cox1 <- cph(fmla1, data = norm.stat.al,surv = T,x = T,y = T)
    fp1 <- predict(cox1)
    fp.al=cbind(fp.al,fp1)
  }
  colnames(fp.al)=c('Nomogram',colnames(clinical_riskscore))
  fp.al=crbind2DataFrame(fp.al)
  fp.al$status=norm.stat.al$status
  mg_plotDCA(fp.al$status
             ,c('Nomogram',colnames(clinical_riskscore))
             ,c('Nomogram',colnames(clinical_riskscore)),fp.al)
  plot(nom)
  options(datadist=NULL)
  return(list(Mod=cox2,Cindex=cindex,CutTime=cut.time))
}

nomogram_buti=function(cox_model,cut.time,title='Nomogram'){
  library(regplot)
  regplot(cox_model#对观测2的六个指标在列线图上进行计分展示
          ,observation=pbc[2,] #也可以不展示
          #预测3年和5年的死亡风险，此处单位是day
          ,title=title
          ,failtime = cut.time
          ,prfail = TRUE #cox回归中需要TRUE
          ,showP = T #是否展示统计学差异
          ,droplines = F#观测2示例计分是否画线
          #,colors = mg_colors[1:3] #用前面自己定义的颜色
          #,rank="decreasing") #根据统计学差异的显著性进行变量的排序
          #,interval="confidence"
          #,rank="decreasing"
          #,clickable=T
          ,points=TRUE)
  
}

pdf('results/nomogram.pdf', width = 20, height = 20)

nom.plot=nomogram_plot(data.frame(Gender=tcga_cox_datas$Gender,
                                  Main.Tumor.Size=tcga_cox_datas$Main.Tumor.Size,
                                  TNM.staging=tcga_cox_datas$TNM.staging,
                                  Multinodular=tcga_cox_datas$Multinodular,
                                  Cirrhosis=tcga_cox_datas$Cirrhosis,
                                  BCLC.staging=tcga_cox_datas$BCLC.staging,
                                  CLIP.staging=tcga_cox_datas$CLIP.staging,
                                  AFP=tcga_cox_datas$AFP,
                                  RiskScore=tcga_cox_datas$Risk),
                       os = tcga_cox_datas$OS.time,
                       status = tcga_cox_datas$OS,
                       mks = c(1,3,5))


dev.off()

nomogram_buti(nom.plot$Mod,cut.time = c(1*365,3*365,5*365))
write.table(tcga_cox_datas,'results/tcga_cox_datas.txt',quote = F,row.names = F,sep='\t')


library(survival)
library(rms)
head(tcga.risk2.cli)
mul_cox <- coxph(Surv(time, status) ~ Gender + Main.Tumor.Size + TNM.staging + Multinodular + Cirrhosis + BCLC.staging + CLIP.staging + AFP + Risk, data = tcga.risk2.cli)

mul_cox 
library(rms)
mul_cox_2 <- cph(Surv(time, status) ~ Gender + Main.Tumor.Size + TNM.staging + Multinodular + Cirrhosis + BCLC.staging + CLIP.staging + AFP + Risk, data = tcga.risk2.cli, x= TRUE, y = TRUE, surv = TRUE)

mul_cox_2

sur <- Survival(mul_cox_2, surv=TRUE)

#定义用于不同时间点的生存函数
sur1<-function(x)sur(365,x) # 1年生存
sur2<-function(x)sur(730,x) # 2年生存
sur3<-function(x)sur(1095,x) # 3年生存

#绘制nomogram，可视化cox比例风险模型的结果
dd <- datadist(tcga.risk2.cli)  # 将 'lung' 替换为您实际使用的数据集
options(datadist = 'dd')  # 设置选项以使用这个datadist对象
nom<-nomogram(mul_cox_2,
              fun=list(sur1, sur2, sur3), # 为每个时间点指定生存函数
              fun.at=c(0.05, seq(0.1, 0.9, by = 0.05), 0.95), # 设置诺莫图上的生存概率点
              funlabel=c("1 year survival", "2 years survival", "3 years survival")) # 设置生存时间标签

plot(nom)
library(regplot)
regplot(mul_cox_2,
        #observation = tcga.risk2.cli[6, ], # 指定进行计分展示的观测数据，这里使用了数据集中的第6行，当然也可以不展示
        points = TRUE, # 在列线图上显示点
        plots = c("density", "no plot"), # 设置要显示的图表类型，这里包括密度图和无图表
        failtime = c(365, 730, 1095), # 指定预测的时间点，这里预测了1年、2年和3年的死亡风险，单位是day
        odds = F, # 是否显示死亡几率
        leftlabel = T, # 是否显示左侧标签
        prfail = TRUE, # 在 Cox 回归中需要设置为 TRUE
        showP = T, # 是否展示统计学差异
        droplines = T, # 示例计分是否画线
        colors = "red", # 可以使用自定义的颜色
        rank = "range", # 根据统计学差异的显著性进行变量的排序
        interval = "confidence", # 使用置信区间进行可视化
        title = "Cox regression") 


# 校准曲线（Calibration Curve)

# time.inc参数用于指定风险比计算的时间单位
mul_cox_2 <- cph(Surv(time, status) ~ Gender + Main.Tumor.Size + TNM.staging + Multinodular + Cirrhosis + BCLC.staging + CLIP.staging + AFP + Risk, data = tcga.risk2.cli, x= TRUE, y = TRUE, surv = TRUE, time.inc = 365)
# 使用calibrate函数创建一个校准对象cal1
cal1 <- calibrate(mul_cox_2, 
                  cmethod = 'KM', # 表示使用Kaplan-Meier（KM）估计方法进行校准
                  method = "boot", # 表示使用自助法（Bootstrap）进行校准，Bootstrap 是一种统计方法，它通过从原始数据中有放回地进行重采样来估计参数的不确定性和分布。在这里，Bootstrap 用于生成多个随机样本来估计校准曲线的分布，以便获得更可靠的校准结果。
                  u = 365, # 设置时间间隔，需要与之前模型中定义的time.inc一致
                  m = 66, # 每次抽样的样本量，根据样本量来确定，标准曲线一般将所有样本分为3组（在图中显示3个点）
                  B = 1000) # 抽样次数

# 绘制校准曲线
par(mar = c(6, 6, 3, 3))
plot(cal1, # 绘制校准曲线的数据
     lwd=1, # 线条宽度为1
     lty=1, # 线条类型为1（实线）
     conf.int=T, # 是否显示置信区间
     errbar.col="blue3", # 直线和曲线的误差线颜色设置为蓝色
     col="red3", # 校准曲线的颜色设置为红色
     xlim=c(0,1), # x轴的限制范围，从0到1
     ylim=c(0,1), # y轴的限制范围，从0到1
     xlab="Nomogram-Predicted Probability of 1-Year OS", # x轴标签
     ylab="Actual 1-Year OS (proportion)", # y轴标签
     subtitles = F) # 不显示副标题
# 在校准曲线上添加一条直线
lines(cal1[,c('mean.predicted',"KM")], 
      type = 'b', # 连线的类型，可以是"p","b","o"
      lwd = 2, # 连线的粗细
      pch = 16, # 点的形状，可以是0-20
      col = c("red3")) # 连线的颜色

# 移除默认的标题
mtext("")

# 绘制图形边框
box(lwd = 1) # 边框粗细

# 添加一条虚线表示对角线
abline(0,1, lty = 3, # 对角线为虚线
       lwd = 2, # 对角线的粗细
       col = c("black")) # 对角线的颜色

#校准曲线
#1-year
cox1 <- cph(Surv(OS.time,OS) ~ Gender + Main.Tumor.Size + TNM.staging + Multinodular + Cirrhosis + BCLC.staging + CLIP.staging + AFP + Risk, 
            surv=T,x=T, y=T,time.inc = 1*365,data=tcga.risk2.cli)
cal <- calibrate(cox1, cmethod="KM", method="boot", u=1*365*1, m= 100, B=1000)

#3-year
cox3 <- cph(Surv(OS.time,OS) ~ Gender + Main.Tumor.Size + TNM.staging + Multinodular + Cirrhosis + BCLC.staging + CLIP.staging + AFP + Risk, 
            surv=T,x=T, y=T,time.inc = 1*365*3,data=tcga.risk2.cli)
ca3 <- calibrate(cox3, cmethod="KM", method="boot", u=1*365*3, m= 100, B=1000)

#5-year
cox5 <- cph(Surv(OS.time,OS) ~ Gender + Main.Tumor.Size + TNM.staging + Multinodular + Cirrhosis + BCLC.staging + CLIP.staging + AFP + Risk,            
            surv=T,x=T, y=T,time.inc = 1*365*5,data=tcga.risk2.cli)
ca5 <- calibrate(cox5, cmethod="KM", method="boot", u=1*365*5, m= 100, B=1000)

plot(cal,lwd=2,lty=1,errbar.col="black",     xlim = c(0.2,1) ,     ylim = c(0.1,1),      
     xlab ="Nomogram-Predicted Probability of 1,3,5Year Survival", ylab="Actual 1,3,5Year Survival",col="blue",sub=F)
plot(ca3,     add = T , lwd=2,lty=1,errbar.col="black",col="red",sub=F)
plot(ca5,     add = T ,     lwd=2,lty=1,errbar.col="black",col="green",sub=F)

3
##DCA
install.packages("dcurves") 
library(dcurves)
library(caret)
#remotes::install_github('yikeshu0611/ggDCA')
library(ggDCA)
pdf('results/DCA.pdf', width = 15, height = 15)
Cox_nomo<-rms::cph(Surv(OS.time, OS) ~ Main.Tumor.Size + TNM.staging+ Multinodular + Risk, data=tcga.risk2.cli)
d_train <- dca(Cox_nomo,  times=1800)

ggplot(d_train)
dev.off()


library("timeROC")
tcga_cox_auc=tcga_cox_datas
tcga_cox_auc$riskscore=as.numeric(tcga_cox_auc$riskscore)
tcga_cox_auc$TNM.staging=as.numeric(as.factor(tcga_cox_auc$TNM.staging))
tcga_cox_auc$Gender=as.numeric(as.factor(tcga_cox_auc$Gender))
tcga_cox_auc$Main.Tumor.Size=as.numeric(as.factor(tcga_cox_auc$Main.Tumor.Size))
tcga_cox_auc$Multinodular=as.numeric(as.factor(tcga_cox_auc$Multinodular))
tcga_cox_auc$Cirrhosis=as.numeric(as.factor(tcga_cox_auc$Cirrhosis))
tcga_cox_auc$BCLC.staging=as.numeric(as.factor(tcga_cox_auc$BCLC.staging))
tcga_cox_auc$CLIP.staging=as.numeric(as.factor(tcga_cox_auc$CLIP.staging))
tcga_cox_auc$AFP=as.numeric(as.factor(tcga_cox_auc$AFP))

ROC.DSST.Risk=timeROC(T=tcga_cox_auc$OS.time/365,
                      delta=tcga_cox_auc$OS,
                      marker=tcga_cox_auc$riskscore,
                      cause=1,weighting="marginal",
                      times=1:5,
                      iid=TRUE)
ROC.DSST.gender <- timeROC(T = tcga_cox_auc$OS.time/365,
                           delta = tcga_cox_auc$OS,
                           marker = tcga_cox_auc$Gender,
                           cause = 1,
                           weighting = "marginal",
                           times = 1:5,
                           iid = TRUE)
ROC.DSST.TNM.staging=timeROC(T=tcga_cox_auc$OS.time/365,
                             delta=tcga_cox_auc$OS,
                             marker=tcga_cox_auc$TNM.staging,
                             cause=1,weighting="marginal",
                             times=1:5,
                             iid=TRUE)

ROC.DSST.Main.Tumor.Size=timeROC(T=tcga_cox_auc$OS.time/365,
                                 delta=tcga_cox_auc$OS,
                                 marker=tcga_cox_auc$Main.Tumor.Size,
                                 cause=1,weighting="marginal",
                                 times=1:5,
                                 iid=TRUE)

ROC.DSST.Multinodular=timeROC(T=tcga_cox_auc$OS.time/365,
                              delta=tcga_cox_auc$OS,
                              marker=tcga_cox_auc$Multinodular,
                              cause=1,weighting="marginal",
                              times=1:5,
                              iid=TRUE)
ROC.DSST.Cirrhosis=timeROC(T=tcga_cox_auc$OS.time/365,
                           delta=tcga_cox_auc$OS,
                           marker=tcga_cox_auc$Cirrhosis,
                           cause=1,weighting="marginal",
                           times=1:5,
                           iid=TRUE)
ROC.DSST.CLIP.staging=timeROC(T=tcga_cox_auc$OS.time/365,
                              delta=tcga_cox_auc$OS,
                              marker=tcga_cox_auc$CLIP.staging,
                              cause=1,weighting="marginal",
                              times=1:5,
                              iid=TRUE)
ROC.DSST.AFP=timeROC(T=tcga_cox_auc$OS.time/365,
                     delta=tcga_cox_auc$OS,
                     marker=tcga_cox_auc$AFP,
                     cause=1,weighting="marginal",
                     times=1:5,
                     iid=TRUE)
ROC.DSST.BCLC.staging=timeROC(T=tcga_cox_auc$OS.time/365,
                              delta=tcga_cox_auc$OS,
                              marker=tcga_cox_auc$BCLC.staging,
                              cause=1,weighting="marginal",
                              times=1:5,
                              iid=TRUE)

ROC.DSST.Nomo<-timeROC(T=tcga_cox_auc$OS.time/365,
                       delta=tcga_cox_auc$OS,
                       marker=tcga_cox_auc$riskscore,
                       other_markers=as.matrix(tcga_cox_auc[,c("Gender",
                                                               "Main.Tumor.Size",
                                                               "TNM.staging",
                                                               "Multinodular",
                                                               "Cirrhosis",
                                                               "BCLC.staging",
                                                               "CLIP.staging",
                                                               "AFP")]),
                       cause=1,
                       weighting="cox",
                       times=1:5,
                       iid=F)

ROC.DSST.gender$AUC
ROC.DSST.T.Stage$AUC
ROC.DSST.N.Stage$AUC
ROC.DSST.M.Stage$AUC
ROC.DSST.Stage$AUC
ROC.DSST.Risk$AUC
ROC.DSST.Nomo$AUC
ROC.DSST.Age$AUC
mg_colors=ggsci::pal_lancet()(9)
pdf('results/AUC.pdf',height = 7,width = 7)
plotAUCcurve(ROC.DSST.Nomo,conf.int=F,col=mg_colors[1])
plotAUCcurve(ROC.DSST.Risk,conf.int=F,col=mg_colors[2],add=TRUE)
plotAUCcurve(ROC.DSST.gender,conf.int=F,col=mg_colors[3],add=TRUE)
plotAUCcurve(ROC.DSST.TNM.staging,conf.int=F,col=mg_colors[4],add=TRUE)
plotAUCcurve(ROC.DSST.Main.Tumor.Size,conf.int=F,col=mg_colors[5],add=TRUE)
plotAUCcurve(ROC.DSST.Multinodular,conf.int=F,col=mg_colors[6],add=TRUE)
plotAUCcurve(ROC.DSST.Cirrhosis,conf.int=F,col=mg_colors[7],add=TRUE)
plotAUCcurve(ROC.DSST.CLIP.staging,conf.int=F,col=mg_colors[8],add=TRUE)
plotAUCcurve(ROC.DSST.AFP,conf.int=F,col=mg_colors[7],add=TRUE)
plotAUCcurve(ROC.DSST.BCLC.staging,conf.int=F,col=mg_colors[8],add=TRUE)


legend("topright",c("Nomogram","RiskScore","Gender",
                    "Main.Tumor.Size",
                    "TNM.staging",
                    "Multinodular",
                    "Cirrhosis",
                    "BCLC.staging",
                    "CLIP.staging",
                    "AFP")
       ,col=mg_colors[c(1:8)],lty=1,lwd=2)

dev.off()

####09.gene.gsva####
setwd('~/LIHC/GEO/09.gene.gsva/')
dir.create('results')
tcga.dat<-read.delim('../03.data.pre/GSE53624/results/geo_T.txt',sep='\t',header = T,check.names = F)
tcga.module<-read.delim('../06.CAF.module/results/tcga.risk.txt',sep='\t',header = T)
tcga.module<- tcga_risk
colnames(tcga.module)
tcga.module.gene.exp=tcga.module[,3:13]
rownames(tcga.module.gene.exp)=tcga.module$Sample
library(GSVA)
library(GSEABase)
gmtFile='c2.cp.kegg.v7.5.1.symbols.gmt'
c2KEGG <- getGmt(gmtFile,
                 collectionType=BroadCollection(category="h"),
                 geneIdType=SymbolIdentifier())
tcga.kegg <- gsva(as.matrix(tcga.dat),
                  c2KEGG,
                  method = 'ssgsea',
                  min.sz = 10,
                  max.sz = 500,
                  verbose = TRUE)
save(tcga.kegg,file = 'tcga.kegg.Rdata')
load('tcga.kegg.Rdata')
tcga.kegg[1:4,1:4]
write.table(tcga.kegg,'results/tcga.kegg.enrichscore.txt',quote = F,row.names = T,sep='\t')
gene.kegg.dat<-cbind.data.frame(tcga.module.gene.exp,
                                t(tcga.kegg[,rownames(tcga.module.gene.exp)]))
gene.kegg.dat[1:4,1:7]
gene.kegg.cor <- Hmisc::rcorr(as.matrix(gene.kegg.dat))
gene.kegg.cor.r=reshape2::melt(gene.kegg.cor$r)
gene.kegg.cor.p=reshape2::melt(gene.kegg.cor$P)
colnames(gene.kegg.cor.r)=c('gene','pathway','cor')
colnames(gene.kegg.cor.p)=c('gene','pathway','p')
gene.kegg.cor.res<-merge(gene.kegg.cor.r,gene.kegg.cor.p,by=c('gene','pathway'))
head(gene.kegg.cor.res)
gene.kegg.cor.res=gene.kegg.cor.res[which(gene.kegg.cor.res$gene %in% colnames(tcga.module.gene.exp) & gene.kegg.cor.res$pathway %in% rownames(tcga.kegg)),]
head(gene.kegg.cor.res)
write.table(gene.kegg.cor.res,'results/gene.kegg.cor.res.txt',quote = F,row.names = F,sep = '\t')
gene.kegg.cor.res.fit=gene.kegg.cor.res[which(abs(gene.kegg.cor.res$cor)>0.6 & gene.kegg.cor.res$p<0.001),]
dim(gene.kegg.cor.res.fit)
sig_pathway<-as.character(unique(gene.kegg.cor.res.fit$pathway))
length(sig_pathway)#37
#相关性分析的热图
gene.kegg.cor1=gene.kegg.cor
rownames(gene.kegg.cor1$r)=gsub('KEGG_','',rownames(gene.kegg.cor1$r))
colnames(gene.kegg.cor1$r)=gsub('KEGG_','',colnames(gene.kegg.cor1$r))

rownames(gene.kegg.cor1$P)=gsub('KEGG_','',rownames(gene.kegg.cor1$P))
colnames(gene.kegg.cor1$P)=gsub('KEGG_','',colnames(gene.kegg.cor1$P))
sig_pathway1=gsub('KEGG_','',sig_pathway)


pdf('results/Fig1A.pdf',height =12,width = 9)
corrplot::corrplot(as.matrix(gene.kegg.cor1$r[colnames(tcga.module.gene.exp),sig_pathway1]), 
                   p.mat = as.matrix(gene.kegg.cor1$P[colnames(tcga.module.gene.exp),sig_pathway1]),
                   mar = c(0,0,1,1),
                   col=colorRampPalette(c('#4DBBD5FF', 'white','#E64B35FF'))(100),
                   tl.srt = 90,
                   tl.cex = 1,
                   tl.col = 'black',
                   tl.offset = 0.5,
                   cl.pos = c("b","r","n")[1], 
                   #cl.lim=c(-0.3,0.3),
                   cl.align.text = 'l',
                   cl.length = 5,
                   cl.ratio = 0.1,
                   cl.cex = 0.8,
                   addgrid.col = 'white',
                   method = 'color',
                   insig = 'label_sig',
                   sig.level=c(0.001,0.01,0.05),
                   pch.cex=1,
                   is.corr=T,
                   xpd=T)

dev.off()
#热图
anno_col=as.data.frame(tcga.module.gene.exp)
colnames(anno_col)
tcga.kegg1=tcga.kegg
tcga.kegg1[1:4,1:4]
rownames(tcga.kegg1)=gsub('KEGG_','',rownames(tcga.kegg1))
cols=ggsci::pal_jco()(9)
library(pheatmap)
pdf('results/Fig1B.pdf',height = 11,width = 12)
pheatmap(mat = as.matrix(tcga.kegg1[sig_pathway1,rownames(anno_col)]),
         scale = 'row',show_colnames = F,show_rownames = T,
         annotation_col = anno_col,
         cluster_cols = T,cluster_rows =T,
         color = colorRampPalette(c("#3C5488FF", "white","#E64B35FF"))(100)
         ,breaks = unique(c(seq(-2, 2, length=100))),
         annotation_colors =data.frame(BTK=colorRampPalette(c("white",cols[1]))(20),
                                       CLEC3B=colorRampPalette(c("white",cols[2]))(20),
                                       ANLN=colorRampPalette(c("white",cols[3]))(20),
                                       CD302=colorRampPalette(c("white",cols[4]))(20),
                                       CYP4B1=colorRampPalette(c("white",cols[5]))(20),
                                       ECT2=colorRampPalette(c("white",cols[6]))(20),
                                       GRIA1=colorRampPalette(c("white",cols[7]))(20),
                                       HMMR=colorRampPalette(c("white",cols[8]))(20)))
dev.off()


####10.gene.immnue####
setwd('~/LIHC/GEO/10.gene.immnue/')
dir.create('results')
#TCGA表达谱
tcga.dat<-read.delim('../03.data.pre/GSE53624/results/geo_T.txt',sep='\t',header = T,row.names = 1,check.names = F)
tcga.dat[1:4,1:4]
tcga.module<-read.delim('../06.CAF.module/results/tcga.risk.txt',sep='\t',header = T)
colnames(tcga.module)
tcga.module.gene.exp=tcga.module[,3:11]
rownames(tcga.module.gene.exp)=tcga.module$Sample
#Estimate
immu_estimate=function(exp,platform='illumina',isTCGA=T){#affymetrix
  library(estimate)
  inf=tempfile()
  ouf=tempfile()
  ouf2=tempfile()
  writeMatrix=function(dat,outpath,row=T,header=T){
    if(row){
      write.table(cbind(Tag=row.names(dat),dat),file=outpath,sep="\t",quote = F,row.names=F,col.names = header)
    }else{
      write.table(dat,file=outpath,sep="\t",quote = F,row.names=F,col.names = header)
    }
  }
  writeMatrix(exp,outpath = inf)
  outputGCT(inf, ouf)
  estimateScore(ouf,ouf2,platform=platform)
  est.score=t(read.csv(ouf2,sep = '\t',row.names = 1,check.names = F,skip = 2))
  est.score=est.score[-1,]
  rnames=row.names(est.score)
  if(isTCGA){
    rnames=gsub('\\.','-',row.names(est.score))
  }
  est.score=apply(est.score, 2, as.numeric)
  row.names(est.score)=rnames
  
  return(est.score)
}
tcga.esti=immu_estimate(exp = tcga.dat,platform = 'illumina',isTCGA = T)
save(tcga.esti,file = 'results/tcga.esti.RData')
cor_point<-function(dat,xlab,ylab,method,leg,cols){
  library(ggpubr)
  colnames(dat)=c('name1','name2','name3')
  p<-ggplot(dat, aes(x = name1, y = name2,color=name3))+
    geom_point()+
    geom_smooth(method = "lm")+
    stat_cor(method = "spearman",aes(x = name1, y = name2,  color = name3))+
    theme_bw()+
    theme(panel.grid.major =element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank())+
    xlab(xlab)+ylab(ylab)+scale_color_manual(values = cols)+
    labs(colour=leg)
  return(p)
}
load('results/tcga.esti.RData')
head(tcga.esti)
tcga.esti1=cbind.data.frame(tcga.esti,tcga.module.gene.exp[rownames(tcga.esti),])
tcga.esti1=reshape2::melt(tcga.esti1,id=colnames(tcga.esti))
head(tcga.esti1)
library(RColorBrewer)
fig1a<-cor_point(dat = tcga.esti1[,c("ImmuneScore","value","variable")],
                 xlab = 'ImmuneScore',ylab = 'Gene Expression',
                 method = 'spearman',leg = 'Gene',
                 cols = brewer.pal(n = 12, name = "Paired"))
fig1a <- fig1a + 
  theme(
    text = element_text(
      family = "",  # 清空自定义字体，使用系统默认字体
      size = 14, 
      color = "black"
    ),
    # 确保所有文本元素均使用系统默认字体
    axis.title = element_text(family = ""),
    axis.text = element_text(family = ""),
    legend.text = element_text(family = ""),
    legend.title = element_text(family = "")
  )
ggsave('results/Fig1a.pdf',height = 14,width = 14)
#
esti.gene=cbind.data.frame(tcga.esti,tcga.module.gene.exp[rownames(tcga.esti),])
esti.gene.cor <- Hmisc::rcorr(as.matrix(esti.gene))
pdf('results/Fig1a1.pdf',height = 5,width = 5)
corrplot::corrplot(as.matrix(esti.gene.cor$r[colnames(tcga.esti),colnames(tcga.module.gene.exp)]), 
                   p.mat = as.matrix(esti.gene.cor$P[colnames(tcga.esti),colnames(tcga.module.gene.exp)]),
                   mar = c(0,0,1,1),
                   col=colorRampPalette(c('#4DBBD5FF', 'white','#E64B35FF'))(100),
                   tl.srt = 60,
                   tl.cex = 1,
                   tl.col = 'black',
                   tl.offset = 0.5,
                   cl.pos = c("b","r","n")[1], 
                   cl.align.text = 'l',
                   cl.length = 5,
                   cl.ratio = 0.1,
                   cl.cex = 0.8,
                   addgrid.col = 'white',
                   method = 'color',
                   insig = 'label_sig',
                   sig.level=c(0.001,0.01,0.05),
                   pch.cex=1,
                   is.corr=T,
                   xpd=T)

dev.off()
#高低表达的比较
library(RColorBrewer)
head(esti.gene)
esti.gene1=esti.gene
esti.gene1 <- na.omit(esti.gene1)
colnames(esti.gene1)
esti.gene1$CD4=ifelse(esti.gene1$CD4>median(esti.gene1$CD4),'High','Low')
esti.gene1$SPP1=ifelse(esti.gene1$SPP1>median(esti.gene1$SPP1),'High','Low')
esti.gene1$TRIM22=ifelse(esti.gene1$TRIM22>median(esti.gene1$TRIM22),'High','Low')
esti.gene1$CD14=ifelse(esti.gene1$CD14>median(esti.gene1$CD14),'High','Low')
esti.gene1$ATP1B3=ifelse(esti.gene1$ATP1B3>median(esti.gene1$ATP1B3),'High','Low')
esti.gene1$COBLL1=ifelse(esti.gene1$COBLL1>median(esti.gene1$COBLL1),'High','Low')
esti.gene1$FABP5=ifelse(esti.gene1$FABP5>median(esti.gene1$FABP5),'High','Low')
esti.gene1$GPLD1=ifelse(esti.gene1$GPLD1>median(esti.gene1$GPLD1),'High','Low')
esti.gene1$LGALS3=ifelse(esti.gene1$LGALS3>median(esti.gene1$LGALS3),'High','Low')
esti.gene1$MMP12=ifelse(esti.gene1$MMP12>median(esti.gene1$MMP12),'High','Low')
esti.gene1$NECAB3=ifelse(esti.gene1$NECAB3>median(esti.gene1$NECAB3),'High','Low')

library(ggplot2)
library(ggsci)
head(esti.gene1)
sig_point<-function(dat,leg,ylab,palette=ggsci::pal_lancet()(12)[3:4]){
  library(ggpubr)
  dat=na.omit(dat)
  colnames(dat)=c('group','gene')
  dat=dat[order(dat$group),]
  all.combn=combn(as.character(unique(dat$group)),2)
  my_comparisons=lapply(seq_len(ncol(all.combn)), function(i) all.combn[,i])
  pp=ggplot(data = dat,aes(x = group, 
                           y = gene, 
                           fill = group))+ 
    scale_fill_manual(values = c("#4DBBD599","#E64B3599")) + #用自定义颜色填充
    geom_violin(alpha = 0.4, position = position_dodge(width = .75), 
                size = 0.8, color="black") +
    geom_boxplot(notch = FALSE, outlier.size = -1, 
                 color="black", lwd=0.8, alpha = 0.7) +
    geom_point(shape = 21, size=2, # 点的性状和大小
               position = position_jitterdodge(), # 让点散开
               color="black", alpha = 1) +
    theme_classic() + 
    ylab("ImmuneScore") +
    xlab("") +
    theme(axis.text.x = element_text(hjust = 1, size = 12,face = "bold.italic"),
          #panel.border = element_rect(colour = "black", fill=NA, size=0.2),
          axis.ticks = element_line(size=0.2, color="black"),
          axis.ticks.length = unit(0.2, "cm"),
          legend.position = "none",
          axis.title = element_text(size = 15,face = "bold.italic"),
          axis.text = element_text(size = 10)) + 
    stat_compare_means(comparisons=my_comparisons,method="wilcox.test",label = "p.signif")
  return(pp)
}

fig1b<-list()

fig1b[[1]] <- sig_point(dat = esti.gene1[, c("CD4", "ImmuneScore")],
                        leg = 'CD4',
                        ylab = 'ImmuneScore',
                        palette = ggsci::pal_lancet()(9)[3:4])

fig1b[[1]]
fig1b[[2]]<-sig_point(dat = esti.gene1[,c("SPP1","ImmuneScore")],
                      leg = 'SPP1',ylab = 'ImmuneScore',
                      palette = ggsci::pal_lancet()(9)[3:4])
fig1b[[2]]
fig1b[[3]]<-sig_point(dat = esti.gene1[,c("TRIM22","ImmuneScore")],
                      leg = 'TRIM22',ylab = 'ImmuneScore',
                      palette = ggsci::pal_lancet()(9)[3:4])
fig1b[[3]]
fig1b[[4]]<-sig_point(dat = esti.gene1[,c("CD14","ImmuneScore")],
                      leg = 'CD14',ylab = 'ImmuneScore',
                      palette = ggsci::pal_lancet()(9)[3:4])
fig1b[[4]]
fig1b[[5]]<-sig_point(dat = esti.gene1[,c("ATP1B3","ImmuneScore")],
                      leg = 'ATP1B3',ylab = 'ImmuneScore',
                      palette = ggsci::pal_lancet()(9)[3:4])
fig1b[[5]]
fig1b[[6]]<-sig_point(dat = esti.gene1[,c("COBLL1","ImmuneScore")],
                      leg = 'COBLL1',ylab = 'ImmuneScore',
                      palette = ggsci::pal_lancet()(9)[3:4])
fig1b[[6]]
fig1b[[7]]<-sig_point(dat = esti.gene1[,c("FABP5","ImmuneScore")],
                      leg = 'FABP5',ylab = 'ImmuneScore',
                      palette = ggsci::pal_lancet()(9)[3:4])
fig1b[[7]]
fig1b[[8]]<-sig_point(dat = esti.gene1[,c("GPLD1","ImmuneScore")],
                      leg = 'GPLD1',ylab = 'ImmuneScore',
                      palette = ggsci::pal_lancet()(9)[3:4])
fig1b[[8]]
fig1b[[9]]<-sig_point(dat = esti.gene1[,c("LGALS3","ImmuneScore")],
                      leg = 'LGALS3',ylab = 'ImmuneScore',
                      palette = ggsci::pal_lancet()(9)[3:4])
fig1b[[9]]
fig1b[[10]]<-sig_point(dat = esti.gene1[,c("MMP12","ImmuneScore")],
                       leg = 'MMP12',ylab = 'ImmuneScore',
                       palette = ggsci::pal_lancet()(9)[3:4])
fig1b[[10]]
fig1b[[11]]<-sig_point(dat = esti.gene1[,c("NECAB3","ImmuneScore")],
                       leg = 'NECAB3',ylab = 'ImmuneScore',
                       palette = ggsci::pal_lancet()(9)[3:4])
fig1b[[11]]





fig1b<-ggpubr::ggarrange(plotlist = fig1b,nrow = 4,ncol = 3)
ggsave('results/Fig1b.pdf',fig1b,height = 18,width = 18)
#CIBERSORT
source('CIBERSORT.R')
tcga.ciber=CIBERSORT("ref.txt", "~/LIHC/GEO/GSE14520-1.txt", perm=100, QN=TRUE)
write.table(tcga.ciber,'results/tcga.ciber.txt',quote = F,sep = '\t',row.names = T)
save(tcga.ciber,file = 'results/tcga.ciber.RData')
load('results/tcga.ciber.RData')
tcga.ciber[1:4,1:4]
#
ciber.gene=cbind.data.frame(tcga.ciber[,1:22],tcga.module.gene.exp[rownames(tcga.ciber),])
ciber.gene.cor <- Hmisc::rcorr(as.matrix(ciber.gene))
pdf('results/Fig1c.pdf',height = 5,width = 10)
corrplot::corrplot(as.matrix(ciber.gene.cor$r[colnames(tcga.ciber)[1:22],colnames(tcga.module.gene.exp)]), 
                   p.mat = as.matrix(ciber.gene.cor$P[colnames(tcga.ciber)[1:22],colnames(tcga.module.gene.exp)]),
                   mar = c(0,0,1,1),
                   col=colorRampPalette(c('blue', 'white','red'))(100),
                   tl.srt = 60,
                   tl.cex = 1,
                   tl.col = 'black',
                   tl.offset = 0.5,
                   cl.pos = c("b","r","n")[1], 
                   cl.align.text = 'l',
                   cl.length = 5,
                   cl.ratio = 0.1,
                   cl.cex = 0.8,
                   addgrid.col = 'white',
                   method = 'color',
                   insig = 'label_sig',
                   sig.level=c(0.001,0.01,0.05),
                   pch.cex=1,
                   is.corr=T,
                   xpd=T)
3
dev.off()
#mcpcount
install.packages("devtools")
library(devtools)
install_github("ebecht/MCPcounter",ref="master", subdir="Source")
library(MCPcounter)
immu_MCPcounter=function(exp,isTCGA=T){
  mcpEstimates=MCPcounter::MCPcounter.estimate(exp,featuresType="HUGO_symbols",
                                               probesets=read.table(curl("http://raw.githubusercontent.com/ebecht/MCPcounter/master/Signatures/probesets.txt"),sep='\t',stringsAsFactors=FALSE,colClasses='character'),
                                               genes=read.table(curl("http://raw.githubusercontent.com/ebecht/MCPcounter/master/Signatures/genes.txt"),sep='\t',stringsAsFactors=FALSE,header=TRUE,colClasses='character',check.names=FALSE))
  mcpEstimates=t(mcpEstimates)
  if(isTCGA){
    rnames=gsub('\\.','-',row.names(mcpEstimates)) 
    row.names(mcpEstimates)=rnames
  }
  return(mcpEstimates)
}
#download.file("https://raw.githubusercontent.com/becherlab/MCPcounter/master/immu_mcp_genes.txt", destfile = "immu_mcp_genes.txt")
tcga.mcp=immu_MCPcounter(exp = tcga.dat,isTCGA = T)

head(tcga.mcp)
library(ggcorrplot) 
mcp.gene=cbind.data.frame(tcga.mcp,tcga.module.gene.exp[rownames(tcga.mcp),])
write.table(mcp.gene,'results/mcp.gene.txt',quote = F,row.names = T,sep='\t')
pmtcars <- cor_pmat(mcp.gene)
cormtcars <- round(cor(mcp.gene), 3)

pdf('results/fig1d.pdf',height = 7,width = 7)
ggcorrplot(cormtcars[colnames(tcga.mcp),colnames(tcga.module.gene.exp)],
           hc.order = F,  #分等级聚类重排矩阵
           ggtheme = ggplot2::theme_void(base_size = 15), #主题修改
           colors = c('blue','white','red'), #自定义颜色，看自己喜欢，或是参考好看的文献Figure用法。
           lab = T,lab_size = 4,    #相关系数文本字体大小
           tl.cex = 10,             #坐标轴字体大小
           p.mat = pmtcars[colnames(tcga.mcp),colnames(tcga.module.gene.exp)],         #添加显著性信息
           sig.level = 0.01,        #显著性水平
           pch = 4,                 #不够显著的色块进行标记，pch表示选择不同的标记方法，可以尝试其他数字表示什么标记方法
           pch.cex = 10)  

dev.off()

####11.PDL####
setwd('~/LIHC/GEO/')
dir.create('results')
tcga.module<-read.delim('~/LIHC/GEO/06.nk.module/results/tcga_risk.txt',sep='\t',header = T)
tcga.module<- tcga_risk
colnames(tcga.module$result)
module.gene=colnames(tcga.module$result)[3:14]
module.gene
library(ggpubr)
library(colorspace)
sig_boxplot<-function(dat,leg,ylab,palette=ggsci::pal_lancet()(10)[3:4]){
  dat=na.omit(dat)
  colnames(dat)=c('group','gene')
  dat=dat[order(dat$group),]
  all.combn=combn(as.character(unique(dat$group)),2)
  my_comparisons=lapply(seq_len(ncol(all.combn)), function(i) all.combn[,i])
  pp=ggplot(dat,aes(x=group,y=gene,fill=group))+
    geom_violin(width =0.8,color='black',size=1)+
    theme_classic() + 
    theme(text = element_text(size=10, colour = "black")) + 
    theme(plot.title = element_text(hjust = 0.5, size = 15),
          axis.text.x = element_text(colour = "black", size = 12),
          axis.text.y = element_text(colour = "black", size = 10),
          axis.title.y = element_text(color = 'black', size = 12),
          axis.line = element_line(size = 1))+ 
    theme(legend.position="none") +  
    stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1), 
                 geom = "pointrange", color = "black", size=1)+
    scale_fill_manual(values = ggsci::pal_npg('nrc',alpha = 0.6)(4))+
    stat_compare_means(comparisons=my_comparisons,method="wilcox.test",label = "p.signif")+
    ylab(ylab)+xlab('')+labs(color=leg)
  return(pp)
}
#计算风险得分
library(survival)
get_riskscore<-function(dat,os,os.time,step=T,direction=c("both", "backward", "forward")[1]){
  crbind2DataFrame=function(dat){
    print(class(dat))
    # 检查 `dat` 的类中是否有 "table"
    if(any(class(dat) == 'table')){
      if(!is.na(ncol(dat))){
        dat=apply(dat, 2, function(x) {
          return(x)
        })
      }
    }
    
    # 使用 is.data.frame() 函数进行检查
    if (!is.data.frame(dat)) {
      dat1 = as.data.frame(as.matrix(dat))
    } else {
      dat1 = dat
    }
    
    for (i in 1:ncol(dat1)) {
      dat1[, i] = as.character(dat1[, i])
      dt = dat1[which(gsub(' ', '', dat1[, i]) != '' & !is.na(dat1[, i])), i]
      dt = dt[which(dt != 'Inf' & dt != 'NaN' & dt != 'NA')]
      if (sum(is.na(as.numeric(dt))) == 0) {
        dat1[, i] = as.numeric(dat1[, i])
      }
    }
    return(dat1)  
  }
  
  tcga_dat1 <- cbind(time=os.time,
                     status=os,
                     dat)
  tcga_dat1=crbind2DataFrame(tcga_dat1)
  colnames(tcga_dat1)=gsub('-','__',colnames(tcga_dat1))
  gene111=gsub('-','__',colnames(dat))
  fmla <- as.formula(paste0("Surv(time, status) ~"
                            ,paste0(gene111,collapse = '+')))
  cox <- coxph(fmla, data =as.data.frame(tcga_dat1))
  if(step==T){
    cox1 <- step(cox,direction =direction)
  }else{
    cox1=cox
  }
  lan <- coef(cox1)
  #round(lan, 3)
  genes <- names(cox1$coefficients)
  mult_results=paste0(round(lan, 3), '*', names(lan),collapse = '+')
  risk.tcga <- as.numeric(lan%*%as.matrix(t(tcga_dat1[,genes])))
  
  data_gene_score_final<-tcga_dat1
  data_gene_score_final$Samples<-rownames(data_gene_score_final)
  data_gene_score_final$riskscore=risk.tcga
  data_gene_score_final$riskscorez=mosaic::zscore(risk.tcga)
  optimalCutoff <- survminer::surv_cutpoint(data.frame(time=data_gene_score_final$time/365,
                                                       event=data_gene_score_final$status,
                                                       risk=data_gene_score_final$riskscore), 
                                            time = "time", event = "event", variables = c("risk"))
  optimalCutoff=optimalCutoff$cutpoint$cutpoint[1]
  #print(optimalCutoff)
  #optimalCutoff=median(data_gene_score_final$riskscore)
  #optimalCutoff=0
  data_gene_score_final$Risk=ifelse(data_gene_score_final$riskscore>optimalCutoff,'High','Low')
  table(data_gene_score_final$Risk)
  data_gene_score_final$cutoff=optimalCutoff
  return(list(result=data_gene_score_final,module.gene=cox1$coefficients,model=mult_results))
}
ggplotKM<-function(time,status,group,labs,palette){
  library(ggplot2)
  library(survival)
  dat1=data.frame(time=time,status=status,group=group)
  colnames(dat1)=c('time','status','groups')
  sf<-survival::survfit(Surv(time,status) ~ groups,data=dat1)
  surv=survminer::ggsurvplot(sf, data = dat1, 
                             palette = c("indianred1","#2E9FDF"), 
                             pval = TRUE,
                             surv.median.line='hv'
                             #,conf.int = T
                             ,conf.int.style ='step'
                             , pval.coord=c(0, 0.2), #Add p-value 
                             risk.table = TRUE, 
                             legend.title = 'Group'
                             ,legend.labs =labs
                             ,conf.int=T
  )
  p1=surv$plot+theme_bw()+
    theme(axis.text.y=element_text(family="Times",face="plain"),
          axis.text.x=element_blank(),
          axis.title.x=element_blank(),
          plot.margin=unit(c(0.2, 0.2, 0, 0.1), "inches"),
          legend.position=c(1,1),
          legend.justification=c(1,1),
          legend.background = element_rect(fill = NA, colour = NA),
          legend.title = element_text(family="Times",face="plain"),
          legend.text = element_text(family="Times",face="plain"))
  p2=surv$table+theme_bw()+
    theme(axis.text.y=element_text(family="Times",face="plain"),
          plot.margin=unit(c(0, 0.2, 0.2, 0.1), "inches"),
          plot.title=element_blank(),
          legend.position=c(1,1), 
          legend.justification=c(1,1),
          legend.title = element_text(family="Times",face="plain"),
          legend.text = element_text(family="Times",face="plain"))
  
  g2=ggpubr::ggarrange(p1,p2, ncol = 1, nrow = 2,heights = c(1,0.3),align = "v")
  return(g2)
}

#1、IMV210数据集####
imv210.exp<-read.delim('~/LIHC/GEO/11.PDL/IMvigor210/IMvigor210_Counts2TPM.txt',sep='\t',header = T)
IMvigor210_anno <- read.delim('~/LIHC/GEO/11.PDL/IMvigor210/IMvigor210_entrez2gene.txt', header = T)
imv210.exp=merge(IMvigor210_anno[,c("entrez_id","symbol")],
                 imv210.exp,by.x='entrez_id',by.y='genes')
imv210.exp[1:4,1:4]
imv210.exp=imv210.exp[,-1]
imv210.exp=aggregate(.~symbol,imv210.exp,mean)
imv210.exp[1:4,1:4]
rownames(imv210.exp)=imv210.exp$symbol
imv210.exp=imv210.exp[-1,-1]
range(imv210.exp)
imv210.exp <- log2(imv210.exp + 1)
IMvigor210_cli <- read.delim('~/LIHC/GEO/11.PDL/IMvigor210/IMvigor210_cli.txt', header = T)
IMvigor210_cli <- IMvigor210_cli[, c("os", "censOS",
                                     "Best.Confirmed.Overall.Response",
                                     "IC.Level", "TC.Level", "Immune.phenotype",
                                     "FMOne.mutation.burden.per.MB",
                                     "Neoantigen.burden.per.MB","TCGA.Subtype")]
colnames(IMvigor210_cli) <- c('OS.time', 'OS', 'Response', 
                              "IC.Level", "TC.Level", "Immune.phenotype",
                              'TMB', 'NEO','Stage')
table(IMvigor210_cli$Response)
IMvigor210_cli=IMvigor210_cli[which(IMvigor210_cli$Response != 'NE'),]
IMvigor210_cli$Response[IMvigor210_cli$Response=='CR'|IMvigor210_cli$Response=='PR']<-'CR/PR'
IMvigor210_cli$Response[IMvigor210_cli$Response=='PD'|IMvigor210_cli$Response=='SD']<-'PD/SD'
imv210.exp=imv210.exp[,rownames(IMvigor210_cli)]

imv.risk<-get_riskscore(dat = t(imv210.exp[intersect(rownames(imv210.exp),module.gene),]),
                        os = IMvigor210_cli$OS,
                        os.time = IMvigor210_cli$OS.time,
                        step = F,direction = 'both')

imv.risk$result$Risk=ifelse(imv.risk$result$riskscorez>0,'High','Low')
library(ggplot2)
library(ggsignif) 
fig1a<-ggplotKM(time = imv.risk$result$time/30,
                status =imv.risk$result$status,
                group = imv.risk$result$Risk,
                labs = c('High','Low'),
                palette = c("indianred1","#2E9FDF") )
fig1a
imv.risk.cli<-merge(imv.risk$result[,c("Samples","riskscorez","Risk")],
                    data.frame(Samples=rownames(IMvigor210_cli),IMvigor210_cli),
                    by='Samples')
write.table(x = imv.risk.cli,'results/imv.risk.cli.txt',quote = F,row.names = F,sep='\t')
fig1b=sig_boxplot(dat = imv.risk.cli[,c("Response","riskscorez")],
                  leg = 'Response',
                  ylab = 'RiskScore',
                  palette = ggsci::pal_nejm()(8))

fig1b <- fig1b +
  theme(axis.text.x = element_text(size = 14),  # X轴字体大小
        axis.text.y = element_text(size = 14),  # Y轴字体大小
        axis.title.x = element_text(size = 16), # X轴标题大小
        axis.title.y = element_text(size = 16), # Y轴标题大小
        panel.grid.major = element_line(size = 0.5),  # 主网格线的粗细
        panel.grid.minor = element_line(size = 0.25))+ # 次网格线的粗细
  geom_signif(comparisons = list(c("CR/PR", "PD/SD")), # 替换为你比较的组
              map_signif_level = TRUE, 
              size = 0.8,           # 星号线的粗细
              textsize = 5)        # 星号的大小

fig1b
response.risk=prop.table(table(imv.risk.cli$Response,imv.risk.cli$Risk),margin=2)

response.risk=reshape2::melt(response.risk)
colnames(response.risk)<-c("type","Risk","Percentage")

response.risk$Percentage<-round(response.risk$Percentage,digits=2)
write.table(response.risk,'results/response.risk.txt',quote = F,row.names = F,sep='\t')

fig1c=ggplot(response.risk,aes(x=Risk,y=Percentage,fill=type))+
  geom_bar(position = "fill",stat="identity")+
  theme_bw()+
  geom_text(aes(label = Percentage),position=position_stack(vjust =0.5),size = 5)+
  ggsci::scale_fill_nejm(alpha = 0.8)
fig1c
#I+II
imv.risk.cli1=imv.risk.cli[which(imv.risk.cli$Stage=='I'|imv.risk.cli$Stage =='II'),]
fig1d<-ggplotKM(time = imv.risk.cli1$OS.time/30,
                status =imv.risk.cli1$OS,
                group = imv.risk.cli1$Risk,
                labs = c('High','Low'),
                palette = c("indianred1","#2E9FDF") )
fig1d
imv.risk.cli2=imv.risk.cli[which(imv.risk.cli$Stage=='III'|imv.risk.cli$Stage=='IV'),]

fig1e<-ggplotKM(time = imv.risk.cli2$OS.time/30,
                status =imv.risk.cli2$OS,
                group = imv.risk.cli2$Risk,
                labs = c('High','Low'),
                palette = c("indianred1","#2E9FDF") )
fig1e
#GSE78220
GSE78220_cli <- read.delim('GSE78220/table_2.xls',sep=',',header = T)

GSE78220_cli1=data.frame(Samples=GSE78220_cli$Accession,
                         Title=GSE78220_cli$Title,
                         Rresponse=GSE78220_cli$anti.pd.1.response,
                         OS.time=GSE78220_cli$overall.survival..days.,
                         OS=GSE78220_cli$vital.status)
GSE78220_cli1 <- na.omit(GSE78220_cli1)
table(GSE78220_cli1$OS)
GSE78220_cli1$OS <- ifelse(GSE78220_cli1$OS == 'Alive', 0, 1)
rownames(GSE78220_cli1) <- GSE78220_cli1$Title

GSE78220_exp <- openxlsx::read.xlsx('GSE78220/GSE78220_PatientFPKM.xlsx',
                                    sheet = 1)
rownames(GSE78220_exp) <- GSE78220_exp$Gene
GSE78220_exp <- GSE78220_exp[, -1]
colnames(GSE78220_exp) <- stringr::str_split_fixed(colnames(GSE78220_exp), '\\.', 3)[, 1]
boxplot(GSE78220_exp[, 1:5])
GSE78220_exp=log2(GSE78220_exp+1)

gse.risk<-get_riskscore(dat = t(GSE78220_exp[intersect(rownames(GSE78220_exp),module.gene),GSE78220_cli1$Title]),
                        os = GSE78220_cli1$OS,
                        os.time = GSE78220_cli1$OS.time,
                        step = F,direction = 'both')
gse.risk$result$Risk=ifelse(gse.risk$result$riskscorez>0,'High','Low')
fig1f<-ggplotKM(time = gse.risk$result$time/30,
                status =gse.risk$result$status,
                group = gse.risk$result$Risk,
                labs = c('High','Low'),
                palette = c("indianred1","#2E9FDF") )
fig1f
gse.risk.cli<-merge(gse.risk$result[,c("Samples","riskscorez","Risk")],
                    data.frame(Samples=rownames(GSE78220_cli1),Response=GSE78220_cli1$Rresponse),
                    by='Samples')
table(gse.risk.cli$Response)
gse.risk.cli$Response[gse.risk.cli$Response=='Complete Response'|gse.risk.cli$Response=='Partial Response']<-'PR/CR'
gse.risk.cli$Response[gse.risk.cli$Response=='Progressive Disease']<-'PD'

write.table(x = gse.risk.cli,'results/gse.risk.cli.txt',quote = F,row.names = F,sep='\t')
gse.risk.cli$Response <- factor(gse.risk.cli$Response,levels = c('PR/CR','PD'))
fig1g=sig_boxplot(dat = gse.risk.cli[,c("Response","riskscorez")],
                  leg = 'Response',
                  ylab = 'RiskScore',
                  palette = ggsci::pal_nejm()(8))+
  theme(axis.text.x = element_text(size = 14),  # X轴字体大小
        axis.text.y = element_text(size = 14),  # Y轴字体大小
        axis.title.x = element_text(size = 16), # X轴标题大小
        axis.title.y = element_text(size = 16), # Y轴标题大小
        panel.grid.major = element_line(size = 0.5),  # 主网格线的粗细
        panel.grid.minor = element_line(size = 0.25))+ # 次网格线的粗细
  geom_signif(comparisons = list(c("CR/PR", "PD/SD")), # 替换为你比较的组
              map_signif_level = TRUE, 
              size = 0.8,           # 星号线的粗细
              textsize = 5)        # 星号的大小
fig1g
response.risk=prop.table(table(gse.risk.cli$Response,gse.risk.cli$Risk),margin=2)

response.risk=reshape2::melt(response.risk)
colnames(response.risk)<-c("type","Risk","Percentage")

response.risk$Percentage<-round(response.risk$Percentage,digits=2)
write.table(response.risk,'results/gse.response.risk.txt',quote = F,row.names = F,sep='\t')

fig1h=ggplot(response.risk,aes(x=Risk,y=Percentage,fill=type))+
  geom_bar(position = "fill",stat="identity")+
  theme_bw()+
  geom_text(aes(label = Percentage),position=position_stack(vjust =0.5),size = 5)+
  ggsci::scale_fill_nejm(alpha = 0.8)
fig1h
fig1<-ggarrange(fig1a,fig1b,fig1c,fig1d,
                fig1e,fig1f,fig1g,fig1h,
                nrow = 2,ncol = 4,labels = '')
fig1
ggsave('results/Fig1.pdf',fig1,height = 12,width = 20)

####单基因####
gene<- c("IRF9")#挑选进行生存分析的目标基因
imv210.exp[gene, ]
imv210.exp1 <- t(imv210.exp)
#在临床矩阵中新增gene列，根据中位数将表达量（连续型变量）分为高表达和低表达两组（分类型变量）：
IMvigor210_cli$gene <- ifelse(imv210.exp[gene,] > median(imv210.exp[gene,]),"High","Low")
IMvigor210_cli$gene[1:6]

imv210.exp[gene, ]

gene_expression <- as.numeric(unlist(imv210.exp[gene, ]))
gene_expression <- as.numeric(as.character(imv210.exp[gene, ]))
IMvigor210_cli$gene <- ifelse(gene_expression > median(gene_expression, na.rm = TRUE), "High", "Low")
library(survminer)
library(survival)
fit2 <- survfit(Surv(OS.time, OS) ~ gene, data = IMvigor210_cli)

ggsurvplot(
  fit2,
  data = IMvigor210_cli,
  censor.shape="|", censor.size = 4,
  conf.int = TRUE,
  conf.int.style = "ribbon",
  conf.int.alpha = 0.2,
  pval= TRUE,
  palette= "lancet",
  #surv.median.line = "hv",
  ggtheme=theme_bw(),
  legend= "top",
  legend.labs = c("High","Low"),
  xlab= "OS_time(days)",
  ylab= "Survival probablity",
  title= "Survival curves",
  break.x.by = 1000,
  break.y.by = 0.2,
  risk.table = TRUE,
  risk.table.col = "strata",
  risk.table.height = 0.2,
  risk.table.y.text = FALSE
)

ggsurvplot(
  fit2,
  data = IMvigor210_cli,
  censor.shape = "|", 
  censor.size = 4,
  conf.int = TRUE,
  conf.int.style = "ribbon",
  conf.int.alpha = 0.2,
  pval = TRUE,
  pval.coord = c(100, 0.25),  # 调整 p 值位置
  pval.size = 4,               # 调整 p 值字体大小
  palette = c("#E69F00", "#56B4E9"),  # 自定义颜色
  size = 1.2,                  # 调整生存曲线线宽
  ggtheme = theme_bw(base_size = 12),  # 调整主题字体大小
  legend = "top",
  legend.title = "Expression Level",  # 设置图例标题
  legend.labs = c("High", "Low"),
  xlab = "OS Time (days)",            # 设置 x 轴标签
  ylab = "Survival Probability",      # 设置 y 轴标签
  title = "Survival Curves by IRF9 Expression",  # 设置主标题
  subtitle = "IMvigor210 Cohort",     # 添加副标题
  break.x.by = 30,                   # 调整 x 轴断点间隔
  break.y.by = 0.1,                   # 调整 y 轴断点间隔
  risk.table = TRUE,
  risk.table.title = "Number at Risk",  # 设置风险表标题
  risk.table.col = "strata",
  risk.table.height = 0.25,            # 调整风险表高度
  risk.table.y.text = FALSE,
  risk.table.y.text.col = TRUE,        # 为风险表 y 轴文本添加颜色
  risk.table.fontsize = 3.5            # 调整风险表字体大小
)


####12.CIBERSORT####
setwd('~/LIHC/GEO/12.CIBERSORT/')
dir.create('results')
install.packages('e1071')

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("preprocessCore")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

library("limma")         #引用包
expFile="uniq.symbol.txt"     #表达输入文件

#读取输入文件，并对输入文件整理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

#数据转换
v=voom(data, plot=F, save.plot=F)
out=v$E
out=rbind(ID=colnames(out), out)
write.table(out,file="uniq.symbol.txt",sep="\t",quote=F,col.names=F)        #输出文件

#运行CIBERSORT，得到免疫细胞浸润的结果
source("IRGPI28.CIBERSORT.R")
results=CIBERSORT("ref.txt", "uniq.symbol.txt", perm=1000, QN=TRUE)

####13.immCor####
setwd('~/LIHC/GEO/13.immCor/')
dir.create('results')

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("reshape2")
#install.packages("ggpubr")


#引用包
library(limma)
library(reshape2)
library(ggpubr)

riskFile="tcga_risk.csv"            #风险文件
immFile="CIBERSORT-Results.txt"     #免疫细胞浸润结果文件
pFilter=0.05                        #免疫细胞浸润结果的过滤条件

#读取免疫细胞结果文件，并对数据进行整理
immune=read.table(immFile, header=T, sep="\t", check.names=F, row.names=1)
immune=immune[immune[,"P-value"]<pFilter,]
data=as.matrix(immune[,1:(ncol(immune)-3)])

#删除正常样品
group=sapply(strsplit(row.names(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
data=data[group==0,]
row.names(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", row.names(data))
data=avereps(data)

#读取风险文件
risk=read.csv(riskFile)
sameSample=intersect(row.names(data), row.names(risk))
rt=cbind(data[sameSample,,drop=F], risk[sameSample,"Risk",drop=F])
rt=rt[order(rt$Risk, decreasing=T),]
conNum=nrow(rt[rt$Risk=="Low",])
treatNum=nrow(rt[rt$Risk=="High",])

##########绘制柱状图##########
data=t(rt[,-ncol(rt)])
pdf("barplot.pdf", height=10, width=18)
col=rainbow(nrow(data), s=0.7, v=0.7)
par(las=1,mar=c(8,5,4,16),mgp=c(3,0.1,0),cex.axis=1.5)
a1=barplot(data,col=col,yaxt="n",ylab="Relative Percent",xaxt="n",cex.lab=1.8)
a2=axis(2,tick=F,labels=F)
axis(2,a2,paste0(a2*100,"%"))
par(srt=0,xpd=T)
rect(xleft = a1[1], ybottom = -0.01, xright = a1[conNum], ytop = -0.06,col="green")
text(a1[conNum]/2,-0.035,"Low risk",cex=2)
rect(xleft = a1[conNum], ybottom = -0.01, xright =a1[length(a1)] , ytop = -0.06,col="red")
text((a1[length(a1)]+a1[conNum])/2,-0.035,"High risk",cex=2)
ytick2 = cumsum(data[,ncol(data)])
ytick1 = c(0,ytick2[-length(ytick2)])
legend(par('usr')[2]*0.98,par('usr')[4],legend=rownames(data),col=col,pch=15,bty="n",cex=1.3)
dev.off()

##################绘制箱线图##################
#把数据转换成ggplot2输入文件
data=rt
data=melt(data, id.vars=c("Risk"))
colnames(data)=c("Risk", "Immune", "Expression")
#绘制箱线图
group=levels(factor(data$Risk))
data$Risk=factor(data$Risk, levels=c("Low","High"))
bioCol=c("#0066FF","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(group)]
boxplot=ggboxplot(data, x="Immune", y="Expression", fill="Risk",
                  xlab="",
                  ylab="Fraction",
                  legend.title="Risk",
                  width=0.8,
                  palette=bioCol)+
  rotate_x_text(50)+
  stat_compare_means(aes(group=Risk),symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1), symbols=c("***", "**", "*", "")), label="p.signif")
#输出图片
pdf(file="immune.diff.pdf", width=8, height=6)
print(boxplot)
dev.off()

####14.immFunction####
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("GSVA")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("GSEABase")

#install.packages("ggpubr")
#install.packages("reshape2")
setwd('~/LIHC/GEO/14.immFunction/')
dir.create('results')

#引用包
library(limma)
library(GSVA)
library(GSEABase)
library(ggpubr)
library(reshape2)

expFile="symbol.txt"           #表达数据文件
gmtFile="immune.gmt" #免疫数据集文件
riskFile <- tcga_risk
riskFile="risk.TCGA.txt"       #风险文件
setwd("C:\\biowolf\\IRGPI\\31.immFunction")   #设置工作目录

#读取表达输入文件，并对输入文件处理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
mat=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
mat=avereps(mat)
mat=mat[rowMeans(mat)>0,]

#读取免疫基因集文件
geneSet=getGmt(gmtFile, geneIdType=SymbolIdentifier())

#ssgsea分析
ssgseaScore=gsva(mat, geneSet, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)
#定义ssGSEA score矫正函数
normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}
#对ssGSEA score进行矫正
data=normalize(ssgseaScore)
ssgseaOut=rbind(id=colnames(data), data)
write.table(ssgseaOut, file="immFunScore.txt", sep="\t", quote=F, col.names=F)

#去除正常样品
group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
data=t(data[,group==0])
rownames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(data))
data=avereps(data)

#读取风险文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
risk <- riskFile
#合并数据
data <- t(data)
sameSample=intersect(row.names(data),row.names(risk))
data=data[sameSample,,drop=F]
risk=risk[sameSample,"Risk",drop=F]
rt1=cbind(data, risk)

#对免疫相关功能绘制箱线图
data=melt(rt1, id.vars=c("Risk"))
colnames(data)=c("Risk","Type","Score")
data$Risk=factor(data$Risk, levels=c("Low","High"))
p=ggboxplot(data, x="Type", y="Score", color = "Risk",
            xlab="",ylab="Score",add = "none",palette = c("blue","red") )
p=p+rotate_x_text(50)
p=p+stat_compare_means(aes(group=Risk),symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "")),label = "p.signif")

#输出图片文件
pdf(file="immFunction.pdf", width=8, height=6)
print(p)
dev.off()

####15.clusters####
setwd('~/LIHC/GEO/15. clusters/')
dir.create('results')
library(limma)
library(ConsensusClusterPlus)

rt=read.csv("tcga_risk.csv")
rownames(rt) <- rt$Sample
data=rt[,(4:15)]
data=t(data)

maxK=9
workDir <- "~/LIHC/GEO/15. clusters/results/"
results=ConsensusClusterPlus(data,
                             maxK=maxK,
                             reps=50,
                             pItem=0.8,
                             pFeature=1,
                             title=workDir,
                             clusterAlg="km",
                             distance="euclidean",
                             seed=123456,
                             plot="png")


clusterNum=3    
Cluster=results[[clusterNum]][["consensusClass"]]
outTab=cbind(rt, Cluster)
outTab[,"Cluster"]=paste0("C", outTab[,"Cluster"])
outTab=rbind(ID=colnames(outTab), outTab)
write.table(outTab, file="cluster.txt", sep="\t", quote=F, col.names=F)

#第二步
library(survival)
library(survminer)
library(tibble)
clusterFile="cluster.txt"
setwd("~/LIHC/GEO/15. clusters/results/")  
rt$time <- rt$time/12

rt=read.table(clusterFile, header=T, sep="\t", check.names=F)
rt <- column_to_rownames(rt,'ID')

length=length(levels(factor(rt$Cluster)))
diff=survdiff(Surv(time, status) ~ Cluster, data = rt)
pValue=1-pchisq(diff$chisq, df=length-1)
if(pValue<0.001){
  pValue="p<0.001"
}else{
  pValue=paste0("p=",sprintf("%.03f",pValue))
}
fit <- survfit(Surv(time, status) ~ Cluster, data = rt)
print(surv_median(fit))

bioCol=c("#0066FF","#FF9900","#00DB00","#FF0000","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[6:8]
library(survival)
library(survminer)
surPlot=ggsurvplot(fit, 
                   data=rt,
                   conf.int=T,
                   pval=pValue,
                   pval.size=10,
                   legend.title="Cluster",
                   legend.labs=levels(factor(rt[,"Cluster"])),
                   legend = c(0.8, 0.8),
                   font.legend=10,
                   xlab="Time(years)",
                   xlab.size = 14,  # Increase x-axis label font size
                   ylab = "Survival Probability",  # Add y-axis label
                   ylab.size = 20,  # Increase y-axis label font size
                   break.time.by = 1,
                   palette = c('#E64B35FF','#377EB8','#FF7F00'),
                   surv.median.line="hv",
                   risk.table=T,
                   cumevents=F,
                   risk.table.height=.25,
                   font.title = 14, # Increase plot title font size
                   font.axis = 14) 

pdf(file="survival.pdf",onefile = FALSE,width=17,height=15.5)
print(surPlot)
dev.off()

# 3.ggalluvial
library(ggalluvial)
library(ggplot2)
library(dplyr)

clusterFile="cluster.txt" 
setwd("C:\\Users\\ZZD\\Desktop\\clusters\\34.ggalluvial") 


rt=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)


rt=rt[,c("Risk", "Cluster")]
colnames(rt)=c("Risk", "Cluster")
corLodes=to_lodes_form(rt, axes = 1:ncol(rt), id = "Cohort")

pdf(file="ggalluvial.pdf", width=7, height=5)
mycol=rep(c("#FF0000","#0000FF","#0066FF","#FF9900","#00DB00","#FF0000","#D8D155","#223D6C","#D20A13","#431A3D","#91612D","#FFD121","#088247","#11AA4D","#58CDD9","#7A142C","#5D90BA","#64495D","#7CC767"),15)
ggplot(corLodes, aes(x = x, stratum = stratum, alluvium = Cohort,fill = stratum, label = stratum)) +
  scale_x_discrete(expand = c(0, 0)) +  
  #??aes.flow??????????ɫ??forward˵??????????ɫ??ǰ??????״ͼһ?£?backward˵???ͺ???????״ͼһ?¡?
  geom_flow(width = 2/10,aes.flow = "backward") + 
  geom_stratum(alpha = .9, width = 3/10) +
  scale_fill_brewer(palette = 'Set2') +
  #size=3??????????С
  geom_text(stat = "stratum", size = 3,color="black") +
  xlab("") + ylab("") + theme_bw() + 
  theme(axis.line = element_blank(),axis.ticks = element_blank(),axis.text.x = element_blank()) + #ȥ????????
  theme(panel.grid =element_blank()) + 
  theme(panel.border = element_blank()) + 
  coord_flip() + ggtitle("") + guides(fill = FALSE)                            
dev.off()

# 4.PCA
library(Rtsne)
library(ggplot2)

clusterFile="cluster.txt"  
setwd("C:\\Users\\ZZD\\Desktop\\clusters\\35.PCA")   


rt=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)
rt = rt[,-1]
data=rt[c(3:(ncol(rt)-6))]

#PCA????
data.pca=prcomp(data, scale. = TRUE)
pcaPredict=predict(data.pca)
PCA = data.frame(PC1=pcaPredict[,1], PC2=pcaPredict[,2], rt[,c("Risk","Cluster")])
#???Ʒ??յ?PCAͼ
pdf(file="PCA.risk.pdf", width=5.5, height=4.5)
p=ggplot(data = PCA, aes(PC1, PC2)) + geom_point(aes(color = Risk)) +
  scale_colour_manual(name="Risk",  values =c("red", "blue"))+
  theme_bw()+
  theme(plot.margin=unit(rep(1.5,4),'lines'))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(p)
dev.off()
#???Ʒ??͵?PCAͼ
bioCol=c("#0066FF","#FF9900","#00DB00","#FF0000","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(unique(rt$Cluster))]
pdf(file="PCA.cluster.pdf", width=5.5, height=4.5)
p=ggplot(data = PCA, aes(PC1, PC2)) + geom_point(aes(color = Cluster)) +
  scale_colour_manual(name="Cluster",  values =bioCol)+
  theme_bw()+
  theme(plot.margin=unit(rep(1.5,4),'lines'))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(p)
dev.off()


#t-SNE????
tsneOut=Rtsne(data, dims=2, perplexity=10, verbose=F, max_iter=500,check_duplicates=F)
tsne=data.frame(tSNE1=tsneOut$Y[,1], tSNE2=tsneOut$Y[,2], rt[,c("Risk","Cluster")])	
#???Ʒ??յ?tSNEͼ
pdf(file="tSNE.risk.pdf", width=5.5, height=4.5)       #???????????ļ?
p=ggplot(data = tsne, aes(tSNE1, tSNE2)) + geom_point(aes(color = Risk)) +
  scale_colour_manual(name="Risk",  values =c("red", "blue"))+
  theme_bw()+
  theme(plot.margin=unit(rep(1.5,4),'lines'))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(p)
dev.off()
#???Ʒ??͵?tSNEͼ
pdf(file="tSNE.Cluster.pdf", width=5.5, height=4.5)       #???????????ļ?
p=ggplot(data = tsne, aes(tSNE1, tSNE2)) + geom_point(aes(color = Cluster)) +
  scale_colour_manual(name="Cluster",  values =bioCol)+
  theme_bw()+
  theme(plot.margin=unit(rep(1.5,4),'lines'))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(p)
dev.off()

# 5.clusterTMEdiff
library(limma)
library(ggpubr)

cluFile="cluster.txt"     
scoreFile="risk.estimate.txt"     
setwd("C:\\Users\\ZZD\\Desktop\\clusters\\36.clusterTMEdiff")   


score=read.table(scoreFile, header=T, sep="\t", check.names=F, row.names=1)
score <- tcga.esti
score=as.matrix(score)
row.names(score)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", row.names(score))
score=avereps(score)


cluster=read.table(cluFile, header=T, sep="\t", check.names=F, row.names=1)


sameSample=intersect(row.names(score), row.names(cluster))
score1=score[sameSample,,drop=F]
cluster1=cluster[sameSample,"Cluster",drop=F]
data=cbind(score1, cluster1)
library(ggrain)
#install.packages("ggrain")

type=levels(factor(data[,"Cluster"]))
data$Cluster=factor(data$Cluster, levels=type)
comp=combn(type, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}


bioCol=c("#4DBBD599","#E64B3599",'#FF7F00')
bioCol=bioCol[1:length(unique(data$Cluster))]


library(ggplot2)
library(ggpubr)  # Ensure you load ggpubr for stat_compare_means

for(i in colnames(data)[1:(ncol(data)-1)]) {
  boxplot <- ggplot(data, aes(x = Cluster, 
                              y = .data[[i]],  # Use .data to refer to variable names
                              fill = Cluster)) + 
    scale_fill_manual(values = c("#4DBBD599","#E64B3599", '#FF7F00')) + 
    geom_violin(alpha = 0.4, position = position_dodge(width = .75), 
                linewidth = 0.8, color = "black") +
    geom_boxplot(notch = FALSE, outlier.size = -1, 
                 color = "black", linewidth = 0.8, alpha = 0.7) +
    geom_point(shape = 21, size = 2, 
               position = position_jitterdodge(), 
               color = "black", alpha = 1) +
    theme_classic() + 
    ylab("") +
    xlab("") +
    theme(axis.text.x = element_text(hjust = 1, size = 12, face = "bold.italic"),
          axis.ticks = element_line(linewidth = 0.2, color = "black"),
          axis.ticks.length = unit(0.2, "cm"),
          legend.position = "none",
          axis.title = element_text(size = 15, face = "bold.italic"),
          axis.text = element_text(size = 10)) +  
    stat_compare_means(comparisons = my_comparisons)  # Ensure 'my_comparisons' is valid
  
  pdf(file = paste0(i, ".pdf"), width = 5, height = 4.5)
  print(boxplot)
  dev.off()
}

# 6.clusterImmCor
library(limma)
library(pheatmap)

clusterFile="cluster.txt"  
immFile="CIBERSORT-Results.txt"   
setwd("C:\\Users\\ZZD\\Desktop\\clusters\\37.clusterImmCor")  

#??ȡ???ͽ????ļ?
cluster=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)


immune=read.table(immFile, header=T, sep="\t", check.names=F, row.names=1)
immune=as.matrix(immune)
rownames(immune)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*)", "\\1\\-\\2\\-\\3", rownames(immune))
immune=avereps(immune)
immune <- immune[,-3]
#???˷??ͺ?????ϸ???ϲ?
sameSample=intersect(row.names(cluster), row.names(immune))
cluster=cluster[sameSample, "Cluster", drop=F]
immune=immune[sameSample, , drop=F]
data=cbind(cluster, immune)
data <- data[,-c(24:25)]
#???͵?????ϸ??????????
outTab=data.frame()
sigCell=c("Cluster")
for(i in colnames(data)[2:ncol(data)]){
  if(sd(data[,i])<0.001){next}
  if(length(levels(factor(data[,"Cluster"])))>2){
    test=kruskal.test(data[,i] ~ data[,"Cluster"])
  }else{
    test=wilcox.test(data[,i] ~ data[,"Cluster"])
  }
  pvalue=test$p.value
  if(pvalue<0.05){
    outTab=rbind(outTab,cbind(immune=i, pvalue))
    sigCell=c(sigCell, i)
  }
}
write.table(file="immuneCor.txt", outTab, sep="\t", quote=F, row.names=F)

#??ͼ????
data=data[,sigCell]
#???͵?ע??
data$Cluster <- factor(data$Cluster,levels = c('C3','C1','C2'))
data=data[order(data[,"Cluster"]),]

annCol=data[,1,drop=F]
annCol[,"Cluster"]=factor(annCol[,"Cluster"], unique(annCol[,"Cluster"]))
data=t(data[,(2:ncol(data))])
#???????͵?ע??
annRow=sapply(strsplit(rownames(data),"_"), '[', 2)
annRow=as.data.frame(annRow)
row.names(annRow)=row.names(data)
colnames(annRow)=c("Methods")
annRow[,"Methods"]=factor(annRow[,"Methods"], unique(annRow[,"Methods"]))
gapCol=as.vector(cumsum(table(annCol[,"Cluster"])))
gapRow=as.vector(cumsum(table(annRow[,"Methods"])))

#??????ͼע?͵???ɫ
bioCol=c("#0066FF","#FF9900","#00DB00","#FF0000","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(unique(annCol[,"Cluster"]))]
bioCol <- c("#D20A13","#FFD121",'#11AA4D')
Cluster=bioCol
names(Cluster)=levels(factor(annCol[,"Cluster"]))
ann_colors=list(Cluster=Cluster)

#??ͼ???ӻ?
pdf("immHeatmap.pdf", width=9, height=6)
pheatmap(data,
         annotation=annCol,
         annotation_colors = ann_colors,
         color = colorRampPalette(c(rep("#4DBBD5FF",5), "white", rep("#E64B35FF",5)))(100),
         cluster_cols =T,
         cluster_rows =F,
         gaps_col=gapCol,
         scale="row",
         show_colnames=F,
         show_rownames=T,
         fontsize=6,
         fontsize_row=5,
         fontsize_col=6)
dev.off()

# 7.clusterCheckpoint
library(limma)
library(reshape2)
library(ggplot2)
library(ggpubr)

expFile="geo_dat.txt"           #?????????ļ?
ClusterFile="cluster.txt"      #???ͽ????ļ?
geneFile="gene.txt"            #???߼??????Ļ????б??ļ?
setwd("C:\\Users\\ZZD\\Desktop\\clusters\\38.clusterCheckpoint")     #???ù???Ŀ¼

rt <- expFile
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(rt),colnames(rt))
data=matrix(as.numeric(as.matrix(rt)),nrow=nrow(rt),dimnames=dimnames)
data=avereps(data)

#??ȡ?????б??ļ?????ȡ???߼????????ػ????ı???��
gene=read.table(geneFile, header=F, sep="\t", check.names=F)
sameGene=intersect(row.names(data), as.vector(gene[,1]))
data=t(data[sameGene,])
data=log2(data+1)

#ɾdata <- data[seq(1,238,by = 2),]????????Ʒ
group=sapply(strsplit(row.names(data),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
data=data[group==0,]
row.names(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", row.names(data))
data=avereps(data)

#?ϲ?????
Cluster=read.table(ClusterFile, sep="\t", header=T, check.names=F, row.names=1)
sameSample=intersect(row.names(data), row.names(Cluster))
rt1=cbind(data[sameSample,], Cluster[sameSample,])
rt1=rt1[,c(sameGene, "Cluster")]

#??ȡ?????????Ļ???
sigGene=c()
for(i in colnames(rt1)[1:(ncol(rt1)-1)]){
  if(sd(rt1[,i])<0.001){next}
  if(length(levels(factor(rt1[,"Cluster"])))>2){
    test=kruskal.test(rt1[,i] ~ rt1[,"Cluster"])
  }else{
    test=wilcox.test(rt1[,i] ~ rt1[,"Cluster"])
  }
  pvalue=test$p.value
  if(pvalue<0.05){
    sigGene=c(sigGene, i)
  }
}
sigGene=c(sigGene, "Cluster")
rt1=rt1[,sigGene]

#??????ת????ggplot2?????ļ?
rt1=melt(rt1,id.vars=c("Cluster"))
colnames(rt1)=c("Cluster","Gene","Expression")

#????????ͼ
bioCol=c("#0066FF","#FF9900","#00DB00","#FF0000","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(unique(rt1$Cluster))]
boxpbioCol <- c("#D20A13","#11AA4D","#FFD121")

lot = ggboxplot(rt1, x="Gene", y="Expression", fill="Cluster",
                xlab="",
                ylab="Gene expression",
                legend.title="Cluster",
                width=0.8,
                palette = bioCol) +
  rotate_x_text(50) +
  stat_compare_means(aes(group=Cluster),
                     symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1), 
                                      symbols=c("***", "**", "*", "")), 
                     label="p.signif")
#????ͼƬ
pdf(file="checkpoint.diff.pdf", width=17.6, height=15)
print(lot)
dev.off()

##8.pRRophetic
BiocManager::install(c("car", "ridge", "preprocessCore", "genefilter", "sva"))
install.packages("pRRophetic_0.5.tar.gz", repos = NULL, type = "source")
library(limma)
library(ggpubr)
library(pRRophetic)
library(ggplot2)
set.seed(12345)

pFilter=0.001                
expFile="data.txt"  
ClusterFile="cluster.txt"   
setwd("C:\\Users\\ZZD\\Desktop\\套路七单细胞+成纤维+免疫治疗\\clusters\\39.pRRophetic")     #???ù???Ŀ¼
alldrugs = c("A.443654","A.770041","ABT.263","ABT.888","AG.014699","AICAR","AKT.inhibitor.VIII","AMG.706","AP.24534","AS601245","ATRA","AUY922","Axitinib","AZ628","AZD.0530","AZD.2281","AZD6244","AZD6482","AZD7762","AZD8055","BAY.61.3606","Bexarotene","BI.2536","BIBW2992","Bicalutamide","BI.D1870","BIRB.0796","Bleomycin","BMS.509744","BMS.536924","BMS.708163","BMS.754807","Bortezomib","Bosutinib","Bryostatin.1","BX.795","Camptothecin","CCT007093","CCT018159","CEP.701","CGP.082996","CGP.60474","CHIR.99021","CI.1040","Cisplatin","CMK","Cyclopamine","Cytarabine","Dasatinib","DMOG","Docetaxel","Doxorubicin","EHT.1864","Elesclomol","Embelin","Epothilone.B","Erlotinib","Etoposide","FH535","FTI.277","GDC.0449","GDC0941","Gefitinib","Gemcitabine","GNF.2","GSK269962A","GSK.650394","GW.441756","GW843682X","Imatinib","IPA.3","JNJ.26854165","JNK.9L","JNK.Inhibitor.VIII","JW.7.52.1","KIN001.135","KU.55933","Lapatinib","Lenalidomide","LFM.A13","Metformin","Methotrexate","MG.132","Midostaurin","Mitomycin.C","MK.2206","MS.275","Nilotinib","NSC.87877","NU.7441","Nutlin.3a","NVP.BEZ235","NVP.TAE684","Obatoclax.Mesylate","OSI.906","PAC.1","Paclitaxel","Parthenolide","Pazopanib","PD.0325901","PD.0332991","PD.173074","PF.02341066","PF.4708671","PF.562271","PHA.665752","PLX4720","Pyrimethamine","QS11","Rapamycin","RDEA119","RO.3306","Roscovitine","Salubrinal","SB.216763","SB590885","Shikonin","SL.0101.1","Sorafenib","S.Trityl.L.cysteine","Sunitinib","Temsirolimus","Thapsigargin","Tipifarnib","TW.37","Vinblastine","Vinorelbine","Vorinostat","VX.680","VX.702","WH.4.023","WO2009093972","WZ.1.84","X17.AAG","X681640","XMD8.85","Z.LLNle.CHO","ZM.447439")
rt = read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
rt=rt[,-1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0.5,]
dimnames=list(rownames(data), colnames(data))
data=matrix(as.numeric(as.matrix(data)), nrow=nrow(data), dimnames=dimnames)
data=data[rowMeans(data)>0.5,]

group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2","1",group)
data=data[,group==0]
data=t(data)
rownames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*)", "\\1\\-\\2\\-\\3", rownames(data))
data=avereps(data)
data=t(data)
#data行名为基因，不是样本名
data  <- as.matrix(data)
ClusterRT=read.table(ClusterFile, header=T, sep="\t", check.names=F, row.names=1)

library(BiocParallel)

# 注册 SerialParam
register(SerialParam())

# 测试并行计算
test_result <- bplapply(1:10, function(x) x^2)

for(drug in allDrugs){
  #Ԥ??ҩ????????
  senstivity=pRRopheticPredict(data, drug, selection=1)
  senstivity=senstivity[senstivity!="NaN"]
  #senstivity[senstivity>quantile(senstivity,0.99)]=quantile(senstivity,0.99)
  
  #?????ͽ????ļ???ҩ???????Խ??????кϲ?
  sameSample=intersect(row.names(ClusterRT), names(senstivity))
  Cluster=ClusterRT[sameSample, "Cluster",drop=F]
  senstivity=senstivity[sameSample]
  rt=cbind(Cluster, senstivity)
  
  #???ñȽ???
  type=levels(factor(rt[,"Cluster"]))
  comp=combn(type, 2)
  my_comparisons=list()
  for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
  
  #??ȡ????֮??????pvalue
  if(length(levels(factor(rt[,"Cluster"])))>2){
    test=kruskal.test(senstivity~Cluster, data=rt)
  }else{
    test=wilcox.test(senstivity~Cluster, data=rt)
  }
  pvalue=test$p.value
  if(pvalue<pFilter){
    #????????ͼ
    boxplot=boxplot=ggplot(data = rt,aes(x = Cluster, 
                                         y = senstivity , 
                                         fill = Cluster))+ 
      scale_fill_manual(values = c('#E64B35FF','#4DBBD5FF','#3C5488FF')) +
      geom_violin(alpha = 0.4, position = position_dodge(width = .75), 
                  size = 0.8, color="black") +
      geom_boxplot(notch = TRUE, outlier.size = -1, 
                   color="black", lwd=0.8, alpha = 0.7) +
      geom_point(shape = 21, size=2, 
                 position = position_jitterdodge(), 
                 color="black", alpha = 1) +
      theme_bw() + 
      ylab(paste0(drug, " senstivity (IC50)")) +
      xlab('Cluster') +
      theme(axis.text.x = element_text(size = 12, color = "black"),
            axis.ticks = element_line(size=0.2, color="black"),
            axis.ticks.length = unit(0.2, "cm"),
            legend.position = "none",
            panel.background = element_blank(),
            panel.grid = element_blank(),
            axis.title = element_text(size = 12),
            axis.text = element_text(size = 12)) + 
      stat_compare_means(comparisons=my_comparisons)
    pdf(file=paste0("durgSenstivity.", drug, ".pdf"), width=5, height=4.5)
    print(boxplot)
    dev.off()
  }
}

library(BiocParallel)

# 注册 SerialParam
register(SerialParam())

# 测试并行计算
test_result <- bplapply(1:10, function(x) x^2)


print(test_result)


trace(calcPhenotype, edit = T)
trace(summarizeGenesByMean, edit = T)
q()
y
n
q()
y

