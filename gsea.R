## GSEA

rm(list = ls())
library(GSEABase)
library(fgsea)
library(org.Hs.eg.db)
library(clusterProfiler)
library(dplyr)
library(ggplot2)
library(readxl)
library(enrichplot)
library(msigdbr)  #提供MSigdb数据库基因集

## 物种设置
organism = 'mmu'    #  人类'hsa' 小鼠'mmu'   
OrgDb = 'org.Mm.eg.db'#人类"org.Hs.eg.db" 小鼠"org.Mm.eg.db"

gene = read_excel("w8-vs-sham-diff-q-val-0.05-FC-2.gene.xlsx")
# 保留gene name和logFC情况
ex_deg <- gene[,c(6,1)]
ex_deg <- data.frame(ex_deg)
need_DEG <- ex_deg

colnames(need_DEG) <- c('log2FoldChange','SYMBOL')
rownames(need_DEG) <- need_DEG$SYMBOL
##### 创建gsea分析的geneList（包含从大到小排列的log2FoldChange和ENTREZID信息）####
#转化id  
df <- bitr(rownames(need_DEG), 
           fromType = "SYMBOL",
           toType =  "ENTREZID",
           OrgDb = OrgDb) #人数据库org.Hs.eg.db 小鼠org.Mm.eg.db
need_DEG <- merge(need_DEG, df, by='SYMBOL')  #按照SYMBOL合并注释信息
geneList <- need_DEG$log2FoldChange
names(geneList) <- need_DEG$ENTREZID
geneList <- sort(geneList, decreasing = T)   #从大到小排序

#-------------------------------------------------
# go
GO_kk_entrez <- gseGO(geneList     = geneList,
                      ont          = "ALL",  # "BP"、"MF"和"CC"或"ALL"
                      OrgDb        = OrgDb,#人类org.Hs.eg.db 鼠org.Mm.eg.db
                      keyType      = "ENTREZID",
                      pvalueCutoff = 1)   #实际为padj阈值可调整
GO_kk <- DOSE::setReadable(GO_kk_entrez, 
                           OrgDb=OrgDb,
                           keyType='ENTREZID')#转化id 

save(GO_kk_entrez, file = "GSEA_go.RData")

class(GO_kk)
a <- GO_kk@result

gseaplot2(GO_kk,
          a$ID[i],#富集的ID编号
          title = a$Description[i],#标题
          color = "red", #GSEA线条颜色
          base_size = 20,#基础字体大小
          rel_heights = c(1.5, 0.5, 1),#副图的相对高度
          subplots = 1:3,   #要显示哪些副图 如subplots=c(1,3) #只要第一和第三个图
          ES_geom = "line", #enrichment score用线还是用点"dot"
          pvalue_table = T) #显示pvalue等信息

#---------------------------------------------------
# kegg
R.utils::setOption("clusterProfiler.download.method",'auto')
KEGG_kk_entrez <- gseKEGG(geneList     = geneList,
                          organism     = organism, #人hsa 鼠mmu
                          pvalueCutoff = 1)  #实际为padj阈值可调整 
KEGG_kk <- DOSE::setReadable(KEGG_kk_entrez, 
                             OrgDb=OrgDb,
                             keyType='ENTREZID')#转化id  
save(KEGG_kk, file = "GSEA_kegg.RData")

class(KEGG_kk)
a <- KEGG_kk@result
rownames(a) <- 1:125

View(KEGG_kk@result)
i=38
gseap1 <- gseaplot2(KEGG_kk,
                    a$ID[i],#富集的ID编号
                    title = a$Description[i],#标题
                    color = "red", #GSEA线条颜色
                    base_size = 20,#基础字体大小
                    rel_heights = c(1.5, 0.5, 1),#副图的相对高度
                    subplots = 1:3,   #要显示哪些副图 如subplots=c(1,3) #只要第一和第三个图
                    ES_geom = "line", #enrichment score用线还是用点"dot"
                    pvalue_table = T) #显示pvalue等信息
gseap1
