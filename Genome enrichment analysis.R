library(data.table)
gse167573<-fread("GSE167573.gene_expression_RNAseq.Normalized_Counts.tsv.gz",
                 header = T)
gse167573_cli<-fread("GSE167573.clinical.tsv.gz",header = T)
gse167573<-as.data.frame(gse167573)
rownames(gse167573)<-gse167573$Gene_Symbol

#SSGSEA
library(GSVA)
library(GSEABase)

fpkm<-gse167573[,-1]
# Read the channel data downloaded from GSEA official website.
cpgmt <- getGmt("D:/DOCUMENTS/c2.cp.v7.4.symbols.gmt")

fpkm<-as.matrix(fpkm)
# gsva 
gs.exp1<-gsva(fpkm,cpgmt,method="ssgsea",kcdf='Gaussian',abs.ranking=TRUE, min.sz = 10)
normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}#Define ssGSEA_Score correction function.
norm_ssGSEA_Score1=normalize(gs.exp1)#Correct ssGSEA_Score.
norm_ssGSEA_Score1=rbind(id=colnames(norm_ssGSEA_Score1),norm_ssGSEA_Score1)

path<-"BIOCARTA_IL3_PATHWAY"
norm_ssGSEA_Score1<-as.data.frame(norm_ssGSEA_Score1)
ssgsea_167573<-norm_ssGSEA_Score1[path,]

gse1<-as.data.frame(ssgsea_167573)
gs1<-as.data.frame(t(gse1))
gs1$Patient<-rownames(gs1)
new_dat1<-as.data.frame(merge(gs1,gse167573_cli))
rownames(new_dat1)<-new_dat1$Patient
str(new_dat1)
new_dat1$BIOCARTA_IL3_PATHWAY<-as.numeric(new_dat1$BIOCARTA_IL3_PATHWAY)
new.cut1<-surv_cutpoint(new_dat1,time = "OS_Time",event = "OS_Status",
                        variables = colnames(new_dat1)[2])
summary(new.cut1)
newF1<-surv_categorize(new.cut1)

gse167573_cli<-cbind(newF1,new_dat1)
colnames(gse167573_cli)[3]<-"il3"

#Extract high and low expression groups respectively
gse167573<-as.data.frame(gse167573[,-1])
IL3_H1<-gse167573_cli[grep("high",gse167573_cli$il3),]
IL3_L1<-gse167573_cli[grep("low",gse167573_cli$il3),]
count_il3h1<-gse167573[,colnames(gse167573)%in%rownames(IL3_H1)]
count_il3l1<-gse167573[,colnames(gse167573)%in%rownames(IL3_L1)]
maa1<-as.data.frame(cbind(count_il3h1,count_il3l1))
#Difference analysis of high expression group of pathway

####IL3####
library(edgeR)

group_list1<-c(IL3_H1$il3,
               IL3_L1$il3)
group_list1<-as.factor(group_list1)

#Data preprocessing

mycounts <- maa1
mycounts=as.matrix(mycounts)  ##The DGEList object requires a matrix format.
exp=mycounts
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)  ##If there are genes with duplicate names, the average value will be automatically taken to facilitate the later analysis.
data=data[rowMeans(data)>1,]  ##This gene is basically not expressed or low expressed in all samples, so it is deleted.

design <- model.matrix(~group_list1)  ##Become an experimental design object
y <- DGEList(counts=data,group=group_list1)  ##To become a data structure that edgr packets can recognize, data needs to be in matrix format.
#Differential expression analysis
y <- calcNormFactors(y)  ##To correct the factors, TMM(trimmed mean of M-values) is used for normalization by default.
y <- estimateCommonDisp(y)  ##Calculate the same dispersion for all genes.
y <- estimateTrendedDisp(y) ##According to the relationship between the mean and variance of different genes, the dispersion is calculated and a trended model is fitted.
y <- estimateTagwiseDisp(y) ##The dispersion of each gene was calculated and compressed to Trend Didspersion by empirical bayes method.
##Personally, I think this item is equivalent to the beta value of each gene in GLM.
plotBCV(y)
et <- exactTest(y,pair = c("high","low")) #Test the coefficient of variation of two groups of time.
## ##Similar paired t test
ordered_tags <- topTags(et, n=100000)##The correction of P value adopts BH method.
#Screening of differentially expressed genes

allDiff=ordered_tags$table  
allDiff=allDiff[is.na(allDiff$FDR)==FALSE,]
diff1=allDiff 
newData=y$pseudo.counts  
###Define two filter conditions
foldChange=1  #The difference is more than 2 times
padj=0.05  ##The corrected p value is less than 0.05.
diff1$change<-ifelse(diff1$PValue<0.05 & abs(diff1$logFC)>=1,
                     ifelse(diff1$logFC>=1,"UP","DOWN"),
                     "NOT")
table(diff1$change)

####GSEA##

#Before analysis, you still need to convert SYMBOL into ENTREZ ID.
diffup1<-diff1[grep("UP",diff1$change),]
diffdown1<-diff1[grep("DOWN",diff1$change),]
diffup1$SYMBOL<-rownames(diffup1)
il3data1<-diffup1[,c(1,6)]
df.id<-bitr(il3data1$SYMBOL, 
            fromType = "SYMBOL",
            toType = "ENTREZID",
            OrgDb = "org.Hs.eg.db")
easy.df1<-merge(il3data1,df.id,by="SYMBOL",all=F)
head(easy.df1)
#In order to do GSEA analysis, all the genes should be sorted first, which is the butterfly diagram that everyone often sees. Here, logFC is used as the quantitative standard to sort the genes from big to small.
sortdf<-easy.df1[order(easy.df1$logFC, decreasing = T),]
head(sortdf)
gene.expr = sortdf$logFC
head(gene.expr)
names(gene.expr) <- sortdf$ENTREZID
head(gene.expr)
kk <- gseKEGG(gene.expr, nPerm = 1000, minGSSize = 10, maxGSSize = 1000, pvalueCutoff=0.05)
head(kk)
dim(kk)
#Sort by enrichment score from high to low, which is convenient for viewing enrichment pathways.
sortkk<-kk[order(kk$enrichmentScore, decreasing = T),]
head(sortkk)
dim(sortkk)
Go <- gseGO(gene.expr, 'org.Hs.eg.db', keyType = "ENTREZID", ont="BP", nPerm = 1000, minGSSize = 10, maxGSSize = 1000, pvalueCutoff=0.05)
head(Go)
sortGo<-Go[order(Go$enrichmentScore, decreasing = T),]
Go_Reactomeresult <- gsePathway(gene.expr, nPerm = 1000, minGSSize = 10, maxGSSize = 1000, pvalueCutoff=0.05)
dim(Go_Reactomeresult )
sortReactome<-Go_Reactomeresult[order(Go_Reactomeresult$enrichmentScore, decreasing = T),]
gsea_result1<-rbind(sortkk,sortReactome)
sortgsea_result1<-gsea_resultn1[order(gsea_result1$enrichmentScore, decreasing = T),]
diffdown1$SYMBOL<-rownames(diffdown1)
il3data2<-diffdown1[,c(1,6)]
df.id<-bitr(il3data2$SYMBOL, 
            fromType = "SYMBOL",
            toType = "ENTREZID",
            OrgDb = "org.Hs.eg.db")
easy.df2<-merge(il3data2,df.id,by="SYMBOL",all=F)
head(easy.df2)

sortdf2<-easy.df2[order(easy.df2$logFC, decreasing = T),]
head(sortdf2)
gene.expr2 = sortdf2$logFC
names(gene.expr2) <- sortdf2$ENTREZID
head(gene.expr2)
kk2 <- gseKEGG(gene.expr2, nPerm = 1000, minGSSize = 10, maxGSSize = 1000, pvalueCutoff=0.05)
head(kk2)
dim(kk2)

sortkk2<-kk2[order(kk2$enrichmentScore, decreasing = T),]
head(sortkk2)
dim(sortkk2)
Go2 <- gseGO(gene.expr2, 'org.Hs.eg.db', keyType = "ENTREZID", ont="BP", nPerm = 1000, minGSSize = 10, maxGSSize = 1000, pvalueCutoff=0.05)
head(Go2)
dim(Go2)
sortGo2<-Go2[order(Go2$enrichmentScore, decreasing = T),]
Go_Reactomeresult2 <- gsePathway(gene.expr2, nPerm = 1000, minGSSize = 10, maxGSSize = 1000, pvalueCutoff=0.05)
dim(Go_Reactomeresult2 )
sortReactome2<-Go_Reactomeresult2[order(Go_Reactomeresult2$enrichmentScore, decreasing = T),]

####H  sortgsea_result1,L sortGo2
result1<-rbind(sortgsea_result1,sortGo2)

####TCGA###
library(readr)
TCGAbiolinks_KIRC_counts <- read_csv("TCGAbiolinks-KIRC-counts.csv")
View(TCGAbiolinks_KIRC_counts)
TCGAbiolinks_KIRC_counts<-as.data.frame(TCGAbiolinks_KIRC_counts)

rownames(TCGAbiolinks_KIRC_counts)<-TCGAbiolinks_KIRC_counts[,1]
TCGAbiolinks_KIRC_counts<-TCGAbiolinks_KIRC_counts[,-1]

TCGAbiolinks_KIRC_counts<-TCGAbiolinks_KIRC_counts[,grep("01A",colnames(TCGAbiolinks_KIRC_counts))]#530个样本
colnames(TCGAbiolinks_KIRC_counts)<-substr(colnames(TCGAbiolinks_KIRC_counts),1,12)

library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
TCGAbiolinks_KIRC_counts$geneid<-rownames(TCGAbiolinks_KIRC_counts)
map_dt <- bitr(rownames(TCGAbiolinks_KIRC_counts), fromType = "ENSEMBL",toType = c( "SYMBOL"),OrgDb = org.Hs.eg.db)
dt_merge <- merge(map_dt,TCGAbiolinks_KIRC_counts, by.y = "geneid", by.x = "ENSEMBL")
dt_merge<-dt_merge[!duplicated(dt_merge$SYMBOL),]
rownames(dt_merge) <- dt_merge$SYMBOL
dt_merge<-dt_merge[,-c(1,2)]


IL3_H<-new4[grep("high",new4$BIOCARTA_IL3_PATHWAY),]
IL3_L<-new4[grep("low",new4$BIOCARTA_IL3_PATHWAY),]
count_il3h<-dt_merge[,colnames(dt_merge)%in%rownames(IL3_H)]
count_il3l<-dt_merge[,colnames(dt_merge)%in%rownames(IL3_L)]
maa<-as.data.frame(cbind(count_il3h,count_il3l))
####IL3####
library(edgeR)

group_list<-c(IL3_H$BIOCARTA_IL3_PATHWAY,
              IL3_L$BIOCARTA_IL3_PATHWAY)
group_list<-as.factor(group_list)

dim(maa)
mycounts <- maa
mycounts=as.matrix(mycounts)  
exp=mycounts
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)  
data=data[rowMeans(data)>1,]  

design <- model.matrix(~group_list)  
y <- DGEList(counts=data,group=group_list)  

y <- calcNormFactors(y)  

y <- estimateCommonDisp(y)  
y <- estimateTrendedDisp(y) 
y <- estimateTagwiseDisp(y) 
plotBCV(y)
et <- exactTest(y,pair = c("high","low"))

ordered_tags <- topTags(et, n=100000)


allDiff=ordered_tags$table  
allDiff=allDiff[is.na(allDiff$FDR)==FALSE,]
diff=allDiff  
newData=y$pseudo.counts  

foldChange=1  
padj=0.05  
diff$change<-ifelse(diff$PValue<0.05 & abs(diff$logFC)>=1,
                    ifelse(diff$logFC>=1,"UP","DOWN"),
                    "NOT")
table(diff$change)

####GSEA##
diffup<-diff[grep("UP",diff$change),]
diffdown<-diff[grep("DOWN",diff$change),]
diffup$SYMBOL<-rownames(diffup)
il3data1<-diffup[,c(1,6)]
df.id<-bitr(il3data1$SYMBOL, 
            fromType = "SYMBOL",
            toType = "ENTREZID",
            OrgDb = "org.Hs.eg.db")
easy.df1<-merge(il3data1,df.id,by="SYMBOL",all=F)
head(easy.df1)
sortdf<-easy.df1[order(easy.df1$logFC, decreasing = T),]
head(sortdf)
gene.expr = sortdf$logFC
head(gene.expr)
names(gene.expr) <- sortdf$ENTREZID
kk <- gseKEGG(gene.expr, nPerm = 1000, minGSSize = 10, maxGSSize = 1000, pvalueCutoff=0.05)
head(kk)
dim(kk)
sortkk<-kk[order(kk$enrichmentScore, decreasing = T),]
head(sortkk)
dim(sortkk)
Go <- gseGO(gene.expr, 'org.Hs.eg.db', keyType = "ENTREZID", ont="BP", nPerm = 1000, minGSSize = 10, maxGSSize = 1000, pvalueCutoff=0.05)
head(Go)
dim(Go)
sortGo<-Go[order(Go$enrichmentScore, decreasing = T),]
Go_Reactomeresult <- gsePathway(gene.expr, nPerm = 1000, minGSSize = 10, maxGSSize = 1000, pvalueCutoff=0.05)
dim(Go_Reactomeresult )
sortReactome<-Go_Reactomeresult[order(Go_Reactomeresult$enrichmentScore, decreasing = T),]
gsea_result1<-rbind(sortkk,sortReactome)
sortgsea_result1<-gsea_result1[order(gsea_result1$enrichmentScore, decreasing = T),]
diffdown$SYMBOL<-rownames(diffdown)
il3data2<-diffdown[,c(1,6)]
df.id<-bitr(il3data2$SYMBOL, 
            fromType = "SYMBOL",
            toType = "ENTREZID",
            OrgDb = "org.Hs.eg.db")
easy.df2<-merge(il3data2,df.id,by="SYMBOL",all=F)
head(easy.df2)

sortdf2<-easy.df2[order(easy.df2$logFC, decreasing = T),]
head(sortdf2)
gene.expr2 = sortdf2$logFC
head(gene.expr2)
names(gene.expr2) <- sortdf2$ENTREZID
head(gene.expr2)
kk2 <- gseKEGG(gene.expr2, nPerm = 1000, minGSSize = 10, maxGSSize = 1000, pvalueCutoff=0.05)
head(kk2)
dim(kk2)

sortkk2<-kk2[order(kk2$enrichmentScore, decreasing = T),]
head(sortkk2)
dim(sortkk2)
Go21 <- gseGO(gene.expr2, 'org.Hs.eg.db', keyType = "ENTREZID", ont="BP", nPerm = 1000, minGSSize = 10, maxGSSize = 1000, pvalueCutoff=0.05)
head(Go21)
dim(Go21)
sortGo2<-Go2[order(Go2$enrichmentScore, decreasing = T),]
Go_Reactomeresult2 <- gsePathway(gene.expr2, nPerm = 1000, minGSSize = 10, maxGSSize = 1000, pvalueCutoff=0.05)
dim(Go_Reactomeresult2 )
sortReactome2<-Go_Reactomeresult2[order(Go_Reactomeresult2$enrichmentScore, decreasing = T),]
gsea_result<-sortGo2
sortgsea_result<-gsea_result[order(gsea_result$enrichmentScore, decreasing = T),]

result<-rbind(sortgsea_result1,sortgsea_result)


result$es_type <- ifelse(result$enrichmentScore< 0, "below", "above")  
result$Description<- factor(result$Description, levels =result$Description)  



id<-intersect(result$ID,result1$ID)
result1<-result1[result1$ID%in%id,]
result<-result[result$ID%in%id,]

gsea1<-result[abs(result$enrichmentScore)>=0.5,]
gsea2<-result1[abs(result1$enrichmentScore)>=0.5,]
id1<-intersect(gsea1$ID,gsea2$ID)

gsea1$es_type <- ifelse(gsea1$enrichmentScore< 0, "below", "above")  
gsea1$Description<- factor(gsea1$Description, levels =gsea1$Description)  
gsea2$es_type <- ifelse(gsea2$enrichmentScore< 0, "below", "above")  
gsea2$Description<- factor(gsea2$Description, levels =gsea2$Description)  

ggplot(gsea1[id1,], aes(x=Description, y=enrichmentScore, label=enrichmentScore)) +  
  geom_bar(stat='identity', aes(fill=es_type), width=.5)  +  
  scale_fill_manual(name="Group", 
                    labels = c("enrichmentScore>0", "enrichmentScore<0"), 
                    values = c("above"="#00ba38", "below"="#f8766d")) +   
  labs(subtitle="Enrichment Score of GSEA,IL3-H versus IL3-L", 
       
       title= "TCGA-KIRC") +  
  coord_flip() + 
  theme_bw()

ggplot(gsea2[id1,], aes(x=Description, y=enrichmentScore, label=enrichmentScore)) +  
  geom_bar(stat='identity', aes(fill=es_type), width=.5)  +  
  scale_fill_manual(name="Group", 
                    labels = c("enrichmentScore>0", "enrichmentScore<0"), 
                    values = c("above"="#00ba38", "below"="#f8766d")) +   
  labs(subtitle="Enrichment Score of GSEA,IL3-H versus IL3-L", 
       
       title= "GSE167573(KIRC)") +  
  coord_flip() + 
  theme_bw()


labels = c("Above Average", "Below Average"), 
values = c("above"="#00ba38", "below"="#f8766d")) +   labs(subtitle="Normalised mileage from 'mtcars'", 
                                                           title= "Diverging Bars") +   coord_flip() + theme_bw()



