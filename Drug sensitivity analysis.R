#Firstly, we need to read the expression matrix and drug treatment information of the training set.
options(stringsAsFactors = F)
library(oncoPredict)
library(data.table)
library(gtools)
library(reshape2)
library(ggpubr)
th=theme(axis.text.x = element_text(angle = 45,vjust = 0.5))
dir="D:/DOCUMENTS/BLCA/blca/kidney_ca/DataFiles/Training Data"
#Import training set
GDSC2_Expr = readRDS(file=file.path(dir,'GDSC2_Expr (RMA Normalized and Log Transformed).rds'))
GDSC2_Res = readRDS(file = file.path(dir,"GDSC2_Res.rds"))
GDSC2_Res <- exp(GDSC2_Res)
write.csv(GDSC2_Res,file ="GDSC2_Res.csv"  )
#Filter KIRC data
library(readxl)
Cell_Lines_Details <- read_excel("DataFiles/GLDS/GDSCv2/Cell_Lines_Details.xlsx")
Cell_Lines_Details$`COSMIC identifier`<-paste0("COSMIC_",Cell_Lines_Details$`COSMIC identifier`)
colnames(Cell_Lines_Details)[10]<-"TCGAtype"
Cell_Lines_Details<-as.data.frame(Cell_Lines_Details)
Cell_Lines_Details<-Cell_Lines_Details[grep("KIRC",Cell_Lines_Details$TCGAtype),]
GDSC2_Res<-GDSC2_Res[rownames(GDSC2_Res)%in%Cell_Lines_Details$`COSMIC identifier`,] 
GDSC2_Expr<-GDSC2_Expr[,colnames(GDSC2_Expr)%in%Cell_Lines_Details$`COSMIC identifier`]

#Import verification set

testExpr<-as.matrix(testExpr)

library(GSVA)
library(GSEABase)
fpkm<-as.matrix(GDSC2_Expr)
str(fpkm)

cpgmt <- getGmt("D:/DOCUMENTS/c2.cp.v7.4.symbols.gmt")


gs.exp1 <- gsva(fpkm,cpgmt,method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE, min.sz = 10)
normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}#定义ssGSEA_Score矫正函数
norm_ssGSEA_Score1=normalize(gs.exp1)#对ssGSEA_Score进行矫正
norm_ssGSEA_Score1=rbind(id=colnames(norm_ssGSEA_Score1),norm_ssGSEA_Score1)


ssGSEA_GDSC_KIRC<-norm_ssGSEA_Score1[-1,]
pathway<-c("BIOCARTA_IL3_PATHWAY")
ssGSEA_GDSC_KIRC<-ssGSEA_GDSC_KIRC[pathway,]
data<-as.data.frame(t(ssGSEA_GDSC_KIRC))
data$id<-rownames(data)
GDSC2_Res<-as.data.frame(GDSC2_Res)
GDSC2_Res$id<-rownames(GDSC2_Res)
data<-merge(data,GDSC2_Res)
rownames(data)<-data$id
data$BIOCARTA_IL3_PATHWAY<-ifelse(data$BIOCARTA_IL3_PATHWAY>median(data$BIOCARTA_IL3_PATHWAY),"high","low")


#Interpretation of drug prediction results
library(ggplot2)
library(ggpubr)

if(T){
  mytheme <- theme(plot.title = element_text(size = 12,color="black",hjust = 0.5),
                   axis.title = element_text(size = 12,color ="black"), 
                   axis.text = element_text(size= 12,color = "black"),
                   panel.grid.minor.y = element_blank(),
                   panel.grid.minor.x = element_blank(),
                   axis.text.x = element_text(angle = 45, hjust = 1 ),
                   panel.grid=element_blank(),
                   legend.position = "top",
                   legend.text = element_text(size= 12),
                   legend.title= element_text(size= 12)
  ) }
#Check the drug IC50 distribution of the sample and make statistical analysis.
for (i in 4:ncol(data)) {
  input = colnames(data)[i]
  box_TME<-ggplot(data, aes(x = BIOCARTA_IL3_PATHWAY, y = data[,input]))+ 
    labs(y=paste0(input,":","IC50(uM)"),x= NULL,title = NULL)+  
    geom_boxplot(aes(fill = BIOCARTA_IL3_PATHWAY),position=position_dodge(0.5),width=0.5,outlier.alpha = 0)+ 
    scale_fill_manual(values = c("#1CB4B8", "#EB7369"))+
    theme_classic() + mytheme + 
    stat_compare_means(aes(group = BIOCARTA_IL3_PATHWAY),
                       label = "p.signif",
                       method = "wilcox.test",
                       hide.ns = F)
  ggsave(paste0("output/drug/",input,"_IC50.png"),box_TME,height=15,width=20,unit="cm")
}

