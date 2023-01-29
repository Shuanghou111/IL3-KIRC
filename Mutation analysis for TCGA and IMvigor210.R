library(maftools)
library(ComplexHeatmap)
#Extract tumor samples:01A
library(readr)
TCGAbiolinks_KIRC_mutect2 <- read_csv("TCGAbiolinks-KIRC-mutect2.csv")
View(TCGAbiolinks_KIRC_mutect2)
TCGAbiolinks_KIRC_mutect2<-as.data.frame(TCGAbiolinks_KIRC_mutect2)
TCGAbiolinks_KIRC_mutect2<-TCGAbiolinks_KIRC_mutect2[grep("01A",TCGAbiolinks_KIRC_mutect2$Tumor_Sample_Barcode),]
TCGAbiolinks_KIRC_mutect2[,17]<-substr(TCGAbiolinks_KIRC_mutect2$Tumor_Sample_Barcode,1,12)
TCGA_mutation<-as.data.frame(TCGAbiolinks_KIRC_mutect2)
#Loading clinical data
TCGA_KIRC_Clinical<-read_csv("TCGA_Clinical.csv")
TCGA_cli<-as.data.frame(TCGA_KIRC_Clinical[,-1])
#Loading TCGA expression data
new4<-read_csv("TCGA_miRNA_data.csv")
#According to IL3 path level grouping
up1<-new4[grep("high",new4$BIOCARTA_IL3_PATHWAY),]
down1<-new4[grep("low",new4$BIOCARTA_IL3_PATHWAY),]
il3up<-TCGA_mutation[TCGA_mutation$Tumor_Sample_Barcode%in%rownames(up1),]
il3down<-TCGA_mutation[TCGA_mutation$Tumor_Sample_Barcode%in%rownames(down1),]
#Genetic mutation analysis using maftools package
library(maftools)
lam_h<-read.maf(maf = il3up,clinicalData =TCGA_cli)
lam_l<-read.maf(maf = il3down,clinicalData =TCGA_cli)
plotmafSummary(maf = lam_h,
               rmOutlier = TRUE, addStat = 'median', 
               dashboard = TRUE, titvRaw = FALSE)
oncoplot(maf=lam_h,top=20,clinicalFeatures
         =c("sex"))
plotmafSummary(maf = lam_l,
               rmOutlier = TRUE, addStat = 'median', 
               dashboard = TRUE, titvRaw = FALSE)
oncoplot(maf=lam_l,top=20,clinicalFeatures
         =c("sex"))
#The maf of two kinds of samples is compared by using the function mafCompare, and the output is the result obtained by fisher test of two kinds of samples.



lam<-read.maf(maf = TCGA_mutation,clinicalData =TCGA_cli)
df<-as.data.frame(lam@data)
gene1<-lam@gene.summary$Hugo_Symbol[1:20]#Top20 statistically significant genes
#IL-3 and TMB Gene Mutation as fisher.test
cg=as.character(lam@clinical.data$Tumor_Sample_Barcode)[lam@clinical.data$il3=='high']
cg#n=294
s1.maf <- read.maf(df[df$Tumor_Sample_Barcode %in% cg, ])

cg=as.character(lam@clinical.data$Tumor_Sample_Barcode)[lam@clinical.data$il3=='low']
cg#n=35
s3.maf <- read.maf(df[df$Tumor_Sample_Barcode %in% cg, ])
#fisher.text
pt.vs.rt <- mafCompare(m1 = s1.maf, m2 = s3.maf, m1Name = 'high',
                       m2Name = 'low', minMut = 1)
View(pt.vs.rt[["results"]])

forestPlot(mafCompareRes = pt.vs.rt, pVal = 0.01)


##construction
#TCGA
result<-pt.vs.rt[["results"]]
result<-result[result$pval<=0.05,]
il3up<-TCGA_mutation[TCGA_mutation$Tumor_Sample_Barcode%in%rownames(up1),]
il3down<-TCGA_mutation[TCGA_mutation$Tumor_Sample_Barcode%in%rownames(down1),]
lam_h<-read.maf(maf = il3up,clinicalData =TCGA_cli)
lam_l<-read.maf(maf = il3down,clinicalData =TCGA_cli)


il3clin_H<-TCGA_cli[grep("high",TCGA_cli$il3),]
il3clin_L<-TCGA_cli[grep("low",TCGA_cli$il3),]


#Obtaining the mutation information matrix of the top 20 mutant genes
oncoplot(maf = lam_h,
         genes = gene1,
         #genes = genelist，
         writeMatrix =T)
mut_matrix <- read.table("onco_matrix.txt",header = T,sep = "\t",as.is = T,stringsAsFactors = F)
colnames(mut_matrix)<-gsub('[.]', '-', colnames(mut_matrix))
mut_matrix[is.na(mut_matrix )] <-""
result<-pt.vs.rt[["results"]]
result<-result[result$pval<=0.05,]
library(gtools)
result$star=stars.pval(result$pval)

nam<-intersect(result$Hugo_Symbol,gene1)


rownames(mut_matrix)[1] <- paste0(rownames(mut_matrix)[1], "**" )
rownames(mut_matrix)[9] <- paste0(rownames(mut_matrix)[9], "*" )
mut_matrix[mut_matrix == "Translation_Start_Site"] <- ""
mut_matrix[mut_matrix == "Nonstop_Mutation"] <- ""
mut_matrix[mut_matrix == "In_Frame_Ins"] <- ""
mut_matrix[mut_matrix == "Multi_Hit"] <- ""
oncoplot(maf = lam_l,
         genes = gene1,
         #genes = genelist，
         writeMatrix =T)
mut_matrix1 <- read.table("onco_matrix.txt",header = T,sep = "\t",as.is = T,stringsAsFactors = F)
colnames(mut_matrix1)<-gsub('[.]', '-', colnames(mut_matrix1))
mut_matrix1[is.na(mut_matrix1 )] <-""
rownames(mut_matrix1)[2] <- paste0(rownames(mut_matrix1)[2], "**" )
rownames(mut_matrix1)[3] <- paste0(rownames(mut_matrix1)[3], "*" )
mut_matrix1[mut_matrix1 == "Translation_Start_Site"] <- ""
mut_matrix1[mut_matrix1 == "Nonstop_Mutation"] <- ""
mut_matrix1[mut_matrix1 == "In_Frame_Ins"] <- ""
mut_matrix1[mut_matrix1 == "Multi_Hit"] <- ""

column_title <- "OncoPrint for TCGA"
heatmap_legend_param <- list(title = "Alternations", at = c("Nonsense_Mutation" , "Missense_Mutation", 
                                                            "Splice_Site" ,"Frame_Shift_Del",
                                                            "Frame_Shift_Ins",
                                                            "In_Frame_Del"), 
                             labels = c("Nonsense_Mutation" , "Missense_Mutation", 
                                        "Splice_Site" ,"Frame_Shift_Del",
                                        "Frame_Shift_Ins",
                                        "In_Frame_Del"))

alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = "#CCCCCC", col = NA))},   
  Nonsense_Mutation = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["Nonsense_Mutation"], col = NA))}
  ,
  Missense_Mutation = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["Missense_Mutation"], col = NA))
  },
  Splice_Site = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, 
              gp = gpar(fill = col["Splice_Site"], col = NA))
  },
  Frame_Shift_Del=function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, 
              gp = gpar(fill = col["Frame_Shift_Del"], col = NA))
  },
  Frame_Shift_Ins=function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, 
              gp = gpar(fill = col["Frame_Shift_Ins"], col = NA))
  },
  In_Frame_Del=function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, 
              gp = gpar(fill = col["In_Frame_Del"], col = NA))
  })

col <-c("Nonsense_Mutation" = "#ff7f00", "Missense_Mutation" = "#984ea3", "Splice_Site"  = "#4daf4a",
        "Frame_Shift_Del"="#E64B35B2","Frame_Shift_Ins"="#4DBBD5B2" ,
        "In_Frame_Del"="#00A087B2")

il3clin_H<-TCGA_cli[grep("high",TCGA_cli$il3),]
il3clin_L<-TCGA_cli[grep("low",TCGA_cli$il3),]
rows<-intersect(colnames(mut_matrix),rownames(il3clin_H))
il3clin_H<-il3clin_H[rows,]
rows1<-intersect(colnames(mut_matrix1),rownames(il3clin_L))
il3clin_L<-il3clin_L[rows1,]
ha =  HeatmapAnnotation (IL3_Status=il3clin_H$il3,
                         Sex=il3clin_H$sex,
                         OS=il3clin_H$os,
                         TNM=il3clin_H$TNM,
                         Race=il3clin_H$race,
                         col = list(IL3_Status = c('high' = '#9afdb4', 'low' = '#0f9a89'),
                                    Sex=c("female"="#1b8f02","male"="#dc8a01")),
                         annotation_legend_param=list(gp = gpar(fontsize =  1 )),  
                         annotation_name_side= "left" ,      
                         annotation_name_gp = gpar(fontsize =  10 ))  
ha1 =  HeatmapAnnotation (IL3_Status=il3clin_L$il3,
                          Sex=il3clin_L$sex,
                          OS=il3clin_L$os,
                          TNM=il3clin_L$TNM,
                          Race=il3clin_L$race,
                          col = list(IL3_Status = c('high' = '#9afdb4', 'low' = '#0f9a89'),
                                     Sex=c("female"="#1b8f02","male"="#dc8a01")),
                          annotation_legend_param=list(gp = gpar(fontsize =  1 )),  
                          show_annotation_name = TRUE,  
                          annotation_name_side= "left" ,      
                          annotation_name_gp = gpar(fontsize =  10 ))  




alter_fun = list(
  background = alter_graphic("rect", fill = "#CCCCCC"),
  Missense_Mutation = alter_graphic("rect", fill = col["Missense_Mutation"]),
  Nonsense_Mutation = alter_graphic("rect", fill = col["Nonsense_Mutation"]),
  Splice_Site = alter_graphic("rect",  fill = col["Splice_Site"]),
  Frame_Shift_Del = alter_graphic("rect",  fill = col["Frame_Shift_Del"]),
  Frame_Shift_Ins = alter_graphic("rect",  fill = col["Frame_Shift_Ins"]),
  In_Frame_Del = alter_graphic("rect",  fill = col["In_Frame_Del"])
)
column_title1 <- "TCGA KIRC IL3H"
column_title2 <- "TCGA KIRC IL3L"
p1<-oncoPrint(mut_matrix[gene1 ,],top_annotation = ha,
              alter_fun = alter_fun, col = col, row_order = order(gene1,decreasing = T),
              column_title = column_title1, heatmap_legend_param = heatmap_legend_param,
              right_annotation = NULL,show_row_names = F)+
  oncoPrint(mut_matrix1[gene1 ,],top_annotation = ha1,
            alter_fun = alter_fun, col = col,row_order = order(gene1,decreasing = T),
            column_title = column_title2, heatmap_legend_param = heatmap_legend_param,
            row_names_side = "left",right_annotation = NULL,
            pct_side = "right")


draw(p1,annotation_legend_side = "bottom")



#Analysis of mutation data in IMvigor210 queue
library(IMvigor210CoreBiologies)
library("biomaRt")
library("circlize")
library("corrplot")
library("dplyr")
library("DT")
library("ggplot2")
library("limma")
library("lsmeans")
library("reshape2")
library("spatstat")
library("survival")
library("plyr")
data(cds)
expr<-counts(cds)
fdata<-fData(cds)
pd<-pData(cds)
View(expr)
View(fdata)
View(pd)
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
pdd<-pd[grep("kidney",pd$Tissue),]#n=67
#Extract mutation data
data(fmone)
ls.str(assayData(fmone))
head(pData(fmone))
head(assayDataElement(fmone, "known_short"))
pfmone<-pData(fmone)
#IL3 path packet data
new2<-read_csv("IMvigor_IL3_path.csv")
new2h<-new2[new2$BIOCARTA_IL3_PATHWAY=="high",]
new2l<-new2[new2$BIOCARTA_IL3_PATHWAY=="low",]
pdd$id<-rownames(pdd)
new2$id<-rownames(new2)
pdd<-merge(pdd,new2[,c(4,5)])
rownames(pdd)<-pdd$id
pddh<-pdd[rownames(pdd)%in%rownames(new2h),]
pddl<-pdd[rownames(pdd)%in%rownames(new2l),]
imv_il3h<-pfmone[pfmone$ANONPT_ID%in%pddh$ANONPT_ID,]
dim(imv_il3h)
imv_il3l<-pfmone[pfmone$ANONPT_ID%in%pddl$ANONPT_ID,]
dim(imv_il3l)
#Extracting integrated mutation data
known_short<-assayDataElement(fmone, "known_short")
amplification<-assayDataElement(fmone, "amplification")
deletion<-assayDataElement(fmone, "deletion")
gain<-assayDataElement(fmone, "gain")
likely_short<-assayDataElement(fmone, "likely_short")
datam1<-known_short
for (i in 1:nrow(datam1)) {
  for (j in 1:ncol(datam1)) {
    ifelse(datam1[i,j]!="",datam1[i,j]<-"known_short",datam1[i,j]<-"")
  }
}
datam1[datam1!=""]
datam2<-amplification
for (i in 1:nrow(datam2)) {
  for (j in 1:ncol(datam2)) {
    ifelse(datam2[i,j]!="",datam2[i,j]<-"amplification",datam2[i,j]<-"")
  }
}
datam2[datam2!=""]
datam3<-deletion
for (i in 1:nrow(datam3)) {
  for (j in 1:ncol(datam3)) {
    ifelse(datam3[i,j]!="",datam3[i,j]<-"deletion",datam3[i,j]<-"")
  }
}
datam3[datam3!=""]
datam4<-gain
for (i in 1:nrow(datam4)) {
  for (j in 1:ncol(datam4)) {
    ifelse(datam4[i,j]!="",datam4[i,j]<-"gain",datam4[i,j]<-"")
  }
}
datam4[datam4!=""]
datam5<-likely_short
for (i in 1:nrow(datam5)) {
  for (j in 1:ncol(datam5)) {
    ifelse(datam5[i,j]!="",datam5[i,j]<-"likely_short",datam5[i,j]<-"")
  }
}
datam5[datam5!=""]
datam<-datam1
m<-which(datam2!="")
datam[m]<-"amplification"
m1<-which(datam3!="")
datam[m1]<-"deletion"
m2<-which(datam4!="")
datam[m2]<-"gain"
m3<-which(datam5!="")
datam[m3]<-"likely_short"
table(datam[datam!=""])

#The following is the analysis of mutation data.
il3up<-datam[,colnames(datam)%in%rownames(imv_il3h)]#392genes，39samples
il3down<-datam[,colnames(datam)%in%rownames(imv_il3l)]#392genes，19samples
pfmone$group<-""
pfmone$group[rownames(pfmone)%in%colnames(il3down)]<-"low"
pfmone$group[rownames(pfmone)%in%colnames(il3up)]<-"high"
pfmone<-pfmone[pfmone$group!="",]
table(pfmone$group)
datam<-datam[,rownames(pfmone)]
name<-colnames(datam)
tmp <- any_mutation(fmone)
tmp<-tmp[,name]
ml<-"FMOne mutation burden per MB"
tmp <- tmp[rowSums(tmp)>0,]#该基因至少在一个病人中发生了突变
pathFMIout<-lapply(rownames(tmp),  function(gene) {
  mutStatus<-tmp[gene,]
  ind <- !is.na(pfmone$group)
  t <- table(group = pfmone$group[ind],
             mutant = factor(mutStatus[ind] > 0,
                             levels=c("TRUE","FALSE")))
  f <- fisher.test( t )
  b <- setNames(c( gene, sum(mutStatus[ind] > 0), "group", 
                   f$estimate, f$p.value ),
                c("Gene","n mutant", "category","estimate","PVal"))
  
})
pathFMIout<- do.call(rbind, pathFMIout)

pathFMIout<- as.data.frame(pathFMIout,stringsAsFactors=FALSE)

pathFMIout$estimate <- signif(as.numeric(pathFMIout$estimate))

pathFMIout$PVal <- signif(as.numeric(pathFMIout$PVal))

## adjust p-values

pathFMIout$AdjP <- NA ##新增加一列矫正后的p值

pathFMIout$AdjP[pathFMIout$category == "group"]<-p.adjust(
  pathFMIout$PVal[pathFMIout$category == "group"],
  method="BH")

pathFMIout$AdjP[pathFMIout$category == "TMB"] <- p.adjust(
  pathFMIout$PVal[pathFMIout$category == "TMB"], method="BH")
library(gtools)
pathFMIout$star=stars.pval(pathFMIout$PVal)

pathFMIout[pathFMIout$PVal<0.05,]                  

gene<-c("TERT","TP53","ARID1A","FGFR3","PIK3CA","CDKN2A",
        "CDKN2B","MLL2","KDM6A","ERBB3","STAG2","FRS2","MDM2",
        "ERBB2","ZNF703","CREBBP","RB1","RICTOR","CCND1","CTNNB1")

mut_matrix_h<-datam[,rownames(imv_il3h)]
mut_matrix_l<-datam[,rownames(imv_il3l)]
mut_matrix_h[is.na(mut_matrix_h )] <-""
mut_matrix_h<-mut_matrix_h[gene,]
mut_matrix_l<-mut_matrix_l[gene,]
rownames(mut_matrix_h)[4] <- paste0(rownames(mut_matrix_h)[4], "*" )
rownames(mut_matrix_h)[6] <- paste0(rownames(mut_matrix_h)[6], "*" )
rownames(mut_matrix_h)[10] <- paste0(rownames(mut_matrix_h)[10], "*" )
rownames(mut_matrix_l)[4] <- paste0(rownames(mut_matrix_l)[4], "*" )
rownames(mut_matrix_l)[6] <- paste0(rownames(mut_matrix_l)[6], "*" )
rownames(mut_matrix_l)[10] <- paste0(rownames(mut_matrix_l)[10], "*" )
column_title <- "OncoPrint for iMVigor210"
heatmap_legend_param <- list(title = "Alternations", at = c("known_short" , "likely_short", 
                                                            "amplification" ,"deletion",
                                                            "gain"), 
                             labels = c("known_short" , "likely_short", 
                                        "amplification" ,"deletion",
                                        "gain"))

alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = "#CCCCCC", col = NA))},   
  known_short = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["known_short"], col = NA))}
  ,
  likely_short = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["likely_short"], col = NA))
  },
  amplification = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, 
              gp = gpar(fill = col["amplification"], col = NA))
  },
  deletion=function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, 
              gp = gpar(fill = col["deletion"], col = NA))
  },
  gain=function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, 
              gp = gpar(fill = col["gain"], col = NA))
  })

col <-c("known_short" = "#ff7f00", "likely_short" = "#984ea3", "amplification"  = "#4daf4a",
        "deletion"="#E64B35B2","gain"="#4DBBD5B2" )

il3clin_H<-pddh[pddh$ANONPT_ID%in%imv_il3h$ANONPT_ID,]
il3clin_L<-pddl[pddl$ANONPT_ID%in%imv_il3l$ANONPT_ID,]
col_fun = colorRamp2(c(0, 50), c("white", "#bb3302"))
ha =  HeatmapAnnotation (IL3_Status=il3clin_H$BIOCARTA_IL3_PATHWAY,
                         Sex=il3clin_H$Sex,
                         OS=il3clin_H$os,
                         TNM=il3clin_H$`TCGA Subtype`,
                         Race=il3clin_H$Race,
                         Response=il3clin_H$binaryResponse,
                         col = list(IL3_Status = c('high' = '#c6dbf0', 'low' = '#072f6a'),
                                    Sex=c("F"="#a3cad4","M"="#d2e3ed"),
                                    TNM=c("I"="#bdd2e7","II"="#8996c1","III"="#8a6ab0","IV"="#8e3c9f"),
                                    Race=c("ASIAN"="#330804","WHITE"="#c5bccd",
                                           "BLACK OR AFRICAN AMERICAN"="#2b7382",
                                           "OTHER"="#b4e2c2"),
                                    Response=c("CR/PR"="#63a4f1","SD/PD"="#ebcb31"),
                                    OS=col_fun),
                         annotation_legend_param=list(gp = gpar(fontsize =  1 )),  
                         annotation_name_side= "left" ,      
                         annotation_name_gp = gpar(fontsize =  10 ))  
ha1 =  HeatmapAnnotation (IL3_Status=il3clin_L$BIOCARTA_IL3_PATHWAY,
                          Sex=il3clin_L$Sex,
                          OS=il3clin_L$os,
                          TNM=il3clin_L$`TCGA Subtype`,
                          Race=il3clin_L$Race,
                          Response=il3clin_L$binaryResponse,
                          col = list(IL3_Status = c('#c6dbf0', 'low' = '#072f6a'),
                                     Sex=c("F"="#a3cad4","M"="#d2e3ed"),
                                     TNM=c("I"="#bdd2e7","II"="#8996c1","III"="#8a6ab0","IV"="#8e3c9f"),
                                     Race=c("ASIAN"="#330804","WHITE"="#c5bccd",
                                            "BLACK OR AFRICAN AMERICAN"="#2b7382",
                                            "OTHER"="#b4e2c2"),
                                     Response=c("CR/PR"="#63a4f1","SD/PD"="#ebcb31"),
                                     OS=col_fun),
                          annotation_legend_param=list(gp = gpar(fontsize =  1 )),  
                          annotation_name_side= "left" ,      
                          annotation_name_gp = gpar(fontsize =  10 ))  





column_title1 <- "IMV KIRC IL3H"
column_title2 <- "IMV KIRC IL3L"
p2<-oncoPrint(mut_matrix_h,top_annotation = ha,
              alter_fun = alter_fun, col = col,,
              column_title = column_title1, heatmap_legend_param = heatmap_legend_param,
              right_annotation = NULL,show_row_names = F)+
  oncoPrint(mut_matrix_l,top_annotation = ha1,
            alter_fun = alter_fun, col = col,
            column_title = column_title2, heatmap_legend_param = heatmap_legend_param,
            row_names_side = "left",right_annotation = NULL,
            pct_side = "right")





