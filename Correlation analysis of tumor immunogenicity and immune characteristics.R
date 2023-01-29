library(IOBR)
library(ggpubr)
#Introducing TCGA data,The file is too large to upload to github after compression.
library(readr)
TCGA_FPKM_1 <- read_csv("TCGA_FPKM_1.csv")
TCGA_FPKM_1 <-as.data.frame(TCGA_FPKM_1 )
rownames(TCGA_FPKM_1 )<-TCGA_FPKM_1 $...1
TCGA_FPKM_1 <-TCGA_FPKM_1 [,-1]
fpkm<-TCGA_FPKM_1
#Convert to log2(FPKM+1) data.
fpkm<-log2(TCGA_FPKM_1+1)

#IOBR collected 255 signatures related to tumor microenvironment and various tumor phenotypes published by CNS.
#Signature collection: signature_ is a list structure, the name represents the name of signature, and the content included is genes related to it; Signature_collection (including all the signatures we collected); We classify some small categories according to the contents to facilitate your specific research: signature_tme (related to tumor microenvironment); Signature_metabolism (metabolism related); Signature_tumor (tumor related);

###Analysis of sigantures related to tumor microenvironment
#In the function of calculate_sig_score, pdata is the phenotypic data of patients, so it is unnecessary to provide it. If it is provided, the parameter: column_of_sample should be further added to clarify the merged column; Eset is the expression matrix, and the name of the behavioral gene (symbol) is listed as the sample name; Sinature is the input signature list (as mentioned above, or you can build one you want to calculate yourself, and check the format_signature function); Method provides four available parameters: pca, ssgsea, zscore, integration (all three methods are calculated and then merged); Mini_gene_count indicates the number of genes intersecting with the target signature gene in the gene expression matrix. If there are few, the reliability of the calculated score is very low, so it needs to be limited. If the total number of genes is less than the minimum number of genes you set, the signature will be skipped.

sig_tme<-calculate_sig_score(pdata= NULL,
                             eset= fpkm,
                             signature= signature_tme,
                             method = "pca",
                             mini_gene_count = 2)



####cibersort
cibersort<-deconvo_tme(eset = fpkm, method = "cibersort", arrays = FALSE, perm = 200 )
View(cibersort)
TME_data<-cibersort
new4$BIOCARTA_IL3_PATHWAY<-as.factor(new4$BIOCARTA_IL3_PATHWAY)
TME_data$group <- new4$BIOCARTA_IL3_PATHWAY
colnames(TME_data)[1] <- "sample"
TME_data<-as.data.frame(TME_data)
# 2.2 Fusion data
variable<-c("T_cells_CD8_CIBERSORT","T_cells_gamma_delta_CIBERSORT",
            "Macrophages_M1_CIBERSORT","Dendritic_cells_resting_CIBERSORT",
            "Macrophages_M0_CIBERSORT")
TME_New = melt(TME_data)
TME_data1<-TME_data[,variable]
TME_data1$Sample<-TME_data$sample
TME_data1$group<-TME_data$group
TME_New1 = melt(TME_data1)
colnames(TME_New1)=c("Sample","group","Celltype","Composition")  #设置行名
# 3.3 plot

mytheme <- theme(plot.title = element_text(size = 12,color="black",hjust = 0.5),
                 axis.title = element_text(size = 12,color ="black"), 
                 axis.text = element_text(size= 12,color = "black"),
                 panel.grid.minor.y = element_blank(),
                 panel.grid.minor.x = element_blank(),
                 axis.text.x = element_text(angle = 45, hjust = 1 ),
                 panel.grid=element_blank(),
                 legend.position = "top",
                 legend.text = element_text(size= 12),
                 legend.title= element_text(size= 12))

TME_New1$Celltype<-as.character(TME_New1$Celltype)
for (i in 1:nrow(TME_New1)) {
  TME_New1$Celltype[i]<-strsplit(TME_New1$Celltype,"_CIBERSORT")[[i]][1] 
}

#"Wilcoxon rank sum test" or called "Mann Whitney U test".#IL3
ggplot(TME_New1, aes(x = Celltype, y = Composition))+ 
  labs(y="signature score",x= NULL,title = "TCGA-KIRC CIBERSORT")+  
  geom_boxplot(aes(fill = group),position=position_dodge(0.5),width=0.5,outlier.alpha = 0)+ 
  scale_fill_manual(values = c("#1CB4B8", "#EB7369"))+
  theme_classic() + mytheme + 
  stat_compare_means(aes(group = group),
                     label = "p.signif",
                     method = "wilcox.test",
                     hide.ns = F)


####Immune cell correlation analysis:xcell
xcell<-deconvo_tme(eset = fpkm, method = "xcell", arrays = FALSE, perm = 200 )
View(xcell)
xcell_data<-xcell
xcell_data$group <- new4$BIOCARTA_IL3_PATHWAY
colnames(xcell_data)[1] <- "sample"
# 2.2 Fusion data
for (i in 2:68) {
  colnames(xcell_data)[i]<-strsplit(colnames(xcell_data),"_xCell")[[i]][1]
}
xcell_New = melt(xcell_data)
variable2<-c("pro_B-cells","Macrophages_M1","Macrophages_M2",
             "CD8+_naive_T-cells","CD8+_T-cells","CD8+_Tcm","CD8+_Tem")
xcell_data1<-xcell_data[,variable2]
xcell_data1$sample<-xcell_data$sample
xcell_data1$group<-xcell_data$group
xcell_data1<-as.data.frame(xcell_data1)
xcell_New1 = melt(xcell_data1)
colnames(xcell_New1)=c("Sample","group","Celltype","Composition")  #设置行名
# 3.3 plot

mytheme <- theme(plot.title = element_text(size = 10,color="black",hjust = 0.5),
                 axis.title = element_text(size = 7,color ="black"), 
                 axis.text = element_text(size= 7,color = "black"),
                 panel.grid.minor.y = element_blank(),
                 panel.grid.minor.x = element_blank(),
                 axis.text.x = element_text(angle = 45, hjust = 1 ),
                 panel.grid=element_blank(),
                 legend.position = "top",
                 legend.text = element_text(size= 8),
                 legend.title= element_text(size= 8))

#Wilcoxon rank sum test或者叫Mann Whitney U test.#IL3
ggplot(xcell_New1, aes(x = Celltype, y = Composition))+ 
  labs(y="signature score",x= NULL,title = "TCGA-KIRC XCell")+  
  geom_boxplot(aes(fill = group),position=position_dodge(0.5),width=0.5,outlier.alpha = 0)+ 
  scale_fill_manual(values = c("#1CB4B8", "#EB7369"))+
  theme_classic() + mytheme + 
  stat_compare_means(aes(group = group),
                     label = "p.signif",
                     method = "wilcox.test",
                     hide.ns = F,size=1.5)


###MCP counter
mcp<-deconvo_tme(eset = fpkm, method = "mcpcounter", arrays = FALSE, perm = 200 )
View(mcp)
mcp_data<-mcp
mcp_data$group <- new4$BIOCARTA_IL3_PATHWAY
colnames(mcp_data)[1] <- "sample"
# 2.2 Fusion data
mcp_New = melt(mcp_data)
colnames(mcp_New)=c("Sample","group","Celltype","Composition")  #设置行名
for (i in 2:11) {
  colnames(mcp_data)[i]<-strsplit(colnames(mcp_data),"_MCPcounter")[[i]][1]
}

variable3<-c("B_lineage","Neutrophils",
             "T_cells","CD8_T_cells")
mcp_data1<-mcp_data[,variable3]
mcp_data1$sample<-mcp_data$sample
mcp_data1$group<-mcp_data$group
mcp_data1<-as.data.frame(mcp_data1)
mcp_New1 = melt(mcp_data1)
colnames(mcp_New1)=c("Sample","group","Celltype","Composition")  #设置行名
# 3.3 plot

mytheme <- theme(plot.title = element_text(size = 10,color="black",hjust = 0.5),
                 axis.title = element_text(size = 9,color ="black"), 
                 axis.text = element_text(size= 9,color = "black"),
                 panel.grid.minor.y = element_blank(),
                 panel.grid.minor.x = element_blank(),
                 axis.text.x = element_text(angle = 45, hjust = 1 ),
                 panel.grid=element_blank(),
                 legend.position = "top",
                 legend.text = element_text(size= 8),
                 legend.title= element_text(size= 8))

#Wilcoxon rank sum test或者叫Mann Whitney U test.#IL3
ggplot(mcp_New1, aes(x = Celltype, y = Composition))+ 
  labs(y="signature score",x= NULL,title = "TCGA-KIRC MCPcounter")+  
  geom_boxplot(aes(fill = group),position=position_dodge(0.5),width=0.5,outlier.alpha = 0)+ 
  scale_fill_manual(values = c("#1CB4B8", "#EB7369"))+
  theme_classic() + mytheme + 
  stat_compare_means(aes(group = group),
                     label = "p.signif",
                     method = "wilcox.test",
                     hide.ns = F)



###TIMER
timer<-deconvo_tme(eset =fpkm, method = "timer", group_list = rep("kirc",dim(fpkm)[2]))
View(timer)
timer_data<-timer
timer_data$group <- new4$BIOCARTA_IL3_PATHWAY
colnames(timer_data)[1] <- "sample"
# 2.2 Fusion data
timer_New = melt(timer_data)
colnames(timer_New)=c("Sample","group","Celltype","Composition")  #设置行名
# 3.3 plot
for (i in 2:7) {
  colnames(timer_data)[i]<-strsplit(colnames(timer_data),"_TIMER")[[i]][1]
}

variable4<-c("B_cell","Neutrophil",
             "Macrophage","DC","T_cell_CD4",
             "T_cell_CD8")
timer_data1<-timer_data[,variable4]
timer_data1$sample<-timer_data$sample
timer_data1$group<-timer_data$group
timer_data1<-as.data.frame(timer_data1)
timer_New1 = melt(timer_data1)
colnames(timer_New1)=c("Sample","group","Celltype","Composition")  #设置行名
mytheme <- theme(plot.title = element_text(size = 10,color="black",hjust = 0.5),
                 axis.title = element_text(size = 10,color ="black"), 
                 axis.text = element_text(size= 10,color = "black"),
                 panel.grid.minor.y = element_blank(),
                 panel.grid.minor.x = element_blank(),
                 axis.text.x = element_text(angle = 45, hjust = 1 ),
                 panel.grid=element_blank(),
                 legend.position = "top",
                 legend.text = element_text(size= 10),
                 legend.title= element_text(size= 10))

#Wilcoxon rank sum test或者叫Mann Whitney U test.#IL3
ggplot(timer_New1, aes(x = Celltype, y = Composition))+ 
  labs(y="signature score",x= NULL,title = "TCGA-KIRC TIMER")+  
  geom_boxplot(aes(fill = group),position=position_dodge(0.5),width=0.5,outlier.alpha = 0)+ 
  scale_fill_manual(values = c("#1CB4B8", "#EB7369"))+
  theme_classic() + mytheme + 
  stat_compare_means(aes(group = group),
                     label = "p.signif",
                     method = "wilcox.test",
                     hide.ns = F)


###quanTIseq
quanTIseq<-deconvo_tme(eset = fpkm, method = "quantiseq", tumor = TRUE, arrays = FALSE, scale_mrna = TRUE)
View(quanTIseq)
quanTIseq_data<-quanTIseq
quanTIseq_data$group <- new4$BIOCARTA_IL3_PATHWAY
colnames(quanTIseq_data)[1] <- "sample"
# 2.2 Fusion data
quanTIseq_New = melt(quanTIseq_data)
colnames(quanTIseq_New)=c("Sample","group","Celltype","Composition")  #设置行名
# 3.3 plot
for (i in 2:12) {
  colnames(quanTIseq_data)[i]<-strsplit(colnames(quanTIseq_data),"_quantiseq")[[i]][1]
}

variable5<-c("B_cells","Neutrophils",
             "Macrophages_M1","Macrophages_M2",
             "Tregs","T_cells_CD8")
quanTIseq_data1<-quanTIseq_data[,variable5]
quanTIseq_data1$sample<-quanTIseq_data$sample
quanTIseq_data1$group<-quanTIseq_data$group
quanTIseq_data1<-as.data.frame(quanTIseq_data1)
quanTIseq_New1 = melt(quanTIseq_data1)
colnames(quanTIseq_New1)=c("Sample","group","Celltype","Composition")  #设置行名
mytheme <- theme(plot.title = element_text(size = 10,color="black",hjust = 0.5),
                 axis.title = element_text(size = 10,color ="black"), 
                 axis.text = element_text(size= 10,color = "black"),
                 panel.grid.minor.y = element_blank(),
                 panel.grid.minor.x = element_blank(),
                 axis.text.x = element_text(angle = 45, hjust = 1 ),
                 panel.grid=element_blank(),
                 legend.position = "top",
                 legend.text = element_text(size= 10),
                 legend.title= element_text(size= 10))


ggplot(quanTIseq_New1, aes(x = Celltype, y = Composition))+ 
  labs(y="signature score",x= NULL,title = "TCGA-KIRC quanTIseq")+  
  geom_boxplot(aes(fill = group),position=position_dodge(0.5),width=0.5,outlier.alpha = 0)+ 
  scale_fill_manual(values = c("#1CB4B8", "#EB7369"))+
  theme_classic() + mytheme + 
  stat_compare_means(aes(group = group),
                     label = "p.signif",
                     method = "wilcox.test",
                     hide.ns = F)


###IPS
ips<-deconvo_tme(eset = fpkm, method = "ips", plot= FALSE)
View(ips)
ips_data<-ips
ips_data$group <- new4$BIOCARTA_IL3_PATHWAY
colnames(ips_data)[1] <- "sample"
# 2.2 Fusion data
ips_New = melt(ips_data)
colnames(ips_New)=c("Sample","group","Celltype","Composition")  #设置行名
# 3.3 plot

mytheme <- theme(plot.title = element_text(size = 10,color="black",hjust = 0.5),
                 axis.title = element_text(size = 10,color ="black"), 
                 axis.text = element_text(size= 10,color = "black"),
                 panel.grid.minor.y = element_blank(),
                 panel.grid.minor.x = element_blank(),
                 axis.text.x = element_text(angle = 45, hjust = 1 ),
                 panel.grid=element_blank(),
                 legend.position = "top",
                 legend.text = element_text(size= 10),
                 legend.title= element_text(size= 10))


ggplot(ips_New, aes(x = Celltype, y = Composition))+ 
  labs(y="signature score",x= NULL,title = "TCGA-KIRC IPS")+  
  geom_boxplot(aes(fill = group),position=position_dodge(0.5),width=0.5,outlier.alpha = 0)+ 
  scale_fill_manual(values = c("#1CB4B8", "#EB7369"))+
  theme_classic() + mytheme + 
  stat_compare_means(aes(group = group),
                     label = "p.signif",
                     method = "wilcox.test",
                     hide.ns = F)



##Immune related gene analysis



signature_score_calculation_methods
sig_tme<-calculate_sig_score(pdata= NULL,
                             eset= fpkm,
                             signature= signature_tme,
                             method= "ssgsea",
                             mini_gene_count = 2)
pdata_group<-new4[,c(3,4)]
pdata_group$ID<-rownames(pdata_group)
res_imucheck<-iobr_cor_plot(pdata_group= pdata_group,
                            id1= "ID",
                            feature_data= fpkm,
                            id2= "ID",
                            target= NULL,
                            group = "BIOCARTA_IL3_PATHWAY",
                            is_target_continuous  = FALSE,
                            padj_cutoff = 1,
                            index  = 1,
                            category= "gene",
                            signature_group= signature_collection[c(4)],
                            palette_box= "paired1",
                            palette_corplot= "pheatmap",
                            palette_heatmap= 4,
                            feature_limit= 26,
                            character_limit  = 30,
                            show_heatmap_col_name = FALSE,
                            show_col = FALSE,
                            show_plot = TRUE)
write.csv(res_imucheck,file="IL3-TCGA-immune_checkpoint.csv")

Chemokines<-iobr_cor_plot(pdata_group= pdata_group,
                          id1= "ID",
                          feature_data= fpkm,
                          id2= "ID",
                          target= NULL,
                          group = "BIOCARTA_IL3_PATHWAY",
                          is_target_continuous  = FALSE,
                          padj_cutoff = 1,
                          index  = 1,
                          category= "gene",
                          signature_group= signature_collection[["Chemokines_Li_et_al"]],
                          palette_box= "paired1",
                          palette_corplot= "pheatmap",
                          palette_heatmap= 4,
                          feature_limit= 26,
                          character_limit  = 30,
                          show_heatmap_col_name = FALSE,
                          show_col = FALSE,
                          show_plot = TRUE)
write.csv(Chemokines,file="IL3-TCGA-immune_Chemokines.csv")

Cytolytic<-iobr_cor_plot(pdata_group= pdata_group,
                         id1= "ID",
                         feature_data= fpkm,
                         id2= "ID",
                         target= NULL,
                         group = "BIOCARTA_IL3_PATHWAY",
                         is_target_continuous  = FALSE,
                         padj_cutoff = 1,
                         index  = 1,
                         category= "gene",
                         signature_group= signature_collection[["Cytolytic_Activity_Rooney_et_al"]],
                         palette_box= "paired1",
                         palette_corplot= "pheatmap",
                         palette_heatmap= 4,
                         feature_limit= 26,
                         character_limit  = 30,
                         show_heatmap_col_name = FALSE,
                         show_col = FALSE,
                         show_plot = TRUE)
write.csv(Cytolytic,file="IL3-TCGA-immune_checkpoint.csv")

#statistic=logfc

res_imucheck1<-iobr_cor_plot(pdata_group= pdata_group,
                             id1= "ID",
                             feature_data= fpkm,
                             id2= "ID",
                             target= NULL,
                             group = "WP_IL4_SIGNALING_PATHWAY",
                             is_target_continuous  = FALSE,
                             padj_cutoff = 1,
                             index  = 1,
                             category= "gene",
                             signature_group= signature_collection[c(4)],
                             palette_box= "paired1",
                             palette_corplot= "pheatmap",
                             palette_heatmap= 4,
                             feature_limit= 26,
                             character_limit  = 30,
                             show_heatmap_col_name = FALSE,
                             show_col = FALSE,
                             show_plot = TRUE)
write.csv(res_imucheck1,file="IL4-TCGA-immune_checkpoint.csv")
Chemokines1<-iobr_cor_plot(pdata_group= pdata_group,
                           id1= "ID",
                           feature_data= fpkm,
                           id2= "ID",
                           target= NULL,
                           group = "WP_IL4_SIGNALING_PATHWAY",
                           is_target_continuous  = FALSE,
                           padj_cutoff = 1,
                           index  = 1,
                           category= "gene",
                           signature_group= signature_collection[["Chemokines_Li_et_al"]],
                           palette_box= "paired1",
                           palette_corplot= "pheatmap",
                           palette_heatmap= 4,
                           feature_limit= 26,
                           character_limit  = 30,
                           show_heatmap_col_name = FALSE,
                           show_col = FALSE,
                           show_plot = TRUE)
write.csv(Chemokines1,file="IL4-TCGA-immune_Chemokines.csv")

Cytolytic1<-iobr_cor_plot(pdata_group= pdata_group,
                          id1= "ID",
                          feature_data= fpkm,
                          id2= "ID",
                          target= NULL,
                          group = "WP_IL4_SIGNALING_PATHWAY",
                          is_target_continuous  = FALSE,
                          padj_cutoff = 1,
                          index  = 1,
                          category= "gene",
                          signature_group= signature_collection[["Cytolytic_Activity_Rooney_et_al"]],
                          palette_box= "paired1",
                          palette_corplot= "pheatmap",
                          palette_heatmap= 4,
                          feature_limit= 26,
                          character_limit  = 30,
                          show_heatmap_col_name = FALSE,
                          show_col = FALSE,
                          show_plot = TRUE)
write.csv(Cytolytic1,file="IL4-TCGA-immune_checkpoint.csv")
##plot

data<-rbind(Chemokines,Cytolytic)
View(data)
data<-as.data.frame(data)
rownames(data)<-data$sig_names
data[,c(5)]<-round(data[,c(5)],2)
colnames(data)[5]<-"logFC"
data<-data[order(data$logFC,decreasing = T),]
data$sig_names<-as.factor(data$sig_names)
ggplot(data,aes(x=group,y=sig_names))+geom_tile(aes(fill=logFC))+
  geom_text(aes(label = logFC),col='black',cex=2)+
  scale_fill_gradient2(limits=c(-2, 2),midpoint = 0,
                       low = "blue",mid="white",high = "red"
  )+
  theme(axis.text.y = element_text(size = 5))
