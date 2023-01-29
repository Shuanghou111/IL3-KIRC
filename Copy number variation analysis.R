#The sample names were modified and divided into high expression and low expression groups.
tumor<-read.table("kirc-cnv.txt",header = T)
tumorC<-tumor[grep("01A",tumor$Sample),]
tumorC$Sample<-substr(tumorC$Sample,1,12)
tumor<-tumorC[tumorC$Sample%in%rownames(new4),]
TCGA_CNV_IL3H1<-tumor[tumor$Sample%in%rownames(new4)[new4$BIOCARTA_IL3_PATHWAY=="high"],]
TCGA_CNV_IL3L1<-tumor[tumor$Sample%in%rownames(new4)[new4$BIOCARTA_IL3_PATHWAY=="low"],]


##Prepare chromosome information

chrom_extract <- function(BSgenome.hg  = NULL) {
  if (is.null(BSgenome.hg )) stop("NULL object !", call. = FALSE)
  obj <- list(species = GenomeInfoDb::organism(BSgenome.hg), genomebuild = BSgenome::providerVersion(BSgenome.hg))
  df <- data.frame(chrom = BSgenome::seqnames(BSgenome.hg), chrN = seq_along(BSgenome::seqnames(BSgenome.hg)), chr.length = GenomeInfoDb::seqlengths(BSgenome.hg), stringsAsFactors = FALSE)
  df <- df[1:24,]
  df$chr.length.sum <- cumsum(as.numeric(df$chr.length))
  df$chr.length.cumsum <- c(0, df$chr.length.sum[-nrow(df)])
  df$middle.chr <- round(diff(c(0, df$chr.length.sum)) /2)
  df$middle.chr.genome <- df$middle.chr + df$chr.length.cumsum
  obj$chromosomes <- df
  obj$chrom2chr <- sapply(obj$chromosomes$chrom, function(k) { obj$chromosomes$chrN[obj$chromosomes$chrom == k]}, simplify = FALSE)
  obj$chr2chrom <- sapply(obj$chromosomes$chrN, function(k) { obj$chromosomes$chrom[obj$chromosomes$chrN == k]}, simplify = FALSE)
  names(obj$chr2chrom) <- obj$chromosomes$chrN
  obj$genome.length <- sum(as.numeric(obj$chromosomes$chr.length), na.rm = TRUE)
  return(obj)
}
# Extract a chromosomes reference loci
BSgenome.hg = "BSgenome.Hsapiens.UCSC.hg19"
BSg.obj <- getExportedValue(BSgenome.hg, BSgenome.hg)
genome.version <- BSgenome::providerVersion(BSg.obj)
chrom <- chrom_extract(BSg.obj)
str(chrom)
#Cnv analysis of IL3-L group,"scores-IL3L.gistic" is derived from online analysis.
scores_L <- read.table("scores-IL3L.gistic", sep="\t",header=T,stringsAsFactors = F)
head(scores_L)
unique(scores_L$Chromosome)
#Change the names of chromosomes from Arabic numerals to the forms of "chr1" and "CHX"
#if is 23，scores[scores$Chromosome==23, "Chromosome"] <- "X"
#if is 24，scores[scores$Chromosome==24, "Chromosome"] <- "Y"
chrID <- unname(unlist(chrom$chrom2chr[as.character(paste0("chr",scores_L$Chromosome))]))
# Important step for accurate length to match back to continual chrom loci
scores_L$Start.geno <- scores_L$Start + chrom$chromosomes$chr.length.cumsum[chrID]
scores_L$End.geno <- scores_L$End + chrom$chromosomes$chr.length.cumsum[chrID]
# Prepare input data for ploting
scores.amp <- scores_L[scores_L$Type=="Amp",]
scores.amp$G.score <- scores.amp$G.score * 1
scores.del <- scores_L[scores_L$Type=="Del",]
scores.del$G.score <- scores.del$G.score * -1
scores <- rbind.data.frame(scores.amp,scores.del)
# seg.col = list(gain = "red", outscale.gain = "darkred", loss = "blue", outscale.red = "midnightblue")
ylim <- c(min(scores$G.score) - 0.1, max(scores$G.score) + 0.1)
table(new4$BIOCARTA_IL3_PATHWAY)

title <- paste0("TCGA KIRC IL3-L copy number gistic score", " ", "n=", 53)

plot(scores.amp$Start.geno, scores.amp$G.score,
     pch = ".", type='h',cex = 2, xaxs = "i", yaxs = "i", 
     xlim = c(0,chrom$genome.length), ylim = ylim,
     main = title, cex.main = 1, ylab = "gistic score", xlab = NA,
     cex.lab = 1, col = adjustcolor("darkred", alpha.f = .8), xaxt = "n", lwd = 2, las=1) # las=1 rotating axis labels in R
lines(scores.del$Start.geno, scores.del$G.score, type='h', lwd = 2, col = adjustcolor("midnightblue", alpha.f = .8))
ink <- chrom$chromosomes$chrN %in% chrID
yrange = abs(diff(ylim))
m.pos <- c(ylim[1]+0.15,ylim[2]-0.15)
m.mod <- -(chrom$chromosomes$chrN[ink] %% 2) +2
try(text(x = chrom$chromosomes$middle.chr.geno[ink], y = m.pos[m.mod], labels = chrom$chromosomes$chrom[ink], cex = 0.6))
abline(h = 0.0, col = 1, lwd = 1, lty = 3)
abline(v = c(0,chrom$chromosomes$chr.length.sum), col = 1, lty = 3, lwd = 1)

col1 <- adjustcolor("darkred", alpha.f = .8)
col2 <- adjustcolor("midnightblue", alpha.f = .8)
# The position of the legend can be specified also using the following keywords : "bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right" and "center".
legend("topleft", c("gain","loss"), cex=0.6, bty="n", fill=c(col1,col2))

#Draw the "percentage/frequency" of the sample
scores.amp <- scores_L[scores_L$Type=="Amp",]
scores.amp$frequency<- scores.amp$frequency * 100
scores.del <- scores_L[scores_L$Type=="Del",]
scores.del$frequency<- scores.del$frequency * -100

# copy number percentage plot
# seg.col = list(gain = "red", outscale.gain = "darkred", loss = "blue", outscale.red = "midnightblue")
ylim<- c(-90,90)
title=paste0("KIRC IL3-L copy number percentage"," ","n=",53)

plot(scores.amp$Start.geno, scores.amp$frequency,
     pch = ".", type='h',cex = 1, xaxs = "i", yaxs = "i", 
     xlim = c(0,chrom$genome.length), ylim = ylim,
     main = title, cex.main = 1, ylab = "gain/loss percentage in cohort", xlab = NA,
     cex.lab = 1, col = adjustcolor("darkred", alpha.f = .8), xaxt = "n", lwd = 2, las=1) # las=1 rotating axis labels in R
lines(scores.del$Start.geno, scores.del$frequency, type='h', lwd = 2, col = adjustcolor("midnightblue", alpha.f = .5))
ink <- chrom$chromosomes$chrN %in% chrID
yrange = abs(diff(ylim))
m.pos <- c(-80,80)
m.mod <- -(chrom$chromosomes$chrN[ink] %% 2) +2
try(text(x = chrom$chromosomes$middle.chr.geno[ink], y = m.pos[m.mod], labels = chrom$chromosomes$chrom[ink], cex = 0.6))
abline(h = 0.0, col = 1, lwd = 1, lty = 3)
abline(v = c(0,chrom$chromosomes$chr.length.sum), col = 1, lty = 3, lwd = 1)

col1 <- adjustcolor("darkred", alpha.f = .8)
col2 <- adjustcolor("midnightblue", alpha.f = .8)
# The position of the legend can be specified also using the following keywords : "bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right" and "center".
legend("topleft", c("gain","loss"), cex=0.6, bty="n", fill=c(col1,col2))





#IL3-H
scores_H <- read.table("scores-IL3H.gistic", sep="\t",header=T,stringsAsFactors = F)
head(scores_H)
unique(scores_H$Chromosome)

chrID <- unname(unlist(chrom$chrom2chr[as.character(paste0("chr",scores_H$Chromosome))]))
# Important step for accurate length to match back to continual chrom loci
scores_H$Start.geno <- scores_H$Start + chrom$chromosomes$chr.length.cumsum[chrID]
scores_H$End.geno <- scores_H$End + chrom$chromosomes$chr.length.cumsum[chrID]
# Prepare input data for ploting
scores.amp <- scores_H[scores_H$Type=="Amp",]
scores.amp$G.score <- scores.amp$G.score * 1
scores.del <- scores_H[scores_H$Type=="Del",]
scores.del$G.score <- scores.del$G.score * -1
scores1 <- rbind.data.frame(scores.amp,scores.del)
# seg.col = list(gain = "red", outscale.gain = "darkred", loss = "blue", outscale.red = "midnightblue")
ylim <- c(min(scores1$G.score) - 0.1, max(scores1$G.score) + 0.1)
title <- paste0("TCGA KIRC IL3-H copy number gistic score", " ", "n=473")

plot(scores.amp$Start.geno, scores.amp$G.score,
     pch = ".", type='h',cex = 1, xaxs = "i", yaxs = "i", 
     xlim = c(0,chrom$genome.length), ylim = ylim,
     main = title, cex.main = 1, ylab = "gistic score", xlab = NA,
     cex.lab = 1, col = adjustcolor("darkred", alpha.f = .8), xaxt = "n", lwd = 2, las=1) # las=1 rotating axis labels in R
lines(scores.del$Start.geno, scores.del$G.score, type='h', lwd = 2, col = adjustcolor("midnightblue", alpha.f = .8))
ink <- chrom$chromosomes$chrN %in% chrID
yrange = abs(diff(ylim))
m.pos <- c(ylim[1]+0.15,ylim[2]-0.15)
m.mod <- -(chrom$chromosomes$chrN[ink] %% 2) +2
try(text(x = chrom$chromosomes$middle.chr.geno[ink], y = m.pos[m.mod], labels = chrom$chromosomes$chrom[ink], cex = 0.6))
abline(h = 0.0, col = 1, lwd = 1, lty = 3)
abline(v = c(0,chrom$chromosomes$chr.length.sum), col = 1, lty = 3, lwd = 1)

col1 <- adjustcolor("darkred", alpha.f = .8)
col2 <- adjustcolor("midnightblue", alpha.f = .8)
# The position of the legend can be specified also using the following keywords : "bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right" and "center".
legend("topleft", c("gain","loss"), cex=0.6, bty="n", fill=c(col1,col2))


scores.amp <- scores_H[scores_H$Type=="Amp",]
scores.amp$frequency<- scores.amp$frequency * 100
scores.del <- scores_H[scores_H$Type=="Del",]
scores.del$frequency<- scores.del$frequency * -100

# copy number percentage plot
# seg.col = list(gain = "red", outscale.gain = "darkred", loss = "blue", outscale.red = "midnightblue")
ylim<- c(-90,90)
title=paste0("KIRC IL3-H copy number percentage"," ","n=",473)

plot(scores.amp$Start.geno, scores.amp$frequency,
     pch = ".", type='h',cex = 1, xaxs = "i", yaxs = "i", 
     xlim = c(0,chrom$genome.length), ylim = ylim,
     main = title, cex.main = 1, ylab = "gain/loss percentage in cohort", xlab = NA,
     cex.lab = 1, col = adjustcolor("darkred", alpha.f = .8), xaxt = "n", lwd = 2, las=1) # las=1 rotating axis labels in R
lines(scores.del$Start.geno, scores.del$frequency, type='h', lwd = 2, col = adjustcolor("midnightblue", alpha.f = .5))
ink <- chrom$chromosomes$chrN %in% chrID
yrange = abs(diff(ylim))
m.pos <- c(-80,80)
m.mod <- -(chrom$chromosomes$chrN[ink] %% 2) +2
try(text(x = chrom$chromosomes$middle.chr.geno[ink], y = m.pos[m.mod], labels = chrom$chromosomes$chrom[ink], cex = 0.6))
abline(h = 0.0, col = 1, lwd = 1, lty = 3)
abline(v = c(0,chrom$chromosomes$chr.length.sum), col = 1, lty = 3, lwd = 1)

col1 <- adjustcolor("darkred", alpha.f = .8)
col2 <- adjustcolor("midnightblue", alpha.f = .8)
# The position of the legend can be specified also using the following keywords : "bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right" and "center".
legend("topleft", c("gain","loss"), cex=0.6, bty="n", fill=c(col1,col2))



