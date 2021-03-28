rm(list = ls()) # Clear variables

library(copynumber)
data(lymphoma)
sub.lymphoma <- subsetData(data=lymphoma,sample=1:3)
sub.lymphoma[1:10,]
lymph.wins <- winsorize(data=sub.lymphoma,verbose=FALSE)
lymph.wins[1:10,]
wins.res <- winsorize(data=sub.lymphoma,return.outliers=TRUE,verbose=FALSE)
wins.res$wins.outliers[1:10,]

single.seg <- pcf(data=lymph.wins,gamma=12,verbose=FALSE)
head(single.seg)
plotGenome(data=sub.lymphoma,segments=single.seg,sample=1,cex=3)
plotSample(data=sub.lymphoma,segments=single.seg,layout=c(5,5),sample=1,cex=3)


lymphoma.res <- pcf(data=lymphoma,gamma=12,verbose=FALSE)

plotFreq(segments=lymphoma.res,thres.gain=0.2,thres.loss=-0.1)

lymphoma.resTTT <- lymphoma.res[,c(-3,-6)]
plotFreq(segments=lymphoma.resTTT,thres.gain=0.2,thres.loss=-0.1)

chr.from <- c(2,12,4)
pos.from <- c(168754669,847879349,121809306)
chr.to <- c(14,21,17)
pos.to <- c(6147539,301955563,12364465)
cl <- c(1,1,2)
arcs <- cbind(chr.from,pos.from,chr.to,pos.to,cl)
plotCircle(segments=lymphoma.res,thres.gain=0.15,arcs=arcs)

plotHeatmap(segments=lymphoma.res,upper.lim=0.3)
plotAberration(segments=lymphoma.res,thres.gain=0.2)

#
# Main file
setwd(getwd()) ## Set current working directory
PathName <- getwd() ## Set output directroy

RVersion = "20210328V1"
dir.create(paste0(PathName,"/",RVersion))

## Load files
CNV_PDAC <- read.table(paste0(PathName,"/SNP6_nocnv_genomicSegmentTTT3.txt"),  # 資料檔名 
                   header=T,          # 資料中的第一列，作為欄位名稱
                   sep="\t")           # 將Tab視為分隔符號來讀取資料
#CNV_PDACTTT <- CNV_PDAC[1:15,]
#colnames(CNV_PDACTTT) <- colnames(lymphoma.resTTT)

CNV_PDAC_CCC <- as.data.frame(CNV_PDAC[,2])
colnames(CNV_PDAC_CCC) <- c("chr")
library(tidyr)
new.dataC <- separate(CNV_PDAC_CCC, chr, c( "chrOri","chr"), "chr")
CNV_PDAC[,2] <- new.dataC[,2]

colnames(CNV_PDAC) <- colnames(lymphoma.resTTT)
#CNV_PDACTTT2 <- gsub(pattern = 'chr',replacement = '',x = CNV_PDACTTT)

plotFreq(segments=CNV_PDAC,thres.gain=0.2,thres.loss=-0.1)


chr.from <- c(2,12,4)
pos.from <- c(168754669,847879349,121809306)
chr.to <- c(14,21,17)
pos.to <- c(6147539,301955563,12364465)
cl <- c(1,1,2)
arcs <- cbind(chr.from,pos.from,chr.to,pos.to,cl)
plotCircle(segments=CNV_PDAC,thres.gain=0.15,arcs=arcs)

#plotHeatmap(segments=CNV_PDAC,upper.lim=0.3)
#plotAberration(segments=CNV_PDAC,thres.gain=0.2)

#-符號沒差
#資料去掉chr: chr1 -> 1
#Hearmap 需要arm和nprobe等資料




############# Load RNA Data #############
FileName <- c("Xena_TCGA_PAAD_RNAGE")
Target_gene_name <- c("TOP2A")


## Import genetic data file
GeneExp_Ori <- read.table(paste0(PathName,"/",FileName),  # 資料檔名
                          header=T,          # 資料中的第一列，作為欄位名稱
                          sep="")
# GeneExp_Ori <- read.table(paste0(PathName,"/Xena_TCGA_LGG_GE"),  # 資料檔名
#                           header=F,          # 資料中的第一列，作為欄位名稱
#                           sep="")
###################################### Note ######################################
# paste0 ==> concatenate strings without any separation/delimiter
# paste("Hello", "World", sep = "-") ==> concatenate strings with seperator "-"
##################################################################################


GeneExp <- GeneExp_Ori 
row.names(GeneExp) <- GeneExp[,1]
#GeneExp <- GeneExp[1:length(GeneExp[,1]), 2:length(GeneExp[1,])]
GeneExp <- GeneExp[, -1]

# load package 'data.table' 
library(data.table)

# Extract data with Target_gene_name
Target_gene <- GeneExp[Target_gene_name,]

# load package 'dplyr'
library(dplyr) # Basic data manupilation tools
Target_gene_Mean <- rowMeans(data.matrix(Target_gene))
Target_gene_SD <- sd(data.matrix(Target_gene))

GeneExp_High <- GeneExp[,GeneExp[Target_gene_name,] >= Target_gene_Mean+Target_gene_SD]
GeneExp_Low <- GeneExp[,GeneExp[Target_gene_name,] <= Target_gene_Mean-Target_gene_SD]
GeneExp_Medium <- GeneExp[,GeneExp[Target_gene_name,] <= Target_gene_Mean+Target_gene_SD & GeneExp[Target_gene_name,] >= Target_gene_Mean-Target_gene_SD]
#####
#library("plyr")
#library(dplyr)

# High
GeneExp_High_ColN <- as.data.frame(colnames(GeneExp_High))
GeneExp_High_ColN[,2] <- GeneExp_High_ColN
colnames(GeneExp_High_ColN) <- c("sampleID","sampleID.Check")
CNV_PDAC_High <- CNV_PDAC
#CNV_PDAC_3 <- CNV_PDAC_2[find(as.character((CNV_PDAC_2[,1]%in%GeneExp_High_ColN)), simple.words = TRUE),]
CNV_PDAC_High_2 <- left_join(CNV_PDAC_High,GeneExp_High_ColN,by="sampleID")
# row.names(CNV_PDAC_2) <- CNV_PDAC_2[,1]
CNV_PDAC_High_3 <- na.omit(CNV_PDAC_High_2)
plotFreq(segments=CNV_PDAC_High_3,thres.gain=0.2,thres.loss=-0.1,ylim=c(-100,100))

# Low
GeneExp_Low_ColN <- as.data.frame(colnames(GeneExp_Low))
GeneExp_Low_ColN[,2] <- GeneExp_Low_ColN
colnames(GeneExp_Low_ColN) <- c("sampleID","sampleID.Check")
CNV_PDAC_Low <- CNV_PDAC
CNV_PDAC_Low_2 <- left_join(CNV_PDAC_Low,GeneExp_Low_ColN,by="sampleID")
CNV_PDAC_Low_3 <- na.omit(CNV_PDAC_Low_2)
plotFreq(segments=CNV_PDAC_Low_3,thres.gain=0.2,thres.loss=-0.1,ylim=c(-100,100))

# Medium
GeneExp_Medium_ColN <- as.data.frame(colnames(GeneExp_Medium))
GeneExp_Medium_ColN[,2] <- GeneExp_Medium_ColN
colnames(GeneExp_Medium_ColN) <- c("sampleID","sampleID.Check")
CNV_PDAC_Medium <- CNV_PDAC
CNV_PDAC_Medium_2 <- left_join(CNV_PDAC_Medium,GeneExp_Medium_ColN,by="sampleID")
CNV_PDAC_Medium_3 <- na.omit(CNV_PDAC_Medium_2)
plotFreq(segments=CNV_PDAC_Medium_3,thres.gain=0.2,thres.loss=-0.1,ylim=c(-100,100))


#All
plotFreq(segments=CNV_PDAC,thres.gain=0.2,thres.loss=-0.1,ylim=c(-100,100))
