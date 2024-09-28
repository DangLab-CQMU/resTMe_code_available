rm(list = ls())

library('dplyr')
library(ggplot2)
library(Cairo)
library(ggrepel)

setwd('~/code_ava/step1-mixData/')

# Ensure output directory exists
output_dir <- file.path(getwd(), "output")
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Read file names with a specific pattern
fileNames <- list.files(pattern = "*sig*")
n <- length(fileNames)

# Read data files
data <- lapply(fileNames, function(x) {
  read.csv(file.path(getwd(), x))
})

geneinfo = read.delim('./geneinfo_beta.txt',sep = '\t')
geneinfo1 = geneinfo[,c(1,3)]

for (i in 1:n){
  temp<-data[[i]]
  tag<-paste0(unlist(lapply(fileNames[i],function(x){strsplit(x, "_")[[1]][1]})),'mixdata.csv')
  tumor = sub(".*\\-(.*)\\_.*", "\\1", fileNames[i])
  ##volcano plot
  options(bitmapType='cairo')
  plot(cars)
  colnames(temp)[1] = 'ensembl_id'
  TCGA_COAD_sigdiff1 = inner_join(temp,geneinfo1,by = 'ensembl_id')
  table(TCGA_COAD_sigdiff1$sig)

  ##build mix matrix
  #read immune matrix
  cancertypes<-paste0(unlist(lapply(fileNames[i],function(x){substr(fileNames[i],6,9)})),'.txt')
  path1 <- "./ImmucellAI"
  filePath1 <- sapply(cancertypes, function(x){paste(path1,x,sep='/')}) 
  COADimmu = read.delim(filePath1,sep = '\t')
  #remove Nomal sample
  tumorSamId = grep('*.0[1-9]A.', rownames(COADimmu),value = T)#选择肿瘤样本
  normalSamId = grep('*.0[1-9]A.', rownames(COADimmu),value = T,invert = T)
  COADimmu = COADimmu[tumorSamId,]
  TCGA_COAD_sigdiff1 = data.frame(TCGA_COAD_sigdiff1)
  table(duplicated(TCGA_COAD_sigdiff1$gene_name))
  TCGA_COAD_sigdiff1 = TCGA_COAD_sigdiff1[!duplicated(TCGA_COAD_sigdiff1$gene_name),]
  rownames(TCGA_COAD_sigdiff1) = TCGA_COAD_sigdiff1$gene_name
  TCGA_COAD_sigdiff2 = data.frame(t(TCGA_COAD_sigdiff1))
  colnames(TCGA_COAD_sigdiff2) = rownames(TCGA_COAD_sigdiff1)

  TCGA_COAD_sigdiff2$sample = rownames(TCGA_COAD_sigdiff2)
  TCGA_COAD_sigdiff2[1:10,1:10]

  #remove immune genes
  immgene = read.table('./ImmucellAI/immcellGeneFiltered.csv',sep = ',')
  TCGA_COAD_sigdiff3 = TCGA_COAD_sigdiff2[,-which(colnames(TCGA_COAD_sigdiff2) %in% immgene[,1])]
  TCGA_COAD_sigdiff4 = TCGA_COAD_sigdiff1[which(TCGA_COAD_sigdiff1$gene_name %in% colnames(TCGA_COAD_sigdiff3)),]
  setdiff(colnames(TCGA_COAD_sigdiff3),TCGA_COAD_sigdiff4$gene_name)
  table(TCGA_COAD_sigdiff4$sig)
  {
  p <- ggplot(data = TCGA_COAD_sigdiff4, aes(x = log2foldChange_median, y = -log10(FDR),color = sig)) +
      geom_point(alpha=0.8,size = 1.2) +  #绘制散点图，透明度，大小
      scale_color_manual(values = c("#D05146", "#567161"), limits = c('up', 'down')) + 
      labs(x = 'Log2 Fold Change', y = '-Log10 FDR', title = tumor, color = '') +  #坐标轴标题
      geom_vline(xintercept = c(-1, 1), lty = 3, color = 'black') +  #添加阈值线
      geom_hline(yintercept = 2, lty = 3, color = 'black') +
      xlim(min(TCGA_COAD_sigdiff4$log2foldChange_median)-0.5, max(TCGA_COAD_sigdiff4$log2foldChange_median)+0.5) + ylim(0, max(-log10(TCGA_COAD_sigdiff4$FDR))+0.5) +  #定义刻度边界 
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5), 
            legend.position="right", 
            legend.title = element_blank())
    TCGA_COAD_sigdiff4$label=ifelse(TCGA_COAD_sigdiff4$FDR < 0.00001 & abs(TCGA_COAD_sigdiff4$log2foldChange_median) >= 5,TCGA_COAD_sigdiff4$gene_name,"")
    VAL <- p+geom_label_repel(data = TCGA_COAD_sigdiff4,aes(label = label), size = 3, max.overlaps = 11)
    VAL
    file0<-paste0(unlist(lapply(fileNames[i],function(x){strsplit(x, "_")[[1]][1]})),'diff.pdf')
    ggsave(VAL, file=paste0("./output/",file0), width=6.5, height=5.5)
  }

  #mix
  COADimmu$sample = rownames(COADimmu)
  mixdata = inner_join(COADimmu,TCGA_COAD_sigdiff3,by = "sample")
  mixdata = mixdata[-grep('sample',colnames(mixdata))]
  write.csv(mixdata,paste0("./output/",tag),row.names = F)

  #print(paste(tumor,'DGE with immune genes (up/down) is',table(TCGA_COAD_sigdiff1$sig)))
  print(paste(tumor,'DGE without immune genes (up/down) is',table(TCGA_COAD_sigdiff4$sig)))
}
