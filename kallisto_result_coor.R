TPM = read.table("kallisto_TPM_merged_table_clean.txt",head=T,row.names = 1)

library(gplots)

b = cor(log2(TPM+0.01),method = "spearman")

#heatmap.2(b,margins = c(10,10))

color.palette = colorRampPalette(c("midnightblue","dodgerblue3","white","goldenrod1","darkorange2"), space="Lab")
#pdf(file="./figs/heatmap.pdf",height=10,width=10)
heatmap.2(b,
          trace="none",density="none",
          #Colv = as.dendrogram(h),
          key=T, #scale="row",
          #Rowv = as.dendrogram(h.r),
          dendrogram="both",
          col=color.palette,
          margins=c(15,15)
          
)

duplicated_colanmes = c("2-cell","2-cell","4-cell","4-cell",
                        "8-cell","8-cell","early_2-cell",
                        "early_2-cell","ICM","ICM","ICM",
                        "mESC_rep1_merge","mESC_rep1","mESC_rep1",
                        "mESC_rep2","mESC_rep2_merge","mESC_rep2",
                        "MII_oocyte","MII_oocyte","zygote","zygote")

####
library(plyr)

TPM.for.aggregate = TPM

colnames(TPM.for.aggregate) = duplicated_colanmes

TPM.after.aggregate = sapply(unique(colnames(TPM.for.aggregate)), 
       function(x) rowMeans(TPM.for.aggregate[, colnames(TPM.for.aggregate) == x, drop = FALSE]))


TPM.big = cbind(TPM.after.aggregate,TPM)
b2 = cor(log2(TPM.big+0.01),method = "spearman")

heatmap.2(b2,
          trace="none",density="none",
          #Colv = as.dendrogram(h),
          key=T, #scale="row",
          #Rowv = as.dendrogram(h.r),
          dendrogram="both",
          col=color.palette,
          margins=c(15,15)
          
)



