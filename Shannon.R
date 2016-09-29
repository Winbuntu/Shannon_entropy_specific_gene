library(gplots)

Stage.fpkm.table = read.table("GSE66582_stage_fpkm.txt",head=T,row.names = 1)



# using Shannon entropy and Q statistics, ideas and formula are come from these papers:
# 1.Promoter features related to tissue specificity as measured by Shannon entropy (Q stst)
# 2.ROKU: a novel method for identification of tissue-specific genes (Shannon)


Pt.g.normalize <- function(x) {
  
  # normalize gene expression level. 
  # this fucntion convert gene expression level into probability,
  # so can be fitted into channon entropy formula
  
  #normalized.values = vector(mode = "numeric",length = length(x))
  normalized.values = x/sum(x)
  return(normalized.values)
  
}


Hg.compute <- function(Pt.g){
  
  #take Pt.g as input
  # compute shannon entropy for each gene
  
  sum(   -(Pt.g)*log2(Pt.g)   )
  
} 

Q.g.t.compute <- function(x){
  
  # compute Q statistics for each gene, 
  # based on shannon entropy and  normalized gene expression level.
  
  #Q.g.t.vector = vector(mode = "numeric",length = length(x))
  
  Pt.g = Pt.g.normalize(x)
  Hg = Hg.compute(Pt.g)
  #print(Hg)
  #print(Pt.g)
  #print(-log2(Pt.g))
  Q.g.t.vector = Hg * (-(log2(Pt.g)))
  return(Q.g.t.vector)
}

#Q.g.t.compute(c(5,5,5,5,50))


########################
# this function compute tissue specific shannon entropy 
# for each gene and each tissue

Q.gt.matrix.compute <- function(FPKM.table){
  
  Q.gt.matrix = matrix(0, nrow = nrow(FPKM.table), ncol = ncol(FPKM.table))
  
  for(i in c(1:nrow(FPKM.table))) {
    
    Q.gt.matrix[i,] = Q.g.t.compute(   as.numeric(FPKM.table[i,] )  )
  }
  #apply(FPKM.table,1,Q.g.t.compute)
  return(Q.gt.matrix)
}

##############################
## following are analysis part.

# we use these method to identify stage specific genes using data from:
# The landscape of accessible chromatin in mammalian preimplantation embryos
###  GSE66582, this is scRNA. This data file is from this record

Q.stat.matrix = Q.gt.matrix.compute(Stage.fpkm.table+0.01)


### find non-maternal stage specific genes
# oocytes (FPKM ≤ 0.5) but that are activated (FPKM>1) after ZGA

not.expressed.in.oocyte = Stage.fpkm.table$MII_oocyte <=0.5
expressed.later.than.oocyte = (apply(Stage.fpkm.table[,c(3: (ncol(Stage.fpkm.table)-1) )],1,max) > 1)

early.2.cell.Q.stst = (Q.stat.matrix[,3]<2)
x2.cell.Q.stat = (Q.stat.matrix[,4]<2)
x4.cell.Q.stat = (Q.stat.matrix[,5]<2)
x8.cell.Q.stat = (Q.stat.matrix[,6]<2)
mESC.cell.Q.stat = (Q.stat.matrix[,7]<2)

early.2.cell.list = not.expressed.in.oocyte & expressed.later.than.oocyte & early.2.cell.Q.stst &
  (Stage.fpkm.table$early_2cell>=3)

x2.cell.list = not.expressed.in.oocyte & expressed.later.than.oocyte & x2.cell.Q.stat &
  (Stage.fpkm.table$X2cell>=3)

x4.cell.list = not.expressed.in.oocyte & expressed.later.than.oocyte & x4.cell.Q.stat &
  (Stage.fpkm.table$X4cell>=3)

x8.cell.list = not.expressed.in.oocyte & expressed.later.than.oocyte & x8.cell.Q.stat &
  (Stage.fpkm.table$X8cell>=3)



color.palette = colorRampPalette(c("midnightblue","dodgerblue3","white","goldenrod1","darkorange2"), space="Lab")
#pdf(file="./figs/heatmap.pdf",height=10,width=10)
heatmap.2(as.matrix((Stage.fpkm.table[ as.vector(na.omit(early.2.cell.list | x2.cell.list|
                                                           x4.cell.list | x8.cell.list)),] )),
          trace="none",density="none",
          Colv = NA,
          key=T,scale="row",
          
          dendrogram="both",
          col=color.palette,
          margins=c(7,6)

)

# using shannon directly, without computing Q statistics


hg.stat = apply(Stage.fpkm.table+0.01,1,Hg.compute)






