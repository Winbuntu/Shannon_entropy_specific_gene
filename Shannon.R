library(gplots)

Stage.fpkm.table = read.table("GSE66582_stage_fpkm.txt",head=T,row.names = 1)



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


Q.stat.matrix = Q.gt.matrix.compute(Stage.fpkm.table+0.01)


