load("TPM.after.aggregate.after.sum.isoform.RData")

bbb = TPM.after.aggregate.after.sum.isoform

#write.table(bbb,file = "a.txt",sep = "\t")

Stage.fpkm.table = TPM.after.aggregate.after.sum.isoform[,c(10,11,4,1,2,3,5,7)]


####

Q.stat.matrix = Q.gt.matrix.compute(Stage.fpkm.table+0.01)



not.expressed.in.oocyte = Stage.fpkm.table$MII_oocyte <=0.5
expressed.later.than.oocyte = (apply(Stage.fpkm.table[,c(3: (ncol(Stage.fpkm.table)-1) )],1,max) > 1)

early.2.cell.Q.stst = (Q.stat.matrix[,3]<2)
x2.cell.Q.stat = (Q.stat.matrix[,4]<2)
x4.cell.Q.stat = (Q.stat.matrix[,5]<2)
x8.cell.Q.stat = (Q.stat.matrix[,6]<2)
mESC.cell.Q.stat = (Q.stat.matrix[,7]<2)

early.2.cell.list = which(not.expressed.in.oocyte & expressed.later.than.oocyte & early.2.cell.Q.stst &
  (Stage.fpkm.table$`early_2-cell`>=1))

x2.cell.list = which(not.expressed.in.oocyte & expressed.later.than.oocyte & x2.cell.Q.stat &
  (Stage.fpkm.table$`2-cell`>=1))

x4.cell.list = which(not.expressed.in.oocyte & expressed.later.than.oocyte & x4.cell.Q.stat &
    (Stage.fpkm.table$`4-cell`>=1))

x8.cell.list = which(not.expressed.in.oocyte & expressed.later.than.oocyte & x8.cell.Q.stat &
  (Stage.fpkm.table$`8-cell`>=1))



color.palette = colorRampPalette(c("midnightblue","dodgerblue3","white","goldenrod1","darkorange2"), space="Lab")
#pdf(file="./figs/heatmap.pdf",height=10,width=10)
heatmap.2(log2(as.matrix((Stage.fpkm.table[ as.vector(na.omit(c(early.2.cell.list,
                                                           x2.cell.list,
                                                           x4.cell.list,
                                                           x8.cell.list))),] ))+0.1),
          trace="none",density="none",
          Colv = NA,
          key=T,scale="row",
          Rowv = NA,
          dendrogram="both",
          col=color.palette,
          margins=c(7,6)
          
)

# using shannon directly, without computing Q statistics


#hg.stat = apply(Stage.fpkm.table+0.01,1,Hg.compute)


g = unique(rownames(Stage.fpkm.table)[na.omit(c(early.2.cell.list,
                                     x2.cell.list,
                                     x4.cell.list,
                                     x8.cell.list))])

match("Xist",g)

length(unique(na.omit(c(early.2.cell.list,
                        x2.cell.list,
                        x4.cell.list,
                        x8.cell.list))))
