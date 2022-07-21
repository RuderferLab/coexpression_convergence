#Calculate coexpression

#input gene expression matrix 
x <- read.table("gene.expression.txt", header=T)

#run pairwise correlation on expression data for pearson 
cor.f <- cor(x, use = "pairwise.complete.obs", method = "pearson")

#write out the coexpression datafile 
write.table(cor.f, "coexpression.txt", quote=F, row.name=F)