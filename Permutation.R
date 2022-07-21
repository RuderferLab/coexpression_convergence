library(tidyverse)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)

#number of permutations
nperm = as.numeric(args[2])

#files
input <- read.table(args[1], header=F)
names(input) <- "gene"

all.genes <- read.table("cmc.genes.txt",header=T)
names(all.genes) = "gene"
input <- as.data.frame(input[which(input$gene %in% all.genes$gene),])
names(input) <- "gene"


# convert to matrix
c4 = fread("coexpression.txt", sep="\t", header=T)
c4 = as.data.frame(c4)

pheno.meta = c4[,colnames(c4) %in% input$gene]
pheno.meta = cbind(c4[,1], pheno.meta)
colnames(pheno.meta)[1] <- "Gene1"

#pheno.meta = c4[,c("Gene1",input$gene)]
pheno.meta$pheno.metaZ = (rowSums(pheno.meta[,-1]))  /sqrt(ncol(pheno.meta)-1)
pheno.meta = pheno.meta[,c(1,ncol(pheno.meta))]

# permutation
for(j in 1:nperm){
   p.pheno = c4[,c(1,sample(2:ncol(c4),nrow(input)))]
   p.pheno$metaZ = (rowSums(p.pheno[,-1]))/sqrt(ncol(p.pheno)-1)
  
  if(j == 1){
cnt = as.integer(as.logical(abs(p.pheno$metaZ) >= abs(pheno.meta$pheno.metaZ))) #test non-abs
    cntNA = as.integer(as.logical(is.na(cnt)))
    cnt[is.na(cnt)] = 0
  }
  if( j > 1){
    cnt.tmp = as.integer(as.logical(abs(p.pheno$metaZ) >= abs(pheno.meta$pheno.metaZ)))
    cntNA = cntNA + as.integer(as.logical(is.na(cnt.tmp)))
    cnt.tmp[is.na(cnt.tmp)] = 0
    cnt = cnt + cnt.tmp
  }
  if( j %% 1000 == 0){
    print(paste("perm",j,sep=": "))
  }
}

pheno.meta$emp.p = cnt / (nperm-cntNA)
pheno.meta$Nperm = nperm - cntNA
pheno.meta$Ncnt <- cnt

pheno.meta$qval <- p.adjust(pheno.meta$emp.p, method = "fdr")
pheno.meta$p.bonf <- p.adjust(pheno.meta$emp.p, method = "bonferroni")

write.table(pheno.meta,paste0(args[3],".perm.txt"),quote=F,row.names=F,col.names=T, sep="\t")


