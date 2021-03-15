#Read in the genes.fpkm_tracking file 
fpkm <- read.table("/projectnb/bf528/users/saxophone/project2/scripts/programmer/P0_1_cufflinks/genes.fpkm_tracking", header = TRUE, sep = "\t")
#Filter out genes where FPKM<500
fpkm <- fpkm[fpkm$FPKM>500,]
#Create histogram of gene name vs. log(FPKM)
par(mar=c(6,5,4,2))
options(scipen=999)
barplot(fpkm$FPKM,names.arg = fpkm$gene_short_name,main="Gene Expression - FPKM",ylab="log(FPKM)",ylim=c(100,10^6),log="y",las=2)
box(which="plot",lty="solid")
