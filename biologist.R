directory1 <- '/project/bf528/project_2/data/samples'
directory2 <- '/projectnb/bf528/users/saxophone/project2'

# load the data file
P0 <- read.csv("P0/genes.fpkm_tracking", sep = "\t", stringsAsFactors = F)
AD_1 <- read.csv("AD_1/genes.fpkm_tracking", sep = "\t", stringsAsFactors = F)
AD_2 <- read.csv("AD_2/genes.fpkm_tracking", sep = "\t", stringsAsFactors = F)

#7_1 plots of the most prominent GO terms
if (!require("tidyverse")) install.packages("tidyverse")
if (!require("reshape2")) install.packages("reshape2")

subset.gene <- function(x1,x2,x3,target.gene){
  result1 <- rep(0,length(target.gene))
  result2 <- rep(0,length(target.gene))
  result3 <- rep(0,length(target.gene))
  for (i in c(1:length(target.gene))){
    if (target.gene[i] %in% x1$gene_short_name) 
      result1[i] = subset(x1,gene_short_name==target.gene[i])$FPKM
    if (target.gene[i] %in% x2$gene_short_name)
      result2[i] = subset(x2,gene_short_name==target.gene[i])$FPKM
    if (target.gene[i] %in% x3$gene_short_name)
      result3[i] = subset(x3,gene_short_name==target.gene[i])$FPKM
  }
  df <- data.frame(gene=target.gene,P0=result1,AD=(result2+result3)/2) %>% 
    melt("gene") %>% mutate(variable=factor(variable,levels = c("P0","AD")))
  return(df)
}
sar_list <- c("Pdlim5","Pygm","Myoz2","Des","Csrp3",
              "Tcap","Cryab")
sar_gene.fpkm <- subset.gene(P0,AD_1,AD_2,sar_list) 
sar_gene.fpkm

mito_list <- c("Prdx3","Acat1","Echs1","Slc25a11","Phyh")
mito_gene.fpkm <- subset.gene(P0,AD_1,AD_2,mito_list)
mito_gene.fpkm

cellc_list <- c("Cdc7","E2f8","Cdk7","Cdc26","Cdc6","Cdc27",
                "E2f1","Cdc45","Rad51","Aurkb","Cdc23")
cellc_gene.fpkm <- subset.gene(P0,AD_1,AD_2,cellc_list)
cellc_gene.fpkm

#7_2 Compare the results obtained from the DAVID analysis in 6.7
up <- read.csv("up_reg.csv",header=F,stringsAsFactors=F)
names(up) <- up[2,]
up <- up %>% filter(Category!="Category") %>% filter(FDR!="")

up_paper <- read.csv("upreg_paper_reference.csv",stringsAsFactors=F) %>% 
  filter(Category %in% c("GOTERM_MF_FAT","GOTERM_BP_FAT","GOTERM_CC_FAT"))
up_select <- up$Term %in% up_paper$Term
up$Overlap <- ifelse(up_select,'Yes','No')  
write.csv(up,"table1.csv")

down <- read.csv("down_reg.csv",header=F,stringsAsFactors=F)
names(down) <- down[2,]
down <- down %>% filter(Category!="Category") %>% filter(FDR!="")

down_paper <- read.csv("downreg_paper_reference.csv",stringsAsFactors=F) %>% 
  filter(Category %in% c("GOTERM_MF_FAT","GOTERM_BP_FAT","GOTERM_CC_FAT"))
down_select <- down$Term %in% down_paper$Term
down$Overlap <- ifelse(down_select,'Yes','No')  
write.csv(down,"table2.csv")


#7_3 FPKM matrix
P4_1 <- read.csv("4_1/genes.fpkm_tracking",sep="\t",header=T,stringsAsFactors = F)
P4_2 <- read.csv("4_2/genes.fpkm_tracking",sep="\t",header=T,stringsAsFactors = F)
P7_1 <- read.csv("7_1/genes.fpkm_tracking",sep="\t",header=T,stringsAsFactors = F)
P7_2 <- read.csv("7_2/genes.fpkm_tracking",sep="\t",header=T,stringsAsFactors = F)
exp_diff <- read.csv("gene_exp.diff",sep="\t",header=T,stringsAsFactors = F)

# Take the average of FPKM for duplicate gene
sP0 <- P0 %>% group_by(gene_short_name) %>% summarise(P0=mean(FPKM))
sP4_1 <- P4_1 %>% group_by(gene_short_name) %>% summarise(P4_1=mean(FPKM))
sP4_2 <- P4_2 %>% group_by(gene_short_name) %>% summarise(P4_2=mean(FPKM))
sP7_1 <- P7_1 %>% group_by(gene_short_name) %>% summarise(P7_1=mean(FPKM))
sP7_2 <- P7_2 %>% group_by(gene_short_name) %>% summarise(P7_2=mean(FPKM))
sAD_1 <- AD_1 %>% group_by(gene_short_name) %>% summarise(AD_1=mean(FPKM))
sAD_2 <- AD_2 %>% group_by(gene_short_name) %>% summarise(AD_2=mean(FPKM))

df_all <- sP0 %>% left_join(sP4_1,"gene_short_name") %>% 
  left_join(sP4_2,"gene_short_name") %>% left_join(sP7_1,"gene_short_name") %>% 
  left_join(sP7_2,"gene_short_name") %>% left_join(sAD_1,"gene_short_name") %>% 
  left_join(sAD_2,"gene_short_name")
exp_diff <- exp_diff %>% arrange(q_value) %>% head(200)
df_all <- df_all %>% filter(gene_short_name %in% exp_diff$gene) %>% filter(complete.cases(.))

## plot
pdf("figure_out.pdf")

sar_gene.fpkm %>% 
  ggplot(aes(x=variable,y=value,col=factor(gene,levels=sar_list),group=factor(gene,levels=sar_list))) + geom_line() +
  geom_point(data=sar_gene.fpkm,aes(x=variable,y=value,fill=factor(gene,levels=sar_list),shape=factor(gene,levels=sar_list)),show.legend=F) +
  labs(x="",y="FPKM",title="Sarcomere",col="") +
  theme_bw()

mito_gene.fpkm %>% 
  ggplot(aes(x=variable,y=value,col=factor(gene,levels=mito_list),group=factor(gene,levels=mito_list))) + geom_line() +
  geom_point(data=mito_gene.fpkm,aes(x=variable,y=value,fill=factor(gene,levels=mito_list),shape=factor(gene,levels=mito_list)),show.legend=F) +
  labs(x="",y="FPKM",title="Mitochondria",col="") +
  theme_bw()

cellc_gene.fpkm %>% 
  ggplot(aes(x=variable,y=value,col=factor(gene,levels=cellc_list),group=factor(gene,levels=cellc_list))) + geom_line() +
  geom_point(data=cellc_gene.fpkm,aes(x=variable,y=value,fill=factor(gene,levels=cellc_list),shape=factor(gene,levels=cellc_list)),show.legend=F) +
  labs(x="",y="FPKM",title="Cell cycle",col="") +
  theme_bw()

heatmap(as.matrix(df_all[2:8]),
        cexCol=1.2,cexRow=0.1,labRow=F,
        ColSideColors=c("cyan",rep(c("deepskyblue",
                            "blue","blue4"),each=2))
)
legend("topright",legend=c("P0","P4","P7","Ad"),
       fill=c("cyan","deepskyblue","blue","blue4"),
       cex=0.5
)
dev.off()
