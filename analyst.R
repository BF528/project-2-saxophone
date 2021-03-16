setwd("C:/Users/jason/Desktop/BU/BF528/project_2")


#Problem 1

d_e<-read.table(file = "gene_exp.diff", header = T)

#Reordering the table by q_value, ascending is default
s_de<- d_e[
  with(d_e,order(q_value)),
]

#removing columns for the top 10 subset
ts_de = subset(s_de, select = -c(1, 2, 4, 5, 6, 7, 11,14))
write.csv(head(ts_de,10), "top10.csv")

#Histogram of all genes by log2.foldchange
hist(s_de[ ,10], breaks = 25, main = paste("Log2 Fold Change of all Genes Differentially Expressed"),col = c("cyan"), xlab = paste("Log2 Fold Change"))

#Subsetting sig to "yes" and into a histogram
sig_de <- subset(s_de, significant == "yes")
hist(sig_de[ ,10],breaks = 25, main = paste("Log2 Fold Change of Significant Genes Differentially Expressed"),col = c("cyan"), xlab = paste("Log2 Fold Change"))

#subsetting sig data into up/down regulated
up_de<- subset(sig_de, log2.fold_change. > 0,select = c(log2.fold_change., gene))
down_de<-subset(sig_de, log2.fold_change. < 0,select = c(log2.fold_change., gene))

#writing up/down regulated genes into two files
write.table(up_de[2], "up_de.txt", row.names = F, col.name = F, quote = F)
write.table(down_de[2], "down_de.txt", row.names = F, col.name = F, quote = F)


