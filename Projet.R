################################################################################
#          E. C O L I   A N T I B I O T I C S   R E S I S T A N C E            #
#                                                                              #
# > September 2021                                                             #
# @ COLAJANNI Antonin                                                          #
# @ ASLOUDJ Yanis                                                              #
################################################################################

# I. Select the genes of interest by using multiple data files :

  # a. the known genes contributing to the antibiotics resistance :
  resi_genes <- read.delim("GOI.txt", header = F)
  # b. the genes producing sRNA :
  srna_genes <- read.delim("sRNA_list_coli.txt", header = T)
  # c. the genes over- or under- expressed in the iron & dipyridil metabolism :
  iron_genes <- read.delim("GOI_FUR_1.txt", header = T)
  expr_iron_genes <- filter_GOI(iron_genes, threshold = 1)
  dipy_genes <- read.delim("GOI_FUR_2.txt", header = T)
  expr_dipy_genes <- filter_GOI(dipy_genes, threshold = 1)
  # d. group the genes to get all the genes of interest :
  inter_genes = unique(c(resi_genes$V1, 
                            srna_genes$sRNA, 
                            expr_iron_genes, 
                            expr_dipy_genes))
  # 
  
# II. Filter the genes expression table using the genes of interest :
  raw_genes_expr <- read.delim("./P_EXPR_antiobios.txt")
  inter_genes_expr <- subset(raw_genes_expr,
                             tolower(raw_genes_expr$IDENTIFIER) %in% 
                               tolower(inter_genes))
  
# III. 
  
# change the format of the dataframe
row.names(Filtered_data) = Filtered_data$IDENTIFIER
# Remove the identifiers columns
Filtered_data = Filtered_data[-c(1,2)]
hist(Filtered_data)


# compute the pearson correlations
Corelations = cor(t(Filtered_data))
Corelations_Spearman = cor(t(Filtered_data), method = "spearman")

#------------------------------------------------
# test Pvalues + Correlations
# Compute the approximate Pvalues
library(WGCNA)
Pval = corAndPvalue(t(Filtered_data))$p

# Change the matrix format
list_pval = melt(Pval)
list_correlation=melt(Corelations)  

table = merge(list_correlation, list_pval, by)
table = cbind(list_correlation, list_pval)
table = table[-c(4,5)]
colnames(table) = c("gene1","gene2","cor","pvalue")

# Filtre pvalue
filtre = as.character(table[,1])<as.character(table[,2])
table = table[filtre,]

# Comparer les valeurs de Pvalues + Corrélations (= Qui est le facteur limitant ?)
table.filtered = subset(table,  table$pvalue <= 0.01 & 
                          vapply(X = table$cor, abs, numeric(1)) >= 0.7 )


#------------------------------------------------





library("reshape2")
list_correlation=melt(Corelations)                    
list_correlation2=melt(Corelations_Spearman) 

# Change the names of the columns
colnames(list_correlation) = c("Gene1","Gene2","values")

# Distribution of the correlation coefficient of Pearson
hist(list_correlation$values)
hist(list_correlation2$value)


# Filtration
Cor.filtered = subset(list_correlation, (!list_correlation$Gene1 == list_correlation$Gene2 ))
Cor.filtered.val1 = subset(Cor.filtered, vapply(X = Cor.filtered$values, abs, numeric(1)) >= 0.8)



