################################################################################
#          E. C O L I   A N T I B I O T I C S   R E S I S T A N C E            #
#                                                                              #
# > September 2021                                                             #                                                
# > Script : projet.R                                                          #                                                        
# > Fonctions "main.R"                                                         #                                                        
# @ COLAJANNI Antonin                                                          #
# @ ASLOUDJ Yanis                                                              #
################################################################################

# I. Select the genes of interest by using multiple data files :

  # a. the known genes contributing to the antibiotics resistance :
  resi_genes <- read.delim("GOI.txt", header = F)
  resi_genes <- as.vector(unique(resi_genes$V1))
  # b. the genes producing sRNA :
  srna_genes <- read.delim("sRNA_list_coli.txt", header = T)
  srna_genes <- as.vector(unique(srna_genes$sRNA))
  # c. the genes that produces transcription factors
  tf_genes <- read.delim("./P_TF_LABEL_COLI.txt")
  tf_genes <- as.vector(unique(tf_genes$TF))
  # d. the genes over- or under- expressed in the iron & dipyridil metabolism :
  iron_genes <- read.delim("GOI_FUR_1.txt", header = T)
  expr_iron_genes <- filter_GOI(iron_genes, threshold = 1)
  dipy_genes <- read.delim("GOI_FUR_2.txt", header = T)
  expr_dipy_genes <- filter_GOI(dipy_genes, threshold = 1)
  # e. group the genes to get all the genes of interest :
  inter_genes = unique(c(resi_genes,
                         srna_genes,
                         tf_genes,
                         expr_iron_genes, 
                         expr_dipy_genes))
  # 
  
  # f. visualize the intersection between the family of genes
  inter_genes_list = list(
    Resistance = resi_genes,
    sRNA = srna_genes,
    TF = tf_genes,
    Iron = expr_iron_genes, 
    Dipyridil = expr_dipy_genes  )
  
  upset(fromList(inter_genes_list), 
        sets.bar.color = "#56B4E9", 
        order.by = "freq", 
        text.scale = 1.5,
        #mb.ratio = c(0.65,0.35),
        #empty.intersections = "on"
        )
  
# II. Filter the genes expression table using the genes of interest :
  raw_genes_expr <- read.delim("./P_EXPR_antiobios.txt")
  inter_genes_expr <- subset(raw_genes_expr,
                             tolower(raw_genes_expr$IDENTIFIER) %in% 
                               tolower(inter_genes))

# III. Find the correlation on expression data
  
  # a. change the format of the dataframe
  row.names(inter_genes_expr) = inter_genes_expr$IDENTIFIER
  # Remove the identifiers columns
  Filtered_data = inter_genes_expr[-c(1,2)]


  # b. compute the pearson correlations to obtain a square matrix
  corelations = cor(t(Filtered_data), method = "spearman")
  
  # c. transform the square matrix into a 3-column dataframe
  library(reshape2)
  list_correlation=melt(corelations)    
  
  # d. change the names of the columns
  colnames(list_correlation) = c("Gene1","Gene2","values")
  
  # e. visualize pearson's rho distribution 
  hist(list_correlation$values)
  
  # f. remove the redundant information
  filter = as.character(list_correlation[,1])<as.character(list_correlation[,2])
  table = list_correlation[filter,]
  hist(table$values)
  
  # g. filtering gene association by absolute value of pearson's rho (-0.8 or 0.8)
  network = filter_correlation(table, threshold = 0.8)
  
  
  
  
  
  
  
  
  
  
  
  

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



