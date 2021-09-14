# Load the raw data
DATA <- read.delim("./P_EXPR_antiobios.txt")

# Load the desired Genes
Desired_genes <- read.delim("GOI.txt", header = F)

# Convert the object into a vector
Desired_genes = as.vector(Desired_genes[[1]])

# Filter the data with the wanted genes
# tolower() function to ensure that there is no case error
Filtered_data = subset(DATA, tolower(DATA$IDENTIFIER) %in% tolower(Desired_genes) )

# Gene ID that are in the Desired_Genes but not in DATA
Not_in_DATA = unique(Desired_genes[!Desired_genes %in% DATA$IDENTIFIER] )

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



