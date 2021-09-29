################################################################################
#          E. C O L I   A N T I B I O T I C S   R E S I S T A N C E            #
#                                                                              #
# > September 2021                                                             #                                                
# > Script : projet.R                                                          #                                                        
# @ COLAJANNI Antonin                                                          #
# @ ASLOUDJ Yanis                                                              #
################################################################################


### I M P O R T S

library(UpSetR)
library(reshape2)


### F U N C T I O N S

#' Filter a data frame by the absolute value of its last column
#'
#' @param df a dataframe.
#' The last column is supposed to contain numerical values.
#' @param threshold a numerical value.
#' It is the absolute value of the threshold used to subset the data by the last column.
#'
#' @return a dataframe
absolute_filter <- function(df, threshold = 1, sign = FALSE) {
  
  df = subset(df, 
              vapply(X = df[, ncol(df)], abs, numeric(1)) >= threshold)
  
  get_sign <- function(x) {
    if (x > 0) {
      return("pos")
    }
    return("neg")
  }
  
  if (sign) {
    df[, ncol(df)] = vapply(X = df[, ncol(df)], FUN = get_sign, character(1))
  }
  
  return(df)
  
}

#' Title
#'
#' @param df a dataframe 
#' The last column is supposed to contain numerical values.
#' @param column_to_keep a character string or a vector of character strings.
#' It depends on the number of columns wanted at the end. 
#' @param threshold a numerical value. 
#' It is the absolute value of the threshold used to subset the data by the last column.
#'
#' @return a dataframe with as mant column as the length of the "column_to_keep" argument.
#' If only one column is passed in argument, a vector is returned.
#' @export
#'
#' @examples
filter_GOI <- function(df, column_to_keep = "Gene", threshold = 1){
  
  GOI = absolute_filter(df, threshold = threshold)
  
  # if there is more than one column to keep, return a data frame.
  if (length(column_to_keep) > 1){
    GOI = GOI[, column_to_keep]
  }
  
  # otherwise, return a vector.
  else if (length(column_to_keep) == 1){
    GOI = as.vector(GOI[[column_to_keep]])
  }
  
  return(GOI)
  
}

#' Returns a 3-column data frame of correlated variables.
#'
#' @param df a dataframe.
#' @param threshold a numerical value.
#' It is the absolute value of of the correlation threshold.
#' @param cor_method a character string.
#' It is the statistical method used to calculate the correlation ("pearson",
#' "spearman", "kendall").
#' @param sign a boolean.
#' If FALSE, the third column indicates the correlation value. If TRUE, it indicates whether
#' the correlation is positive or negative.
#'
#' @return a dataframe.
#' The first two columns are the pairs of correlated variables, and the third column
#' is an information about the correlation (its value or its sign).
get_associations <- function(df, threshold = 0.8, cor_method = "pearson", sign = TRUE) {
  
  # a. compute the correlations to obtain a square matrix
  cor_matrix = cor(t(df), method = cor_method)
  
  # b. transform the square matrix into a 3-column data frame
  cor_table = melt(cor_matrix)
  
  # c. remove the redundant information
  filter = as.character(cor_table[, 1]) < as.character(cor_table[, 2])
  cor_table = cor_table[filter, ]
  
  # filter the associations using the absolute value of Pearson's rho
  network = absolute_filter(cor_table, threshold, sign)
  
  return(network)
}


################################################################################
################################################################################
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
  expr_iron_genes <- filter_GOI(iron_genes)
  dipy_genes <- read.delim("GOI_FUR_2.txt", header = T)
  expr_dipy_genes <- filter_GOI(dipy_genes)
  # e. group the genes to get all the genes of interest :
  inter_genes = unique(c(resi_genes,
                         srna_genes,
                         tf_genes,
                         expr_iron_genes, 
                         expr_dipy_genes))
  
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

  # a. use the identifier as row names
  row.names(inter_genes_expr) = inter_genes_expr$IDENTIFIER
  inter_genes_expr <- inter_genes_expr[-c(1,2)]
  
  # b. compute a list of edges (co-expressing genes)
  edges <- get_associations(inter_genes_expr)
  
  # c. pre-format the data for a use on Cytoscape
  tmp <- edges[, 3]
  edges[, 3] = edges[, 2]
  edges[, 2] = tmp
  colnames(edges) = NULL
  write.table(edges, "./network.sif", sep = " ", quote = F, row.names = F)
  