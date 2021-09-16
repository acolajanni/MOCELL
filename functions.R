#' Filtering Gene of Interest in a Dataframe
#'
#' @param df Dataframe that summarise
#' with a column named 'Gene' and the last column is supposed to be the log2FoldChange
#' @param threshold numerical value. 
#' It is the absolute value of threshold used to subset the data by the log2foldChange
#'
#' @return A vector of the 'Gene' column
#' @export
#'
#' @examples
filter_GOI<- function(df, column_to_keep = "Gene", threshold = 1) {

  # Using the log2 Fold change as the last column of the dataframe
  GOI = subset(df[column_to_keep], 
               vapply(X = df[,ncol(df)], abs, numeric(1)) >= threshold )
  #Returning the gene inside a vector
  return(as.vector(GOI))
}

GOI_FUR_1[,c('Gene','bnum')]
test=filter_GOI(GOI_FUR_1,c('Gene','bnum'))

#Regarder upset plot pour voir les intersections entre gene impliques méto fer / resistance / dypiridil / snRNA
