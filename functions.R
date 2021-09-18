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
  #Returning the genes inside a vector
  return(as.vector(unlist(GOI[[column_to_keep]])))

}

GOI_FUR_1 = read.delim("./GOI_FUR_1.txt")
GOI_FUR_2 = read.delim("./GOI_FUR_2.txt")
sRNA = read.delim("./sRNA_list_coli.txt")


GOI1=filter_GOI(GOI_FUR_1)
GOI2=filter_GOI(GOI_FUR_2)
GOI3=as.vector(unique(sRNA$sRNA))

library(UpSetR)


#Regarder upset plot pour voir les intersections entre gene impliques métabo fer / resistance / dypiridil / snRNA / !!!!!!!!!! Facteur de transcription !!!!!!!!!!
# Fonction filtration coexpression : en reutilisant filter_GOI


########## Pas besoin de fonction pour deux lignes de code ?
GOI_list <- list(fur1 = GOI1, 
                  fur2 = GOI2, 
                  sRNA = GOI3)

upset(fromList(GOI_list), 
      sets.bar.color = "#56B4E9", 
      order.by = "freq", 
      text.scale = 1.5,
      #mb.ratio = c(0.65,0.35),
      empty.intersections = "on")



#### Ou en une seule ligne :

upset(fromList(list(fur1 = GOI1, 
                    fur2 = GOI2, 
                    sRNA = GOI3)), 
      
      sets.bar.color = "#56B4E9", 
      order.by = "freq", 
      text.scale = 1.5,
      empty.intersections = "on")
