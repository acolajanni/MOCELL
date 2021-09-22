################################################################################
#          E. C O L I   A N T I B I O T I C S   R E S I S T A N C E            #
#                                                                              #
# > September 2021                                                             #                                                
# > Script : functions.R                                                       #                                                        
# > Fonction : Fonctions permettant la réalisation du projet                   #                                                        
# @ COLAJANNI Antonin                                                          #
# @ ASLOUDJ Yanis                                                              #
################################################################################


#' Filtering Gene of Interest in a Dataframe
#'
#' @param df Dataframe that summarise
#' with a column named 'Gene' and the last column is supposed to be the log2FoldChange
#' @param threshold numerical value. 
#' It is the absolute value of threshold used to subset the data by the log2foldChange
#' @param column_to_keep Character string.
#' It should be name of the column that contains the gene of interest.
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

filter_correlation <-function(df, threshold = 0.8){
  
  # Filter both gene columns on the absolute value of correlation coefficient
  GOI1 = filter_GOI(df, column_to_keep = "Gene1", threshold = threshold)
  GOI2 = filter_GOI(df, column_to_keep = "Gene2", threshold = threshold)
  
  # Returning the result into a dataframe
  GOI = data.frame("Gene1" = GOI1, "Gene2" = GOI2)
  return(GOI)
  
}

#### Question
## Faire une fonction qui renvoit toutes les colonnes après filtration
## Faire une seconde fonction qui supprime les colonnes non voulues ?
# Qu'en penses-tu ?
# Réutilisation +++++++

# exemple :

#' Filtering by absolute value of the last column of the dataframe
#'
#' @param df Dataframe 
#' The last column is supposed to contain numerical values.
#' @param threshold numerical value. 
#' It is the absolute value of threshold used to subset the data by the last column.
#'
#' @return a dataframe
#' @export
#'
#' @examples
absolute_filter<- function(df, threshold = 1) {
  
  # Using the log2 Fold change as the last column of the dataframe
  GOI = subset(df, 
               vapply(X = df[,ncol(df)], abs, numeric(1)) >= threshold )
  #Returning the genes inside a vector
  return(GOI)
  
}

#' Title
#'
#' @param df Dataframe 
#' The last column is supposed to contain numerical values.
#' @param column_to_keep Character string or a vector of character string depending on the number of column wanted at the end. 
#' @param threshold numerical value. 
#' It is the absolute value of threshold used to subset the data by the last column.
#'
#' @return a dataframe with as much column as the length of the "column_to_keep" argument.
#' If only one column is passed in argument, a vector is returned.
#' @export
#'
#' @examples
filter_GOI2 <- function(df, column_to_keep, threshold = 1){
  
  # On utilise la fonction précédente qui renvoit le tableau complet
  GOI = absolute_filter(df, threshold = threshold)
  
  # On dit à la fonction quelle colonne on veut garder
  
  # Si on a plus d'une colonne on renvoit un dataframe
  if (length(column_to_keep) > 1){
    GOI = GOI[,column_to_keep]
    
  }
  
  # Si on a qu'une colonne on renvoit un vecteur
  else if (length(column_to_keep) == 1){
    GOI = as.vector(GOI[[column_to_keep]])
  }
  
  return(GOI)
  
}
# Exemple :
  
GOI_FUR_1 = read.delim("./GOI_FUR_1.txt")
TEST = filter_GOI2(GOI_FUR_1, column_to_keep = c("bnum",'Gene'))
TEST2 = filter_GOI2(GOI_FUR_1, "Gene")
class(TEST2)






###################################################


GOI_FUR_1 = read.delim("./GOI_FUR_1.txt")
GOI_FUR_2 = read.delim("./GOI_FUR_2.txt")
TF = read.delim("./P_TF_LABEL_COLI.txt")
GOI = read.delim("./GOI.txt")
sRNA = read.delim("./sRNA_list_coli.txt")


GOI1=filter_GOI(GOI_FUR_1)
GOI2=filter_GOI(GOI_FUR_2)
GOI_rna=as.vector(unique(sRNA$sRNA))
GOI_TF = as.vector(unique(TF$TF))

library(UpSetR)


#Regarder upset plot pour voir les intersections entre gene impliques métabo fer / resistance / dypiridil / snRNA / !!!!!!!!!! Facteur de transcription !!!!!!!!!!
# Fonction filtration coexpression : en reutilisant filter_GOI


########## Pas besoin de fonction pour deux lignes de code ?
GOI_list <- list(fur1 = GOI1, 
                  fur2 = GOI2, 
                  sRNA = GOI_rna,
                  TF = GOI_TF)

upset(fromList(GOI_list), 
      sets.bar.color = "#56B4E9", 
      order.by = "freq", 
      text.scale = 1.5,
      #mb.ratio = c(0.65,0.35),
      empty.intersections = "on")



#### Ou en une seule ligne :

upset(fromList(list(fur1 = GOI1, 
                    fur2 = GOI2, 
                    sRNA = GOI_rna,
                    TF = GOI_TF)), 
      sets.bar.color = "#56B4E9", 
      order.by = "freq", 
      text.scale = 1.5,
      empty.intersections = "on")
