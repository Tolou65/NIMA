PlotInput <- function(x, index_nominal,
                      index_numeric,
                      color = NULL,
                      verbose = FALSE) {
  # Plot the given matrix with heatmap plot.
  #
  # Args:
  #   x: Data Frame containing the missing values. 3 columns (at least 2 nominal 
  #   columns)
  #   index_nominal: a vector with two values. the indexes for the nominal columns
  #   (the first value indicating the rows objects and the second value indicating
  #   the column object)
  #   index_numeric: index for numeric values. (this is value for selecting the 
  #   column which contains our numeric values and we change it to
  #   the matrix for missing value investigation and imputation) 
  #   color: color plate use for the heatmap, default is NULL. (Find the option 
  #   for the future used function later and include it as options).
  #   verbose: If TRUE, the plot is saved as the .png file in the working directory.
  #   Default is FALSE.
  #
  # Returns:
  #   Plot the heatmap of the input matrix.
  
  
  # show the axis values by appropriate rate.
  library(reshape2)
  library(tidyverse)
  library(RColorBrewer)
  library(purrr)
  
  input_type <- c("matrix", "data.frame")
  if(!(class(x)[1] %in% input_type))
    stop("The input should be matrix or data frame")
  
  if(is_empty(colnames(x)))
    stop("The input data frame should include column names")
  
  if(!(index_nominal %in% c(1:dim(x)[2])))
    stop("The nominal indecies should be two columns of the input data frame")
  
  if(!(index_numeric %in% c(1:dim(x)[2])))
    stop("The numerical index should be a column of the input data frame")
  
  cn <- colnames(x)
  
  x_mat <- acast( x, eval(parse(text = cn[index_nominal[2]])) ~ eval(parse(text = cn[index_nominal[1]]))
                   , value.var = cn[index_numeric])
  if (nrow(x_mat)>ncol(x_mat)) {
    par(pin = c(3,5))} else {
      par(pin = c(5,3))}
  if (is.null(color)){
    col_palette = brewer.pal(5,"Blues")
  } else {
    col_palette = color 
  }
  
  if ( verbose == TRUE){ 
    png(filename = "./Plot_missing.png")
    image(1:ncol(x_mat), 1:nrow(x_mat), t(apply(x_mat, 2, rev)),
        ylab=cn[index_nominal[2]], xlab=cn[index_nominal[1]], col = col_palette)
    dev.off()} else{
      image(1:ncol(x_mat), 1:nrow(x_mat), t(apply(x_mat, 2, rev)),
            ylab=cn[index_nominal[2]], xlab=cn[index_nominal[1]], col = col_palette)
    }
  return(x_mat)
}



dt <- readRDS("~/NIMA/NIMA/beatamldrugsInddf.rds")
dt2 <- readRDS("~/NIMA/NIMA/drugcombdrugsInddf (1).rds") 


index_nominal <- c(1,2)
index_numeric <- 3
