ExtractNonMissing <- function(x, shape = "Square",
                              verbose = FALSE,
                              color = NULL,
                              row.vars = NULL,
                              col.vars = NULL) {
  # Extract the non-missing islands in a given matrix.
  #
  # Args:
  #   x: matrix containing the missing values.
  #   shape: string array indicating the desire shape of the output, other options are "Square", "Rectangular_row",
  #          "Rectangular_col", "Rectangular_element_max", and "All".
  #   verbose: If TRUE, plot the heatmap of the extracted non-missing matrix; if not, not. Default is FALSE.
  #
  # Returns:
  #   The bigest non-missing matrix.
  library(dplyr)
  library(stringr)
  library(reshape2)
  library(tidyverse)
  library(RColorBrewer)
  input_type <- c("matrix", "data.frame")
  if(!(class(x)[1] %in% input_type))
    stop("The input should be matrix or data frame")
  #making new col and row names in order using mutate properly and replace back the original names at the end.
  R.names <- row.names(x)
  C.names <- c(colnames(x))
  
  if(is.null(R.names) | is.null(C.names))
    warning("There is no row names or column names for given input")
  
  
  if(is.null(row.vars) | is.null(col.vars))
    warning("There is no name specified for samples on x-axis lable or y-axis lable")
  
  row_size <- dim(x)[1]
  col_size <- dim(x)[2]
  
  R.Names.mod <- paste0("x", seq(1:row_size))
  C.Names.mod <- paste0("y", seq(1:col_size))
  row.names.Htable <- as.data.frame(cbind(R.names, R.Names.mod))
  col.names.Htable <- as.data.frame(cbind(C.names, C.Names.mod))
  x <- as.data.frame(x, row.names = R.Names.mod)
  colnames(x) <- C.Names.mod
  ###### copy the original matrix in order to start the arranging rows and columns
  x_1 <- x
  x_1[!is.na(x_1)]=1
  # 3 rearrenging of the row and columns
  #1 col arranging
  x_arrage1 <- as.data.frame(x_1, col.names =C.Names.mod, row.names = R.Names.mod)%>%
    rownames_to_column(var = "lab_id") %>%
    mutate_if(is.numeric, replace_na, replace = 0) %>%
    mutate(sumVar = rowSums(.[2:col_size+1]))%>%
    column_to_rownames(var = "lab_id") %>%
    arrange(desc(sumVar)) %>%
    select(-sumVar)%>%
    t %>% data.frame()
  #2 row arranging
  x_arrage2 <- x_arrage1 %>% 
    data.frame() %>%
    rownames_to_column(var = "inh_id") %>%
    mutate(sumVar = rowSums(.[2:row_size+1]))%>%
    column_to_rownames(var = "inh_id") %>%
    arrange(desc(sumVar)) %>% 
    select(-sumVar) %>% 
    t(.) %>% data.frame()
  #3 total rearranging after replacing back NAs as the NAs would be ignor and result in a new arrangment
  x_arrage2 <- na_if(x_arrage2, 0)
  
  x_arrage3  <- x_arrage2 %>% 
    t %>%
    data.frame() %>%
    arrange_if(is_double)
  
  ## Giving the xlab and ylab if they are given by user
  
  if (is.null(row.vars) | is.null(col.vars)) {
    c_lable= "columns"
    r_lable= "rows"
    } else {
      c_lable = as.character(col.vars)
      r_lable = as.character(row.vars)
      }
  
  ## Plotting after rearranging the row and columns
  if (nrow(x_arrage3)>ncol(x_arrage3)) {
    par(pin = c(3,5))} else {
      par(pin = c(5,3))}
  
  if (is.null(color)){
    col_palette = brewer.pal(5,"Blues")
  } else {
    col_palette = color 
  }
  if ( verbose == TRUE){
    png(filename = "./arranged_missing.png")
    image(1:ncol(x_mat), 1:nrow(x_mat), t(apply(x_mat, 2, rev)),
          ylab=cn[index_nominal[2]], xlab=cn[index_nominal[1]], col = col_palette)
    dev.off()} else{
    image(1:ncol(x_mat), 1:nrow(x_mat), t(apply(x_mat, 2, rev)),
          ylab=cn[index_nominal[2]], xlab=cn[index_nominal[1]], col = col_palette)
  }
  
  # replacing back the original row and columns
  col_name_rep <- as.data.frame(colnames(x_arrage3))  
  colnames(col_name_rep) <- "v1"
  col_name_rep <- left_join(col_name_rep, row.names.Htable, by=c("v1"="R.Names.mod"))
  
  row_name_rep <- as.data.frame(rownames(x_arrage3))  
  colnames(row_name_rep) <- "v1"
  row_name_rep <- left_join(row_name_rep, col.names.Htable, by=c("v1"="C.Names.mod"))
  
  row.names(x_arrage3) <- row_name_rep$C.names
  colnames(x_arrage3) <- col_name_rep$R.names
  
  
  ### selecting the largest square
  i <- 1
  j <- 1
  while(sum(x_arrage3[1:i,1:j], na.rm = TRUE) == i*j){
    i <- i+1
    j <- j+1
    max_sq_nomiss <- x_arrage3[1:i-1,1:j-1]
    #print(max_sq_nomiss)
  }
  
  
  ### row, col, and element wise maximum non-missing
  
  ##finding the first appearance of the missing value
  #in col
  col.na <- apply(is.na(x_arrage3), 2, which)
  col.1st.na <- unlist(lapply(col.na, min))
  
  #in row
  row.na <- apply(is.na(x_arrage3), 1, which)
  row.1st.na <- unlist(lapply(row.na, min))
  
  ####row-wise biggest module:
  ind <- max(which(col.1st.na==max(col.1st.na[2:length(col.1st.na)])))
  row_ind <- col.1st.na[ind]
  
  row_max_rec <- x_arrage3[1:row_ind-1, 1:ind]
  
  
  ####col-wise biggest module:
  ind <- max(which(row.1st.na==max(row.1st.na[2:length(row.1st.na)])))
  col_ind <- row.1st.na[ind]
  
  col_max_rec <- x_arrage3[1:ind, 1:col_ind-1]
  
  
  # element wise selection
  max_row <- col.1st.na*seq_along(col.1st.na)
  max_col <- row.1st.na*seq_along(row.1st.na)
  
  if (max(max_row)>max(max_col)){
    ind <- which.max(max_row)
    max_rec <- x_arrage3[1:col.1st.na[ind]-1, 1:ind]
  } else {
    ind <- which.max(max_col)
    max_rec <- x_arrage3[1:ind, 1:row.1st.na[ind]-1]
  }
  
  req_shape <- c("Square", "Rectangular_row", "Rectangular_col", "Rectangular_element_max", "All")
  if (!(shape %in% req_shape)) {
    stop("Please specify the desire shape of the output")
  } else if ( shape == "Square") {
    return(max_sq_nomiss)
  } else if ( shape == "Rectangular_col"){
    return(col_max_rec)
  } else if ( shape == "Rectangular_row"){
    return(row_max_rec)
  } else if ( shape == "Rectangular_element_max"){
    return(max_rec)
  } else if ( shape == "All"){
    return(list(max_sq_nomiss, col_max_rec, row_max_rec, max_rec))
  }

}






