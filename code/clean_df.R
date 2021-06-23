# This function makes comparable the cell lines of two dfs (columns) 

clean_dfs <- function(df1, df2){
  cl1 <- colnames(df1)
  cl2 <- colnames(df2)
  if(length(cl1)!=length(cl2)){
    cl_common <- intersect(cl1, cl2)
  }else{
    cl_common <- cl1
  }
  df1 <- df1[ , cl_common]
  df2 <- df2[ , cl_common]
  list(df1,df2)
  return(list(df1,df2))
}








