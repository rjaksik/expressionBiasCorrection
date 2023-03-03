#' Rescale microarray values using a factor variable
#'
#' @param expr_matrix expression matrix 
#' @param factor_values factor values used for expression level scaling (contunous or categorial variable)
#'
#' @return
#' Corrected expression level values
#' 
#' @export
#'
#' @examples
#' rescale_values_factor_simple(expr_matrix,factor_values)
rescale_values_factor_simple = function(expr_matrix,factor_values) {
  ModelsDF = data.frame()
  for (i in 1:ncol(expr_matrix)) {
    input_df = data.frame(exprs = expr_matrix[,i],Feature=factor_values)
    loess_model = loess(log2(exprs)~Feature, input_df)
    ModelsDF = rbind(ModelsDF,unique(data.frame(sample_id = colnames(expr_matrix)[i], x=loess_model$x[,1],fitted=loess_model$fitted)))
  }
  
  ModelsDF2 = ModelsDF %>%
    tibble() %>%
    group_by(x) %>%
    dplyr::mutate(scaling=mean(fitted)/fitted) %>%
    dplyr::mutate(scaled=fitted*scaling) %>%
    select(x,sample_id,scaling)
  
  #calculate scaled values
  expr_matrix_scaled = expr_matrix %>%
    as_tibble() %>%
    mutate(id=paste0('N',1:nrow(expr_matrix)),factor_values = factor_values) %>%
    pivot_longer(colnames(expr_matrix)) %>%
    left_join(ModelsDF2,by=c('factor_values'='x','name'='sample_id')) %>%
    mutate(scaled_value=ifelse(is.na(scaling),value,value*scaling)) %>%
    select(id,name,scaled_value) %>%
    pivot_wider(names_from =name,values_from =scaled_value) %>%
    select(-id) %>%
    as.matrix()
  
  return(expr_matrix_scaled)
}
