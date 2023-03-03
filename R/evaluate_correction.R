#' Evaluate expression level correction
#'
#' @param before expression table obtained before correction
#' @param after expression table obtained after correction
#' @param data_cols columns containing expression data
#'
#' @return
#' Median variance of expression levels for a single sample, calculated within exons and median variance between samples.
#' Both statistics are calculated independently for dataset before and after correction.
#' 
#' @export
#' 
#' @examples
#' evaluate_correction(exprs_original,exprs_correctted,data_cols = SampleIDs)
evaluate_correction = function(before,after,data_cols) {
  
  median_var_exon = function(data) {
    data %>%
      as_tibble() %>%
      pivot_longer(data_cols) %>%
      dplyr::group_by(Exon.exon_id,name) %>%
      dplyr::summarize(var = sqrt(var(value))) %>%
      ungroup() %>%
      dplyr::summarize(median_var = median(var,na.rm=T))
  }
  
  median_var_rep = function(data) {
    data %>%
      as_tibble() %>%
      pivot_longer(data_cols) %>%
      mutate(cells = strsplit2(name,'_')[,1]) %>%
      group_by(ProbeID,cells) %>%
      dplyr::summarize(var = sqrt(var(value))) %>%
      ungroup() %>%
      dplyr::summarize(median_var = median(var,na.rm=T))
  }
  
  #median variance between probes targetting specific exon
  exon_var_before = before %>% median_var_exon()
  exon_var_after  = after  %>% median_var_exon()
  
  rep_var_before = before %>% median_var_rep()
  rep_var_after  = after  %>% median_var_rep()
  
  return(data.frame(dataset = c('before','after'),
                    exon_var = as.numeric(c(exon_var_before,exon_var_after)),
                    rep_var = as.numeric(c(rep_var_before,rep_var_after))))
}