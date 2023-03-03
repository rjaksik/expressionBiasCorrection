
<!-- README.md is generated from README.Rmd. Please edit that file -->

# microarrayBiasCorrection

<!-- badges: start -->
<!-- badges: end -->

Functions used to correct expression level data using factor values.

## Instalation:

The package requires R 3.5.0 or later

``` r
install.packages("devtools")  
devtools::install_github("rjaksik/microarrayBiasCorrection")
```

## Example use:

``` r
library('microarrayBiasCorrection')

path='/RawData'
targets_file = "sample_desc.txt"

targets <- readTargets(targets_file,path=path)

exprs_original = standardize_data(targets,FALSE)
exprs_correctted = standardize_data(targets,TRUE)
SampleIDs = colnames(exprs_original)

evaluate_correction(exprs_original,
                    exprs_correctted,
                    data_cols = SampleIDs)




```
