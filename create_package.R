library(Rcpp)
library(RcppArmadillo)
Rcpp.package.skeleton("CompRiskRE", path = "../package/",
                      code_files = c("source_code/CompRiskRE_CP.R",
                                     "source_code/CompRiskRE_FT.R",
                                     "source_code/CompRiskRE_XGBoost_CP.R",
                                     "source_code/CompRiskRE_XGBoost_FT.R",
                                     "source_code/sim_CP.R",
                                     "source_code/utils.R"),
                      cpp_files = c("source_code/utils.cpp",
                                    "source_code/Makevars",
                                    "source_code/Makevars.win"),
                      example_code = FALSE)


## After run "Rcpp.package.skeleton" line, we need to manually add ", RcppArmadillo" in "DESCRIPTION" file,
## and replace "Rcpp.h" by "RcppArmadillo.h" in "src/RcppExports.cpp" file

install.packages("/pkg-name", repos = NULL, type = "source")
library(TmpLasso)