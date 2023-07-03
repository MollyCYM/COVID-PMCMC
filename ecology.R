install.packages("devtools")
library(devtools)
install_github("sbfnk/rbi",ref="master")
install_github("sbfnk/rbi.helpers",ref="master")

rm(list = ls(all.names=TRUE))
unlink(".RData")

library('rbi')
try(detach(package:rbi, unload = TRUE), silent = TRUE)
library(rbi, quietly = TRUE)

library('rbi.helpers')

library('ggplot2', quietly = TRUE)
library('gridExtra', quietly = TRUE)
endTime <- 50

PP <- bi_model("ecology.bi")
synthetic_dataset_PP <- rbi::bi_generate_dataset(endtime=endTime,
                                            model=PP,
                                            seed="42",
                                            verbose=TRUE,
                                            add_options = list(
                                              noutputs=500),output_every = 1)

rdata_PP <- bi_read(synthetic_dataset_PP)