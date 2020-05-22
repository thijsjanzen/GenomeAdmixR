library("covr")


setwd("/Users/janzen/MEGasync/GitHub/GenomeAdmixR/")
library(lintr)
lintr::lint_package()

covr::codecov(type = "tests", token = "d40809f7-5a82-4f09-a357-bbdbf618ea62", quiet = FALSE)

library(goodpractice)
goodpractice::gp(quiet = FALSE)
