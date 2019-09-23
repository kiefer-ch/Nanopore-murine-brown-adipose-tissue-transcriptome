#!/usr/bin/Rscript --no-restore --no-environ --no-save

rm(list = ls())

###############################################################################
#
# author: Christoph Kiefer
# email: christophak@bmb.sdu.dk
#
################################################################################

# Set repositories to include bioconductor. This allows packrat to use
# install.packages() to install bioconductor packages.
setRepositories(ind = c(1, 2, 3, 4))

source("packrat/init.R")
packrat::restore()
