#!/usr/bin/Rscript --no-restore --no-environ --no-save

rm(list = ls())

###############################################################################
#
# author: Christoph Kiefer
# email: christophak@bmb.sdu.dk
#
################################################################################

source("packrat/init.R")
packrat::restore()
