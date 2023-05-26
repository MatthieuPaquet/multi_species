#this script is meant to draw a figure comparing precision of density dependencies with and without random time variation
#the output imported were obtained by running the script "diagnostics_noddinter_centered.R"
library(viridis)
setwd("/home/matpaquet")
load("/home/matpaquet/list_precis_nostoch_nodd.R")
load("/home/matpaquet/list_precis_stoch_nodd.R")
