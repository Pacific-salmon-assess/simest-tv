#!/usr/local/bin/Rscript --slave
args <- commandArgs(trailingOnly=TRUE)

packages <- c("ggplot2", "remotes", "dplyr", "here", "rstan", "TMB")
install.packages(setdiff(packages, rownames(installed.packages()))) 

#install samsim 
remotes::install_github("Pacific-salmon-assess/samSim", ref="timevar", force=TRUE, upgrade="always" )

#install samest
remotes::install_git('https://github.com/Pacific-salmon-assess/samEst', force=TRUE)


library(samEst)
library(samSim)
library(ggplot2)
library(dplyr)
library(here)
library(rstan)
source("code/utils.R")

print(args)
