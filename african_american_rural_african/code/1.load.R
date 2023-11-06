# Load libraries
library(microbiome)
library(tidyverse)
library(ggpubr)
library("SmartEDA")
library(gridExtra)
library(vegan)
library(rstatix)
library(ANCOMBC)

# Functions
get_upper_tri = function(cormat){
  cormat[lower.tri(cormat)] = NA
  diag(cormat) = NA
  return(cormat)
}
