#Two Stage


rm(list=ls())
gc()


# Library -----------------------------------------------------------------

# Install pacman if not already installed
# if (!require("pacman")) {
#   install.packages("pacman")
#   library(pacman)
# }

# List of packages to load
packages <- c("tidyverse", "asreml", "lme4", "car", "snpReady", "AGHmatrix",
              "ASRgenomics", "redxl", "data.table", "gridExtra", "PupillometryR",
              "readxl")

library(data.table)
library(tidyr)
library(dplyr)
library(readr)


#install.packages("StageWise")

# Use pacman to load and install packages if needed
pacman::p_load(char = packages)

# Load functions 
#source("ic_reml_asr.R")

# Phenotypic data ---------------------------------------------------------

pheno <- read_excel("C:/Users/LENOVO/OneDrive/POS-DOC/TPC/data/new_data/dados_TPC.xlsx", 
                    sheet = "phenotype")
head(pheno)
str(pheno)
dim(pheno)


pheno$REP   = as.factor(pheno$Rep)
pheno$IBLC  = as.factor(pheno$Bloco)
pheno$GEN   = as.factor(pheno$Subject_id)
pheno$CLONE = as.factor(pheno$Tratamento)
pheno$SITE  = as.factor(pheno$Local)
pheno$TRIAL = as.factor(pheno$EXP)
pheno$YEAR  = as.factor(pheno$Ano)
pheno$Linha  = as.factor(pheno$Linha)
pheno$CHECK = as.factor(ifelse(pheno$Testemunha == "SIM", 1, 0))
str(pheno)

pheno$VOL_d3  = as.numeric(pheno$Volume * 1000)  

pheno$Parcela  = as.factor(pheno$Parcela)
pheno$ROW  = as.factor(pheno$ROW)
pheno$COL  = as.factor(pheno$COL)

pheno = pheno[order(pheno$COL, pheno$ROW), ]

# Replace the Trials names

replacements =c("T_1" = "1", "T_2"= "2", "T_3" = "3", 
                "T_4" = "4", "T_5" = "5", "T_6" = "6", 
                "T_7" = "7", "T_8" = "8",  "T_9" = "9", 
                "T_10" ="10", "T_11" ="11", "T_12" ="12" )

# Use ifelse to replace strings
pheno$TRIAL  <- as.factor(ifelse(pheno$TRIAL  %in% names(replacements),
                                 replacements[as.factor(pheno$TRIAL)],
                                 pheno$TRIAL ))

replacements2 =c("Matão" = "1", "Palmito"= "2")

# Use ifelse to replace strings
pheno$SITE  <- as.factor(ifelse(pheno$SITE  %in% names(replacements2),
                                replacements2[as.factor(pheno$SITE)],
                                pheno$SITE ))
str(pheno)

pheno = pheno[ , -c(5:6)]

pheno$group = as.factor(paste(pheno$TRIAL, pheno$SITE, sep= "_"))
levels(pheno$group)


write_csv(pheno, "C:/Users/LENOVO/OneDrive/POS-DOC/TPC/data/new_pheno.csv", col_names = TRUE)


# Take only the last measurement ------------------------------------------

levels(pheno$YEAR)
pheno_3y = droplevels(subset(pheno, YEAR ==3)) 


str(pheno_3y)

# Standardization of the data  ------------------------------------------

###  y_std =(y-siteMean_y)/siteSD_y   ####

pheno_3y=  pheno_3y %>%
  group_by(SITE) %>%
  mutate(Y_std = scale(VOL_d3))

str(pheno_3y)

write.csv(pheno_3y, "data/pheno_3y.csv")



# Matao ---------------------------------------------------------------------

levels(pheno_3y$SITE)
pheno_mat = droplevels(subset(pheno_3y, SITE =="1")) 


source("utilities/spats_fill.R") # Home-made function to fill spatial gaps


#  Fill the spatial gaps --------------------------------------------------


str(pheno_mat)

# Get the levels of the TRIAL column
trial_levels <- levels(pheno_mat$TRIAL)
trial_levels

mat_list <- vector("list", length(trial_levels))     # Initialize mat_list with the number of TRIAL levels
names(mat_list) <- trial_levels                      # Name the list by TRIAL levels

for (i in seq_along(trial_levels)) {
  
  temp_pheno <- droplevels(subset(pheno_mat, TRIAL == trial_levels[i]))
  
  # Remove duplicates based on the "Parcela" column
  temp_pheno <- temp_pheno %>%
    distinct(Parcela, .keep_all = TRUE)
  
  # Fill the matrix using spats_fill
  mat_list[[i]] <- spats_fill(temp_pheno, "COL", "ROW", "GENx", "TRIAL") # Home-made function
}


saveRDS(mat_list, "data/mat_list.RDS")


# Palmito ---------------------------------------------------------------------

levels(pheno_3y$SITE)
pheno_pal = droplevels(subset(pheno_3y, SITE =="2")) 

#  Fill the spatial gaps --------------------------------------------------

str(pheno_pal)

# Get the levels of the TRIAL column

trial_levels2 <- levels(pheno_pal$TRIAL)
trial_levels2

pal_list <- vector("list", length(trial_levels2))     # Initialize pal_list with the number of TRIAL levels
names(pal_list) <- trial_levels2                      # Name the list by TRIAL levels

for (i in seq_along(trial_levels2)) {
  
  temp_pheno <- droplevels(subset(pheno_pal, TRIAL == trial_levels2[i]))
  
  # Remove duplicates based on the "Parcela" column
  temp_pheno <- temp_pheno %>%
    distinct(Parcela, .keep_all = TRUE)
  
  # Fill the matrix using spats_fill
  pal_list[[i]] <- spats_fill(temp_pheno, "COL", "ROW", "GENx", "TRIAL") # Home-made function
}

saveRDS(pal_list, "data/pal_list.RDS")



