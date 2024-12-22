

# Descriptive Analyses ----------------------------------------------------
# TPC: MET (random Ibloc) -------------------------------------------------------------------

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

# Use pacman to load and install packages if needed
pacman::p_load(char = packages)

# Load functions 
#source("ic_reml_asr.R")

# Phenotypic data ---------------------------------------------------------

pheno <- readxl::read_excel("data_raw/TPC-MS_info_v2_unesp.xlsx", 
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
pheno$CHECK = as.factor(ifelse(pheno$Testemunha == "SIM", 1, 0))
str(pheno)

pheno$VOL_d3  = pheno$Volume * 1000


# Replace the Trials names

replacements =c("T_1" = "T_01", "T_2"= "T_02", "T_3" = "T_03", 
                "T_4" = "T_04", "T_5" = "T_05", "T_6" = "T_06", 
                "T_7" = "T_07", "T_8" = "T_08",  "T_9" = "T_09", 
                "T_10" ="T_10", "T_11" ="T_11", "T_12" ="T_12" )

# Use ifelse to replace strings
pheno$TRIAL  <- as.factor(ifelse(pheno$TRIAL  %in% names(replacements),
                                 replacements[as.character(pheno$TRIAL)],
                                 pheno$TRIAL ))
str(pheno)

#write_csv(pheno, "data/new_pheno.csv", col_names = TRUE)


# Take only the last measurement ------------------------------------------

levels(pheno$YEAR)
pheno_3y = droplevels(subset(pheno, YEAR ==3)) 

pheno_3y <- pheno_3y[order(pheno_3y$SITE, pheno_3y$TRIAL, pheno_3y$REP, pheno_3y$IBLC),]
str(pheno_3y)


# Descriptive analysis ----------------------------------------------------

# pheno = setDT(pheno)
# ans <- pheno[TRIAL == "T_01", .(menor = min(Parcela, na.rm = T), maior = max(Parcela, na.rm = T)), by = .(Bloco)]

# Função
source("utilities/reps_joint_ADGBLUP.R") # nossa funcao para considerar tudo junto


# Vol (d3) ------------------------------------------------------------------

G.inv = readRDS("data/Ga_inv.rds")
Gd.inv = readRDS("data/D_inv.rds")

pheno_3y$GENA = as.factor(pheno_3y$GEN)
pheno_3y$GEND = as.factor(pheno_3y$GEN)

# Combinações -------------------------------------------------------------

source("utilities/reps_joint_ADGBLUP.R")

# 1 REP -------------------------------------------------------------------

one_rep = reps_joint_ADGBLUP(df= pheno_3y, form_fixed = VOL_d3 ~ SITE + SITE:TRIAL + SITE:TRIAL:REP,
                           form_random = ~  vm(GENA, G.inv) +      #additive
                                            vm(GEND, Gd.inv) +               #dominance
                                            SITE:TRIAL:REP:IBLC,             #IBLC
                           form_residual = ~dsum(~units|SITE),
                           comb_size = 1,  
                           trial = "TRIAL",
                           rep= "REP",
                           ibloc = "IBLC",
                           site = "SITE", 
                           year = "YEAR",
                           pattern = c("vm\\(GENA", "vm\\(GEND") )



#Heritability 

# Use lapply to apply the condition to each third-level list

h2a <- lapply(one_rep[["VarComp"]], function(age_list) {
  lapply(age_list, function(info_list) {
    # Access the second element in the info_list and calculate its mean
    mean(as.numeric(info_list[[2]]))
  })
})

h2a = data.frame(unlist(h2a))
h2a$Num_REP= rep("1_REP", nrow(h2a))

d2 <- lapply(one_rep[["VarComp"]], function(age_list) {
  lapply(age_list, function(info_list) {
    # Access the second element in the info_list and calculate its mean
    mean(as.numeric(info_list[[3]]))
  })
})

d2 = data.frame(unlist(d2))
d2$Num_REP= rep("1_REP", nrow(d2))


h2g <- lapply(one_rep[["VarComp"]], function(age_list) {
  lapply(age_list, function(info_list) {
    # Access the second element in the info_list and calculate its mean
    mean(as.numeric(info_list[[4]]))
  })
})

h2g = data.frame(unlist(h2g))
h2g$Num_REP= rep("1_REP", nrow(h2g))

#Accuracy PEV Adiitive effect

acc_PEV_a <- lapply(one_rep[["VarComp"]], function(age_list) {
  lapply(age_list, function(info_list) {
    # Access the second element in the info_list and calculate its mean
    mean(as.numeric(info_list[[5]]))
  })
})

acc_PEV_a = data.frame(unlist(acc_PEV_a))
acc_PEV_a$Num_REP= rep("1_REP", nrow(acc_PEV_a))


#Accuracy PEV Dominance effect

acc_PEV_d <- lapply(one_rep[["VarComp"]], function(age_list) {
  lapply(age_list, function(info_list) {
    # Access the second element in the info_list and calculate its mean
    mean(as.numeric(info_list[[6]]))
  })
})

acc_PEV_d = data.frame(unlist(acc_PEV_d))
acc_PEV_d$Num_REP= rep("1_REP", nrow(acc_PEV_d))


#Results

res1 = data.frame(cbind(h2a, d2, h2g, acc_PEV_a, acc_PEV_d))

res1 = res1 %>% dplyr::select(c(1,3,5,7,9,10)) %>% 
  dplyr::rename(h2a = unlist.h2a., d2= unlist.d2.,
                h2g = unlist.h2g., Acc_a = unlist.acc_PEV_a.,
                Acc_d = unlist.acc_PEV_d., Num_REP = Num_REP.4)
res1



# 2 REPS ------------------------------------------------------------------

two_rep = reps_joint_ADGBLUP(df= pheno_3y, form_fixed = VOL_d3 ~ SITE + SITE:TRIAL + SITE:TRIAL:REP,
                             form_random = ~  vm(GENA, G.inv) +      #additive
                                              vm(GEND, Gd.inv) +               #dominance
                                              SITE:TRIAL:REP:IBLC,             #IBLC
                             form_residual = ~dsum(~units|SITE),
                             comb_size = 2,  
                             trial = "TRIAL",
                             rep= "REP",
                             ibloc = "IBLC",
                             site = "SITE", 
                             year = "YEAR",
                             pattern = c("vm\\(GENA", "vm\\(GEND") )



#Heritability 

# Use lapply to apply the condition to each third-level list

h2a <- lapply(two_rep[["VarComp"]], function(age_list) {
  lapply(age_list, function(info_list) {
    # Access the second element in the info_list and calculate its mean
    mean(as.numeric(info_list[[2]]))
  })
})

h2a = data.frame(unlist(h2a))
h2a$Num_REP= rep("2_REP", nrow(h2a))

d2 <- lapply(two_rep[["VarComp"]], function(age_list) {
  lapply(age_list, function(info_list) {
    # Access the second element in the info_list and calculate its mean
    mean(as.numeric(info_list[[3]]))
  })
})

d2 = data.frame(unlist(d2))
d2$Num_REP= rep("2_REP", nrow(d2))


h2g <- lapply(two_rep[["VarComp"]], function(age_list) {
  lapply(age_list, function(info_list) {
    # Access the second element in the info_list and calculate its mean
    mean(as.numeric(info_list[[4]]))
  })
})

h2g = data.frame(unlist(h2g))
h2g$Num_REP= rep("2_REP", nrow(h2g))

#Accuracy PEV Adiitive effect

acc_PEV_a <- lapply(two_rep[["VarComp"]], function(age_list) {
  lapply(age_list, function(info_list) {
    # Access the second element in the info_list and calculate its mean
    mean(as.numeric(info_list[[5]]))
  })
})

acc_PEV_a = data.frame(unlist(acc_PEV_a))
acc_PEV_a$Num_REP= rep("2_REP", nrow(acc_PEV_a))


#Accuracy PEV Dominance effect

acc_PEV_d <- lapply(two_rep[["VarComp"]], function(age_list) {
  lapply(age_list, function(info_list) {
    # Access the second element in the info_list and calculate its mean
    mean(as.numeric(info_list[[6]]))
  })
})

acc_PEV_d = data.frame(unlist(acc_PEV_d))
acc_PEV_d$Num_REP= rep("2_REP", nrow(acc_PEV_d))


#Results

res2 = data.frame(cbind(h2a, d2, h2g, acc_PEV_a, acc_PEV_d))

res2 = res2 %>% dplyr::select(c(1,3,5,7,9,10)) %>% 
  dplyr::rename(h2a = unlist.h2a., d2= unlist.d2.,
                h2g = unlist.h2g., Acc_a = unlist.acc_PEV_a.,
                Acc_d = unlist.acc_PEV_d., Num_REP = Num_REP.4)
res2


# 3 REPS ------------------------------------------------------------------

three_rep = reps_joint_ADGBLUP(df= pheno_3y, form_fixed = VOL_d3 ~ SITE + SITE:TRIAL + SITE:TRIAL:REP,
                             form_random = ~  vm(GENA, G.inv) +      #additive
                                              vm(GEND, Gd.inv) +               #dominance
                                              SITE:TRIAL:REP:IBLC,             #IBLC
                             form_residual = ~dsum(~units|SITE),
                             comb_size = 3,  
                             trial = "TRIAL",
                             rep= "REP",
                             ibloc = "IBLC",
                             site = "SITE", 
                             year = "YEAR",
                             pattern = c("vm\\(GENA", "vm\\(GEND") )



#Heritability 

# Use lapply to apply the condition to each third-level list

h2a <- lapply(three_rep[["VarComp"]], function(age_list) {
  lapply(age_list, function(info_list) {
    # Access the second element in the info_list and calculate its mean
    mean(as.numeric(info_list[[2]]))
  })
})

h2a = data.frame(unlist(h2a))
h2a$Num_REP= rep("3_REP", nrow(h2a))

d2 <- lapply(three_rep[["VarComp"]], function(age_list) {
  lapply(age_list, function(info_list) {
    # Access the second element in the info_list and calculate its mean
    mean(as.numeric(info_list[[3]]))
  })
})

d2 = data.frame(unlist(d2))
d2$Num_REP= rep("3_REP", nrow(d2))


h2g <- lapply(three_rep[["VarComp"]], function(age_list) {
  lapply(age_list, function(info_list) {
    # Access the second element in the info_list and calculate its mean
    mean(as.numeric(info_list[[4]]))
  })
})

h2g = data.frame(unlist(h2g))
h2g$Num_REP= rep("3_REP", nrow(h2g))

#Accuracy PEV Adiitive effect

acc_PEV_a <- lapply(three_rep[["VarComp"]], function(age_list) {
  lapply(age_list, function(info_list) {
    # Access the second element in the info_list and calculate its mean
    mean(as.numeric(info_list[[5]]))
  })
})

acc_PEV_a = data.frame(unlist(acc_PEV_a))
acc_PEV_a$Num_REP= rep("3_REP", nrow(acc_PEV_a))


#Accuracy PEV Dominance effect

acc_PEV_d <- lapply(three_rep[["VarComp"]], function(age_list) {
  lapply(age_list, function(info_list) {
    # Access the second element in the info_list and calculate its mean
    mean(as.numeric(info_list[[6]]))
  })
})

acc_PEV_d = data.frame(unlist(acc_PEV_d))
acc_PEV_d$Num_REP= rep("3_REP", nrow(acc_PEV_d))


#Results

res3 = data.frame(cbind(h2a, d2, h2g, acc_PEV_a, acc_PEV_d))

res3 = res3 %>% dplyr::select(c(1,3,5,7,9,10)) %>% 
  dplyr::rename(h2a = unlist.h2a., d2= unlist.d2.,
                h2g = unlist.h2g., Acc_a = unlist.acc_PEV_a.,
                Acc_d = unlist.acc_PEV_d., Num_REP = Num_REP.4)
res3



# 4 REPS ------------------------------------------------------------------

four_rep = reps_joint_ADGBLUP(df= pheno_3y, form_fixed = VOL_d3 ~ SITE + SITE:TRIAL + SITE:TRIAL:REP,
                             form_random = ~  vm(GENA, G.inv) +      #additive
                                              vm(GEND, Gd.inv) +               #dominance
                                              SITE:TRIAL:REP:IBLC,             #IBLC
                             form_residual = ~dsum(~units|SITE),
                             comb_size = 4,  
                             trial = "TRIAL",
                             rep= "REP",
                             ibloc = "IBLC",
                             site = "SITE", 
                             year = "YEAR",
                             pattern = c("vm\\(GENA", "vm\\(GEND") )

#Heritability 

# Use lapply to apply the condition to each third-level list

h2a <- lapply(four_rep[["VarComp"]], function(age_list) {
  lapply(age_list, function(info_list) {
    # Access the second element in the info_list and calculate its mean
    mean(as.numeric(info_list[[2]]))
  })
})

h2a = data.frame(unlist(h2a))
h2a$Num_REP= rep("4_REP", nrow(h2a))

d2 <- lapply(four_rep[["VarComp"]], function(age_list) {
  lapply(age_list, function(info_list) {
    # Access the second element in the info_list and calculate its mean
    mean(as.numeric(info_list[[3]]))
  })
})

d2 = data.frame(unlist(d2))
d2$Num_REP= rep("4_REP", nrow(d2))


h2g <- lapply(four_rep[["VarComp"]], function(age_list) {
  lapply(age_list, function(info_list) {
    # Access the second element in the info_list and calculate its mean
    mean(as.numeric(info_list[[4]]))
  })
})

h2g = data.frame(unlist(h2g))
h2g$Num_REP= rep("4_REP", nrow(h2g))

#Accuracy PEV Adiitive effect

acc_PEV_a <- lapply(four_rep[["VarComp"]], function(age_list) {
  lapply(age_list, function(info_list) {
    # Access the second element in the info_list and calculate its mean
    mean(as.numeric(info_list[[5]]))
  })
})

acc_PEV_a = data.frame(unlist(acc_PEV_a))
acc_PEV_a$Num_REP= rep("4_REP", nrow(acc_PEV_a))


#Accuracy PEV Dominance effect

acc_PEV_d <- lapply(four_rep[["VarComp"]], function(age_list) {
  lapply(age_list, function(info_list) {
    # Access the second element in the info_list and calculate its mean
    mean(as.numeric(info_list[[6]]))
  })
})

acc_PEV_d = data.frame(unlist(acc_PEV_d))
acc_PEV_d$Num_REP= rep("4_REP", nrow(acc_PEV_d))


#Results

res4 = data.frame(cbind(h2a, d2, h2g, acc_PEV_a, acc_PEV_d))

res4 = res4 %>% dplyr::select(c(1,3,5,7,9,10)) %>% 
  dplyr::rename(h2a = unlist.h2a., d2= unlist.d2.,
                h2g = unlist.h2g., Acc_a = unlist.acc_PEV_a.,
                Acc_d = unlist.acc_PEV_d., Num_REP = Num_REP.4)
res4



# 5 REPS ------------------------------------------------------------------
  
five_rep = reps_joint_ADGBLUP(df= pheno_3y, form_fixed = VOL_d3 ~ SITE + SITE:TRIAL + SITE:TRIAL:REP,
                              form_random = ~  vm(GENA, G.inv) +      #additive
                                vm(GEND, Gd.inv) +               #dominance
                                SITE:TRIAL:REP:IBLC,             #IBLC
                              form_residual = ~dsum(~units|SITE),
                              comb_size = 5,  
                              trial = "TRIAL",
                              rep= "REP",
                              ibloc = "IBLC",
                              site = "SITE", 
                              year = "YEAR",
                              pattern = c("vm\\(GENA", "vm\\(GEND") )

#Heritability 

# Use lapply to apply the condition to each third-level list

h2a <- lapply(five_rep[["VarComp"]], function(age_list) {
  lapply(age_list, function(info_list) {
    # Access the second element in the info_list and calculate its mean
    mean(as.numeric(info_list[[2]]))
  })
})

h2a = data.frame(unlist(h2a))
h2a$Num_REP= rep("5_REP", nrow(h2a))

d2 <- lapply(five_rep[["VarComp"]], function(age_list) {
  lapply(age_list, function(info_list) {
    # Access the second element in the info_list and calculate its mean
    mean(as.numeric(info_list[[3]]))
  })
})

d2 = data.frame(unlist(d2))
d2$Num_REP= rep("5_REP", nrow(d2))


h2g <- lapply(five_rep[["VarComp"]], function(age_list) {
  lapply(age_list, function(info_list) {
    # Access the second element in the info_list and calculate its mean
    mean(as.numeric(info_list[[4]]))
  })
})

h2g = data.frame(unlist(h2g))
h2g$Num_REP= rep("5_REP", nrow(h2g))

#Accuracy PEV Adiitive effect

acc_PEV_a <- lapply(five_rep[["VarComp"]], function(age_list) {
  lapply(age_list, function(info_list) {
    # Access the second element in the info_list and calculate its mean
    mean(as.numeric(info_list[[5]]))
  })
})

acc_PEV_a = data.frame(unlist(acc_PEV_a))
acc_PEV_a$Num_REP= rep("5_REP", nrow(acc_PEV_a))


#Accuracy PEV Dominance effect

acc_PEV_d <- lapply(five_rep[["VarComp"]], function(age_list) {
  lapply(age_list, function(info_list) {
    # Access the second element in the info_list and calculate its mean
    mean(as.numeric(info_list[[6]]))
  })
})

acc_PEV_d = data.frame(unlist(acc_PEV_d))
acc_PEV_d$Num_REP= rep("5_REP", nrow(acc_PEV_d))


#Results

res5 = data.frame(cbind(h2a, d2, h2g, acc_PEV_a, acc_PEV_d))

res5 = res5 %>% dplyr::select(c(1,3,5,7,9,10)) %>% 
  dplyr::rename(h2a = unlist.h2a., d2= unlist.d2.,
                h2g = unlist.h2g., Acc_a = unlist.acc_PEV_a.,
                Acc_d = unlist.acc_PEV_d., Num_REP = Num_REP.4)
res5

save.image("images/image_H2_AC_AD-GBLUP.RData")  
#load("image_H2_AC_AD-GBLUP.RData")  

   

# All combinations --------------------------------------------------------

all_h2_acc = data.frame(rbind(res1, res2, res3, res4, res5))

all_h2_acc= all_h2_acc[order(all_h2_acc$Num_REP), ]
head(all_h2_acc)


#dir.create("tables")
write.table("results/tab_new_AD_GBLUP_all_h2_acc.txt")

# Plots -------------------------------------------------------------------


#cbPalette <- c( "#115f9a", "#48b5c4", "#76c68f","#d0ee11")

cbPalette <- c( "#b3bfa1",   "#b3bfd1",  "#a4a2a8",  "#c86558",  "#991f17")

# Butterfly plot

library(PupillometryR)
library(ggthemes) #theme clean are needed
source("utilities/theme_jace.R")


p1 = ggplot(all_h2_acc) +
  aes(x = Num_REP, 
      y = h2a, 
      fill = Num_REP) + # Code changed here
  geom_flat_violin(trim = FALSE,
                   alpha = .5) +
  geom_boxplot(width = .25, 
               outlier.shape = NA,
               alpha = 1) +
  scale_fill_manual(values = cbPalette) +
  scale_color_manual(values = cbPalette)  +
  labs(x = "Number of Repetition",
       y = "Heritability",
       #title = "Raincloud plot of heritability from m combinations of r repetitions in Matão"
  )+
  theme_jace +
  theme(legend.position = "bottom",
        legend.margin = margin(-5, 0, 0, 0));p1


p2 = ggplot(all_h2_acc) +
  aes(x = Num_REP, 
      y = d2, 
      fill = Num_REP) + # Code changed here
  geom_flat_violin(trim = FALSE,
                   alpha = .5) +
  geom_boxplot(width = .25, 
               outlier.shape = NA,
               alpha = 1) +
  scale_fill_manual(values = cbPalette) +
  scale_color_manual(values = cbPalette)  +
  labs(x = "Number of Repetition",
       y = "Accuracy"#,
       #title = "Raincloud plot of Accuracy from m combinations of r repetitions in Matão"
  )  +
  theme_jace +
  theme(legend.position = "bottom",
        legend.margin = margin(-5, 0, 0, 0)); p2


p3 = ggplot(all_h2_acc) +
  aes(x = Num_REP, 
      y = h2g, 
      fill = Num_REP) + # Code changed here
  geom_flat_violin(trim = FALSE,
                   alpha = .5) +
  geom_boxplot(width = .25, 
               outlier.shape = NA,
               alpha = 1) +
  scale_fill_manual(values = cbPalette) +
  scale_color_manual(values = cbPalette)  +
  labs(x = "Number of Repetition",
       y = "Heritability",
       #title = "Raincloud plot of heritability from m combinations of r repetitions in Matão"
  )+
  theme_jace +
  theme(legend.position = "bottom",
        legend.margin = margin(-5, 0, 0, 0));p3


p4 = ggplot(all_h2_acc) +
  aes(x = Num_REP, 
      y = Acc_a, 
      fill = Num_REP) + # Code changed here
  geom_flat_violin(trim = FALSE,
                   alpha = .5) +
  geom_boxplot(width = .25, 
               outlier.shape = NA,
               alpha = 1) +
  scale_fill_manual(values = cbPalette) +
  scale_color_manual(values = cbPalette)  +
  labs(x = "Number of Repetition",
       y = "Accuracy"#,
       #title = "Raincloud plot of Accuracy from m combinations of r repetitions in Matão"
  )  +
  theme_jace +
  theme(legend.position = "bottom",
        legend.margin = margin(-5, 0, 0, 0)); p4

p5 = ggplot(all_h2_acc) +
  aes(x = Num_REP, 
      y = Acc_d, 
      fill = Num_REP) + # Code changed here
  geom_flat_violin(trim = FALSE,
                   alpha = .5) +
  geom_boxplot(width = .25, 
               outlier.shape = NA,
               alpha = 1) +
  scale_fill_manual(values = cbPalette) +
  scale_color_manual(values = cbPalette)  +
  labs(x = "Number of Repetition",
       y = "Accuracy"#,
       #title = "Raincloud plot of Accuracy from m combinations of r repetitions in Matão"
  )  +
  theme_jace +
  theme(legend.position = "bottom",
        legend.margin = margin(-5, 0, 0, 0)); p5


library(gridExtra)
svg("results/new_fig_reps_AD_GBLUP.svg", width = 7, height = 7)
grid.arrange(p1, p2, nrow=2, ncol=1)
dev.off()


  

# Average value
mean_h2a = all_h2_acc %>%
  group_by(Num_REP) %>%
  summarise_at(vars(h2a), list(name = ~mean(., na.rm = TRUE)))

mean_d2 = all_h2_acc %>%
  group_by(Num_REP) %>%
  summarise_at(vars(d2), list(name = ~mean(., na.rm = TRUE)))

mean_h2g = all_h2_acc %>%
  group_by(Num_REP) %>%
  summarise_at(vars(h2g), list(name = ~mean(., na.rm = TRUE)))

mean_acc_a = all_h2_acc %>%
  group_by(Num_REP) %>%
  summarise_at(vars(Acc_a), list(name = ~mean(., na.rm = TRUE)))

mean_acc_d = all_h2_acc %>%
  group_by(Num_REP) %>%
  summarise_at(vars(Acc_d), list(name = ~mean(., na.rm = TRUE)))

mean_acc_d = all_h2_acc %>%
  group_by(Num_REP) %>%
  summarise_at(vars(Acc_d), list(name = ~mean(., na.rm = TRUE)))

media = cbind(mean_h2, mean_acc_PEV, mean_acc_H2)[ ,c(1,2,4,6)]
colnames(media)= c("Num_REP", "H2", "ACC_PEV", "ACC_H2")

write.table(media, "results/new_A_BLUP_means_vol_d3.txt")


# Save Image --------------------------------------------------------------

#save.image("images/image_H2_AC_AD-GBLUP.RData")  
load("C:/Users/LENOVO/OneDrive/POS-DOC/TPC/image_H2_AC_AD-GBLUP.RData") 


# Heritabilities and SE ---------------------------------------------------


h1 = asreml::vpredict(one_rep[["Summary"]][["3"]][[1]][[1]], 
                      h2a~((V2+V3)/(V1+V2+V3+((V4+V5)/2)))); h1

h2 = asreml::vpredict(two_rep[["Summary"]][["3"]][[1]][[1]], 
                      h2a~((V2+V3)/(V1+V2+V3+((V4+V5)/2)))); h2

h3 = asreml::vpredict(three_rep[["Summary"]][["3"]][[1]][[1]], 
                      h2a~((V2+V3)/(V1+V2+V3+((V4+V5)/2)))); h3

h4 = asreml::vpredict(four_rep[["Summary"]][["3"]][[1]][[1]], 
                      h2a~((V2+V3)/(V1+V2+V3+((V4+V5)/2)))); h4

h5 = asreml::vpredict(five_rep[["Summary"]][["3"]][[1]][[1]], 
                      h2a~((V2+V3)/(V1+V2+V3+((V4+V5)/2)))); h5

hs= rbind(h1, h2, h3, h4, h5)

write.csv(hs, "results/h2_SE_combinations.csv")

# Variance Components - One Rep
v1.1 = summary(one_rep[["Summary"]][["3"]][[1]][[1]])$varcomp[ ,c(1:2)]
v1.2 = summary(one_rep[["Summary"]][["3"]][[2]][[1]])$varcomp[ ,c(1:2)]
v1.3 = summary(one_rep[["Summary"]][["3"]][[3]][[1]])$varcomp[ ,c(1:2)]
v1.4 = summary(one_rep[["Summary"]][["3"]][[4]][[1]])$varcomp[ ,c(1:2)]
v1.5 = summary(one_rep[["Summary"]][["3"]][[5]][[1]])$varcomp[ ,c(1:2)]

v1 = cbind(v1.1, v1.2, v1.3, v1.4, v1.5 )

v1$mean_component = rowMeans(v1[, c(1, 3, 5, 7, 9)])
v1$mean_std_error = rowMeans(v1[, c(2, 4, 6, 8, 10)])


# Variance Components - Two Rep
v2.1 = summary(two_rep[["Summary"]][["3"]][[1]][[1]])$varcomp[ ,c(1:2)]
v2.2 = summary(two_rep[["Summary"]][["3"]][[2]][[1]])$varcomp[ ,c(1:2)]
v2.3 = summary(two_rep[["Summary"]][["3"]][[3]][[1]])$varcomp[ ,c(1:2)]
v2.4 = summary(two_rep[["Summary"]][["3"]][[4]][[1]])$varcomp[ ,c(1:2)]
v2.5 = summary(two_rep[["Summary"]][["3"]][[5]][[1]])$varcomp[ ,c(1:2)]
v2.6 = summary(two_rep[["Summary"]][["3"]][[6]][[1]])$varcomp[ ,c(1:2)]
v2.7 = summary(two_rep[["Summary"]][["3"]][[7]][[1]])$varcomp[ ,c(1:2)]
v2.8 = summary(two_rep[["Summary"]][["3"]][[8]][[1]])$varcomp[ ,c(1:2)]
v2.9 = summary(two_rep[["Summary"]][["3"]][[9]][[1]])$varcomp[ ,c(1:2)]
v2.10 = summary(two_rep[["Summary"]][["3"]][[10]][[1]])$varcomp[ ,c(1:2)]

v2 = cbind(v2.1, v2.2, v2.3, v2.4, v2.5, v2.6, v2.7, v2.8, v2.9, v2.10)

v2$mean_component = rowMeans(v2[, c(1, 3, 5, 7, 9, 11, 13, 15, 17, 19)])
v2$mean_std_error = rowMeans(v2[, c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20)])


# Variance Components - Two Rep
v3.1 = summary(three_rep[["Summary"]][["3"]][[1]][[1]])$varcomp[ ,c(1:2)]
v3.2 = summary(three_rep[["Summary"]][["3"]][[2]][[1]])$varcomp[ ,c(1:2)]
v3.3 = summary(three_rep[["Summary"]][["3"]][[3]][[1]])$varcomp[ ,c(1:2)]
v3.4 = summary(three_rep[["Summary"]][["3"]][[4]][[1]])$varcomp[ ,c(1:2)]
v3.5 = summary(three_rep[["Summary"]][["3"]][[5]][[1]])$varcomp[ ,c(1:2)]
v3.6 = summary(three_rep[["Summary"]][["3"]][[6]][[1]])$varcomp[ ,c(1:2)]
v3.7 = summary(three_rep[["Summary"]][["3"]][[7]][[1]])$varcomp[ ,c(1:2)]
v3.8 = summary(three_rep[["Summary"]][["3"]][[8]][[1]])$varcomp[ ,c(1:2)]
v3.9 = summary(three_rep[["Summary"]][["3"]][[9]][[1]])$varcomp[ ,c(1:2)]
v3.10 = summary(three_rep[["Summary"]][["3"]][[10]][[1]])$varcomp[ ,c(1:2)]

v3 = cbind(v3.1, v3.2, v3.3, v3.4, v3.5, v3.6, v3.7, v3.8, v3.9, v3.10)

v3$mean_component = rowMeans(v3[, c(1, 3, 5, 7, 9, 11, 13, 15, 17, 19)])
v3$mean_std_error = rowMeans(v3[, c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20)])



# Variance Components - Four Rep
v4.1 = summary(four_rep[["Summary"]][["3"]][[1]][[1]])$varcomp[ ,c(1:2)]
v4.2 = summary(four_rep[["Summary"]][["3"]][[2]][[1]])$varcomp[ ,c(1:2)]
v4.3 = summary(four_rep[["Summary"]][["3"]][[3]][[1]])$varcomp[ ,c(1:2)]
v4.4 = summary(four_rep[["Summary"]][["3"]][[4]][[1]])$varcomp[ ,c(1:2)]
v4.5 = summary(four_rep[["Summary"]][["3"]][[5]][[1]])$varcomp[ ,c(1:2)]

v4 = cbind(v4.1, v4.2, v4.3, v4.4, v4.5 )

v4$mean_component = rowMeans(v4[, c(1, 3, 5, 7, 9)])
v4$mean_std_error = rowMeans(v4[, c(2, 4, 6, 8, 10)])


# Variance Components - Five Rep
v5 = summary(five_rep[["Summary"]][["3"]][[1]][[1]])$varcomp[ ,c(1:2)]
v5

comp = data.frame(cbind(v1$mean_component, v1$mean_std_error, 
                        v2$mean_component, v2$mean_std_error,
                        v3$mean_component, v3$mean_std_error,
                        v4$mean_component, v4$mean_std_error,
                        v5$component, v5$std.error))

rownames(comp)= c("Ibloco", "a", "d", "res_1", "res2")


write.csv(vc_single, "results/varcomp_single_model_combinations.csv")



