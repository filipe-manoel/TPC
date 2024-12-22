#05_fit_stage

rm(list = ls())

library(tidyverse)
library(readr)
library(readxl)
library(asreml)
library(ASRgenomics)
library(AGHmatrix)

# Matão -------------------------------------------------------------------
mat_list = read_rds("data/mat_list.RDS")

# Fit Stage One Model -----------------------------------------------------

# Get the levels of the TRIAL column
trial_levels <- seq_along(mat_list)
trial_levels

asr_list_matao <- vector("list", length(trial_levels))     # Initialize mat_list with the number of TRIAL levels
names(asr_list_matao) <- trial_levels 

blues_list <- vector("list", length(trial_levels))     # Initialize mat_list with the number of TRIAL levels
names(blues_list) <- trial_levels

vcov_list <- vector("list", length(trial_levels))     # Initialize mat_list with the number of TRIAL levels
names(vcov_list) <- trial_levels

names1_list <- vector("list", length(trial_levels))     # Initialize mat_list with the number of TRIAL levels
names(names1_list) <- trial_levels


for (i in seq_along(mat_list)) {
  
  sdat = mat_list[[i]] 
  
  sdat = sdat[order(sdat$COL, sdat$ROW), ]
  
  m1 = asreml(
    fixed = VOL_d3 ~ GEN,
    random = ~  REP + REP:IBLC,
    residual = ~ ar1(COL):ar1(ROW),
    maxit = 50, #trace = F,
    workspace = "1gb",
    na.action = na.method(x = "include", y = "include"),
    data = sdat)
  
  m1 = update.asreml(m1)
  m1 = update.asreml(m1)
  m1 = update.asreml(m1)
  m1 = update.asreml(m1)
  m1 = update.asreml(m1)
  m1 = update.asreml(m1)
  
  asr_list_matao[[i]] = m1
  
  blues_temp = m1$coefficients$fixed
  
  temp1 = data.frame(blues_temp[grep("GEN_", rownames(blues_temp)),1])
  blues_list[[i]] = temp1
  attr(blues_list[[i]], "names") <- "effects"
  
  temp2 = predict.asreml(m1,classify = "GEN",vcov = T)$vcov
  rownames(temp2) = rownames(temp1)
  colnames(temp2) = rownames(temp1)
  vcov_list[[i]] = temp2
  
  names1_list[[i]] = rownames(temp1)
  
}

# Save the BLUES -----------------------------------------------------
# Save on the format used on stagewise package

b1_list <- vector("list", length(trial_levels))     # Initialize mat_list with the number of TRIAL levels
names(b1_list) <- trial_levels

for (i in seq_along(blues_list)) {
  
  b1 = data.frame(blues_list[[i]])
  b1$env = i
  b1$id = rownames(b1)
  
  b1 = dplyr::select(b1, env, id, effects)
  colnames(b1) = c("env", "id", "BLUE")
  b1_list[[i]] = b1
  
}


b1 = data.frame(rbind(b1_list[[1]],b1_list[[2]], b1_list[[3]],
                      b1_list[[4]],b1_list[[5]], b1_list[[6]],
                      b1_list[[7]],b1_list[[8]], b1_list[[9]],
                      b1_list[[10]],b1_list[[11]]
))

head(b1)
dim(b1)

# Palmito ---------------------------------------------------------------------

pal_list = read_rds("data/pal_list.RDS")

# Fit Stage One Model -----------------------------------------------------

# Get the levels of the TRIAL column
trial_levels <- seq_along(pal_list)
trial_levels

asr_list_palmito <- vector("list", length(trial_levels))     # Initialize pal_list with the number of TRIAL levels
names(asr_list_palmito) <- trial_levels 

blues_list2 <- vector("list", length(trial_levels))     # Initialize pal_list with the number of TRIAL levels
names(blues_list2) <- trial_levels

vcov_list2 <- vector("list", length(trial_levels))     # Initialize pal_list with the number of TRIAL levels
names(vcov_list2) <- trial_levels

names2_list <- vector("list", length(trial_levels))     # Initialize pal_list with the number of TRIAL levels
names(names2_list) <- trial_levels

for (i in seq_along(pal_list)) {
  
  sdat = pal_list[[i]] 
  
  sdat = sdat[order(sdat$COL, sdat$ROW), ]
  
  m1 = asreml(
    fixed = VOL_d3 ~ GEN,
    random = ~  REP + REP:IBLC,
    residual = ~ ar1(COL):ar1(ROW),
    maxit = 50, #trace = F,
    workspace = "1gb",
    na.action = na.method(x = "include", y = "include"),
    data = sdat)
  
  m1 = update.asreml(m1)
  m1 = update.asreml(m1)
  m1 = update.asreml(m1)
  m1 = update.asreml(m1)
  m1 = update.asreml(m1)
  m1 = update.asreml(m1)
  
  asr_list_palmito[[i]] = m1
  
  blues_temp = m1$coefficients$fixed
  
  temp1 = data.frame(blues_temp[grep("GEN_", rownames(blues_temp)),1])
  blues_list2[[i]] = temp1
  attr(blues_list2[[i]], "names") <- "effects"
  
  temp2 = predict.asreml(m1,classify = "GEN",vcov = T)$vcov
  rownames(temp2) = rownames(temp1)
  colnames(temp2) = rownames(temp1)
  vcov_list2[[i]] = temp2
  
  names2_list[[i]] = rownames(temp1)
  
}


# Save the BLUES -----------------------------------------------------
# Saved on the format used on stagewise package

b2_list <- vector("list", length(trial_levels))     # Initialize mat_list with the number of TRIAL levels
names(b2_list) <- trial_levels

for (i in seq_along(blues_list2)) {
  
  b2 = data.frame(blues_list2[[i]])
  b2$env = i
  b2$id = rownames(b2)
  
  b2 = dplyr::select(b2, env, id, effects)
  colnames(b2) = c("env", "id", "BLUE")
  b2_list[[i]] = b2
  
}


b2 = data.frame(rbind(b2_list[[1]],b2_list[[2]], b2_list[[3]],
                      b2_list[[4]],b2_list[[5]], b2_list[[6]],
                      b2_list[[7]],b2_list[[8]], b2_list[[9]],
                      b2_list[[10]],b2_list[[11]], b2_list[[12]]
))

head(b2)
dim(b2)

# BLUES -------------------------------------------------------------------

b1$loc = "1"
b2$loc = "2"

blues <- rbind(b1, b2)
blues$GEN_TRIAL_SITE = paste(blues$id, blues$env, blues$loc, sep="_")

dim(blues)
head(blues); tail(blues)
str(blues)

# Save the WEIGHTS -----------------------------------------------------------------

matrices1 <- lapply(vcov_list, function(x) as.matrix(x)) # Convert each list into a matrix

diag1_list = list()

for (i in seq_along(matrices1) ) {
  
  m1 = as.matrix(matrices1[[i]])
  
  # Specify the row and column names to exclude
  row_to_exclude <- "GEN_GENx"
  col_to_exclude <- "GEN_GENx"
  
  # Get the indices of the rows and columns to exclude
  row_index <- which(rownames(m1) == row_to_exclude)
  col_index <- which(colnames(m1) == col_to_exclude)
  
  # Exclude the specified row and column
  m1.red <- m1[-row_index, -col_index]
  
  
  diag1_list[[i]] = diag(solve(m1.red))
  
}

matrices2 <- lapply(vcov_list2, function(x) as.matrix(x)) # Convert each list into a matrix 2

diag2_list = list()

for (i in seq_along(matrices2) ) {
  
  m1 = as.matrix(matrices2[[i]])
  
  # Specify the row and column names to exclude
  row_to_exclude <- "GEN_GENx"
  col_to_exclude <- "GEN_GENx"
  
  # Get the indices of the rows and columns to exclude
  row_index <- which(rownames(m1) == row_to_exclude)
  col_index <- which(colnames(m1) == col_to_exclude)
  
  # Exclude the specified row and column
  m1.red <- m1[-row_index, -col_index]
  
  diag2_list[[i]] = diag(solve(m1.red))
  
}


list_df <- list()

for(i in seq_along(diag1_list) ) {
  
  w1 = data.frame(diag1_list[[i]])
  w1= cbind(rownames(w1), w1)
  w1$TRIAL = i
  w1$SITE = "1"
  names(w1) = c("GEN", "W", "TRIAL", "SITE")
  w1$GEN_TRIAL_SITE = paste(w1$GEN, w1$TRIAL, w1$SITE, sep = "_")
  list_df[[i]] <- data.frame(w1) 
}

weights1 <- do.call(rbind, list_df)


list_df <- list()

for(i in seq_along(diag2_list) ) {
  
  w1 = data.frame(diag2_list[[i]])
  w1= cbind(rownames(w1), w1)
  w1$TRIAL = i
  w1$SITE = "2"
  names(w1) = c("GEN", "W", "TRIAL", "SITE")
  w1$GEN_TRIAL_SITE = paste(w1$GEN, w1$TRIAL, w1$SITE, sep = "_")
  list_df[[i]] <- data.frame(w1) 
}

weights2 <- do.call(rbind, list_df)


# WHEIGHTS

W_df = rbind(weights1, weights2)



# STAGE 2 -----------------------------------------------------

head(blues)
head(W_df)
dim(W_df)

different_rows_blues <- data.frame(setdiff(blues$GEN_TRIAL_SITE, W_df$GEN_TRIAL_SITE))
different_rows_weights <- data.frame(setdiff(W_df$GEN_TRIAL_SITE, blues$GEN_TRIAL_SITE))


W_df_filt <- W_df[W_df$GEN_TRIAL_SITE %in% blues$GEN_TRIAL_SITE, ]
dim(W_df_filt)

blues_filt <- blues[blues$GEN_TRIAL_SITE %in% W_df_filt$GEN_TRIAL_SITE, ]
dim(blues_filt)


data_stage2 =  merge(blues_filt, W_df_filt, by = "GEN_TRIAL_SITE")
dim(data_stage2)

data_stage2 = data_stage2 %>% dplyr::select(GEN_TRIAL_SITE, GEN, TRIAL, SITE, BLUE, W)
dim(data_stage2)
head(data_stage2)

data_stage2$GEN_TRIAL_SITE = as.factor(data_stage2$GEN_TRIAL_SITE)
data_stage2$GEN = gsub("GEN_", "", data_stage2$GEN)
data_stage2$GEN = as.factor(data_stage2$GEN)

data_stage2$TRIAL = as.factor(data_stage2$TRIAL)
data_stage2$SITE = as.factor(data_stage2$SITE)


# GENO DATA -----------------------------------------------------------------

### A e D matrix


Ac = readRDS("C:/Users/LENOVO/OneDrive/POS-DOC/TPC/output/Ac.rds")
G = readRDS("C:/Users/LENOVO/OneDrive/POS-DOC/TPC/output/G.rds")

library(ASRgenomics)

Aclean <- ASRgenomics::match.G2A(A = as.matrix(Ac), G = as.matrix(G), clean = TRUE, ord = TRUE, mism = TRUE)$Aclean
Aclean[1:5, 1:5]

Ga.blend <- ASRgenomics::G.tuneup(G = G, A= Aclean, blend = TRUE, pblend = 0.02)$Gb
dim(Ga.blend)

Ga.blend_filt = Ga.blend[rownames(Ga.blend) %in% data_stage2$GEN, colnames(Ga.blend) %in% data_stage2$GEN]
dim(Ga.blend_filt)

pheno_stage2 = data_stage2[data_stage2$GEN %in% rownames(Ga.blend_filt), ]
dim(pheno_stage2)

Ga.inv <- ASRgenomics::G.inverse(G = Ga.blend_filt, sparseform = TRUE)$Ginv  # sparseform = TRUE

attr(Ga.inv,"INVERSE")<-TRUE
attr(Ga.inv, "rowNames") # irá atribuir o nome da Ga.blend_filt a Ga.inv


# Modelo Aditivo Dominante ----------------------------------------------------------

D.inv = readRDS("C:/Users/LENOVO/OneDrive/POS-DOC/TPC/output/Gd_inv.rds")

# A and D
pheno_stage2$GENA = as.factor(pheno_stage2$GEN)
pheno_stage2$GEND = as.factor(pheno_stage2$GEN)

pheno_stage2 <- pheno_stage2[order(pheno_stage2$GEN),]
pheno_stage2 <- pheno_stage2[order(pheno_stage2$TRIAL),]
pheno_stage2 <- pheno_stage2[order(pheno_stage2$SITE),]
str(pheno_stage2)

write.csv(pheno_stage2, "data/pheno_stage2.csv")



# idv

asreml.options(ai.sing = TRUE, threads= -1)

mod.ad1 <- asreml(
  fixed = BLUE ~ SITE,
   random = ~  vm(GENA, Ga.inv) +            #aditivo
               vm(GEND, D.inv) +             #dominante
               SITE:TRIAL +                  #trial (possivel pq a testemunha repete em cada trial)
               vm(GENA, Ga.inv):idv(SITE) +  #int axs
               vm(GEND, D.inv):idv(SITE),    #int dxs
  #weights = "W", 
  na.action = na.method(x = "include", y = "include"), 
  workspace = "4gb",
  data = pheno_stage2
)

vc.ad1 = summary(mod.ad1)$varcomp; vc.ad1
aic.ad1 = summary(mod.ad1)$aic; aic.ad1


# Results Stage One -------------------------------------------------------

## Matão

# Variance Components Stage one - Matao
vc1_matao =
  data.frame(round(rbind(asr_list_matao[["1"]][["vparameters"]],asr_list_matao[["2"]][["vparameters"]],
                         asr_list_matao[["3"]][["vparameters"]], asr_list_matao[["4"]][["vparameters"]],
                         asr_list_matao[["5"]][["vparameters"]], asr_list_matao[["6"]][["vparameters"]],
                         asr_list_matao[["7"]][["vparameters"]], asr_list_matao[["8"]][["vparameters"]],
                         asr_list_matao[["9"]][["vparameters"]], asr_list_matao[["10"]][["vparameters"]],
                         asr_list_matao[["11"]][["vparameters"]]),2))


vc1_matao$TRIAL = seq(1:nrow(vc1_matao))
vc1_matao$SITE = "MAT"
vc1_matao


# Variance Components Stage one - Palmito

vc1_palmito =
  data.frame(round(rbind(asr_list_palmito[["1"]][["vparameters"]], asr_list_palmito[["2"]][["vparameters"]],
                         asr_list_palmito[["3"]][["vparameters"]], asr_list_palmito[["4"]][["vparameters"]],
                         asr_list_palmito[["5"]][["vparameters"]], asr_list_palmito[["6"]][["vparameters"]],
                         asr_list_palmito[["7"]][["vparameters"]], asr_list_palmito[["8"]][["vparameters"]],
                         asr_list_palmito[["9"]][["vparameters"]], asr_list_palmito[["10"]][["vparameters"]],
                         asr_list_palmito[["11"]][["vparameters"]],asr_list_palmito[["12"]][["vparameters"]]),2))


vc1_palmito$TRIAL = seq(1:nrow(vc1_palmito))
vc1_palmito$SITE = "PAL"
vc1_palmito

# Results Stage Two -------------------------------------------------------


# Variance parameters and their index (number).
vc.ad = cbind.data.frame(
  summary(mod.ad1)$varcomp,
  number = seq_along(mod.ad1$vparameters))
vc.ad

#BLUEs
BLUE <- data.frame(summary(mod.ad1, coef=TRUE)$coef.fixed) ; BLUE

#BLUPs
BLUP <- data.frame(summary(mod.ad1, coef=TRUE)$coef.random); head(BLUP)

# create BLUP wide
pattern <- "vm\\(GENA"
indices_a <- grep(pattern, rownames(BLUP))
length(indices_a)
BLUP_wide <- BLUP[indices_a, ]

# ADITIVO: Usar grep, grepl, strsplit(), etc.

# Usar gsub para remover o prefixo 'vm(GENA, Ga.inv)_' e obter o restante da string
rownames(BLUP_wide) = gsub("^vm\\(GENA, Ga\\.inv\\)_", "", rownames(BLUP_wide))

# Separar linhas com nomes com menos de 6 caracteres 

short_names = nchar(rownames(BLUP_wide)) < 6   # ex: "01.15" (5 characteres)

# Criar o data frame com essas linhas
a =  BLUP_wide[short_names, -3 ]

# Separar linhas com nomes com mais de 6 caracteres 

long_names = nchar(rownames(BLUP_wide)) > 6   # ex: "01.15" (5 characteres)

# Criar o data frame com essas linhas
as =  BLUP_wide[long_names, ]


# grep para selecionar as linhas cujo nome termina com ":SITE_1"
as1 = as[grep(":SITE_1$", rownames(as)), -3]
rownames(as1) = gsub(":SITE_1", "", rownames(as1))

# grep para selecionar as linhas cujo nome termina com ":SITE_2"
as2 = as[grep(":SITE_2$", rownames(as)), -3]
rownames(as2) = gsub(":SITE_2", "", rownames(as2))

aditivo = cbind(a,as1,as2)
colnames(aditivo) = c("a", "se_a", "as1", "se_as1", "as2", "se_as2")
head(aditivo)

# DOMINANCIA: Usar grep, grepl, strsplit(), etc.

#BLUPs
BLUP <- data.frame(summary(mod.ad1, coef=TRUE)$coef.random); head(BLUP)

# create BLUP wide
pattern <- "vm\\(GEND"
indices_d <- grep(pattern, rownames(BLUP))
length(indices_d)
BLUP_wide <- BLUP[indices_d, ]

# Usar gsub para remover o prefixo 'vm(GEND, D.inv)_' e obter o restante da string
rownames(BLUP_wide) = gsub("^vm\\(GEND, D\\.inv\\)_", "", rownames(BLUP_wide))

# Separar linhas com nomes com menos de 6 caracteres 

short_names = nchar(rownames(BLUP_wide)) < 6   # ex: "01.15" (5 characteres)

# Criar o data frame com essas linhas
d =  BLUP_wide[short_names, -3 ]
d = d[rownames(aditivo), ]
dim(d)


# Separar linhas com nomes com mais de 6 caracteres 

long_names = nchar(rownames(BLUP_wide)) > 6   # ex: "01.15" (5 characteres)

# Criar o data frame com essas linhas
ds =  BLUP_wide[long_names, ]
dim(ds)

# grep para selecionar as linhas cujo nome termina com ":SITE_1"
ds1 = ds[grep(":SITE_1$", rownames(ds)), -3]
rownames(ds1) = gsub(":SITE_1", "", rownames(ds1))
ds1 = ds1[rownames(aditivo), ]
dim(ds1)

# grep para selecionar as linhas cujo nome termina com ":SITE_2"
ds2 = ds[grep(":SITE_2$", rownames(ds)), -3]
rownames(ds2) = gsub(":SITE_2", "", rownames(ds2))
ds2 = ds2[rownames(aditivo), ]
dim(ds2)

dominancia = cbind(d,ds1,ds2)
colnames(dominancia) = c("d", "se_d", "ds1", "se_ds1", "ds2", "se_ds2")
head(dominancia)


# BLUP oficial ------------------------------------------------------------

eBLUP = cbind(aditivo, dominancia)
head(eBLUP)

#Accuracy and Reliability - Additive
Va <- vc.ad1[2,1]; Va
eBLUP$PEV_a <- eBLUP$se_a^2
eBLUP$acc_a <- sqrt(1 - (eBLUP$PEV_a)/Va)
eBLUP$rel_a <- 1 - eBLUP$PEV_a/Va

mean_acc.a = mean(eBLUP$acc_a, na.rm = TRUE); mean_acc.a
mean_real.a = mean(eBLUP$rel_a, na.rm = TRUE); mean_real.a

#Accuracy and Reliability - Dominance
Vd <- vc.ad1[3,1]; Vd
eBLUP$PEV_d <- eBLUP$se_d^2
eBLUP$acc_d <- sqrt(1 - eBLUP$PEV_d/Vd)
eBLUP$rel_d <- 1 - eBLUP$PEV_d/Vd

mean_acc.d = mean(eBLUP$acc_d, na.rm = TRUE); mean_acc.d 
mean_real.d = mean(eBLUP$rel_d, na.rm = TRUE); mean_real.d


# Herdabilidade-------------------------------------------------------------------
vc.ad1
h2a = asreml::vpredict(mod.ad1, h2a~((V2)/(V1+V2+V3+V4+V5+V6))) ; h2a
d2 = asreml::vpredict(mod.ad1, d2~((V3)/(V1+V2+V3+V4+V5+V6))) ; d2
h2g = asreml::vpredict(mod.ad1, h2g~((V2 + V3)/(V1+V2+V3+V4+V5+V6))) ; h2g

herit.mod.ad1 = rbind(h2a, d2, h2g); herit.mod.ad1


# Ganho -------------------------------------------------------------------
str(eBLUP)

eBLUP$g = eBLUP$a + eBLUP$d
sel_g = eBLUP[order(eBLUP$g, decreasing = T),][1:30, ] # 1% dos selecionados equivale a aproximadamente 10% dos individuos no experimento
sel_g$Rank = seq(1:nrow(sel_g))

# Selected Clones

pheno_3y = read.csv("data/pheno_3y.csv")

# Encontra a familia, parcela, etc com base no pheno_3y
sel_gen = pheno_3y[pheno_3y$GEN  %in% rownames(sel_g) , ]
dim(sel_gen)
head(sel_g)

ordering <- factor(sel_gen$GEN, levels = rownames(sel_g))

# Order df1 based on this factor

selected_gen <- data.frame(sel_gen[order(ordering), , drop = FALSE])

mean_sel_fen = mean(selected_gen$VOL_d3,na.rm = TRUE); mean_sel_fen

mean_zero = mean(pheno_3y$VOL_d3); mean_zero

gain.g = (mean(sel_g$g)/mean(pheno_3y$VOL_d3))*100 ; gain.g


write_excel_csv(selected_gen, "C:/Users/LENOVO/OneDrive/POS-DOC/TPC/output2/selected_genotypes.csv")

# Save Image --------------------------------------------------------------

save.image("images/image_04.10.24.RData")

# load("images/image_01.10.24")



