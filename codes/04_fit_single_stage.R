# fit_models

rm(list = ls())

zdata = read.csv("data/zdata.csv")
names(zdata)

# Transfomando em fatores
colnames_string = "COL,ROW,REP,IBLC,GEN,CLONE,SITE,TRIAL,YEAR,CHECK,group"  #no space after the comma

colnames_list = unlist(strsplit(colnames_string, ","))
class(colnames_list)

zdata[colnames_list] = lapply(zdata[colnames_list], as.factor)
str(zdata)

# Ordering data : the data need to be ordered as specified by the residual term 

zdata <- zdata[order(zdata$GEN),]
zdata <- zdata[order(zdata$TRIAL),]
zdata <- zdata[order(zdata$SITE),]
str(zdata)

# Modelo Aditivo Dominante ----------------------------------------------------------

A.inv = readRDS("data/Ga_inv.rds")
D.inv = readRDS("data/D_inv.rds")

# A and D
zdata$GENA = as.factor(zdata$GEN)
zdata$GEND = as.factor(zdata$GEN)

asreml.options(sing.ai=TRUE, threads=-1)
m9 <- asreml(
  fixed = VOL_d3 ~ SITE + group + SITE:TRIAL:REP,
  random = ~ vm(GENA, A.inv) +     # aditivo
    vm(GEND, D.inv) +               # dominante
    SITE:TRIAL:REP:IBLC +   
    at(group):ar1(ROW):ar1(COL),    # residual in the random part
  residual = ~dsum(~units|SITE),    # residual
  #residual = ~dsum(~units|SITE),   
  na.action = na.method(x = "exclude", y = "exclude"), 
  workspace = "4gb",
  data = zdata,
  sing.ai=TRUE
)

write_rds(m9, "outputs/outputs_m9.rds")


# Results -----------------------------------------------------------------

# BLUES e BLUP ------------------------------------------------------------

#BLUEs
BLUE = data.frame(summary(m9, coef=TRUE)$coef.fixed)
dim(BLUE)
#BLUPs
BLUP = data.frame(summary(m9, coef=TRUE)$coef.random); dim(BLUP)

# create BLUP wide
pattern <- "vm\\(GENA"
indices_a <- grep(pattern, rownames(BLUP))
length(indices_a)
BLUP_wide <- BLUP[indices_a, ]
dim(BLUP_wide)

# ADITIVO: Usar grep, grepl, strsplit(), etc.

# Usar gsub para remover o prefixo 'vm(GENA, Ga.inv)_' e obter o restante da string
rownames(BLUP_wide) = gsub("^vm\\(GENA, A\\.inv\\)_", "", rownames(BLUP_wide))

# Separar linhas com nomes com menos de 6 caracteres 

short_names = nchar(rownames(BLUP_wide)) < 6   # ex: "01.15" (5 characteres)

# Criar o data frame com essas linhas
a =  BLUP_wide[short_names, -3 ]
dim(a)
head(a)


# DOMINANCIA: Usar grep, grepl, strsplit(), etc.

#BLUPs
BLUP = data.frame(summary(m9, coef=TRUE)$coef.random); dim(BLUP)

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
d = d[rownames(a), ]


# eBLUP ------------------------------------------------------------
head(a)
head(d)
eBLUP = data.frame(cbind(a, d))
head(eBLUP)
colnames(eBLUP) = c("a", "se_a", "d", "se_d")


# Variance Components -----------------------------------------------------

vc.m9 = cbind.data.frame(
  summary(m9)$varcomp)[ ,c(1,2)]
vc.m9

vc.m9$component = round(vc.m9$component,2)
vc.m9$std.error= round(vc.m9$std.error,2)

write.csv(vc.m9, "results/varcomp_m9.csv")

varcomp = vc.m9[!grepl("cor$", rownames(vc.m9)), ]
varcomp
dim(varcomp)

varcomp[nrow(varcomp)+1, ] = cbind(mean(varcomp$component[1:23]), 
                                   mean(varcomp$std.error[1:23]))


varcomp.m9 = varcomp[c(24:29), ]
rownames(varcomp.m9) = c("IBLOCO", "A", "D", "AR1_site1", "AR1_site2", "residual_nugget")
varcomp.m9$num = seq(1:nrow(varcomp.m9))
varcomp.m9

# Herdabilidade -----------------------------------------------------------

vc.m9$num =seq(1:nrow(vc.m9))
vc.m9
vc.m9[!grepl("cor$", rownames(vc.m9)), ]

h2a = asreml::vpredict(m9,h2a~((V71)/(V70+V71+ V72 + ((V73 + V74)/2) +
                     ((V1+V4+V7+V10+V13+V16+V19+V22+V25+V28+
                      V31+V34+V37+V40+V43+V46+V49+V52+V55+V58+
                      V61+V64+V67)/23)))); h2a


d2 = asreml::vpredict(m9,d2~((V72)/(V70+V71+ V72 + ((V73 + V74)/2) +
                                        ((V1+V4+V7+V10+V13+V16+V19+V22+V25+V28+
                                            V31+V34+V37+V40+V43+V46+V49+V52+V55+V58+
                                            V61+V64+V67)/23)))); d2

H2 = asreml::vpredict(m9,H2~((V71 + V72)/(V70+V71+ V72 + ((V73 + V74)/2) +
                                       ((V1+V4+V7+V10+V13+V16+V19+V22+V25+V28+
                                           V31+V34+V37+V40+V43+V46+V49+V52+V55+V58+
                                           V61+V64+V67)/23)))); H2

# Additive: Accuracy and Reliability ------------------------------

Va <- varcomp.m9[2,1]; Va
eBLUP$PEV_a <- eBLUP$se_a^2
eBLUP$acc_a <- sqrt(1 - (eBLUP$PEV_a)/Va)
eBLUP$rel_a <- 1 - eBLUP$PEV_a/Va

mean_acc.a = mean(eBLUP$acc_a, na.rm = TRUE); mean_acc.a
mean_real.a = mean(eBLUP$rel_a, na.rm = TRUE); mean_real.a

#Accuracy and Reliability - Dominance
Vd <- varcomp.m9[3,1]; Vd
eBLUP$PEV_d <- eBLUP$se_d^2
eBLUP$acc_d <- sqrt(1 - eBLUP$PEV_d/Vd)
eBLUP$rel_d <- 1 - eBLUP$PEV_d/Vd

mean_acc.d = mean(eBLUP$acc_d, na.rm = TRUE); mean_acc.d
mean_real.d = mean(eBLUP$rel_d, na.rm = TRUE); mean_real.d


# Ganho -------------------------------------------------------------------
str(eBLUP)

eBLUP$g = eBLUP$a + eBLUP$d
sel_g = eBLUP[order(eBLUP$g, decreasing = T),][1:30, ] # 1% dos selecionados equivale a aproximadamente 10% dos individuos no experimento
sel_g$Rank = seq(1:nrow(sel_g))

gain.g = (mean(sel_g$g))/
  mean(zdata$VOL_d3, na.rm=T) * 100 ; gain.g


# Encontra a familia, parcela, etc com base no pheno_3y
sel_gen = zdata[zdata$GENA  %in% rownames(sel_g) , ]
dim(sel_gen)
head(sel_gen)

ordering <- factor(sel_gen$GEN, levels = rownames(sel_g))

# Order df1 based on this factor

selected_gen <- sel_gen[order(ordering), , drop = FALSE]

mean_sel = mean(selected_gen$VOL_d3,na.rm = TRUE); mean_sel
mean_zero = mean(as.numeric(zdata$VOL_d3), na.rm=TRUE); mean_zero

# ROW COL Correlations -----------------------------------------------------------

cor_RC = vc.m9[grepl("cor$", rownames(vc.m9)), ]
cor_RC

# Extract text inside quotes for each row in the 'text' column
cor_RC$extracted = str_extract(rownames(cor_RC), "'([^']+)'")
cor_RC$extracted = gsub("'", "", cor_RC$extracted)
cor_RC$direction = rep(c("ROW", "COL"), 23)
cor_RC$direction = paste(cor_RC$extracted, cor_RC$direction, sep = "_")
cor_RC= cor_RC %>% select(direction, component, std.error) 
cor_RC

write.csv(cor_RC, "results/cor_ROW_COL_AR1.csv")

write.csv(selected_gen, "results/selected_genotypes.csv")


# save image --------------------------------------------------------------

save.image("images/image_04_fit_models_27.09.R")


