# Tidy data
rm(list=ls()); gc()

# Increase memory size in R

# library(usethis)  
# usethis::edit_r_environ()  
# #Paste in the pop-up tab:  R_MAX_VSIZE=8Gb
# #Re-start R

library(tidyverse)
library(asreml)
library(lme4)
library(car)
library(snpReady)
library(AGHmatrix)
library(MCMCglmm)
library(ASRgenomics)
library(readxl)


# Phenotypic data ---------------------------------------------------------

pheno <- read_csv("data/new_pheno.csv")
head(pheno)
str(pheno)
dim(pheno)


# Transfomando em fatores
colnames_string = "COL,ROW,REP,IBLC,GEN,CLONE,SITE,TRIAL,YEAR,CHECK,group"  #no space after the comma

colnames_list = unlist(strsplit(colnames_string, ","))
class(colnames_list)

pheno[colnames_list] = lapply(pheno[colnames_list], as.factor)
str(pheno)

pheno$FAM2 = pheno$CLONE
dim(pheno)

# Pedigree ----------------------------------------------------------------

# Get pedigree file from phenotype file.
ped = read_excel("data_raw/TPC-MS_info_v2_unesp.xlsx", 
                          sheet = "pedigree")
head(ped)
ped[ped==0] = NA



# Pedigree Corrigido----------------------------------------------------------------

# Get pedigree file from phenotype file.
ped_c = read_excel("data_raw/TPC-MS_info_v2_unesp.xlsx", 
                 sheet = "pedigreeC")
head(ped_c)
ped_c[ped_c==0] = NA

# Get A.inv -----------------------------------------------------

# For ABLUP - get the inverse of the A-matrix
#Obtaining the inverse of the A-matrix (package asreml)

A.inv = ainverse(ped) #Returns the inverse (package ASReml).
head(A.inv)
attr(A.inv, "rowNames")
saveRDS(A.inv, "data/A_inv.rds")

# Get A.inv Corrigida -----------------------------------------------------

A.inv.C = ainverse(ped_c) #Returns the inverse (package ASReml).
head(A.inv.C)
attr(A.inv.C, "rowNames")

saveRDS(A.inv.C, "data/A_inv.C.rds")


# Molecular Data ----------------------------------------------------------

# FAMotypes must be coded as follow:
# 0 = 0 minor allele (homozygote for the major allele).
# 1 = 1 minor allele (heterozygote).
# 2 = 2 minor alleles (homozygote for the minor allele).

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("VariantAnnotation")

library(VariantAnnotation)
vcfFile <- "data_raw//TPC_MS_geno_20230713.vcf"
vcf <- VariantAnnotation::readVcf(vcfFile)

# Extract genotypes
genotypes <- VariantAnnotation::geno(vcf)
sample_names <- colnames(genotypes)

# Convert genotypes to a data frame
genotype_df <- as.data.frame(genotypes)
dim(genotype_df)
class(genotype_df)

geno = genotype_df[ , -c(1:2)]
geno = t(geno)
geno[1:5, 1:10]


# Match the names on pheno e geno data ------------------------------------

subject_id = read_xlsx("data_raw/marker-data/TPC-MS_geno_info.xlsx")

subject_id = subject_id[subject_id$VCF_IDGEN %in% rownames(geno), ]
dim(subject_id)

geno = as.data.frame(geno)

geno = add_column(geno, VCF_IDGEN = rownames(geno), .before = 1)
geno[1:5, 1:5]

geno_ok = data.frame(left_join(subject_id, geno, by= "VCF_IDGEN"))

rownames(geno_ok) = geno_ok$Subject_id
geno_ok[1:5, 1:7]

# Retire as colunas extras
geno = data.frame(geno_ok[ , -c(1:4)])

geno[1:5, 1:5]
dim(geno)
class(geno)

geno[geno=="0/0"] <- 0
geno[geno=="0/1"] <- 1
geno[geno=="1/1"] <- 2
geno[geno=="./."] <- NA
dim(geno)


geno012 <- apply(geno, 2, function(x) as.numeric(as.character(x)))
geno012[1:5, 1:5]
dim(geno012)
rownames(geno012)=rownames(geno)
colnames(geno012)=colnames(geno)


M <- snpReady::raw.data(
  geno012,
  frame = "wide",
  base = FALSE,
  sweep.sample = 0.95,
  call.rate = 0.90,
  maf = 0.01,
  imput = T,
  imput.type = "wright",
  outfile = "012"      # the allel with minor frequency receives 2
)


M <- M$M.clean |>
  round()

dim(M)
M[1:5, 1:5]
class(M)

saveRDS(M, "data/Marker_Matrix_012_TPC.rds")
#M = readRDS("output/Marker_Matrix_012_TPC.rds")
dim(M)


# Gblend matrix ---------------------------------------------------------------

## Calculating the additive genomic relationship matrix (Ga) and its inverse

# Get Gblend matrix.

M_mat = as.matrix(M)
rownames(M_mat)= rownames(M)
colnames(M_mat)= colnames(M)

M_mat[1:5, 1:5]  #confere os nomes 

# Do a little cleaning to speed it up
remove(vcfFile)
remove(vcf)
remove(geno)
remove(geno_ok)
remove(geno012)
gc()

 
G <- ASRgenomics::G.matrix(M = M_mat, method = "VanRaden")$G
G[1:5, 1:5]

# Match G and A.
attr(G, "rowNames") <- as.character(rownames(M_mat))
attr(G, "colNames") <- as.character(colnames(M_mat))

# Get A matrix
## Obtaining the pedigree-based additive relationship matrix A for blending
# We obtain the matrix A (not the inverse) using the package AGHmatrix:

ped_c= data.frame(ped_c)
ped_c[is.na(ped_c)] = 0

prep.A = nadiv::prepPed(ped_c)
Ac = as.matrix(nadiv::makeA(prep.A))

Ac = Ac[order(rownames(Ac)), order(colnames(Ac))]
G = G[order(rownames(G)), order(colnames(G))]

#Subsetting the matrix to only individuals present in the matrix Ga (i.e., excluding the parents).
Ac <- Ac[rownames(Ac) %in% rownames(G),  colnames(Ac) %in% colnames(G)]  #droplevels in rows of missing genotypes
G <- G[rownames(G) %in% rownames(Ac),  colnames(G) %in% colnames(Ac)]  #droplevels in rows of missing genotypes

#Making sure that the matrices are in the same order before blending:
dim(ped)
dim(Ac)
dim(G)

all(rownames(A) == rownames(G)) #order and name of individuals

saveRDS(Ac, "data/Ac.rds")
saveRDS(G, "data/G.rds")

Aclean <- match.G2A(A = as.matrix(A), G = as.matrix(G), clean = TRUE, ord = TRUE, mism = TRUE)$Aclean
Aclean[1:5, 1:5]

Ga.blend <- G.tuneup(G = G, A= Aclean, blend = TRUE, pblend = 0.02)$Gb
Ga.inv <- G.inverse(G = Ga.blend, sparseform = TRUE)$Ginv  # sparseform = TRUE

attr(Ga.inv,"INVERSE")<- TRUE
attr(Ga.inv, "rowNames")
saveRDS(Ga.inv, "data/Ga_inv.rds")

# H matrix ----------------------------------------------------------------
M = readRDS("data/Marker_Matrix_012_TPC.rds")

library(AGHmatrix)
#Computing Ac matrix (Munoz)
H_Munoz <- Hmatrix(A=Ac, G=G, markers = M,
                      ploidy=2, method="Munoz",
                      roundVar=2,
                      maf=0.05)
# ‘Munoz‘ shrinks the G matrix towards the A matrix scaling the molecular relatadness
#by each relationship classes; 

H_Munoz.inv = solve(H_Munoz) #Get inverse of matrix Ga_blended by solving.

## Getting ready for ASReml-r V4 package

#Transform to sparse matrix (i.e., sorting in a different way, only records were there are no zeros).
H_Munoz.inv <-as(H_Munoz.inv ,"sparseMatrix")
#Get lower diagonal of the matrix and get ready for use in ASReml (package MCMCglmm).
H_Munoz.inv  <-MCMCglmm::sm2asreml(H_Munoz.inv)
#Set "INVERSE" as an attribute. This is needed for ASReml-R V4.
attr(H_Munoz.inv, "INVERSE") <- TRUE
attr(H_Munoz.inv, "rowNames")
dim(H_Munoz.inv)

saveRDS(H_Munoz.inv, "data/Ac_inv.rds")


# H matrix ----------------------------------------------------------------
#Computing H matrix (Martini)/ Legarra)2009)

Hmat_Legarra <- Hmatrix(A=Ac, G=G, method="Martini",
                        ploidy=2,
                        maf=0.05,
                        tau = 1,
                        omega = 1)


# ‘Martini‘ is a modified version from Legarra et al. (2009) where
# combines A and G matrix using scaling factors. When method is equal ‘Martini‘ and ‘tau=1‘ and
# ‘omega=1‘ you have the same H matrix as in Legarra et al. (2009).

H_Legarra.inv = solve(Hmat_Legarra) #Get inverse of matrix Ga_blended by solving.

## Getting ready for ASReml-r V4 package
#Transform to sparse matrix (i.e., sorting in a different way, only records were there are no zeros).

H_Legarra.inv <-as(H_Legarra.inv,"sparseMatrix")

#Get lower diagonal of the matrix and get ready for use in ASReml (package MCMCglmm).

H_Legarra.inv <-MCMCglmm::sm2asreml(H_Legarra.inv)

#Set "INVERSE" as an attribute. This is needed for ASReml-R V4.
attr(H_Legarra.inv, "INVERSE") <- TRUE
attr(H_Legarra.inv, "rowNames")

saveRDS(H_Legarra.inv, "outputs/H_inv.rds")


# Dominance Matrix --------------------------------------------------------

## Calculating the dominant genomic relationship matrix (Gd) and its inverse

#Obtain the matrix Gd (package AGHmatrix).
Gd = Gmatrix(SNPmatrix = as.matrix(M), method = "Vitezica")
Gd[1:5, 1:5]
# It is not possible to obtain the pedigree-based dominant relationship matrix (D)
# for blending, since we are working with a half-sif families

attr(Gd, "rowNames") <- as.character(rownames(M))
attr(Gd, "colNames") <- as.character(colnames(M))

D.inv <- G.inverse(G = Gd, sparseform = TRUE)$Ginv
D.inv[1:5, ]

# Gyinv.dom  <- G.inverse(G = Gyb.dom, sparseform = TRUE)$Ginv.sparse
# Gyinv.dom[1:5, 1:3]

attr(D.inv,"INVERSE")<-TRUE
attr(D.inv, "rowNames")
saveRDS(D.inv, "data/D_inv.rds")


# Ga ######## -------------------------------------------------------------

G <- ASRgenomics::G.matrix(M = M_mat, method = "VanRaden")$G
G[1:5, 1:5]

# Match G and A.
attr(G, "rowNames") <- as.character(rownames(M_mat))
attr(G, "colNames") <- as.character(colnames(M_mat))

# Get A matrix
## Obtaining the pedigree-based additive relationship matrix A for blending
# We obtain the matrix A (not the inverse) using the package AGHmatrix:

ped_c= data.frame(ped_c)
ped_c[is.na(ped_c)] = 0


prep.A = nadiv::prepPed(ped_c)
Ac = as.matrix(nadiv::makeA(prep.A))


Ac = Ac[order(rownames(Ac)), order(colnames(Ac))]
G = G[order(rownames(G)), order(colnames(G))]

#Subsetting the matrix to only individuals present in the matrix Ga (i.e., excluding the parents).
Ac <- Ac[rownames(Ac) %in% rownames(G),  colnames(Ac) %in% colnames(G)]  #droplevels in rows of missing genotypes
G <- G[rownames(G) %in% rownames(Ac),  colnames(G) %in% colnames(Ac)]  #droplevels in rows of missing genotypes

#Making sure that the matrices are in the same order before blending:
dim(ped)
dim(Ac)
dim(G)

all(rownames(Ac) == rownames(G)) #order and name of individuals

saveRDS(Ac, "data/Ac.rds")
saveRDS(G, "data/G.rds")


Aclean <- match.G2A(A = as.matrix(Ac), G = as.matrix(G), clean = TRUE, ord = TRUE, mism = TRUE)$Aclean
Aclean[1:5, 1:5]

Ga.blend <- G.tuneup(G = G, A= Aclean, blend = TRUE, pblend = 0.02)$Gb

Ga.inv <- G.inverse(G = Ga.blend, sparseform = TRUE)$Ginv  # sparseform = TRUE

attr(Ga.inv,"INVERSE")<-TRUE
attr(Ga.inv, "rowNames")
saveRDS(Ga.inv, "data/Ga_inv.rds")


