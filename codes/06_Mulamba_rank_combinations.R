
# Mulamba Ranking among the combinations for each scenario

# Load required objects from previously analysis

load("images/image_H2_AC_AD-GBLUP.RData") 

# Pacotes
library(tidyverse)


BLUP_list = list()

y = list(one_rep, two_rep, three_rep, four_rep, five_rep)

for (j in 1:length(y)) {

  results_list1 <- list()
  
  x = combn((length(y[[j]][["Summary"]][["3"]])),1)
  
  for (i in 1:ncol(x)){
    
    BLUP =   y[[j]][["Summary"]][["3"]][[i]][[1]][["coefficients"]][["random"]]
    
    pattern <- "vm\\(GENA"
    indices_a <- grep(pattern, rownames(BLUP))
    BLUP_wide <- data.frame(BLUP[indices_a, ])
    
    pattern <- "vm\\(GEND"
    indices_d <- grep(pattern, rownames(BLUP))
    d <- data.frame(BLUP[indices_d, ])
    BLUP_wide = data.frame(cbind(BLUP_wide, d))
    colnames(BLUP_wide) = c("a", "d")
    
    BLUP_wide$g = BLUP_wide$a + BLUP_wide$d
    
    # Change name
    rownames(BLUP_wide) = rownames(BLUP_wide) %>% str_replace_all( "[^0-9.]", "") %>% str_sub(2)
    
    BLUP_wide = BLUP_wide[order(BLUP_wide$g, decreasing = T),] 
    BLUP_wide$Rank = seq(1:nrow(BLUP_wide))
    
    results_list1[[i]] = BLUP_wide
    
  } # fecha o i
 
  BLUP_list[[j]] = results_list1
  
} # fecha o j


# Extract the "Rank" columns, align them by row names, and sum them row-wise

values_aligned <- lapply(BLUP_list, function(rep_list) {
  lapply(rep_list, function(comb_list) {
    # Access the second element in the info_list and calculate its mean
    comb_list[order(rownames(comb_list)), "Rank"]
  })
})


Mulamba_list_all = list()
Mulamba_list_sel = list()


for (k in 1:5) {
  
values_aligned2 = values_aligned[[k]]

# Create a matrix from the aligned values
values_matrix <- do.call(cbind, values_aligned2)
rownames(values_matrix) <- sort(rownames(BLUP_list[[1]][[1]]))

# Sum the values row-wise
summed_values <- rowSums(values_matrix)

# Create a new data frame with the summed values
summed_df <- data.frame(Sum_Rank = summed_values)
rownames(summed_df) <- rownames(values_matrix)
summed_df$IND = rownames(summed_df)

# order summed data frame (from smallest to highest) - Mulamba Rank and select the 30 best-ranked

Mulamba_list_all[[k]] = summed_df[order(summed_df$Sum_Rank), ]

Mulamba_list_sel[[k]] = summed_df[order(summed_df$Sum_Rank), ][1:30, ] # 1% dos selecionados equivale a aproximadamente 10% dos individuos no experimento

}

# How many clone selections were required for n=1,2,3, and 4 to ensure that 
# the best-10 clones (for n=5) were included in the selections?

best_5rep = as.vector(rownames(Mulamba_list_all[[5]])[1:10])

# Find the positions of names in the dataframe - one rep
positions = which(rownames(Mulamba_list_all[[1]]) %in% best_5rep)
positions1 = data.frame(Name = rownames(Mulamba_list_all[[1]])[positions], Position = positions)
positions1

# Find the positions of names in the dataframe - two rep
positions = which(rownames(Mulamba_list_all[[2]]) %in% best_5rep)
positions2 = data.frame(Name = rownames(Mulamba_list_all[[2]])[positions], Position = positions)
positions2

# Find the positions of names in the dataframe - three rep
positions = which(rownames(Mulamba_list_all[[3]]) %in% best_5rep)
positions3 = data.frame(Name = rownames(Mulamba_list_all[[3]])[positions], Position = positions)
positions3

# Find the positions of names in the dataframe - four rep
positions = which(rownames(Mulamba_list_all[[4]]) %in% best_5rep)
positions4 = data.frame(Name = rownames(Mulamba_list_all[[4]])[positions], Position = positions)
positions4

positions5 = data.frame(cbind(best_5rep, positon = seq(1:10)))

howmany_10 = data.frame(cbind(position5, positions4, positions3, positions2, positions1))
howmany_10
# Make the changes in ranking figure --------------------------------------

data_rank =  cbind(Mulamba_list[[1]], Mulamba_list[[2]], Mulamba_list[[3]], Mulamba_list[[4]], Mulamba_list[[5]])
head(data_rank)

data_rank = data_rank[ , c(2, 4, 6, 8, 10)]
colnames(data_rank) = c("Rank1", "Rank2", "Rank3", "Rank4", "Rank5")

# Convert character columns to numeric
data_rank2 = data.frame(lapply(data_rank, function(x) {
  if (is.character(x)) {
    as.numeric(x)
  } else {
    x
  }
}))


# Correspondence ----------------------------------------------------------

# Correspondence 10 
rank_df = data_rank2[1:10, ]

# Function to find the number of common elements between two vectors
count_common_elements <- function(vec1, vec2) {
  intersect_length <- length(intersect(vec1, vec2))
  return(intersect_length)
}

# Initialize a matrix to store the number of common elements between pairwise columns
common_elements_matrix <- matrix(NA, ncol = ncol(rank_df), nrow = ncol(rank_df), 
                                 dimnames = list(colnames(rank_df), colnames(rank_df)))

# Find the number of common elements between pairwise columns
for (i in 1:(ncol(rank_df) - 1)) {
  for (j in (i + 1):ncol(rank_df)) {
    common_elements_matrix[i, j] <- count_common_elements(rank_df[[i]], rank_df[[j]])
    common_elements_matrix[j, i] <- common_elements_matrix[i, j]  # Reflect across the diagonal
  }
}

correspondence_10 =  ggcorrplot::ggcorrplot(common_elements_matrix/10, hc.order = FALSE, type = "upper",
                       lab = TRUE, colors = c("red", "white", "red"))

# Correspondence 20 
rank_df = data_rank2[1:20, ]

# Function to find the number of common elements between two vectors
count_common_elements <- function(vec1, vec2) {
  intersect_length <- length(intersect(vec1, vec2))
  return(intersect_length)
}

# Initialize a matrix to store the number of common elements between pairwise columns
common_elements_matrix <- matrix(NA, ncol = ncol(rank_df), nrow = ncol(rank_df), 
                                 dimnames = list(colnames(rank_df), colnames(rank_df)))

# Find the number of common elements between pairwise columns
for (i in 1:(ncol(rank_df) - 1)) {
  for (j in (i + 1):ncol(rank_df)) {
    common_elements_matrix[i, j] <- count_common_elements(rank_df[[i]], rank_df[[j]])
    common_elements_matrix[j, i] <- common_elements_matrix[i, j]  # Reflect across the diagonal
  }
}

correspondence_20 =  ggcorrplot::ggcorrplot(common_elements_matrix/20, hc.order = FALSE, type = "upper",
                                            lab = TRUE, colors = c("red", "white", "red"))

# Correspondence 30 
rank_df = data_rank2[1:30, ]

# Function to find the number of common elements between two vectors
count_common_elements <- function(vec1, vec2) {
  intersect_length <- length(intersect(vec1, vec2))
  return(intersect_length)
}

# Initialize a matrix to store the number of common elements between pairwise columns
common_elements_matrix <- matrix(NA, ncol = ncol(rank_df), nrow = ncol(rank_df), 
                                 dimnames = list(colnames(rank_df), colnames(rank_df)))

# Find the number of common elements between pairwise columns
for (i in 1:(ncol(rank_df) - 1)) {
  for (j in (i + 1):ncol(rank_df)) {
    common_elements_matrix[i, j] <- count_common_elements(rank_df[[i]], rank_df[[j]])
    common_elements_matrix[j, i] <- common_elements_matrix[i, j]  # Reflect across the diagonal
  }
}

correspondence_30 =  ggcorrplot::ggcorrplot(common_elements_matrix/30, hc.order = FALSE, type = "upper",
                                            lab = TRUE, colors = c("red",  "white",  "red"))


library(gridExtra)
svg("figures/correspondence_AD_GBLUP.svg", width = 7, height = 7)
grid.arrange(correspondence_10, correspondence_20, correspondence_30, nrow=1, ncol =3)
dev.off()


data_rank2$Position = rownames(data_rank2)
head(data_rank2)


data_rank3 = data_rank2[1:10, ]

library(tidyr)

rank_long = data.frame(data_rank3 %>% 
  tidyr::pivot_longer(
    cols = `Rank1`:`Rank5`, 
    names_to = "Ranking",
    values_to = "Gen"
  )
)

rank_long$Gen = as.character(rank_long$Gen)  # sua coluna IND
rank_long$Ranking = as.numeric(gsub("Rank", "", rank_long$Ranking)) #sua coluna MODEL
rank_long$Position = as.numeric(rank_long$Position) #sua coluna Position


str(rank_long)

head(rank_long)

rank_long


colours1 <- c("#7F8C8D", "#000000", "#56B4E9", "#76448A", "#F0E442",
           "#0072B2", "#D55E00", "#CC79A7", "#48647B", "#CA6F1E")

library(ggbump)

ggplot(rank_long , aes(x = Ranking, y = Position, color = Gen))+
  geom_point(size = 5) +
  geom_text(data = rank_long %>% filter(Ranking == min(Ranking)),
            aes(x = Ranking - .1, label = Gen), size = 5, hjust = 1) +
  geom_text(data = rank_long %>% filter(Ranking == max(Ranking)),
            aes(x = Ranking + .1, label = Gen), size = 5, hjust = 0) +
  geom_bump(size = 1.5, geom = "line") + scale_y_reverse() + scale_x_reverse() +
  cowplot::theme_minimal_grid(font_size = 14, line_size = 0) +
  theme(legend.position = "none",
        panel.grid.major = element_blank()) +
       scale_colour_discrete(type = colours1)


