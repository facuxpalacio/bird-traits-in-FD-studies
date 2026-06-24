# Load packages
library(dplyr)
library(tidyr)
library(iNEXT)

# Load data
data <- read.delim2("C:/RD/Database bird funct traits.txt", comment.char="#")
data$trait_study <- iconv(data$trait_study, from = "", to = "UTF-8") # proper encoding
data$trait_study <- paste(data$studyID, data$trait_synonym, sep = "_")

# Collapse duplicated rows
data_clean <- data %>%
  mutate(
    studyID = iconv(studyID, from = "", to = "UTF-8"),
    trait_synonym = iconv(trait_synonym, from = "", to = "UTF-8")
  ) %>%
  filter(
    !is.na(studyID), studyID != "",
    !is.na(trait_synonym), trait_synonym != ""
  ) %>%
  distinct(studyID, trait_synonym) %>%
  mutate(presence = 1)

# Matrix of studies x traits
mat <- data_clean %>%
  pivot_wider(
    names_from = studyID,
    values_from = presence,
    values_fill = 0
  )
mat_sp_site <- as.matrix(mat[,-1])
rownames(mat_sp_site) <- mat$trait_synonym
mat_site_sp <- t(mat_sp_site)

# Rarefaction curves
# number of studies
n_sites <- nrow(mat_site_sp)

# incidence per trait
incidence <- colSums(mat_site_sp)
incidence_vec <- c(n_sites, incidence)

out <- iNEXT(
  incidence_vec,
  datatype = "incidence_freq",
  q = 0   # "species" richness (traits)
)

# plot
ggiNEXT(out, type = 3)

# estimate "richness" (number of studies)
D_est <- estimateD(
  incidence_vec,
  datatype = "incidence_freq",
  base = "size",
  nboot = 1000,
  level = NULL,
  conf = 0.95
)

D_est
