######################################################
# Lepidoptera Habitat Specialization Analysis
# Moths
# Absolute niche breadth calculations (ANB)
######################################################



#Load Packages ----

library(tidyverse)
library(ade4)
library(lme4)
library(performance)
library(openxlsx)
library(readxl)
library(corrplot)
library(here)


#Data preparation ----

# Load moth community data
matcom_moths <- read_excel(here("Data/Raw/Lepidoptera_data.xlsx"), sheet = 1) %>%
  filter(Locality != "BF") %>%
  mutate(Plot = as.character(Plot))

# Basic dataset summary
n_individuals <- sum(matcom_moths$Abundance, na.rm = TRUE)
n_species <- n_distinct(matcom_moths$Species_complete)

# Load environmental data
matenv <- read_excel(here("Data/Raw/Lepidoptera_data.xlsx"), sheet = 3) %>%
  filter(Locality != "BF") %>%
  mutate(Plot = as.character(Plot))


#Environmental data processing ----

env_global <- matenv %>%
  select(
    Locality, Plot,
    Number_of_tree_individuals, Mean_DBH, Mean_tree_height,
    Stand_wood_volume, Dead_trees, Stem_slenderness_index,
    `Mean_%_Cnpy_Open___`
  ) %>%
  mutate(across(where(is.numeric), ~ scale(.x))) %>%
  na.omit()

# Principal Component Analysis of environmental variables
dudi_global <- dudi.pca(
  env_global %>% select(-Locality, -Plot),
  scale = FALSE,
  scannf = FALSE,
  nf = 3
)

#Niche Breadth Calculation Function ----

calculate_niche_metrics <- function(species_matrix, dudi_obj, group_name = "All") {
  species_df <- species_matrix %>%
    select(-any_of(c("Locality", "Plot", "Elevation", "Real_Elevation"))) %>%
    as.data.frame()
  
  if (nrow(species_df) != nrow(dudi_obj$li)) {
    stop("Row mismatch between species data (", nrow(species_df),
         ") and PCA (", nrow(dudi_obj$li), ")")
  }
  
  niche_result <- niche(dudi_obj, log1p(species_df), scannf = FALSE)
  tolerance <- niche.param(niche_result)[, 3]
  
  tolerance_matrix <- species_df
  for (sp in colnames(tolerance_matrix)) {
    tolerance_matrix[tolerance_matrix[[sp]] > 0, sp] <- tolerance[sp]
  }
  
  species_matrix %>%
    mutate(
      CWM_tolerance = sapply(1:nrow(tolerance_matrix), function(i) {
        x <- tolerance_matrix[i, ]
        abundances <- species_df[i, ]
        present <- x > 0
        if (sum(present) == 0) return(NA)
        weighted.mean(x[present], as.numeric(abundances[present]))
      }),
      Group = group_name
    ) %>%
    select(Locality, Plot, Elevation, CWM_tolerance, Group)
}

#Community-Level Analysis ----

community_matrix <- matcom_moths %>%
  group_by(Locality, Plot, Elevation, Species_complete) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  pivot_wider(
    names_from = Species_complete,
    values_from = Abundance,
    values_fill = list(Abundance = 0)
  ) %>%
  mutate(across(-c(Locality, Plot, Elevation), as.numeric))

community_results <- calculate_niche_metrics(
  species_matrix = community_matrix,
  dudi_obj = dudi_global,
  group_name = "All_moths"
)

#Family-Level Analysis ----

target_families <- c("Geometridae", "Erebidae", "Noctuidae")
family_results <- list()

# Elevation lookup table
elevation_lookup <- matenv %>%
  select(Locality, Plot, Real_Elevation) %>%
  distinct(Locality, Plot, .keep_all = TRUE)

for (family in target_families) {
  family_data <- matcom_moths %>% filter(Family == family)
  
  family_matrix <- env_global %>%
    select(Locality, Plot) %>%
    left_join(
      family_data %>%
        group_by(Locality, Plot, Species_complete) %>%
        summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop"),
      by = c("Locality", "Plot")
    ) %>%
    pivot_wider(
      names_from = Species_complete,
      values_from = Abundance,
      values_fill = list(Abundance = 0)
    ) %>%
    left_join(
      matcom_moths %>% select(Locality, Plot, Elevation) %>% distinct(),
      by = c("Locality", "Plot")
    ) %>%
    mutate(across(-c(Locality, Plot, Elevation), as.numeric)) %>%
    select(where(~ is.numeric(.) && sum(., na.rm = TRUE) > 0) | c(Locality, Plot, Elevation))
  
  if (ncol(family_matrix) > 3) {
    family_matrix <- family_matrix %>%
      left_join(elevation_lookup, by = c("Locality", "Plot")) %>%
      select(-any_of("Real_Elevation"))
    
    family_niche <- calculate_niche_metrics(
      species_matrix = family_matrix,
      dudi_obj = dudi_global,
      group_name = family
    )
    family_results[[family]] <- family_niche
  }
}

#Combine Results and Export ----

community_results <- community_results %>%
  select(-any_of("Real_Elevation")) %>%
  left_join(elevation_lookup, by = c("Locality", "Plot"))

family_results <- lapply(family_results, function(df) {
  df %>%
    select(-any_of("Real_Elevation")) %>%
    left_join(elevation_lookup, by = c("Locality", "Plot"))
})

all_results_moths <- bind_rows(
  community_results,
  bind_rows(family_results)
)

# Load existing workbook from butterfly analysis
wb <- loadWorkbook(here("Data/Raw/ANB.xlsx"))

# Add moth results as a new sheet
addWorksheet(wb, "Moths")
writeData(wb, "Moths", all_results_moths)

# Save workbook
saveWorkbook(wb, file = here("Data/Raw/ANB.xlsx"), overwrite = TRUE)

