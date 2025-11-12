######################################################
# Lepidoptera Habitat Specialization Analysis
# Moths
# Within-Elevation niche breadth calculations (WENB)
######################################################


# Load Packages ----
library(readxl)
library(ade4)
library(reshape2)
library(dplyr)
library(tidyr)
library(here)
library(openxlsx)
library(tibble)   

# Data Preparation ----
matcom_moths <- read_excel("Data/Raw/Lepidoptera_data.xlsx", sheet = 1) %>%
  mutate(Plot = as.character(Plot))

matenv <- read_excel(here("Data/Raw/Lepidoptera_data.xlsx"), sheet = 3) %>%
  filter(Locality != "BF") %>%
  mutate(Plot = as.character(Plot))

# Environmental Data Processing ----
env_global <- matenv %>%
  select(
    Locality, Plot,
    Number_of_tree_individuals, Mean_DBH, Mean_tree_height,
    Stand_wood_volume, Dead_trees, Stem_slenderness_index,
    `Mean_%_Cnpy_Open___`
  ) %>%
  mutate(across(where(is.numeric), ~ scale(.x))) %>%
  na.omit()

# Community Matrix ----
moths_loc <- dcast(
  matcom_moths, Locality + Plot ~ Species_complete,
  value.var = "Abundance", fun.aggregate = sum
)

# Analyses per Locality ----
BCm <- moths_loc[(moths_loc$Locality == "BC"), -c(1, 2)]
CLm <- moths_loc[(moths_loc$Locality == "CL"), -c(1, 2)]
DGm <- moths_loc[(moths_loc$Locality == "DG"), -c(1, 2)]
ECm <- moths_loc[(moths_loc$Locality == "EC"), -c(1, 2)]
MSm <- moths_loc[(moths_loc$Locality == "MS"), -c(1, 2)]
PCm <- moths_loc[(moths_loc$Locality == "PC"), -c(1, 2)]

# Exclude Singletons ----
BCm1 <- BCm[, colSums(BCm > 0) > 0, drop = FALSE]
CLm1 <- CLm[, colSums(CLm > 0) > 0, drop = FALSE]
DGm1 <- DGm[, colSums(DGm > 0) > 0, drop = FALSE]
ECm1 <- ECm[, colSums(ECm > 0) > 0, drop = FALSE]
MSm1 <- MSm[, colSums(MSm > 0) > 0, drop = FALSE]
PCm1 <- PCm[, colSums(PCm > 0) > 0, drop = FALSE]

# Copies for OMI/Tolerance Outcomes ----
BCm1t <- BCm1; CLm1t <- CLm1; DGm1t <- DGm1
ECm1t <- ECm1; MSm1t <- MSm1; PCm1t <- PCm1

# Environmental Descriptors by Locality ----
BCe <- env_global[(env_global$Locality == "BC"), ]
CLe <- env_global[(env_global$Locality == "CL"), ]
DGe <- env_global[(env_global$Locality == "DG"), ]
ECe <- env_global[(env_global$Locality == "EC"), ]
MSe <- env_global[(env_global$Locality == "MS"), ]
PCe <- env_global[(env_global$Locality == "PC"), ]

# OMI Method (scale = TRUE) ----
dudiBC <- dudi.pca(BCe[, -c(1:4)], scale = TRUE, scan = FALSE, nf = 3)
dudiCL <- dudi.pca(CLe[, -c(1:4)], scale = TRUE, scan = FALSE, nf = 3)
dudiDG <- dudi.pca(DGe[, -c(1:4)], scale = TRUE, scan = FALSE, nf = 3)
dudiEC <- dudi.pca(ECe[, -c(1:4)], scale = TRUE, scan = FALSE, nf = 3)
dudiMS <- dudi.pca(MSe[, -c(1:4)], scale = TRUE, scan = FALSE, nf = 3)
dudiPC <- dudi.pca(PCe[, -c(1:4)], scale = TRUE, scan = FALSE, nf = 3)

# Community Data on OMI Space (scale = TRUE) ----
nicBCm <- niche(dudiBC, log(BCm1 + 1), scann = FALSE)
nicCLm <- niche(dudiCL, log(CLm1 + 1), scann = FALSE)
nicDGm <- niche(dudiDG, log(DGm1 + 1), scann = FALSE)
nicECm <- niche(dudiEC, log(ECm1 + 1), scann = FALSE)
nicMSm <- niche(dudiMS, log(MSm1 + 1), scann = FALSE)
nicPCm <- niche(dudiPC, log(PCm1 + 1), scann = FALSE)

# OMI Method (scale = FALSE) ----
dudiBC <- dudi.pca(BCe[, -c(1:4)], scale = FALSE, scan = FALSE, nf = 3)
dudiCL <- dudi.pca(CLe[, -c(1:4)], scale = FALSE, scan = FALSE, nf = 3)
dudiDG <- dudi.pca(DGe[, -c(1:4)], scale = FALSE, scan = FALSE, nf = 3)
dudiEC <- dudi.pca(ECe[, -c(1:4)], scale = FALSE, scan = FALSE, nf = 3)
dudiMS <- dudi.pca(MSe[, -c(1:4)], scale = FALSE, scan = FALSE, nf = 3)
dudiPC <- dudi.pca(PCe[, -c(1:4)], scale = FALSE, scan = FALSE, nf = 3)

total_inertia <- c(
  BC = sum(dudiBC$eig),
  CL = sum(dudiCL$eig),
  DG = sum(dudiDG$eig),
  EC = sum(dudiEC$eig),
  MS = sum(dudiMS$eig),
  PC = sum(dudiPC$eig)
)

# Community Data on OMI Space (scale = FALSE) ----
nicBCm <- niche(dudiBC, log(BCm1 + 1), scann = FALSE)
nicCLm <- niche(dudiCL, log(CLm1 + 1), scann = FALSE)
nicDGm <- niche(dudiDG, log(DGm1 + 1), scann = FALSE)
nicECm <- niche(dudiEC, log(ECm1 + 1), scann = FALSE)
nicMSm <- niche(dudiMS, log(MSm1 + 1), scann = FALSE)
nicPCm <- niche(dudiPC, log(PCm1 + 1), scann = FALSE)

# Tolerance and OMI Results (use column 6) ----
aBCm <- niche.param(nicBCm)
aCLm <- niche.param(nicCLm)
aDGm <- niche.param(nicDGm)
aECm <- niche.param(nicECm)
aMSm <- niche.param(nicMSm)
aPCm <- niche.param(nicPCm)

get_tol_df <- function(tol_matrix, elev) {
  tibble(
    Species   = rownames(tol_matrix),
    Tolerance = tol_matrix[, 6],
    Locality  = elev
  )
}

all_tol <- bind_rows(
  get_tol_df(aBCm, "BC"),
  get_tol_df(aCLm, "CL"),
  get_tol_df(aDGm, "DG"),
  get_tol_df(aECm, "EC"),
  get_tol_df(aMSm, "MS"),
  get_tol_df(aPCm, "PC")
)

# Abundance and Elevation Lookups ----
abund <- matcom_moths %>%
  group_by(Locality, Plot, Species_complete, Family) %>%
  summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop")

elev_lookup <- matenv %>%
  select(Locality, Plot, Real_Elevation) %>%
  distinct()

# Combine Abundances with Tolerances and Elevation ----
abund_tol <- abund %>%
  left_join(all_tol,     by = c("Species_complete" = "Species", "Locality")) %>%
  left_join(elev_lookup, by = c("Locality", "Plot")) %>%
  filter(!is.na(Tolerance))

# CWM for All Moths ----
cwm_all <- abund_tol %>%
  group_by(Locality, Plot, Real_Elevation) %>%
  summarise(
    CWM_All = weighted.mean(Tolerance, Abundance, na.rm = TRUE),
    .groups = "drop"
  )

# CWM for Target Families ----
target_families <- c("Geometridae", "Erebidae", "Noctuidae")

cwm_fams <- abund_tol %>%
  filter(Family %in% target_families) %>%
  group_by(Locality, Plot, Real_Elevation, Family) %>%
  summarise(
    CWM = weighted.mean(Tolerance, Abundance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from  = Family,
    values_from = CWM,
    names_prefix = "CWM_"
  )

# Combine and Export ----
moth_nb <- cwm_all %>%
  left_join(cwm_fams, by = c("Locality", "Plot", "Real_Elevation"))

wb <- loadWorkbook(here("Data/Processed/WENB.xlsx"))
addWorksheet(wb, "Moths_WENB")
writeData(wb, "Moths_WENB", moth_nb)
saveWorkbook(wb, file = here("Data/Processed/WENB.xlsx"), overwrite = TRUE)
