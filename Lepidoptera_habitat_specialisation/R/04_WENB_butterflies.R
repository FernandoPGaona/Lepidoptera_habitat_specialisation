######################################################
# Lepidoptera Habitat Specialization Analysis
# Butterflies
# Within-Elevation niche breadth (WENB)
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
butt_mat <- read_excel("Data/Raw/Lepidoptera_data.xlsx", sheet = 2) %>%
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
butt_loc <- dcast(
  butt_mat, Locality + Plot ~ Species,
  value.var = "Abundance", fun.aggregate = sum
)

# Analyses per Locality ----
BCb <- butt_loc[(butt_loc$Locality == "BC"), -c(1, 2)]
CLb <- butt_loc[(butt_loc$Locality == "CL"), -c(1, 2)]
DGb <- butt_loc[(butt_loc$Locality == "DG"), -c(1, 2)]
ECb <- butt_loc[(butt_loc$Locality == "EC"), -c(1, 2)]
MSb <- butt_loc[(butt_loc$Locality == "MS"), -c(1, 2)]
PCb <- butt_loc[(butt_loc$Locality == "PC"), -c(1, 2)]

BCb1 <- BCb[, colSums(BCb > 0) > 0, drop = FALSE]
CLb1 <- CLb[, colSums(CLb > 0) > 0, drop = FALSE]
DGb1 <- DGb[, colSums(DGb > 0) > 0, drop = FALSE]
ECb1 <- ECb[, colSums(ECb > 0) > 0, drop = FALSE]
MSb1 <- MSb[, colSums(MSb > 0) > 0, drop = FALSE]
PCb1 <- PCb[, colSums(PCb > 0) > 0, drop = FALSE]

BCb1t <- BCb1; CLb1t <- CLb1; DGb1t <- DGb1
ECb1t <- ECb1; MSb1t <- MSb1; PCb1t <- PCb1

# Environmental Descriptors by Locality ----
BCe <- env_global[(env_global$Locality == "BC"), ]
CLe <- env_global[(env_global$Locality == "CL"), ]
DGe <- env_global[(env_global$Locality == "DG"), ]
ECe <- env_global[(env_global$Locality == "EC"), ]
MSe <- env_global[(env_global$Locality == "MS"), ]
PCe <- env_global[(env_global$Locality == "PC"), ]

# OMI Method ----
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

# Community Data on OMI Space ----
nicBCb <- niche(dudiBC, log(BCb1 + 1), scann = FALSE)
nicCLb <- niche(dudiCL, log(CLb1 + 1), scann = FALSE)
nicDGb <- niche(dudiDG, log(DGb1 + 1), scann = FALSE)
nicECb <- niche(dudiEC, log(ECb1 + 1), scann = FALSE)
nicMSb <- niche(dudiMS, log(MSb1 + 1), scann = FALSE)
nicPCb <- niche(dudiPC, log(PCb1 + 1), scann = FALSE)

# Tolerance Extraction (column 6) ----
aBCb <- niche.param(nicBCb)
aCLb <- niche.param(nicCLb)
aDGb <- niche.param(nicDGb)
aECb <- niche.param(nicECb)
aMSb <- niche.param(nicMSb)
aPCb <- niche.param(nicPCb)

get_tol_df <- function(tol_matrix, elev) {
  tibble(
    Species   = rownames(tol_matrix),
    Tolerance = tol_matrix[, 6],
    Locality  = elev
  )
}

all_tol <- bind_rows(
  get_tol_df(aBCb, "BC"),
  get_tol_df(aCLb, "CL"),
  get_tol_df(aDGb, "DG"),
  get_tol_df(aECb, "EC"),
  get_tol_df(aMSb, "MS"),
  get_tol_df(aPCb, "PC")
)

# Join Tolerance and Elevation ----
env_lookup <- matenv %>%
  select(Locality, Plot, Real_Elevation) %>%
  distinct()

butt_all <- butt_mat %>%
  filter(Abundance > 0) %>%
  left_join(all_tol,    by = c("Species", "Locality")) %>%
  left_join(env_lookup, by = c("Locality", "Plot")) %>%
  filter(!is.na(Tolerance))

# CWM Aggregation ----
cwm_all <- butt_all %>%
  group_by(Locality, Plot, Real_Elevation) %>%
  summarise(CWM_All = weighted.mean(Tolerance, Abundance), .groups = "drop")

target_subfams <- c("Satyrinae", "Limenitidinae", "Charaxinae")

cwm_subfamilies <- butt_all %>%
  filter(Subfamily %in% target_subfams) %>%
  group_by(Locality, Plot, Real_Elevation, Subfamily) %>%
  summarise(CWM = weighted.mean(Tolerance, Abundance), .groups = "drop") %>%
  pivot_wider(names_from = Subfamily, values_from = CWM, names_prefix = "CWM_")

cwm_combined <- left_join(
  cwm_all, cwm_subfamilies,
  by = c("Locality", "Plot", "Real_Elevation")
)

# Export ----
wb <- createWorkbook()
addWorksheet(wb, "Butterflies_WENB")
writeData(wb, "Butterflies_WENB", cwm_combined)
saveWorkbook(wb, here("Data/Processed/WENB.xlsx"), overwrite = TRUE)
