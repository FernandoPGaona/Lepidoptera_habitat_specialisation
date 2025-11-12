######################################################
# Lepidoptera Habitat Specialization Analysis
# GAM Modeling & Plotting
######################################################


#Load Packages ----
library(tidyverse)
library(mgcv)
library(patchwork)
library(here)
library(ggplot2)
library(svglite)
library(readxl)

#Load Packages ----
DATA_FILE <- here("Data/Processed/ANB.xlsx")
SHEET_BUTTS <- 1
SHEET_MOTHS <- 2

#Plot Theme ----

BioR_theme <- theme(
  panel.background = element_blank(),
  panel.border = element_rect(colour = "black", fill = NA, size = 1),
  panel.grid.major = element_line(color = "gray80", linewidth = 0.3),
  panel.grid.minor = element_line(color = "gray99", linewidth = 0.1),
  axis.line = element_line("black"),
  text = element_text(size = 15),
  axis.text = element_text(size = 20, colour = "black"),
  axis.title = element_text(size = 20, colour = "black"),
  legend.title = element_text(size = 16),
  legend.text = element_text(size = 16),
  legend.key = element_blank(),
  plot.title = element_text(size = 20, hjust = 0, margin = margin(b = 5)),
  axis.ticks.x = element_line(color = "black", linewidth = 1),
  axis.ticks.length.x = unit(0.3, "cm"),
  axis.ticks.y = element_line(color = "black", linewidth = 1),
  axis.ticks.length.y = unit(0.3, "cm"),
  plot.margin = margin(t = -1, r = -1, b = -1, l = -1, unit = "pt")
)

#Data Import ----

# Columns expected: Locality, Plot, Real_Elevation, Group, CWM_tolerance
load_cwm_sheet <- function(file, sheet) {
  read_excel(file, sheet = sheet) %>%
    select(Locality, Plot, Real_Elevation, Group, Value = CWM_tolerance) %>%
    mutate(
      Plot = factor(as.character(Plot)),
      Group = factor(as.character(Group))
    ) %>%
    filter(!is.na(Value), !is.na(Real_Elevation))
}

moths_df <- load_cwm_sheet(DATA_FILE, SHEET_MOTHS)
butts_df <- load_cwm_sheet(DATA_FILE, SHEET_BUTTS)
cwm_data <- bind_rows(moths_df, butts_df)


#Gam modeling ----

gamma_groups <- c("All_moths", "Erebidae", "Noctuidae", "Geometridae", "All_butts")

fit_gam_models <- function(data, k_default = 10) {
  data %>%
    group_by(Group) %>%
    filter(n() >= 5) %>%
    nest() %>%
    mutate(
      k_value = if_else(Group == "Satyrinae", 10L, k_default),
      fam     = if_else(Group %in% gamma_groups, "gamma", "gaussian"),
      model   = pmap(list(data, k_value, fam), function(df, k, fam) {
        if (fam == "gamma") {
          gam(Value ~ s(Real_Elevation, k = k),
              data = df,
              family = Gamma(link = "log"),
              method = "REML",
              select = TRUE)
        } else {
          gam(Value ~ s(Real_Elevation, k = k),
              data = df,
              family = gaussian(),
              method = "REML",
              select = TRUE)
        }
      }),
      predictions = map2(model, data, ~ {
        pr  <- predict(.x, newdata = .y, type = "link", se.fit = TRUE)
        eta <- as.numeric(pr$fit)
        se  <- as.numeric(pr$se.fit)
        inv <- .x$family$linkinv
        bind_cols(
          .y,
          tibble(
            fit = inv(eta),
            lwr = inv(eta - 1.96 * se),
            upr = inv(eta + 1.96 * se)
          )
        )
      })
    ) %>%
    unnest(predictions)
}

model_results <- fit_gam_models(cwm_data)


#Plottig ----

color_palette <- list(
  moths = list(point = "#7280bbff", line = "#384690ff"),
  butterflies = list(point = "#ac5353ff", line = "#8b2323ff")
)

create_group_plot <- function(data, group_name, point_color, line_color, title, y_label = "CWM Absolute niche breadth") {
  plot_data <- data %>% filter(Group == group_name)
  if (nrow(plot_data) == 0) return(ggplot() + theme_void() + ggtitle(paste0(title, " (no data)")))
  
  plot_data <- plot_data %>% arrange(Real_Elevation)
  
  ggplot(plot_data, aes(x = Real_Elevation)) +
    geom_point(aes(y = Value), shape = 19, size = 3.5, color = point_color) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.15) +
    geom_line(aes(y = fit), linewidth = 2.2, color = line_color) +
    labs(x = "", y = y_label, title = title) +
    scale_x_continuous(
      limits = c(300, 2300),
      breaks = c(500, 1000, 1500, 2000),
      labels = c("500", "1000", "1500", "2000")
    ) +
    coord_cartesian(xlim = c(300, 2300)) +
    BioR_theme
}

make_panel <- function(df, pal) {
  wrap_plots(
    create_group_plot(df, "All_moths",   pal$moths$point, pal$moths$line, "a. All moths"),
    create_group_plot(df, "Erebidae",    pal$moths$point, pal$moths$line, "b. Erebidae"),
    create_group_plot(df, "Noctuidae",   pal$moths$point, pal$moths$line, "c. Noctuidae"),
    create_group_plot(df, "Geometridae", pal$moths$point, pal$moths$line, "d. Geometridae"),
    create_group_plot(df, "All_butts",   pal$butterflies$point, pal$butterflies$line, "e. All butterflies"),
    create_group_plot(df, "Satyrinae",   pal$butterflies$point, pal$butterflies$line, "f. Satyrinae"),
    create_group_plot(df, "Limenitidinae", pal$butterflies$point, pal$butterflies$line, "g. Limenitidinae"),
    create_group_plot(df, "Charaxinae",  pal$butterflies$point, pal$butterflies$line, "h. Charaxinae"),
    nrow = 2, ncol = 4
  ) + plot_annotation(caption = "Elevation (m a.s.l.)") &
    theme(plot.caption = element_text(size = 25, hjust = 0.5, margin = margin(t = 0, b = -2)))
}

panel_plot <- make_panel(model_results, color_palette)
panel_plot

#Extract Model Statistics ----

extract_gam_stats <- function(gam_model) {
  s <- summary(gam_model)
  data.frame(
    Intercept = unname(s$p.coeff["(Intercept)"]),
    edf_elev  = s$s.table["s(Real_Elevation)", "edf"],
    F_elev    = s$s.table["s(Real_Elevation)", "F"],
    p_elev    = s$s.table["s(Real_Elevation)", "p-value"],
    R2_adj    = s$r.sq,
    Dev_exp   = s$dev.expl * 100,
    n         = s$n
  )
}

build_model_list <- function(data, k_default = 10) {
  tbl <- data %>%
    group_by(Group) %>%
    filter(n() >= 5) %>%
    nest() %>%
    mutate(
      k_value = if_else(Group == "Satyrinae", 7L, k_default),
      fam     = if_else(Group %in% gamma_groups, "gamma", "gaussian"),
      model   = pmap(list(data, k_value, fam), function(df, k, fam) {
        if (fam == "gamma") {
          gam(Value ~ s(Real_Elevation, k = k),
              data = df,
              family = Gamma(link = "log"),
              method = "REML",
              select = TRUE)
        } else {
          gam(Value ~ s(Real_Elevation, k = k),
              data = df,
              family = gaussian(),
              method = "REML",
              select = TRUE)
        }
      })
    )
  mdl <- tbl$model; names(mdl) <- tbl$Group; mdl
}

mdl <- build_model_list(cwm_data)

do.call(rbind, lapply(names(mdl), function(g) {
  transform(extract_gam_stats(mdl[[g]]), Group = g)
})) %>% dplyr::relocate(Group)

#######################################