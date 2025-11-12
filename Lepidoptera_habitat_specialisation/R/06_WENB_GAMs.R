######################################################
# Lepidoptera Habitat Specialization Analysis
# GAM Modeling & Plotting
######################################################

# Load Packages ----
library(tidyverse)   # dplyr, tidyr, ggplot2, purrr, stringr, tibble
library(readxl)
library(mgcv)
library(patchwork)
library(broom)
library(here)

# Plot Theme ----
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

# Data Paths ----
xlsx <- here("Data/Processed/WENB.xlsx")


# Read CWMs ----
moths <- read_excel(xlsx, sheet = "Moths_WENB") %>%
  mutate(Locality = as.character(Locality),
         Plot     = as.character(Plot)) %>%
  rename(
    All_moths   = CWM_All,
    Erebidae    = CWM_Erebidae,
    Noctuidae   = CWM_Noctuidae,
    Geometridae = CWM_Geometridae
  ) %>%
  select(Locality, Plot, Real_Elevation, All_moths, Erebidae, Noctuidae, Geometridae)

butts <- read_excel(xlsx, sheet = "Butterflies_WENB") %>%
  mutate(Locality = as.character(Locality),
         Plot     = as.character(Plot)) %>%
  rename(
    All_butts     = CWM_All,
    Charaxinae    = CWM_Charaxinae,
    Limenitidinae = CWM_Limenitidinae,
    Satyrinae     = CWM_Satyrinae
  ) %>%
  select(Locality, Plot, Real_Elevation, All_butts, Charaxinae, Limenitidinae, Satyrinae)

wenb_merged <- full_join(
  moths, butts,
  by = c("Locality", "Plot", "Real_Elevation")
) %>%
  arrange(Locality, Plot)

wenb_simple <- wenb_merged %>%
  select(
    Plot, Real_Elevation,
    All_moths, Erebidae, Noctuidae, Geometridae,
    All_butts, Charaxinae, Limenitidinae, Satyrinae
  )

cwm_data <- wenb_simple %>%
  pivot_longer(
    cols = c(All_moths, Erebidae, Noctuidae, Geometridae,
             All_butts, Charaxinae, Limenitidinae, Satyrinae),
    names_to = "Group",
    values_to = "Value"
  )

# Modeling Functions ----
fit_gam_models <- function(data, k_value = 8) {
  data %>%
    group_by(Group) %>%
    filter(n() >= 5) %>%
    nest() %>%
    mutate(
      cleaned = map(data, ~ drop_na(.x, Value, Real_Elevation)),
      model   = map(cleaned, ~ gam(
        Value ~ s(Real_Elevation, k = k_value),
        data = .x,
        method = "REML",
        select = TRUE
      )),
      preds = map2(model, cleaned, ~ {
        pr <- predict(.x, newdata = .y, se.fit = TRUE, type = "link")
        inv <- .x$family$linkinv
        tibble(
          Real_Elevation = .y$Real_Elevation,
          Value          = .y$Value,
          fit            = inv(as.numeric(pr$fit)),
          se             = as.numeric(pr$se.fit)
        ) %>%
          mutate(
            lwr = inv(qlogis(pmin(pmax(plogis(as.numeric(pr$fit)) - 1.96 * se, .Machine$double.eps), 1 - .Machine$double.eps))),
            upr = inv(qlogis(pmin(pmax(plogis(as.numeric(pr$fit)) + 1.96 * se, .Machine$double.eps), 1 - .Machine$double.eps)))
          )
      })
    ) %>%
    select(Group, preds) %>%
    unnest(preds) %>%
    arrange(Group, Real_Elevation)
}

# Plotting ----
create_group_plot <- function(data, group_name, point_color, line_color,
                              title, y_label = "") {
  plot_data <- data %>% filter(Group == group_name)
  if (nrow(plot_data) == 0) return(ggplot() + theme_void())
  
  ggplot(plot_data, aes(x = Real_Elevation)) +
    geom_point(aes(y = Value), shape = 19, size = 4, color = point_color) +
    geom_ribbon(aes(ymin = fit - 1.96 * se, ymax = fit + 1.96 * se),
                alpha = 0.15) +
    geom_line(aes(y = fit), linewidth = 3, color = line_color) +
    labs(x = "", y = y_label, title = title) +
    scale_x_continuous(
      breaks = c(500, 1000, 1500, 2000),
      labels = c("500", "1000", "1500", "2000")
    ) +
    coord_cartesian(xlim = c(300, 2300)) +
    BioR_theme
}

color_palette <- list(
  moths = list(point = "#7280bbff", line = "#384690ff"),
  butterflies = list(point = "#ac5353ff", line = "#8b2323ff")
)

make_panel <- function(df, pal) {
  wrap_plots(
    create_group_plot(df, "All_moths",   pal$moths$point,       pal$moths$line,       "i.  All moths",        "CWM niche breadth"),
    create_group_plot(df, "Erebidae",    pal$moths$point,       pal$moths$line,       "j.  Erebidae"),
    create_group_plot(df, "Noctuidae",   pal$moths$point,       pal$moths$line,       "k.  Noctuidae"),
    create_group_plot(df, "Geometridae", pal$moths$point,       pal$moths$line,       "l.  Geometridae"),
    create_group_plot(df, "All_butts",   pal$butterflies$point, pal$butterflies$line, "m.  All butterflies",  "CWM niche breadth"),
    create_group_plot(df, "Satyrinae",   pal$butterflies$point, pal$butterflies$line, "n.  Satyrinae"),
    create_group_plot(df, "Limenitidinae", pal$butterflies$point, pal$butterflies$line, "o.  Limenitidinae"),
    create_group_plot(df, "Charaxinae",  pal$butterflies$point, pal$butterflies$line, "p.  Charaxinae"),
    nrow = 2, ncol = 4
  ) + plot_annotation(caption = "Elevation (m a.s.l.)") &
    theme(plot.caption = element_text(size = 25, hjust = 0.5, margin = margin(t = 0, b = -2)))
}

# Run Analysis ----
model_results <- fit_gam_models(cwm_data, k_value = 8)
panel_plot    <- make_panel(model_results, color_palette)
print(panel_plot)

ggsave(filename = "LOCAL_FINAL.svg", plot = panel_plot, width = 16.7, height = 9)

# Optional: model list and diagnostics ----
build_model_list <- function(data, k_value = 8) {
  tbl <- data %>%
    group_by(Group) %>%
    filter(n() >= 5) %>%
    nest() %>%
    mutate(
      model_data = map(data, ~ drop_na(.x, Value, Real_Elevation)),
      model = map(model_data, ~ gam(
        Value ~ s(Real_Elevation, k = k_value),
        data = .x,
        method = "REML",
        select = TRUE
      ))
    )
  mdl <- tbl$model
  names(mdl) <- tbl$Group
  mdl
}

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

model_list_cwm <- build_model_list(cwm_data, k_value = 8)

cwm_summary <- do.call(
  rbind,
  lapply(names(model_list_cwm), function(g) {
    out <- try(extract_gam_stats(model_list_cwm[[g]]), silent = TRUE)
    if (inherits(out, "try-error")) return(NULL)
    transform(out, Group = g)
  })
) %>% tibble::as_tibble() %>% relocate(Group)

print(cwm_summary)

# Optional checks
# gam.check(model_list_cwm[["All_moths"]])
# gam.check(model_list_cwm[["Erebidae"]])
# gam.check(model_list_cwm[["Noctuidae"]])
# gam.check(model_list_cwm[["Geometridae"]])
# gam.check(model_list_cwm[["All_butts"]])
# gam.check(model_list_cwm[["Satyrinae"]])
# gam.check(model_list_cwm[["Limenitidinae"]])
# gam.check(model_list_cwm[["Charaxinae"]])
