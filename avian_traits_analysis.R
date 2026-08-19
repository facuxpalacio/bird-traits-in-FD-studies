###### This code accompany the manuscript titled "Trait selection and use in avian functional diversity research"
###### code by Facundo Palacio and Marta Jarzyna
library(tidyr)
library(dplyr)
library(ggplot2)
library(viridis)
library(stringr)
library(tidyverse)
library(scales)
library(circlize)
library(iNEXT)

set.seed(357)

###### Load data
dat <- read.csv("./data/avian_traits_database.csv")


###### Rarefaction curves (Fig. 1)

#dat$trait_study <- iconv(dat$trait_study, from = "", to = "UTF-8") # proper encoding
dat$trait_study <- paste(dat$studyID, dat$trait_synonym, sep = "_")

# Collapse duplicated rows
dat_clean <- dat %>%
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
mat <- dat_clean %>%
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
p <- ggiNEXT(out, type = 1) +
  scale_color_viridis_d(option = "D") +
  scale_fill_viridis_d(option = "D") +
  labs(x = "Number of sampling units", y = "Trait richness",
       color = "Diversity type", fill = "Diversity type", shape = "Diversity type") +
  theme_minimal(base_size = 18) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "grey85", linewidth = 0.5),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        panel.background = element_rect(fill = "white", color = NA),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.position = "inside",
        legend.position.inside = c(0.8, 0.2),
        legend.background = element_rect(color = "black", fill = "white", linewidth = 0.5),
        legend.key.width = unit(1.5, "cm"),   # widen the key so lines are visible
        legend.key.height = unit(0.6, "cm"))
p

ggsave(plot=p, filename = "rarefaction_1.jpeg", path = "./output",
       width = 8, height = 6,  units = "in", dpi = 600)

p <- ggiNEXT(out, type = 3) +
  scale_color_viridis_d(option = "D") +
  scale_fill_viridis_d(option = "D") +
  labs(x = "Sample coverage", y = "Trait richness",
       color = "Diversity type", fill = "Diversity type", shape = "Diversity type") +
  theme_minimal(base_size = 18) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "grey85", linewidth = 0.5),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        panel.background = element_rect(fill = "white", color = NA),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.position = "inside",
        legend.position.inside = c(0.2, 0.8),
        legend.background = element_rect(color = "black", fill = "white", linewidth = 0.5),
        legend.key.width = unit(1.5, "cm"),   # widen the key so lines are visible
        legend.key.height = unit(0.6, "cm"))
p

ggsave(plot=p, filename = "rarefaction_3.jpeg", path = "./output",
       width = 8, height = 6,  units = "in", dpi = 600)


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



###### Temporal distribution of studies (Fig. 2D)

dat <- read.csv("./data/avian_traits_database.csv")
length(unique(dat$studyID))
dat_study <- dat %>%
  distinct(studyID, year, country, continent)

# there are studies that encompass many continents
# to account for this, we assign weights to each continent for such studies
# for example, if a study is only done on 1 continent, the weight is 1
# if a study spans 3 continents, each continent gets a weight of 1/3
dat_weighted <- dat_study %>%
  # Split continent text at commas
  mutate(continent = str_split(continent, ",\\s*")) %>%
  # Create one row per continent
  unnest(continent) %>%
  # Assign fractional weight: 1 / number of continents for that study
  group_by(studyID) %>%
  mutate(weight = 1 / n()) %>%
  ungroup()

unique(dat_weighted$continent)

dat_weighted$continent[dat_weighted$continent == "All"] <- "Global"

dat_weighted$continent <- factor(dat_weighted$continent,
                                 levels = c("Africa", "Asia", "Europe", "Central America",
                                            "North America", "South America", "Oceania",
                                            "Island systems", "Global"))
p <- ggplot(dat_weighted,
            aes(x = factor(year), y = weight, fill = continent)) +
  geom_col() +
  scale_fill_viridis_d(option = "D") +
  #scale_x_discrete(breaks = levels(factor(dat_weighted$year))[c(TRUE, FALSE)]) +
  scale_x_discrete(breaks = c(2003, 2008, 2010, 2012, 2014, 2016, 2018, 2020, 2022, 2024)) +
  labs(x = "Year", y = "Weighted number of studies",
       fill = "Continent") +
  theme_minimal(base_size = 18) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "grey85", linewidth = 0.5),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        panel.background = element_rect(fill = "white", color = NA),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.position = "inside",
        legend.position.inside = c(0.2, 0.62),
        legend.background = element_rect(color = "black", fill = "white", linewidth = 0.5))
p

ggsave(plot=p, filename = "studiesperyear.jpeg", path = "./output",
       width = 9, height = 6,  units = "in", dpi = 600)



###### Spatial distribution of studies, proportional contributions (Fig. 2C)

dat_pie <- dat_weighted %>%
  group_by(continent) %>%
  summarise(total_weight = sum(weight), .groups = "drop") %>%
  mutate(pct = total_weight / sum(total_weight))

p <- ggplot(dat_pie, aes(x = 2, y = total_weight, fill = continent)) +
  geom_col(width = 1) +
  coord_polar(theta = "y") +
  xlim(0.5, 2.5) +
  scale_fill_viridis_d(option = "D") +
  labs(fill = "Continent") +
  theme_void(base_size = 18) +
  theme(legend.text = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.margin = margin(8, 8, 8, 8),
        legend.background = element_rect(color = "black", fill = "white", linewidth = 0.5))
p

ggsave(plot=p, filename = "studiespercontinent.jpeg", path = "./output",
       width = 8, height = 6,  units = "in", dpi = 600)

# plot with percent labels
p <- ggplot(dat_pie, aes(x = 2, y = total_weight, fill = continent)) +
  geom_col(width = 1) +
  coord_polar(theta = "y") +
  xlim(0.5, 2.5) +
  geom_text(aes(label = scales::percent(pct, accuracy = 0.1)),
            position = position_stack(vjust = 0.5),
            color = "white", size = 18 / .pt) +
  scale_fill_viridis_d(option = "D") +
  labs(fill = "Continent") +
  theme_void(base_size = 18) +
  theme(legend.text = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.background = element_rect(color = "black", fill = "white", linewidth = 0.5))
p



###### Spatial distribution of studies, maps (Figs. 2A,B)

dat <- read.csv("./data/avian_traits_database.csv")

# there are studies that encompass many countries
# to account for this, we assign weights to each country for such studies
dat_study <- dat %>%
  distinct(studyID, year, country, continent)

# clean up
dat_clean <- dat_study[!is.na(dat_study$country), ]

# split multi-country studies and assign weights
dat_weighted <- dat_clean %>%
  # Split countries
  mutate(country = str_split(country, ",\\s*")) %>%
  unnest(country) %>%
  # Assign weight per study
  group_by(studyID) %>%
  mutate(weight = 1 / n()) %>%
  ungroup()

# calculate total weighted studies per country
country_summary <- dat_weighted %>%
  group_by(country) %>%
  summarise(weighted_studies = sum(weight), .groups = "drop")

country_summary <- country_summary %>%
  mutate(country = case_when(
    country %in% c("USA", "United States") ~ "United States of America",
    TRUE ~ country
  ))

country_summary <- country_summary %>% arrange(desc(weighted_studies))

# load a world map and join data
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

world <- ne_countries(scale = "medium", returnclass = "sf")

world_data <- world %>%
  left_join(country_summary, by = c("name" = "country"))

graticule <- sf::st_graticule(
  lat = seq(-90, 90, by = 30),
  lon = seq(-180, 180, by = 60),
  crs = "+proj=moll"
)

p <- ggplot(world_data) +
  geom_sf(data = graticule, color = "grey85", linewidth = 0.4) +
  geom_sf(aes(fill = weighted_studies), color = NA) +
  scale_fill_viridis(
    option = "plasma",
    na.value = "grey90",
    name = "Weighted number of studies",
    breaks = c(10, 30, 50, 70)) +
  labs(title = element_blank()) +
  coord_sf(crs = "+proj=moll") +
  guides(fill = guide_colorbar(direction = "horizontal",
                               title.position = "top",
                               barwidth = 15, barheight = 0.8)) +
  theme_minimal(base_size = 18) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.margin = margin(8, 8, 8, 8),
        legend.background = element_rect(color = "black", fill = "white", linewidth = 0.5),
        panel.grid = element_blank())
p
ggsave(plot=p, filename = "studiespercountry_map.jpeg", path = "./output",
       width = 8, height = 6,  units = "in", dpi = 600)


# plot localities for local-level studies
dat <- read.csv("./data/avian_traits_database.csv")
dat_study <- dat %>%
  distinct(studyID, longitude, latitude)

dat_clean <- dat_study[!is.na(dat_study$longitude), ]
dat_clean <- dat_clean[!is.na(dat_clean$latitude), ]

dat_clean$latitude  <- as.numeric(dat_clean$latitude)
dat_clean$longitude <- as.numeric(dat_clean$longitude)

dat_clean <- dat_clean[!is.na(dat_clean$longitude), ]
dat_clean <- dat_clean[!is.na(dat_clean$latitude), ]

summary(dat_clean[, c("longitude", "latitude")])

points_sf <- st_as_sf(dat_clean, coords = c("longitude", "latitude"), crs = 4326)

world <- ne_countries(scale = "medium", returnclass = "sf")

graticule <- sf::st_graticule(
  lat = seq(-90, 90, by = 30),
  lon = seq(-180, 180, by = 60),
  crs = "+proj=moll"
)

p <- ggplot() +
  geom_sf(data = graticule, color = "grey85", linewidth = 0.4) +
  geom_sf(data = world, fill = "grey90", color = "grey70", linewidth = 0.2) +
  geom_sf(data = points_sf, color = "black", size = 1.5, alpha = 0.6, shape = 16) +
  labs(title = element_blank()) +
  coord_sf(crs = "+proj=moll") +
  theme_minimal(base_size = 18) +
  theme(panel.grid = element_blank())
p


ggsave(plot=p, filename = "studies-localities_map.jpeg", path = "./output",
       width = 8, height = 6,  units = "in", dpi = 600)



###### Spatial extent of studies (Fig. 3A)

library(scales)
dat <- read.csv("./data/avian_traits_database.csv")

dat_study <- dat %>%
  distinct(studyID, spatial_extent)

nrow(dat_study)
study_count <- length(unique(dat_study$studyID))

# Check if number of rows equals number of unique studyIDs
nrow(dat_study) == length(unique(dat_study$studyID)) #there are no studies that cover two spatial scales, so there is no need for weights

extent_summary <- dat_study %>%
  count(spatial_extent) %>%
  mutate(pct = n / study_count,
         label = paste0(spatial_extent, " (", percent(pct, accuracy = 0.1), ")"))

extent_summary <- extent_summary %>%
  mutate(spatial_extent = recode(spatial_extent,
                                 "local"       = "Local",
                                 "global"      = "Global",
                                 "continental" = "Continental",
                                 "biome"       = "Biome",
                                 "regional"    = "Regional"),
         spatial_extent = factor(spatial_extent, levels = c("Local", "Regional", "Biome", "Continental", "Global")))

p <- ggplot(extent_summary, aes(x = 2, y = pct, fill = spatial_extent)) +
  geom_col(width = 1) +
  coord_polar(theta = "y") +
  xlim(0.5, 2.5) +
  scale_fill_viridis_d(option = "D") +
  labs(fill = "Spatial extent") +
  theme_void(base_size = 18) +
  theme(legend.text = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.margin = margin(8, 8, 8, 8),
        legend.background = element_rect(color = "black", fill = "white", linewidth = 0.5),
        legend.position = "none")
p

ggsave(plot=p, filename = "studiesperspatialextent.jpeg", path = "./output",
       width = 6, height = 6,  units = "in", dpi = 600)


p <- ggplot(extent_summary, aes(x = 2, y = pct, fill = label)) +
  geom_col(width = 1) +
  coord_polar(theta = "y") +
  xlim(0.5, 2.5) +
  scale_fill_viridis_d(option = "D") +
  labs(fill = "Spatial extent") +
  theme_void(base_size = 18) +
  theme(legend.text = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.margin = margin(8, 8, 8, 8),
        legend.background = element_rect(color = "black", fill = "white", linewidth = 0.5))
p


###### Temporal extent of studies (Fig. 3B)

dat <- read.csv("./data/avian_traits_database.csv")

dat_study
dat_study <- dat %>%
  distinct(studyID, temporal_scale_category)

nrow(dat_study)
study_count <- length(unique(dat_study$studyID))

# Check if number of rows equals number of unique studyIDs
nrow(dat_study) == length(unique(dat_study$studyID)) #there are some (4) studies that cover >1 temporal scale

# because this number is different, we will employ weights
temp_summary <- dat_study %>%
  # Calculate weight for each study based on number of entries
  group_by(studyID) %>%
  mutate(weight = 1 / n()) %>%
  ungroup() %>%
  # Sum weighted contributions per temporal_scale_category
  group_by(temporal_scale_category) %>%
  summarise(weighted_count = sum(weight)) %>%
  ungroup() %>%
  # Calculate percentage (weighted counts sum to total number of unique studies)
  mutate(percentage = weighted_count / sum(weighted_count) * 100) %>%
  mutate(label = paste0(temporal_scale_category, " (", round(percentage, 1), "%)"))

# if a given temporal extent comprises <3% of all studies, they become a category "other"
temp2_summary <- temp_summary %>%
  mutate(temporal_scale_category2 = if_else(percentage < 3, "Other", temporal_scale_category)) |>
  summarise(weighted_count = sum(weighted_count),
            percentage = sum(percentage),
            .by = temporal_scale_category2)

temp2_summary <- temp2_summary %>%
  mutate(temporal_scale_category2 = recode(temporal_scale_category2,
                                           "unknown" = "Unknown",
                                           "average of multiple years" = "Multiple-year average"),
         temporal_scale_category2 = factor(temporal_scale_category2, levels = c("1 season", "2-10 seasons", "1 year", "Multiple-year average", "Other", "Unknown")))

temp2_summary

p <- ggplot(temp2_summary, aes(x = 2, y = percentage, fill = temporal_scale_category2)) +
  geom_col(width = 1) +
  coord_polar(theta = "y") +
  xlim(0.5, 2.5) +
  scale_fill_viridis_d(option = "D") +
  labs(fill = "Temporal extent") +
  theme_void(base_size = 18) +
  theme(legend.text = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.margin = margin(8, 8, 8, 8),
        legend.background = element_rect(color = "black", fill = "white", linewidth = 0.5),
        legend.position = "none")
p

ggsave(plot=p, filename = "studiespertemporalextent.jpeg", path = "./output",
       width = 6, height = 6,  units = "in", dpi = 600)



###### Correlations between spatial and temporal extent of studies (Fig. 3C)

dat <- read.csv("./data/avian_traits_database.csv")
dat_study <- dat %>%
  distinct(studyID, spatial_extent, temporal_scale_category)

# if a given temporal extent comprises <3% of all studies, they become a category "other"
categs <- unique(temp2_summary$temporal_scale_category2)
dat_study <- dat_study %>%
  mutate(temporal_scale_category = recode(temporal_scale_category,
                                          "unknown" = "Unknown",
                                          "average of multiple years" = "Multiple-year average"))
dat_study2 <- dat_study |>
  mutate(temporal_scale_category = ifelse(temporal_scale_category %in% categs, temporal_scale_category, "Other"))


library(circlize)
# clean and prep
dat_chord <- dat_study2 %>%
  filter(
    !is.na(spatial_extent),         !is.na(temporal_scale_category),
    spatial_extent != "",            temporal_scale_category != ""
  ) %>%
  select(studyID, spatial_extent, temporal_scale_category) %>%
  distinct()

# count co-occurrences between spatial and temporal extent categories
pair_counts <- dat_chord %>%
  count(spatial_extent, temporal_scale_category, sort = TRUE)

# pivot to matrix (rows = spatial, cols = temporal)
mat <- pair_counts %>%
  pivot_wider(
    names_from  = temporal_scale_category,
    values_from = n,
    values_fill = 0
  ) %>%
  column_to_rownames("spatial_extent") %>%
  as.matrix()

# test if the associations are random (null H) or not (alternative H)
chisq_result <- chisq.test(mat) #chi-square
print(chisq_result)

fisher_result <- fisher.test(mat, simulate.p.value = TRUE, B = 9999) #fisher's
print(fisher_result)

cramers_v <- function(x) {
  ct <- chisq.test(x)
  sqrt(ct$statistic / (sum(x) * (min(dim(x)) - 1)))
}

obs_v <- cramers_v(mat) # Cramér's V (effect size for chi-square)
cat("Observed Cramér's V:", round(obs_v, 4), "\n")

# permutation distribution
set.seed(42)
n_perm <- 9999

perm_v <- replicate(n_perm, {
  perm_dat <- dat_chord %>%
    mutate(temporal_scale_category = sample(temporal_scale_category))  # shuffle labels

  perm_mat <- perm_dat %>%
    count(spatial_extent, temporal_scale_category) %>%
    pivot_wider(names_from = temporal_scale_category, values_from = n, values_fill = 0) %>%
    column_to_rownames("spatial_extent") %>%
    as.matrix()

  # Align columns to original matrix (some may be missing after shuffle)
  perm_mat <- perm_mat[rownames(mat),
                       intersect(colnames(mat), colnames(perm_mat)), drop = FALSE]
  cramers_v(perm_mat)
})

p_perm <- mean(perm_v >= obs_v)
cat("Permutation p-value:", p_perm, "\n")

spatial_cats  <- rownames(mat)
temporal_cats <- colnames(mat)
all_cats      <- c(spatial_cats, temporal_cats)

# spatial = warm tones, temporal = cool tones
grid_colors <- setNames(
  c(rev(RColorBrewer::brewer.pal(max(3, length(spatial_cats)),  "Oranges")[1:length(spatial_cats)]),
    rev(RColorBrewer::brewer.pal(max(3, length(temporal_cats)), "Blues")[1:length(temporal_cats)])),
  all_cats
)

# plot the chord diagram
png("./output/spatial_temporal_chord.png", width = 1400, height = 1400, res = 300)

circos.clear()
circos.par(gap.after = c(rep(3, length(spatial_cats) - 1),  15,   # gap after last spatial
                         rep(3, length(temporal_cats) - 1), 15))  # gap after last temporal

chordDiagram(
  mat,
  grid.col        = grid_colors,
  transparency    = 0.3,
  annotationTrack = "grid",
  preAllocateTracks = list(track.height = 0.12)
)

# Add sector labels
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  sector <- get.cell.meta.data("sector.index")
  circos.text(
    x         = get.cell.meta.data("xcenter"),
    y         = get.cell.meta.data("ycenter"),
    labels    = sector,
    facing    = "clockwise",
    niceFacing = TRUE,
    cex       = 0.75,
    font      = 2
  )
}, bg.border = NA)

# Add a legend to distinguish spatial vs temporal sectors
legend("bottomleft",
       legend = c(spatial_cats, temporal_cats),
       fill   = grid_colors,
       border = NA,
       title  = NULL,
       ncol   = 2,
       bty    = "n",
       cex    = 0.65)

dev.off()

# plot with no labels
png("./output/spatial_temporal_chord_TBU.png", width = 1400, height = 1400, res = 300)
circos.clear()
circos.par(gap.after = c(rep(3, length(spatial_cats) - 1),  15,
                         rep(3, length(temporal_cats) - 1), 15))
chordDiagram(
  mat,
  grid.col        = grid_colors,
  transparency    = 0.3,
  annotationTrack = "grid"
)
dev.off()



###### Research topic (Fig. 3D)

dat <- read.csv("./data/avian_traits_database.csv")
dat_study <- dat %>%
  distinct(studyID, year, country, continent, topic_covered)

unique(dat_study$topic_covered)

topic_summary <- dat_study %>%
  count(topic_covered) %>%
  mutate(pct = n / sum(n),
         label = paste0(topic_covered, " (", percent(pct, accuracy = 0.1), ")"))

p <- ggplot(topic_summary, aes(x = "", y = pct, fill = label)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  scale_fill_viridis_d(option = "D", name = "Topic") +
  labs(title = "Percentage of studies by topic covered") +
  theme_void()


# there are many research foci
# we will lump some of them together for better visualization
dat_study <- dat %>%
  distinct(studyID, year, country, continent, topic_covered)

unique(dat_study$topic_covered)

keep_categories <- c("land-use change", "climate change", "management",
                     "temporal drivers", "environmental drivers", "community assembly",
                     "ecosystem functioning", "conservation", "methodological")
dat_study <- dat_study %>%
  mutate(topic_covered_lumped = if_else(topic_covered %in% keep_categories,
                                        topic_covered,
                                        "other"))
topic_summary <- dat_study %>%
  count(topic_covered_lumped) %>%
  mutate(pct = n / sum(n),
         label = paste0(topic_covered_lumped, " (", percent(pct, accuracy = 0.1), ")"))

# data prep
topic_summary_plot <- topic_summary %>%
  mutate(topic_covered_lumped = reorder(topic_covered_lumped, -pct))

topic_summary_plot <- topic_summary_plot %>%
  mutate(topic_covered_lumped = recode(topic_covered_lumped,
                                       "climate change"       = "Climate change",
                                       "community assembly"      = "Community assembly",
                                       "conservation" = "Conservation",
                                       "ecosystem functioning"       = "Ecosystem functioning",
                                       "environmental drivers"    = "Environmental drivers",
                                       "land-use change" = "Land use change",
                                       "management" = "Management",
                                       "methodological" = "Methodology",
                                       "other" = "Other",
                                       "temporal drivers" = "Temporal drivers"))


# plot
p <- ggplot(topic_summary_plot, aes(x = pct, y = topic_covered_lumped)) +
  geom_col(fill = "cadetblue4", color = "black") +
  geom_text(aes(label = paste0(n, " (", round(pct * 100, 1), "%)")),
            hjust = -0.1, size = 4.5) +
  scale_x_continuous(limits = c(0, max(topic_summary_plot$pct) * 1.4),
                     expand = c(0, 0),
                     labels = scales::percent) +
  labs(x = "Percentage of studies (%)", y = "Research topic",
       title = element_blank()) +
  theme_minimal(base_size = 18) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "grey85", linewidth = 0.5),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        panel.background = element_rect(fill = "white", color = NA),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.position = "inside",
        legend.position.inside = c(0.2, 0.62),
        legend.background = element_rect(color = "black", fill = "white", linewidth = 0.5))
p

ggsave(plot=p, filename = "studies-per-topic.jpeg", path = "./output",
       width = 8, height = 6,  units = "in", dpi = 600)



###### Taxonomic grouping (Fig. S3)

dat <- read.csv("./data/avian_traits_database.csv")
dat_study <- dat %>%
  distinct(studyID, group)

group_summary <- dat_study %>%
  count(group) %>%
  mutate(pct = n / sum(n),
         label = paste0(group, " (", percent(pct, accuracy = 0.001), ")"))

# there are too many groups to get anything informative
# so we will group them into the following categories (see paper for details):
# phylogeny (e.g. hummingbirds, raptors),
# morphology (e.g. small birds),
# function (e.g. cavity nesters, pollinators),
# habitat affinity (e.g. waterbirds, forest-dependent species),
# abundance (e.g. common species),
# distribution (e.g. island birds),
# phenology (e.g. breeding birds),
# conservation status (e.g. reintroduced birds, traded species),
# behavior (e.g. mixed-species flocks)

groups_uni <- unique(dat_study$group)

dat_study <- dat_study %>%
  mutate(group_cat = recode(group,
                            "all birds" = "All birds",
                            "frugivores" = "Function",
                            "waterbirds" = "Habitat affinity",
                            "Passeriformes" = "Phylogeny",
                            "breeding birds" = "Phenology",
                            "threatened birds" = "Conservation status",
                            "Passeriformes and Picidae" = "Phylogeny",
                            "seabirds" = "Habitat affinity",
                            "landbirds" = "Habitat affinity",
                            "riparian birds" = "Habitat affinity",
                            "frugivores and insectivores" = "Function",
                            "Psittaciformes" = "Phylogeny",
                            "seed predators" = "Function",
                            "insectivores" = "Function",
                            "Anatidae" = "Phylogeny",
                            "mixed-species flocks" = "Behavior",
                            "breeding resident birds" = "Phenology",
                            "Tyrannida and Furnariida" = "Phylogeny",
                            "forest-breeding birds" = "Habitat affinity, Phenology",
                            "most frequent birds" = "Abundance",
                            "ground-nesting birds" = "Function",
                            "scavengers" = "Function",
                            "breeding forest birds" = "Habitat affinity, Phenology",
                            "forest-dependent birds" = "Habitat affinity",
                            "forest Passeriformes" = "Habitat affinity, Phylogeny",
                            "nonbreeding birds" = "Phenology",
                            "understory birds" = "Function",
                            "traded species" = "Conservation status",
                            "cavity nesters" = "Function",
                            "native breeding passerines" = "Phenology, Phylogeny",
                            "dominant species" = "Abundance",
                            "Falconiformes" = "Phylogeny",
                            "proxy and represented species" = "Other",
                            "Trochilidae" = "Phylogeny",
                            "Tyrannidae, Vireonidae, Parulidae "= "Phylogeny",
                            "island birds" = "Distribution",
                            "waders" = "Habitat affinity",
                            "pollinators" = "Function",
                            "breeding insectivores" = "Phenology, Function",
                            "four selected species" = "Other",
                            "resident birds" = "Phenology",
                            "lowland forest birds" = "Habitat affinity",
                            "wintering birds" = "Phenology",
                            "small birds" = "Morphology",
                            "Passeriformes, Columbiformes and Psittaciformes" = "Phylogeny",
                            "granivores" = "Function",
                            "farmland birds" = "Habitat affinity",
                            "reintroduced birds" = "Conservation status",
                            "Galliformes" = "Phylogeny",
                            "endangered birds" = "Conservation status",
                            "wetland birds" = "Habitat affinity",
                            "Phasianidae" = "Phylogeny",
                            "breeding Anatidae" =  "Phenology, Phylogeny",
                            "resident Passeriformes" = "Phenology, Phylogeny",
                            "Timaliidae, Pycnonotidae, Sylviidae" = "Phylogeny",
                            "Tyrannidae, Vireonidae, Parulidae" = "Phylogeny"
  ))

dat_study_expanded <- dat_study %>%
  mutate(group_cat_sep = strsplit(group_cat, ",\\s*")) %>%
  unnest(group_cat_sep)

group_summary <- dat_study_expanded %>%
  count(group_cat_sep) %>%
  mutate(pct = n / sum(n),
         label = paste0(group_cat_sep, " (", percent(pct, accuracy = 0.001), ")"))

# plot
p <- ggplot(group_summary, aes(x = pct, y = fct_reorder(group_cat_sep, -pct))) +
  geom_col(fill = "cadetblue4", color = "black") +
  geom_text(aes(label = paste0(n, " (", round(pct * 100, 1), "%)")),
            hjust = -0.1, size = 4.5) +
  scale_x_continuous(limits = c(0, max(group_summary$pct) * 1.4),
                     expand = c(0, 0),
                     labels = scales::percent) +
  labs(x = "Percentage of studies (%)", y = "Taxonomic grouping",
       title = element_blank()) +
  theme_minimal(base_size = 18) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "grey85", linewidth = 0.5),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        panel.background = element_rect(fill = "white", color = NA),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.position = "inside",
        legend.position.inside = c(0.2, 0.62),
        legend.background = element_rect(color = "black", fill = "white", linewidth = 0.5))
p

ggsave(plot=p, filename = "taxonomicgroups_categories.jpeg", path = "./output",
       width = 8, height = 8,  units = "in", dpi = 600)



###### Species numbers (Fig. S2)

dat <- read.csv("./data/avian_traits_database.csv")
dat_study <- dat %>%
  distinct(studyID, year, country, continent, nspecies)

mean(dat_study$nspecies, na.rm=TRUE)
median(dat_study$nspecies, na.rm=TRUE)
quantile(dat_study$nspecies, probs = c(0.25, 0.5, 0.75), na.rm=TRUE)
sd(dat_study$nspecies, na.rm=TRUE)

# plot on a log scale
p <- ggplot(dat_study, aes(x = nspecies)) +
  geom_histogram(bins = 30, fill = "mediumorchid4", color = "black") +
  scale_x_log10() +
  labs(
    x = "Number of species",
    y = "Number of studies",
    title = element_blank()) +
  theme_minimal(base_size = 18) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "grey85", linewidth = 0.5),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        panel.background = element_rect(fill = "white", color = NA),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.position = "inside",
        legend.position.inside = c(0.2, 0.62),
        legend.background = element_rect(color = "black", fill = "white", linewidth = 0.5))
p

ggsave(plot=p, filename = "nspecies_logscale.jpeg", path = "./output",
       width = 8, height = 6,  units = "in", dpi = 600)



###### Analysis focused on traits (Fig. 4A-F, Fig. S4)

dat <- read.csv("./data/avian_traits_database.csv")

### number of traits used per study (Fig. 4A)
study_trait_counts <- dat %>%
  group_by(studyID) %>%
  summarise(n_traits = n_distinct(trait_synonym), .groups = "drop")

trait_summary <- study_trait_counts %>%
  count(n_traits)
median(study_trait_counts$n_traits)
mean(study_trait_counts$n_traits)


# binned
dat <- read.csv("./data/avian_traits_database.csv")
study_trait_counts <- dat %>%
  group_by(studyID) %>%
  summarise(n_traits = n_distinct(trait_synonym), .groups = "drop")

trait_summary <- study_trait_counts %>%
  count(n_traits)
median(study_trait_counts$n_traits)
mean(study_trait_counts$n_traits)
quantile(study_trait_counts$n_traits, probs = c(0.25, 0.75))

study_trait_counts <- study_trait_counts %>%
  mutate(trait_bin = cut(
    n_traits,
    breaks = c(0, 1,2,3,4,5,6,7,8,9,10,15,35),
    labels = c("1", "2", "3", "4", "5", "6", "7", "8","9","10","11-15",">15"),
    right = TRUE
  ))

trait_summary <- study_trait_counts %>%
  count(trait_bin)

# plot
p <- ggplot(trait_summary, aes(x = trait_bin, y = n)) +
  geom_col(fill = "mediumorchid4", color = "black") +
  labs(
    x = "Number of traits used",
    y = "Number of studies",
    title = element_blank()) +
  theme_minimal(base_size = 18) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "grey85", linewidth = 0.5),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        panel.background = element_rect(fill = "white", color = NA),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.position = "inside",
        legend.position.inside = c(0.2, 0.62),
        legend.background = element_rect(color = "black", fill = "white", linewidth = 0.5))
p

ggsave(plot=p, filename = "number-traits-per-study_bins.jpeg", path = "./output",
       width = 8, height = 6,  units = "in", dpi = 600)


### the relationship between no of traits and no of species (Fig. S4)
dat <- read.csv("./data/avian_traits_database.csv")
dat_study <- dat %>%
  distinct(studyID, nspecies, trait)

dat_study <- dat_study %>%
  group_by(studyID) %>%
  mutate(ntraits = n()) %>%
  ungroup()

dat_study <- dat_study %>%
  distinct(studyID, nspecies, ntraits)
dat_study <- dat_study %>%
  filter(!is.na(nspecies), !is.na(ntraits))

cor(dat_study$nspecies, dat_study$ntraits)#r = -0.031

# plot
p <- ggplot(dat_study, aes(x = nspecies, y = ntraits)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1.2,
              fill = "grey90", alpha = 0.3) +
  labs(x = "Number of species", y = "Number of traits") +
  guides(linetype = "none") +
  theme_minimal(base_size = 18) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "grey85", linewidth = 0.5),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        panel.background = element_rect(fill = "white", color = NA),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.position = "inside",
        legend.position.inside = c(0.2, 0.62),
        legend.background = element_rect(color = "black", fill = "white", linewidth = 0.5))
p

ggsave(plot=p, filename = "nspecies-ntraits.jpeg", path = "./output",
       width = 8, height = 6,  units = "in", dpi = 600)


### trait source (Fig. 4B)
dat <- read.csv("./data/avian_traits_database.csv")
dat_study <- dat %>%
  distinct(studyID, trait_source)

dim(dat_study)
study_count <- length(unique(dat_study$studyID))

# in some cases, more than one trait source is listed per study; let's make separate rows for those cases
dat_study_expanded <- dat_study %>%
  separate_rows(trait_source, sep = ",\\s*") %>%
  mutate(trait_source = trimws(trait_source))
dim(dat_study_expanded)

studies_no <- length(unique(dat_study_expanded$studyID))

dat_study_expanded <- dat_study_expanded %>%
  mutate(trait_source = if_else(trait_source == "", "Unknown", trait_source)) %>%
  mutate(trait_source = if_else(is.na(trait_source), "Unknown", trait_source))

study_counts <- dat_study_expanded %>%
  group_by(studyID) %>%
  summarise(n_traits = n())
distribution <- study_counts %>%
  count(n_traits, name = "n_studies") %>%
  rename(n_trait_sources = n_traits) %>%
  mutate(percent = round(100 * n_studies / sum(n_studies), 1))
distribution

trait_source_summary <- dat_study_expanded %>%
  select(studyID, trait_source) %>%
  distinct() %>%
  group_by(trait_source) %>%
  summarise(n_studies = n()) %>%
  ungroup() %>%
  mutate(percentage = 100 * n_studies / studies_no)

# if trait sources occurs in <3% of studies, it becomes "other"
trait_source_summary2 <- trait_source_summary %>%
  mutate(trait_source2 = if_else(percentage < 3, "Other", trait_source)) |>
  summarise(n_studies = sum(n_studies),
            percentage = sum(percentage),
            .by = trait_source2)

trait_source_summary2 <- trait_source_summary2 %>%
  mutate(trait_source2 = recode(trait_source2,
                                "database"       = "Database",
                                "field measurements"      = "Field measurements",
                                "field observations" = "Field observations",
                                "literature"       = "Literature",
                                "museum measurements"    = "Museum measurements",
                                "website" = "Website",
                                "habitat" = "Habitat"))


# plot
p <- ggplot(trait_source_summary2, aes(x = reorder(trait_source2, -percentage), y = percentage)) +
  geom_col(fill = "cadetblue4", color = "black") +
  geom_text(aes(label = paste0(n_studies, " (", round(percentage, 1), "%)")),
            hjust = -0.1, size = 4.5) +
  scale_y_continuous(limits = c(0, 120),
                     breaks = c(0, 25, 50, 75, 100),
                     expand = c(0, 0),
                     labels = scales::label_percent(scale = 1)) +
  labs(y = "Percentage of studies (%)", x = "Trait source",
       title = element_blank()) +
  coord_flip() +
  theme_minimal(base_size = 18) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "grey85", linewidth = 0.5),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        panel.background = element_rect(fill = "white", color = NA),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.position = "inside",
        legend.position.inside = c(0.2, 0.62),
        legend.background = element_rect(color = "black", fill = "white", linewidth = 0.5))
p

ggsave(plot=p, filename = "perc-studies-traits-source.jpeg", path = "./output",
       width = 8, height = 6,  units = "in", dpi = 600)



### trait source - spatial extent chord diagram (Fig. S5A)
dat <- read.csv("./data/avian_traits_database.csv")
dat_study <- dat %>%
  distinct(studyID, trait_source, spatial_extent)

# in some cases, more than one trait source is listed per study; let's make separate rows for those cases
dat_study_expanded <- dat_study %>%
  separate_rows(trait_source, sep = ",\\s*") %>%
  mutate(trait_source = trimws(trait_source))
dim(dat_study_expanded)

studies_no <- length(unique(dat_study_expanded$studyID))

dat_study_expanded <- dat_study_expanded %>%
  mutate(trait_source = if_else(trait_source == "", "Unknown", trait_source)) %>%
  mutate(trait_source = if_else(is.na(trait_source), "Unknown", trait_source))

# some trait course comprises less than 3% - make those Other category
dat_study_expanded <- dat_study_expanded %>%
  mutate(trait_source = recode(trait_source,
                               "expert assessment" = "Other",
                               "mobile application" = "Other",
                               "DVD-ROM" = "Other",
                               "data analysis"= "Other"))

# clean and prep
dat_chord2 <- dat_study_expanded %>%
  filter(
    !is.na(trait_source),    !is.na(spatial_extent),
    trait_source != "",       spatial_extent != ""
  ) %>%
  select(studyID, trait_source, spatial_extent) %>%
  distinct()

# count co-occurrences
pair_counts2 <- dat_chord2 %>%
  count(spatial_extent, trait_source, sort = TRUE)

# pivot to matrix (rows = spatial, cols = trait_source)
mat2 <- pair_counts2 %>%
  pivot_wider(
    names_from  = trait_source,
    values_from = n,
    values_fill = 0
  ) %>%
  column_to_rownames("spatial_extent") %>%
  as.matrix()

# test for random associations
chisq_result <- chisq.test(mat2)
print(chisq_result) #chi-suqre
fisher_result <- fisher.test(mat2, simulate.p.value = TRUE, B = 9999) #fisher's
print(fisher_result)
cramers_v <- function(x) {
  ct <- chisq.test(x)
  sqrt(ct$statistic / (sum(x) * (min(dim(x)) - 1)))
}
obs_v <- cramers_v(mat2) # Caremr's V
cat("Observed Cramér's V:", round(obs_v, 4), "\n")

spatial_cats2      <- rownames(mat2)
trait_source_cats  <- colnames(mat2)
all_cats2          <- c(spatial_cats2, trait_source_cats)

# spatial = warm tones, trait_source = cool tones
grid_colors2 <- setNames(
  c(rev(colorRampPalette(RColorBrewer::brewer.pal(9, "Oranges"))(length(spatial_cats2))),
    rev(colorRampPalette(RColorBrewer::brewer.pal(9, "Blues"))(length(trait_source_cats)))),
  all_cats2
)

# plot with labels
png("./output/spatial_traitsource_chord.png", width = 1400, height = 1400, res = 300)
circos.clear()
circos.par(gap.after = c(rep(3, length(spatial_cats2) - 1),     15,
                         rep(3, length(trait_source_cats) - 1), 15))
chordDiagram(
  mat2,
  grid.col          = grid_colors2,
  transparency      = 0.3,
  annotationTrack   = "grid",
  preAllocateTracks = list(track.height = 0.12)
)

circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  sector <- get.cell.meta.data("sector.index")
  circos.text(
    x          = get.cell.meta.data("xcenter"),
    y          = get.cell.meta.data("ycenter"),
    labels     = sector,
    facing     = "clockwise",
    niceFacing = TRUE,
    cex        = 0.75,
    font       = 2
  )
}, bg.border = NA)

legend("bottomleft",
       legend = c(spatial_cats2, trait_source_cats),
       fill   = grid_colors2,
       border = NA,
       title  = NULL,
       ncol   = 2,
       bty    = "n",
       cex    = 0.65)
dev.off()

# plot without labels
png("./output/spatial_traitsource_chord.png", width = 1400, height = 1400, res = 300)
circos.clear()
circos.par(gap.after = c(rep(3, length(spatial_cats2) - 1),     15,
                         rep(3, length(trait_source_cats) - 1), 15))
chordDiagram(
  mat2,
  grid.col        = grid_colors2,
  transparency    = 0.3,
  annotationTrack = "grid"
)
dev.off()


### trait source - temporal extent chord diagram (Fig. S5B)
dat_study <- dat %>%
  distinct(studyID, spatial_extent, temporal_scale_category)

# for temporal scale if there are <3% of studies, they become 'other'
dat_study <- dat_study %>%
  mutate(temporal_scale_category = recode(temporal_scale_category,
                                          "unknown" = "Unknown_temp",
                                          "average of multiple years" = "Multiple-year average"))
categs <- c("1 season", "1 year", "Unknown_temp", "Multiple-year average", "2-10 seasons")

dat_study2 <- dat_study |>
  mutate(temporal_scale_category = ifelse(temporal_scale_category %in% categs, temporal_scale_category, "Other_temp"))

# in some cases, more than one trait source is listed per study; let's make separate rows for those cases
dat_study_expanded <- dat_study2 %>%
  separate_rows(trait_source, sep = ",\\s*") %>%
  mutate(trait_source = trimws(trait_source))
dim(dat_study_expanded)

studies_no <- length(unique(dat_study_expanded$studyID))

dat_study_expanded <- dat_study_expanded %>%
  mutate(trait_source = if_else(trait_source == "", "Unknown", trait_source)) %>%
  mutate(trait_source = if_else(is.na(trait_source), "Unknown", trait_source))

# some trait sources comprises less than 3% - make those Other category
dat_study_expanded <- dat_study_expanded %>%
  mutate(trait_source = recode(trait_source,
                               "expert assessment" = "Other",
                               "mobile application" = "Other",
                               "DVD-ROM" = "Other",
                               "data analysis"= "Other"))

# clean and prep
dat_chord2 <- dat_study_expanded %>%
  filter(
    !is.na(trait_source),    !is.na(temporal_scale_category),
    trait_source != "",       temporal_scale_category != ""
  ) %>%
  select(studyID, trait_source, temporal_scale_category) %>%
  distinct()

# count co-occurrences
pair_counts2 <- dat_chord2 %>%
  count(temporal_scale_category, trait_source, sort = TRUE)

# pivot to matrix (rows = spatial, cols = trait_source)
mat2 <- pair_counts2 %>%
  pivot_wider(
    names_from  = trait_source,
    values_from = n,
    values_fill = 0
  ) %>%
  column_to_rownames("temporal_scale_category") %>%
  as.matrix()

# test for random associations
chisq_result <- chisq.test(mat2)
print(chisq_result) #chi-square
fisher_result <- fisher.test(mat2, simulate.p.value = TRUE, B = 9999) #fisher's
print(fisher_result)
cramers_v <- function(x) {
  ct <- chisq.test(x)
  sqrt(ct$statistic / (sum(x) * (min(dim(x)) - 1)))
}
obs_v <- cramers_v(mat2) #Cramer's V
cat("Observed Cramér's V:", round(obs_v, 4), "\n")

temporal_cats2      <- rownames(mat2)
trait_source_cats  <- colnames(mat2)
all_cats2          <- c(temporal_cats2, trait_source_cats)

# Spatial = warm tones, trait_source = cool tones
grid_colors2 <- setNames(
  c(rev(colorRampPalette(RColorBrewer::brewer.pal(9, "Oranges"))(length(temporal_cats2))),
    rev(colorRampPalette(RColorBrewer::brewer.pal(9, "Blues"))(length(trait_source_cats)))),
  all_cats2
)

# plot with labels
png("./output/temporal_traitsource_chord.png", width = 1400, height = 1400, res = 300)
circos.clear()
circos.par(gap.after = c(rep(3, length(temporal_cats2) - 1),     15,
                         rep(3, length(trait_source_cats) - 1), 15))
chordDiagram(
  mat2,
  grid.col          = grid_colors2,
  transparency      = 0.3,
  annotationTrack   = "grid",
  preAllocateTracks = list(track.height = 0.12)
)

circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  sector <- get.cell.meta.data("sector.index")
  circos.text(
    x          = get.cell.meta.data("xcenter"),
    y          = get.cell.meta.data("ycenter"),
    labels     = sector,
    facing     = "clockwise",
    niceFacing = TRUE,
    cex        = 0.75,
    font       = 2
  )
}, bg.border = NA)

legend("bottomleft",
       legend = c(temporal_cats2, trait_source_cats),
       fill   = grid_colors2,
       border = NA,
       title  = NULL,
       ncol   = 2,
       bty    = "n",
       cex    = 0.65)
dev.off()

# plot without labels
png("./output/temporal_traitsource_chord.png", width = 1400, height = 1400, res = 300)
circos.clear()
circos.par(gap.after = c(rep(3, length(temporal_cats2) - 1),     15,
                         rep(3, length(trait_source_cats) - 1), 15))
chordDiagram(
  mat2,
  grid.col        = grid_colors2,
  transparency    = 0.3,
  annotationTrack = "grid"
)
dev.off()


### which traits are used (Fig. 4C,D)
dat <- read.csv("./data/avian_traits_database.csv")

# percent studies using a given trait
trait_summary <- dat %>%
  select(studyID, trait_synonym) %>%
  distinct() %>%
  group_by(trait_synonym) %>%
  summarise(n_studies = n()) %>%
  ungroup() %>%
  mutate(percentage = 100 * n_studies / sum(n_studies))

library(forcats)

# trait_level 1
# there are too many traits in trait_synonym field to be informative
# Instead look at trait_level1
dat1 <- dat %>%
  mutate(trait_level1 = if_else(trait_level1 == "", "unknown", trait_level1)) %>%
  mutate(trait_level1 = if_else(is.na(trait_level1), "unknown", trait_level1))

studies_no <- length(unique(dat1$studyID))

trait1_summary <- dat1 %>%
  select(studyID, trait_level1) %>%
  distinct() %>%
  group_by(trait_level1) %>%
  summarise(n_studies = n()) %>%
  ungroup() %>%
  mutate(percentage = 100 * n_studies / studies_no)

sum(trait1_summary$percentage)

trait1_summary <- trait1_summary %>%
  mutate(trait_level1 = recode(trait_level1,
                               "anatomy"       = "Anatomy",
                               "physiology"      = "Physiology",
                               "unknown" = "Unknown",
                               "behavior_general"       = "Behavior",
                               "activity"    = "Activity",
                               "dispersal" = "Dispersal",
                               "habitat" = "Habitat",
                               "life_history" = "Life history",
                               "trophic" = "Trophic",
                               "morphology" = "Morphology"))

p <- ggplot(trait1_summary, aes(x = reorder(trait_level1, -percentage), y = percentage)) +
  geom_col(fill = "cadetblue4", color = "black") +
  geom_text(aes(label = paste0(n_studies, " (", round(percentage, 1), "%)")),
            hjust = -0.1, size = 4.5) +
  scale_y_continuous(limits = c(0, 120),
                     breaks = c(0, 25, 50, 75, 100),
                     expand = c(0, 0),
                     labels = scales::label_percent(scale = 1)) +
  labs(y = "Percentage of studies (%)", x = "Trait category: Level 1",
       title = element_blank()) +
  coord_flip() +
  theme_minimal(base_size = 18) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "grey85", linewidth = 0.5),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        panel.background = element_rect(fill = "white", color = NA),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.position = "inside",
        legend.position.inside = c(0.2, 0.62),
        legend.background = element_rect(color = "black", fill = "white", linewidth = 0.5))
p

ggsave(plot=p, filename = "perc-studies-traits-level1.jpeg", path = "./output",
       width = 8, height = 6,  units = "in", dpi = 600)


# trait_level 2
dat1 <- dat %>%
  mutate(trait_level2 = if_else(trait_level2 == "", "unknown", trait_level2)) %>%
  mutate(trait_level2 = if_else(is.na(trait_level2), "unknown", trait_level2))

trait2_summary <- dat1 %>%
  select(studyID, trait_level2) %>%
  distinct() %>%
  group_by(trait_level2) %>%
  summarise(n_studies = n()) %>%
  ungroup() %>%
  mutate(percentage = 100 * n_studies / studies_no)

# let's group some categories (if they cover less than 5% of studies)
trait2_summary <- trait2_summary %>%
  mutate(trait_level2 = if_else(percentage < 5, "Other", trait_level2)) |>
  summarise(n_studies = sum(n_studies),
            percentage = sum(percentage),
            .by = trait_level2)

trait2_summary <- trait2_summary %>%
  mutate(trait_level2 = recode(trait_level2,
                               "lifespan"       = "Lifespan",
                               "sociality"      = "Sociality",
                               "unknown" = "Unknown",
                               "habitat_breadth"       = "Habitat breadth",
                               "activity_time"    = "Activity time",
                               "Other" = "Other",
                               "tail" = "Tail",
                               "tarsus" = "Tarsus",
                               "habitat_type" = "Habitat type",
                               "migration" = "Migration",
                               "bill" = "Bill",
                               "wing" = "Wing",
                               "reproduction" = "Reproduction",
                               "foraging" = "Foraging",
                               "diet" = "Diet",
                               "body" = "Body mass"))

p <- ggplot(trait2_summary, aes(x = reorder(trait_level2, -percentage), y = percentage)) +
  geom_col(fill = "cadetblue4", color = "black") +
  geom_text(aes(label = paste0(n_studies, " (", round(percentage, 1), "%)")),
            hjust = -0.1, size = 4.5) +
  scale_y_continuous(limits = c(0, 120),
                     breaks = c(0, 25, 50, 75, 100),
                     expand = c(0, 0),
                     labels = scales::label_percent(scale = 1)) +
  labs(y = "Percentage of studies (%)", x = "Trait category: Level 2",
       title = element_blank()) +
  coord_flip() +
  theme_minimal(base_size = 18) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "grey85", linewidth = 0.5),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        panel.background = element_rect(fill = "white", color = NA),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.position = "inside",
        legend.position.inside = c(0.2, 0.62),
        legend.background = element_rect(color = "black", fill = "white", linewidth = 0.5))
p

ggsave(plot=p, filename = "perc-studies-traits-level2.jpeg", path = "./output",
       width = 8, height = 6,  units = "in", dpi = 600)


### Understand associations between traits within each level (i.e., which traits are often studies together) (Fig. 4E,F)
library(purrr)

dat <- read.csv("./data/avian_traits_database.csv")

# basic chord diagram
# start with trait_level2
# group traits_level2 that constitute less than 5% of studies
dat1 <- dat %>%
  mutate(trait_level2 = if_else(trait_level2 == "", "unknown", trait_level2)) %>%
  mutate(trait_level2 = if_else(is.na(trait_level2), "unknown", trait_level2))

dat_t2 <- dat1 %>%
  select(studyID, trait_level2) %>%
  distinct()

unique(dat_t2$trait_level2)

trait2_summary <- dat1 %>%
  select(studyID, trait_level2) %>%
  distinct() %>%
  group_by(trait_level2) %>%
  summarise(n_studies = n()) %>%
  ungroup() %>%
  mutate(percentage = 100 * n_studies / studies_no)

trait2_summary <- trait2_summary %>%
  mutate(trait_level2 = if_else(percentage < 5, "Other", trait_level2)) |>
  summarise(n_studies = sum(n_studies),
            percentage = sum(percentage),
            .by = trait_level2)

traits <- unique(trait2_summary$trait_level2)

dat_t2 <- dat_t2 |>
  mutate(trait_level2 = if_else(trait_level2 %in% traits, trait_level2, "Other"))

dat_t2 <- dat_t2 %>%
  mutate(trait_level2 = recode(trait_level2,
                               "lifespan"       = "Lifespan",
                               "sociality"      = "Sociality",
                               "unknown" = "Unknown",
                               "habitat_breadth"       = "Habitat breadth",
                               "activity_time"    = "Activity time",
                               "Other" = "Other",
                               "tail" = "Tail",
                               "tarsus" = "Tarsus",
                               "habitat_type" = "Habitat type",
                               "migration" = "Migration",
                               "bill" = "Bill",
                               "wing" = "Wing",
                               "reproduction" = "Reproduction",
                               "foraging" = "Foraging",
                               "diet" = "Diet",
                               "body" = "Body mass"))


# remove studies with only 1 trait or missing values
dat_t2_multi <- dat_t2 %>%
  #filter(!is.na(trait_level2)) %>%   # remove rows with missing trait_level2
  group_by(studyID) %>%
  filter(n() > 1) %>%                # keep only studies with >1 trait_level2
  ungroup()

# Get all pairs of traits within the same study
co_occur <- dat_t2_multi |>
  inner_join(dat_t2_multi, by = "studyID") |>        # join table with itself on studyID
  filter(trait_level2.x < trait_level2.y) |>          # avoid duplicates and self-pairs
  group_by(trait_level2.x, trait_level2.y) |>
  summarise(n = n_distinct(studyID), .groups = "drop") # count how many studies share each pair

# convert to matrix
# Get all unique traits
all_traits <- unique(c(co_occur$trait_level2.x, co_occur$trait_level2.y))

# Build symmetric matrix
mat <- matrix(0, nrow = length(all_traits), ncol = length(all_traits),
              dimnames = list(all_traits, all_traits))

for (i in seq_len(nrow(co_occur))) {
  t1 <- co_occur$trait_level2.x[i]
  t2 <- co_occur$trait_level2.y[i]
  mat[t1, t2] <- co_occur$n[i]
  mat[t2, t1] <- co_occur$n[i]
}

#test for random associations
chisq_result <- chisq.test(mat)
print(chisq_result)
fisher_result <- fisher.test(mat, simulate.p.value = TRUE, B = 9999)
print(fisher_result)
cramers_v <- function(x) {
  ct <- chisq.test(x)
  sqrt(ct$statistic / (sum(x) * (min(dim(x)) - 1)))
}
obs_v <- cramers_v(mat) #this is the value we will report
cat("Observed Cramér's V:", round(obs_v, 4), "\n")

# plot
trait_colors <- setNames(viridis(nrow(mat)), rownames(mat))
trait_colors <- setNames(viridis(nrow(mat), option = "turbo"), rownames(mat))

png("./output/trait_level2_chord.png", width = 7000, height = 7000, res = 600)
chordDiagram(mat,
             grid.col = trait_colors,
             transparency = 0.4,
             annotationTrack = "grid",
             preAllocateTracks = list(track.height = 0.1))

# Add trait labels
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  sector.name = get.cell.meta.data("sector.index")
  circos.text(CELL_META$xcenter, CELL_META$ylim[1] + 0.5,
              sector.name, facing = "clockwise",
              niceFacing = TRUE, adj = c(0, 0.5), cex = 1.4)
}, bg.border = NA)

dev.off()

png("./output/trait_level2_chord.png", width = 7000, height = 7000, res = 600)
chordDiagram(mat,
             grid.col = trait_colors,
             transparency = 0.4,
             annotationTrack = "grid",
             annotationTrackHeight = 0.05)
dev.off()


# do the same for trait_level1
dat <- read.csv("./data/avian_traits_database.csv")

dat1 <- dat %>%
  mutate(trait_level1 = if_else(trait_level1 == "", "unknown", trait_level1)) %>%
  mutate(trait_level1 = if_else(is.na(trait_level1), "unknown", trait_level1))

dat1 <- dat1 %>%
  mutate(trait_level1 = recode(trait_level1,
                               "anatomy"       = "Anatomy",
                               "physiology"      = "Physiology",
                               "unknown" = "Unknown",
                               "behavior_general"       = "Behavior",
                               "activity"    = "Activity",
                               "dispersal" = "Dispersal",
                               "habitat" = "Habitat",
                               "life_history" = "Life history",
                               "trophic" = "Trophic",
                               "morphology" = "Morphology"))

dat_t1 <- dat1 %>%
  select(studyID, trait_level1) %>%
  distinct()

# remove studies with only 1 trait or missing values
dat_t1_multi <- dat_t1 %>%
  #filter(!is.na(trait_level2)) %>%   # remove rows with missing trait_level2
  group_by(studyID) %>%
  filter(n() > 1) %>%                # keep only studies with >1 trait_level2
  ungroup()


# Get all pairs of traits within the same study
co_occur <- dat_t1_multi |>
  inner_join(dat_t1_multi, by = "studyID") |>        # join table with itself on studyID
  filter(trait_level1.x < trait_level1.y) |>          # avoid duplicates and self-pairs
  group_by(trait_level1.x, trait_level1.y) |>
  summarise(n = n_distinct(studyID), .groups = "drop") # count how many studies share each pair

# convert to matrix
# Get all unique traits
all_traits <- unique(c(co_occur$trait_level1.x, co_occur$trait_level1.y))

# Build symmetric matrix
mat <- matrix(0, nrow = length(all_traits), ncol = length(all_traits),
              dimnames = list(all_traits, all_traits))

for (i in seq_len(nrow(co_occur))) {
  t1 <- co_occur$trait_level1.x[i]
  t2 <- co_occur$trait_level1.y[i]
  mat[t1, t2] <- co_occur$n[i]
  mat[t2, t1] <- co_occur$n[i]
}

# test for random associations
chisq_result <- chisq.test(mat)
print(chisq_result)
fisher_result <- fisher.test(mat, simulate.p.value = TRUE, B = 9999)
print(fisher_result)
cramers_v <- function(x) {
  ct <- chisq.test(x)
  sqrt(ct$statistic / (sum(x) * (min(dim(x)) - 1)))
}
obs_v <- cramers_v(mat) #this is the value we will report
cat("Observed Cramér's V:", round(obs_v, 4), "\n")

# plot
trait_colors <- setNames(viridis(nrow(mat), option = "turbo"), rownames(mat))

png("./output/trait_level1_chord.png", width = 7000, height = 7000, res = 600)
chordDiagram(mat,
             grid.col = trait_colors,
             transparency = 0.4,
             annotationTrack = "grid",
             preAllocateTracks = list(track.height = 0.1))

# Add trait labels
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  sector.name = get.cell.meta.data("sector.index")
  circos.text(CELL_META$xcenter, CELL_META$ylim[1] + 0.5,
              sector.name, facing = "clockwise",
              niceFacing = TRUE, adj = c(0, 0.5), cex = 1.4)
}, bg.border = NA)

dev.off()

# plot with no labels
png("./output/trait_level1_chord.png", width = 7000, height = 7000, res = 600)
chordDiagram(mat,
             grid.col = trait_colors,
             transparency = 0.4,
             annotationTrack = "grid",
             annotationTrackHeight = 0.05)
dev.off()




###### Methodlogical transparency (Fig. 5A-D))
dat <- read.csv("./data/avian_traits_database.csv")

# intra-specific trait variation (Fig. 5A)
dat_itv <- dat %>%
  distinct(studyID, ITV)

dat_itv_count <- dat_itv %>%
  count(ITV) %>%
  mutate(percentage = n / sum(n) * 100,
         ITV = case_when(
           ITV == "No"  ~ "ITV not included",
           ITV == "Yes" ~ "ITV included",
           is.na(ITV)   ~ "No information",
           TRUE ~ ITV
         ))

# Bar plot
p <- ggplot(dat_itv_count, aes(x = factor(ITV, levels = c("ITV included", "ITV not included", "No information")),
                               y = percentage, fill = ITV)) +
  geom_col() +
  geom_text(aes(label = paste0(n, " (", round(percentage, 1), "%)")),
            vjust = -0.5, size = 5) +
  scale_y_continuous(limits = c(0, max(dat_itv_count$percentage) * 1.15),
                     expand = c(0, 0)) +
  scale_fill_manual(values = c("ITV included"     = "darkgoldenrod1",
                               "ITV not included" = "dodgerblue1",
                               "No information"   = "grey"),
                    breaks = c("ITV included", "ITV not included", "No information")) +
  labs(
    x = "ITV status",
    y = "Percentage of studies (%)",
    title = element_blank()) +
  theme_minimal(base_size = 18) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "grey85", linewidth = 0.5),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        panel.background = element_rect(fill = "white", color = NA),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.position = "inside",
        legend.position.inside = c(0.2, 0.8),
        legend.title = element_blank(),
        legend.background = element_rect(color = "black", fill = "white", linewidth = 0.5))
p

ggsave(plot=p, filename = "number-studies-itv.jpeg", path = "./output",
       width = 8, height = 6,  units = "in", dpi = 600)


# accounting for probability of species detection (Fig. 5B)
dat_det <- dat %>%
  distinct(studyID, detection)

dat_det_count <- dat_det %>%
  count(detection) %>%
  mutate(percentage = n / sum(n) * 100,
         detection = case_when(
           detection == "No"  ~ "Detection not included",
           detection == "Yes" ~ "Detection included",
           is.na(detection)   ~ "No information",
           TRUE ~ detection
         ))

# Bar plot
p <- ggplot(dat_det_count, aes(x = factor(detection, levels = c("Detection included", "Detection not included", "No information")),
                               y = percentage, fill = detection)) +
  geom_col() +
  geom_text(aes(label = paste0(n, " (", round(percentage, 1), "%)")),
            vjust = -0.5, size = 5) +
  scale_y_continuous(limits = c(0, max(dat_det_count$percentage) * 1.15),
                     expand = c(0, 0)) +
  scale_fill_manual(values = c("Detection included"     = "darkgoldenrod1",
                               "Detection not included" = "dodgerblue1",
                               "No information"         = "grey"),
                    breaks = c("Detection included", "Detection not included", "No information")) +
  labs(
    x = "Imperfect detection status",
    y = "Percentage of studies (%)",
    title = element_blank()) +
  theme_minimal(base_size = 18) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "grey85", linewidth = 0.5),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        panel.background = element_rect(fill = "white", color = NA),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.position = "inside",
        legend.position.inside = c(0.25, 0.8),
        legend.title = element_blank(),
        legend.background = element_rect(color = "black", fill = "white", linewidth = 0.5))
p

ggsave(plot=p, filename = "prop-studies-detection.jpeg", path = "./output",
       width = 8, height = 6,  units = "in", dpi = 600)


# trait completeness (Fig. 5C)
dat_tcompl <- dat %>%
  distinct(studyID, trait_completeness)

dat_tcompl <- dat_tcompl |>
  mutate(trait_completeness = case_when(
    trait_completeness %in% c("Yes", "Yes ")  ~ "Trait completeness",
    trait_completeness %in% c("No", "No ")    ~ "No trait completeness",
    is.na(trait_completeness) | trait_completeness == "" ~ "No information"
  ))

dat_tcompl_count <- dat_tcompl %>%
  count(trait_completeness) %>%
  mutate(percentage = n / sum(n) * 100)

# Bar plot
p <- ggplot(dat_tcompl_count, aes(x = factor(trait_completeness, levels = c("Trait completeness", "No trait completeness", "No information")),
                                  y = percentage, fill = trait_completeness)) +
  geom_col() +
  geom_text(aes(label = paste0(n, " (", round(percentage, 1), "%)")),
            vjust = -0.5, size = 5) +
  scale_y_continuous(limits = c(0, max(dat_tcompl_count$percentage) * 1.15),
                     expand = c(0, 0)) +
  scale_fill_manual(values = c("Trait completeness"     = "darkgoldenrod1",
                               "No trait completeness" = "dodgerblue1",
                               "No information"         = "grey"),
                    breaks = c("Trait completeness", "No trait completeness", "No information")) +
  labs(
    x = "Trait completeness status",
    y = "Percentage of studies (%)",
    title = element_blank()) +
  theme_minimal(base_size = 18) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "grey85", linewidth = 0.5),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        panel.background = element_rect(fill = "white", color = NA),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.position = "inside",
        legend.position.inside = c(0.25, 0.8),
        legend.title = element_blank(),
        legend.background = element_rect(color = "black", fill = "white", linewidth = 0.5))
p

ggsave(plot=p, filename = "number-studies-tcomplt.jpeg", path = "./output",
       width = 8, height = 6,  units = "in", dpi = 600)


# trait imputation (Fig. 5D)
dat_timp <- dat %>%
  distinct(studyID, trait_imputation)

dat_timp <- dat_timp |>
  mutate(trait_imputation = case_when(
    trait_imputation == "Yes"  ~ "Trait imputation",
    trait_imputation == "No"  ~ "No trait imputation",
    is.na(trait_imputation) | trait_imputation == "" ~ "No information"
  ))

dat_timp_count <- dat_timp %>%
  count(trait_imputation) %>%
  mutate(percentage = n / sum(n) * 100)


# Bar plot
p <- ggplot(dat_timp_count, aes(x = factor(trait_imputation, levels = c("Trait imputation", "No trait imputation", "No information")),
                                y = percentage, fill = trait_imputation)) +
  geom_col() +
  geom_text(aes(label = paste0(n, " (", round(percentage, 1), "%)")),
            vjust = -0.5, size = 5) +
  scale_y_continuous(limits = c(0, max(dat_timp_count$percentage) * 1.15),
                     expand = c(0, 0)) +
  scale_fill_manual(values = c("Trait imputation"     = "darkgoldenrod1",
                               "No trait imputation" = "dodgerblue1",
                               "No information"         = "grey"),
                    breaks = c("Trait imputation", "No trait imputation", "No information")) +
  labs(
    x = "Trait imputation status",
    y = "Percentage of studies (%)",
    title = element_blank()) +
  theme_minimal(base_size = 18) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "grey85", linewidth = 0.5),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        panel.background = element_rect(fill = "white", color = NA),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.position = "inside",
        legend.position.inside = c(0.25, 0.8),
        legend.title = element_blank(),
        legend.background = element_rect(color = "black", fill = "white", linewidth = 0.5))
p

ggsave(plot=p, filename = "number-studies-timp.jpeg", path = "./output",
       width = 8, height = 6,  units = "in", dpi = 600)



###### Methodological transparency  - continued
dat <- read.csv("./data/avian_traits_database.csv")

### the proportion of studies that provide any functional rationale for trait selection
# (i.e., the average proportion of traits per study for which the ecological significance is explicitly described; functional_significance field)
dat_sub <- dat %>%
  group_by(studyID) %>%
  summarise(
    has_justification = any(!is.na(functional_significance) &
                              functional_significance != ""),
    .groups = "drop"
  ) %>%
  summarise(
    n_studies = n(),
    n_with = sum(has_justification),
    n_without = sum(!has_justification),
    prop_with = mean(has_justification)
  )

dat_function <- dat %>%
  distinct(studyID, trait_synonym, functional_significance)
study_proportions <- dat_function %>%
  group_by(studyID) %>%
  summarise(
    n_traits = n(),
    n_described = sum(!is.na(functional_significance)),
    prop_described = n_described / n_traits
  )
mean(study_proportions$prop_described)

# of studies that provide at ecol significance for at least 1 trait, how many traits are assigned signifcance
study_proportions_filtered <- study_proportions %>%
  filter(prop_described > 0)
mean(study_proportions_filtered$prop_described)


### proportion of studies that provide info on effect versus response trait
dat_sub <- dat %>%
  group_by(studyID) %>%
  summarise(
    has_justification = any(!is.na(trait_type) &
                              trait_type != ""),
    .groups = "drop"
  ) %>%
  summarise(
    n_studies = n(),
    n_with = sum(has_justification),
    n_without = sum(!has_justification),
    prop_with = mean(has_justification)
  )

dat_sub



###### Methodological decisions: fucntional diversity metrics (Fig. 6)
dat <- read.csv("./data/avian_traits_database.csv")
dat <- dat %>%
  mutate(main_metric_group = if_else(main_metric_group == "", "Unknown", main_metric_group)) %>%
  mutate(main_metric_group = if_else(is.na(main_metric_group), "Unknown", main_metric_group))
dat <- dat %>%
  mutate(main_metric_group = recode(main_metric_group,
                                    "functional diversity"       = "Functional\ndiversity",
                                    "functional composition"      = "Functional\ncomposition",
                                    "unknown" = "Unknown",
                                    "functional rarity"       = "Functional\nrarity",
                                    "functional turnover"    = "Functional\nturnover",
                                    "functional nestedness" = "Functional\nnestedness",
                                    "functional group diversity" = "Functional\ngroup diversity",
                                    "functional structure" = "Functional\nstructure"))

# main metric (Fig. 6A)
unique(dat$main_metric_group)

studies_no <- length(unique(dat$studyID))

dat_summary <- dat %>%
  group_by(main_metric_group) %>%
  summarise(n_studies = n_distinct(studyID), .groups = "drop") %>%
  mutate(percentage = 100 * n_studies / studies_no) %>%
  arrange(desc(n_studies))  # optional: sort descending

# bar plot
p <- ggplot(dat_summary, aes(x = reorder(main_metric_group, -percentage), y = percentage,
                             fill = reorder(main_metric_group, -percentage))) +
  geom_col(color = "black") +
  geom_text(aes(label = paste0(n_studies, " (", round(percentage, 1), "%)")),
            hjust = -0.1, size = 4.5) +
  scale_y_continuous(limits = c(0, 130),
                     breaks = c(0, 25, 50, 75, 100),
                     expand = c(0, 0),
                     labels = scales::label_percent(scale = 1)) +
  scale_fill_viridis_d(option="turbo", direction = -1) +  # change option to "viridis", "magma", "plasma" etc.
  labs(y = "Percentage of studies (%)", x = "Main functional metric",
       title = element_blank()) +
  coord_flip() +
  theme_minimal(base_size = 18) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "grey85", linewidth = 0.5),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        panel.background = element_rect(fill = "white", color = NA),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.position = "none")  # legend not needed since x axis labels the bars
p

ggsave(plot=p, filename = "main_metric.jpeg", path = "./output",
       width = 8, height = 6,  units = "in", dpi = 600)


### functional metrics within the main metric category (Fig. 6B)
# subset to one main metric group, e.g., "Functional Diversity"
# we will only look at f diversity ,rarity, composition, turnover, and nestedness
cat <- unique(dat$main_metric_group)

### go with i=1, functional diversity
i=1
dat_fd <- dat %>%
  filter(main_metric_group == cat[[i]]) %>%
  filter(!is.na(functional_metric))  # remove NAs

unique(dat_fd$functional_metric)

# summarize counts per functional_metric
studies_no <- length(unique(dat_fd$studyID))

dat_fd_summary <- dat_fd %>%
  group_by(functional_metric) %>%
  summarise(n_studies = n_distinct(studyID), .groups = "drop") %>%
  mutate(percentage = 100 * n_studies / studies_no)

dat_fd_summary2 <- dat_fd_summary %>%
  mutate(functional_metric = if_else(percentage < 5, "Other", functional_metric)) |>
  summarise(n_studies = sum(n_studies),
            percentage = sum(percentage),
            .by = functional_metric)

dat_fd_summary2 <- dat_fd_summary2 %>%
  mutate(functional_metric = recode(functional_metric,
                                    "Other"       = "Other",
                                    "functional dispersion"       = "Functional\ndispersion",
                                    "functional beta diversity"       = "Functional\nbeta diversity",
                                    "functional divergence"    = "Functional\ndivergence",
                                    "functional diversity" = "Functional\ndiversity",
                                    "functional richness" = "Functional\nrichness",
                                    "functional evenness" = "Functional\nevenness",
                                    "mean pairwise trait dissimilarity" = "Trait\ndissimilarity",
                                    "Rao's quadratic entropy" = "Rao's quadratic\nentropy"))

# plot proportions as a bar chart
p <- ggplot(dat_fd_summary2, aes(x = reorder(functional_metric, -percentage), y = percentage,
                                 fill = reorder(functional_metric, -percentage))) +
  geom_col(fill="brown4", color = "black") +
  geom_text(aes(label = paste0(n_studies, " (", round(percentage, 1), "%)")),
            hjust = -0.1, size = 4.5) +
  scale_y_continuous(limits = c(0, 80),
                     breaks = c(0, 25, 50, 75, 100),
                     expand = c(0, 0),
                     labels = scales::label_percent(scale = 1)) +
  #scale_fill_viridis_d(option="turbo", direction = -1) +  # change option to "viridis", "magma", "plasma" etc.
  labs(y = "Percentage of studies (%)", x = "Functional diversity metric",
       title = element_blank()) +
  coord_flip() +
  theme_minimal(base_size = 18) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "grey85", linewidth = 0.5),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        panel.background = element_rect(fill = "white", color = NA),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.position = "none")  # legend not needed since x axis labels the bars
p

ggsave(plot=p, filename = paste0("div_metrics_fd.jpeg"), path = "./output",
       width = 8, height = 6,  units = "in", dpi = 600)

### go with i=2, functional rarity
i=2
dat_fd <- dat %>%
  filter(main_metric_group == cat[[i]]) %>%
  filter(!is.na(functional_metric))  # remove NAs

unique(dat_fd$functional_metric)

# summarize counts per functional_metric
studies_no <- length(unique(dat_fd$studyID))

dat_fd_summary <- dat_fd %>%
  group_by(functional_metric) %>%
  summarise(n_studies = n_distinct(studyID), .groups = "drop") %>%
  mutate(percentage = 100 * n_studies / studies_no)

dat_fd_summary2 <- dat_fd_summary %>%
  mutate(functional_metric = if_else(percentage < 5, "Other", functional_metric)) |>
  summarise(n_studies = sum(n_studies),
            percentage = sum(percentage),
            .by = functional_metric)

dat_fd_summary2 <- dat_fd_summary2 %>%
  mutate(functional_metric = recode(functional_metric,
                                    "Other"       = "Other",
                                    "functional distinctiveness"       = "Functional\ndistinctinctiveness",
                                    "functional uniqueness"    = "Functional\nuniqueness",
                                    "mean functional distinctiveness" = "Functional\ndistinctinctiveness\n(mean)",
                                    "mean functional uniqueness" = "Functional\nuniqueness\n(mean)"))

# plot proportions as a bar chart
p <- ggplot(dat_fd_summary2, aes(x = reorder(functional_metric, -percentage), y = percentage,
                                 fill = reorder(functional_metric, -percentage))) +
  geom_col(fill="orangered3", color = "black") +
  geom_text(aes(label = paste0(n_studies, " (", round(percentage, 1), "%)")),
            hjust = -0.1, size = 4.5) +
  scale_y_continuous(limits = c(0, 80),
                     breaks = c(0, 25, 50, 75, 100),
                     expand = c(0, 0),
                     labels = scales::label_percent(scale = 1)) +
  #scale_fill_viridis_d(option="turbo", direction = -1) +  # change option to "viridis", "magma", "plasma" etc.
  labs(y = "Percentage of studies (%)", x = "Functional rarity metric",
       title = element_blank()) +
  coord_flip() +
  theme_minimal(base_size = 18) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "grey85", linewidth = 0.5),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        panel.background = element_rect(fill = "white", color = NA),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.position = "none")  # legend not needed since x axis labels the bars
p

ggsave(plot=p, filename = paste0("div_metrics_frar.jpeg"), path = "./output",
       width = 8, height = 6,  units = "in", dpi = 600)


### go with i=3, functional composition
i=3
dat_fd <- dat %>%
  filter(main_metric_group == cat[[i]]) %>%
  filter(!is.na(functional_metric))  # remove NAs

unique(dat_fd$functional_metric)
dat_fd <- dat_fd %>%
  mutate(functional_metric = recode(functional_metric,
                                    "Community-weighted mean" = "community-weighted mean"))

# summarize counts per functional_metric
studies_no <- length(unique(dat_fd$studyID))

dat_fd_summary <- dat_fd %>%
  group_by(functional_metric) %>%
  summarise(n_studies = n_distinct(studyID), .groups = "drop") %>%
  mutate(percentage = 100 * n_studies / studies_no)

dat_fd_summary2 <- dat_fd_summary %>%
  mutate(functional_metric = if_else(percentage < 5, "Other", functional_metric)) |>
  summarise(n_studies = sum(n_studies),
            percentage = sum(percentage),
            .by = functional_metric)

dat_fd_summary2 <- dat_fd_summary2 %>%
  mutate(functional_metric = recode(functional_metric,
                                    "Other"       = "Other",
                                    "community-weighted mean" = "Community-\nweighted\nmean"))

# plot proportions as a bar chart

p <- ggplot(dat_fd_summary2, aes(x = reorder(functional_metric, -percentage), y = percentage,
                                 fill = reorder(functional_metric, -percentage))) +
  geom_col(fill="darkorange1", color = "black") +
  geom_text(aes(label = paste0(n_studies, " (", round(percentage, 1), "%)")),
            hjust = -0.1, size = 4.5) +
  scale_y_continuous(limits = c(0, 125),
                     breaks = c(0, 25, 50, 75, 100),
                     expand = c(0, 0),
                     labels = scales::label_percent(scale = 1)) +
  #scale_fill_viridis_d(option="turbo", direction = -1) +  # change option to "viridis", "magma", "plasma" etc.
  labs(y = "Percentage of studies (%)", x = "Functional composition metric",
       title = element_blank()) +
  coord_flip() +
  theme_minimal(base_size = 18) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "grey85", linewidth = 0.5),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        panel.background = element_rect(fill = "white", color = NA),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.position = "none")  # legend not needed since x axis labels the bars
p

ggsave(plot=p, filename = paste0("div_metrics_fcomp.jpeg"), path = "./output",
       width = 8, height = 6,  units = "in", dpi = 600)


### go with i=5, functional turnover
i=5
dat_fd <- dat %>%
  filter(main_metric_group == cat[[i]]) %>%
  filter(!is.na(functional_metric))  # remove NAs

unique(dat_fd$functional_metric)

# summarize counts per functional_metric
studies_no <- length(unique(dat_fd$studyID))

dat_fd_summary <- dat_fd %>%
  group_by(functional_metric) %>%
  summarise(n_studies = n_distinct(studyID), .groups = "drop") %>%
  mutate(percentage = 100 * n_studies / studies_no)

dat_fd_summary2 <- dat_fd_summary %>%
  mutate(functional_metric = recode(functional_metric,
                                    "functional turnover"       = "Functional\nturnover",
                                    "overlap-based dissimilarity" = "Overlap-based\ndissimilarity",
                                    "unique fraction of functional space occupied" = "Unique fraction\nof functional\nspace"))

# plot proportions as a bar chart

p <- ggplot(dat_fd_summary2, aes(x = reorder(functional_metric, -percentage), y = percentage,
                                 fill = reorder(functional_metric, -percentage))) +
  geom_col(fill="green1", color = "black") +
  geom_text(aes(label = paste0(n_studies, " (", round(percentage, 1), "%)")),
            hjust = -0.1, size = 4.5) +
  scale_y_continuous(limits = c(0, 115),
                     breaks = c(0, 25, 50, 75, 100),
                     expand = c(0, 0),
                     labels = scales::label_percent(scale = 1)) +
  #scale_fill_viridis_d(option="turbo", direction = -1) +  # change option to "viridis", "magma", "plasma" etc.
  labs(y = "Percentage of studies (%)", x = "Functional turnover metric",
       title = element_blank()) +
  coord_flip() +
  theme_minimal(base_size = 18) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "grey85", linewidth = 0.5),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        panel.background = element_rect(fill = "white", color = NA),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.position = "none")  # legend not needed since x axis labels the bars
p

ggsave(plot=p, filename = paste0("div_metrics_fturn.jpeg"), path = "./output",
       width = 8, height = 6,  units = "in", dpi = 600)


### go with i=6, functional nestedness
i=6
dat_fd <- dat %>%
  filter(main_metric_group == cat[[i]]) %>%
  filter(!is.na(functional_metric))  # remove NAs

unique(dat_fd$functional_metric)

# summarize counts per functional_metric
studies_no <- length(unique(dat_fd$studyID))

dat_fd_summary <- dat_fd %>%
  group_by(functional_metric) %>%
  summarise(n_studies = n_distinct(studyID), .groups = "drop") %>%
  mutate(percentage = 100 * n_studies / studies_no)

dat_fd_summary2 <- dat_fd_summary %>%
  mutate(functional_metric = if_else(percentage < 5, "Other", functional_metric)) |>
  summarise(n_studies = sum(n_studies),
            percentage = sum(percentage),
            .by = functional_metric)

dat_fd_summary2 <- dat_fd_summary2 %>%
  mutate(functional_metric = recode(functional_metric,
                                    "functional nestedness" = "Functional\nnestedness"))

# 3. Plot proportions as a bar chart

p <- ggplot(dat_fd_summary2, aes(x = reorder(functional_metric, -percentage), y = percentage,
                                 fill = reorder(functional_metric, -percentage))) +
  geom_col(fill="turquoise3", color = "black") +
  geom_text(aes(label = paste0(n_studies, " (", round(percentage, 1), "%)")),
            hjust = -0.1, size = 4.5) +
  scale_y_continuous(limits = c(0, 110),
                     breaks = c(0, 25, 50, 75, 100),
                     expand = c(0, 0),
                     labels = scales::label_percent(scale = 1)) +
  #scale_fill_viridis_d(option="turbo", direction = -1) +  # change option to "viridis", "magma", "plasma" etc.
  labs(y = "Percentage of studies (%)", x = "Functional nestedness metric",
       title = element_blank()) +
  coord_flip() +
  theme_minimal(base_size = 18) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "grey85", linewidth = 0.5),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        panel.background = element_rect(fill = "white", color = NA),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.position = "none")  # legend not needed since x axis labels the bars
p

ggsave(plot=p, filename = paste0("div_metrics_fnest.jpeg"), path = "./output",
       width = 8, height = 6,  units = "in", dpi = 600)


### trait space estimation (Fig. 6C)
dat <- read.csv("./data/avian_traits_database.csv")

unique(dat$trait_space_estimation)
dat <- dat %>%
  mutate(trait_space_estimation = if_else(trait_space_estimation == "", "Unknown", trait_space_estimation)) %>%
  mutate(trait_space_estimation = if_else(is.na(trait_space_estimation), "Unknown", trait_space_estimation))

unique(dat$trait_space_estimation)
studies_no <- length(unique(dat$studyID))

dat_summary <- dat %>%
  group_by(trait_space_estimation) %>%
  summarise(n_studies = n_distinct(studyID), .groups = "drop") %>%
  mutate(percentage = 100 * n_studies / studies_no) %>%
  arrange(desc(n_studies))  # optional: sort descending

dat_summary2 <- dat_summary %>%
  mutate(trait_space_estimation = if_else(percentage < 5, "Other", trait_space_estimation)) |>
  summarise(n_studies = sum(n_studies),
            percentage = sum(percentage),
            .by = trait_space_estimation)

dat_summary2 <- dat_summary2 %>%
  mutate(trait_space_estimation = recode(trait_space_estimation,
                                         "none"       = "None",
                                         "convex hull"      = "Convex hull",
                                         "minimum spanning tree" = "Minimum\nspanning tree",
                                         "functional dendrogram"       = "Functional\ndendrogram",
                                         "Other"    = "Other"))
# bar plot
p <- ggplot(dat_summary2, aes(x = reorder(trait_space_estimation, -percentage), y = percentage,
                              fill = reorder(trait_space_estimation, -percentage))) +
  geom_col(color = "black") +
  geom_text(aes(label = paste0(n_studies, " (", round(percentage, 1), "%)")),
            hjust = -0.1, size = 4.5) +
  scale_y_continuous(limits = c(0, 90),
                     breaks = c(0, 25, 50, 75),
                     expand = c(0, 0),
                     labels = scales::label_percent(scale = 1)) +
  scale_fill_viridis_d(option="turbo", direction = -1) +  # change option to "viridis", "magma", "plasma" etc.
  labs(y = "Percentage of studies (%)", x = "Main trait space metric",
       title = element_blank()) +
  coord_flip() +
  theme_minimal(base_size = 18) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "grey85", linewidth = 0.5),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        panel.background = element_rect(fill = "white", color = NA),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.position = "none")  # legend not needed since x axis labels the bars
p

ggsave(plot=p, filename = "trait_space_metric.jpeg", path = "./output",
       width = 8, height = 6,  units = "in", dpi = 600)



