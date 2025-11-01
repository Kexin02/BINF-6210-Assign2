# ====================================================
# Comparing the BIN Composition of Subfamily Sciurinae 
# between North America and Eurasia
# Author: Kexin Gong
# ====================================================

# 0. Load required packages ----
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(vegan)
  library(readr)
  library(stringr)

# 1. Read data file ----
dat <- read.delim("../data/result.tsv", sep = "\t", header = TRUE, quote = "", check.names = FALSE)
# Replace the inconvenient column name used in the following code
dat <- dat %>% rename(country = `country/ocean`)

# 2. Map countries to continents ----
# Define lists of countries belonging to North America and Eurasia
# The following used country names based on the data file
north_america <- c(
  "Canada", "United States", "Mexico", "Panama", "Costa Rica"
)
eurasia <- c(
  "Austria", "China", "Denmark", "Georgia", "Germany", "India", "Italy",
  "Laos", "Lebanon", "Mongolia", "Norway", "Russia", "Slovakia", "South Korea",
  "Switzerland", "Taiwan", "Turkiye", "Vietnam"
)

# Assign each record to a continent
dat <- dat %>%
  mutate(continent = case_when(
    country %in% north_america ~ "North America",
    country %in% eurasia       ~ "Eurasia",
    TRUE ~ NA_character_
  ))

# Keep only the two continents of interest
dat2 <- dat %>% filter(continent %in% c("North America", "Eurasia"))

# 3. Build BIN-by-continent matrix ----
bin_counts <- dat2 %>%
  filter(!is.na(bin_uri)) %>%
  group_by(continent, bin_uri) %>%
  summarise(n = n(), .groups = "drop")

# Convert to wide format: one row per BIN, one column per continent
bin_wide <- bin_counts %>%
  tidyr::pivot_wider(names_from = continent, values_from = n, values_fill = 0)

# Presence/absence matrix (1 = present, 0 = absent)
pa <- bin_wide %>%
  transmute(bin_uri,
            `North America` = as.integer(`North America` > 0),
            `Eurasia`       = as.integer(`Eurasia` > 0))

# 4. Compute similarity ----
# Convert to 2×N matrix (rows = continents, columns = BINs)
mat <- as.matrix(t(pa[, -1]))
rownames(mat) <- colnames(pa[, -1])

# Compute binary Jaccard distance
dist_jac <- vegdist(mat, method = "jaccard", binary = TRUE)
print(dist_jac)

# Compute Bray–Curtis distance
mat_bray <- as.matrix(t(bin_wide %>% select(`North America`,`Eurasia`)))
dist_bray <- vegdist(mat_bray, method = "bray")

# Count shared BINs
shared_bins <- pa %>% filter(`North America` == 1 & `Eurasia` == 1) %>% pull(bin_uri)

message(sprintf("Unique BINs NA=%d, EUAS=%d, Shared=%d",
                sum(pa$`North America`), sum(pa$Eurasia), length(shared_bins)))

# 5. Create figures ----
# Figure 1: BIN richness per continent (bar plot)
bin_rich <- pa |>
  summarise(`North America` = sum(`North America`),
            `Eurasia` = sum(Eurasia)) |>
  pivot_longer(everything(), names_to = "continent", values_to = "unique_bins")

p1 <- ggplot(bin_rich, aes(continent, unique_bins, fill = continent)) +
  geom_col(width = 0.7) +
  labs(title = "Sciurinae BIN richness by continent",
       x = NULL, y = "Number of unique BINs") +
  scale_fill_manual(values = c("Eurasia" = "#FFC3C3", "North America" = "#FFE3A9"))+
  theme_minimal(base_size = 12) +
  theme(legend.position = "none",
        panel.grid = element_blank()) +
  ylim(0,15)

ggsave("../figs/Fig1_BIN_richness_by_continent.png", p1, width = 6, height = 4, dpi = 300)

# Figure 2: Presence/absence heatmap (BIN × continent)
pa_long <- pa %>%
  pivot_longer(cols = c(`North America`, `Eurasia`),
               names_to = "continent", values_to = "present")

# Sort BINs by total occurrences to make the heatmap cleaner
keep_ids <- pa %>% mutate(tot = `North America` + Eurasia) %>%
  arrange(desc(tot)) %>% pull(bin_uri)
pa_long$bin_uri <- factor(pa_long$bin_uri, levels = keep_ids)

p2 <- ggplot(pa_long, aes(continent, bin_uri, fill = factor(present))) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c("0" = "#eeeeee", "1" = "#B0CE88"), name = "Present") +
  labs(title = "Presence/absence of Sciurinae BINs across continents",
       x = NULL, y = "BIN") +
  theme_minimal(base_size = 10) +
  theme(legend.position = "right",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank())

ggsave("../figs/Fig2_BIN_presence_heatmap.png", p2, width = 6, height = 6, dpi = 300)

# Figure 3: Number of public records by continent
rec_counts <- dat2 %>% count(continent)
p3 <- ggplot(rec_counts, aes(continent, n, fill = continent)) +
  geom_col(width = 0.7) +
  labs(title = "Number of public records by continent",
       x = NULL, y = "Records") +
  scale_fill_manual(values = c("Eurasia" = "#FFC3C3", "North America" = "#FFE3A9"))+
  theme_minimal(base_size = 12) +
  theme(legend.position = "none",
        panel.grid = element_blank())

ggsave("../figs/Fig3_public_records_by_continent.png", p3, width = 6, height = 4, dpi = 300)


# 6. Summary ----
cat("\n=== SUMMARY ===\n")
cat("File:", "../data/result.tsv", "\n")
cat("Records used:", nrow(dat2), "\n")
cat("Continents:\n"); print(table(dat2$continent))
cat("Unique BINs per continent:\n"); print(bin_rich)
cat("Jaccard distance (0 = identical, 1 = completely different):\n"); print(dist_jac)
cat("Bray–Curtis distance:\n"); print(as.matrix(dist_bray))
cat("Shared BIN IDs (first 10):\n"); print(head(shared_bins, 10))