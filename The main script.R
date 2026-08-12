###Linear regression###

#Upload the OTU table and creating versions we will need downstream
library("readxl")
OTU.tab<-read_excel("OTU_tab.xlsx")
library("purrr")
library("tidyr") ######- PRIDANO
library("dplyr") ######- PRIDANO 
OTU.tab1<-OTU.tab[,-2]
OTU.tab1<-map_df(OTU.tab1, as.numeric)
OTU.tab1.1<-OTU.tab1[,-1]
OTU.tab2<-as.data.frame(t(OTU.tab1.1))
otu_long<-pivot_longer(OTU.tab1.1, cols = starts_with("0"), names_to = "sample",values_to = "count")
otu_data<-OTU.tab2

OTU.tab<-read.table(OTUtab)

#_____________________________________________________________________________________________________
###Richness###

richness<-read.table("richness_hotovo.txt", header = TRUE, stringsAsFactors = TRUE)

#separating samples based on the altitudinal zonation
richness<- richness %>%
  mutate(zone = case_when (
    Altitude <= 1700 ~ "Forest",
    Altitude > 1700 & Altitude <= 2200 ~ "Transient",
    Altitude > 2200 ~ "Alpine"
  ))

library(ggplot2)

#Comparison of species richness of specific samples across altitudinal zones
ggplot(richness, aes(x = sample, y = Richness, fill = zone)) + 
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = c(
    "Forest" = "forestgreen",
    "Transient" = "deeppink",
    "Alpine" = "dodgerblue"
  )) + 
  labs(title = "Species richness across altitudinal zones",
       x = "Sample",
       y = "Species richness",
       fill = "Zone") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

#ordering of the zones
richness<- richness %>%
  mutate(zone = factor(zone, levels=c("Forest", "Transient", "Alpine")))

#visualization of richness in boxplots
ggplot(richness, aes(x = zone, y = Richness, fill = zone)) + 
  geom_boxplot() + 
  scale_fill_manual(values = c(
    "Forest" = "forestgreen",
    "Transient" = "deeppink",
    "Alpine" = "dodgerblue"
  )) + 
  labs(title = "Species richness across altitudinal zones",
       x = "Zone",
       y = "Species richness",
       fill = "Zone") +
  theme_minimal()

#analysis to find out if there is a significant difference in richness between the zones
anova_richness<-aov(Richness ~ zone, data = richness)
summary(anova_richness)

#comparison of the zones among each other
TukeyHSD(anova_richness)

#__________________________________________________________________________________________

###Shannon-Wiener Index###

library(vegan)

#uploading metadata including altitude, elevation, and Shannon-Wiener index
metadata<-read_excel("metadata_complete.xlsx", sheet = 2)

#transfering to numeric formate if needed
metadata$Shannon<-as.numeric(as.character(metadata$Shannon))
metadata$elevation<-as.numeric(as.character(metadata$elevation))
summary(metadata$Shannon)
summary(metadata$elevation)

#corellate index with the altitude    ####### - VYMAZAT???
cor.test(metadata$elevation, metadata$Shannon, method= "spearman")

#linear model creation
lm_model<-lm(Shannon~elevation, data = metadata)
summary(lm_model)


#visualization
ggplot(metadata, aes(x=elevation, y=Shannon)) + 
  geom_point(color = "turquoise4") + 
  geom_smooth(method = "lm", color = "blueviolet") + 
  theme_minimal() + 
  labs(x = "Altitude (m)", y = "Shannon-Wiener Index",
       title = "Linear reggression")

#_______________________________________________________________________

###Eurytopic X Stenotopic OTUs
#Connecting sample names with altitude#
elevation_data <- data.frame(
  sample = c("030", "029", "028", "009","010", "027", "008", "026", "016", "007","024","003","004","006","001","023","014","002","013","012","022","011","021","020","005","019","018","017"),
  elevation = c(850,	900,	950,	980,	1100,	1190,	1460,	1630, 1670,	1920,	1920,	1960,	1970,	1980,	2040,	2050,	2150,	2190,	2230,	2250,	2360,	2400,	2460,	2490,	2500,	2590,	2620,	2655)
)
elevation_data$zone <- cut(elevation_data$elevation, 
                           breaks = c(800, 1700, 2200, 2700), 
                           labels = c("forest", "transient", "alpine"), 
                           include.lowest = TRUE)

#adding OTU table to the created zones
library(tidyr)
#transfering the OTU table to the long formate
otu_long<-pivot_longer(OTU.tab1, cols=starts_with("0"),names_to="sample", values_to="count")
otu_long<-merge(otu_long, elevation_data, by="sample")

#now we need to filter if the OTU is present
library(dplyr)
presence_absence<-otu_long %>%
  group_by(No., zone) %>%
  summarise(present = any(count > 0), .groups = "drop")


#now count in how many zones is each OTU present
otu_per_zone <- presence_absence %>% 
  group_by(zone) %>% 
  summarise(unique_otu = sum(present))

print(otu_per_zone)

#OTUs present in all zones
otu_all_zones<-presence_absence %>%
  group_by(No.) %>%
  summarise(zones_present = sum(present)) %>%
  filter(zones_present ==3)

print(otu_all_zones)

#visualization of unique OTUs
library(ggplot2)
ggplot(otu_per_zone, aes(x = zone, y = unique_otu, fill = zone)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Number of unique OTUs per altitudinal zone",
       x = "Altitudina zone", y = "Number of unique OTUs")


#transfering the presence/absence of OTU to binary (1=present, 0=absent)
otu_binary<-OTU.tab1
otu_binary<-otu_binary[,-1]
otu_binary<-map_df(otu_binary, as.numeric)

otu_binary[otu_binary > 0] <- 1
sample_names<-read_excel("sample_names.xlsx")
otu_binary <- cbind(Sample_ID = sample_names, otu_binary)
otu_binary <- pivot_longer(otu_binary, cols = -sample, names_to = "sample", values_to = "abundance")
otu_binary <- otu_binary %>%
  pivot_longer(cols = -species, names_to = "sample", values_to = "abundance")

otu_binary<-merge(otu_binary, elevation_data, by="sample")

#merge by zones, so each OTU has the information about presence in each zone
otu_binary <- otu_binary %>%
  left_join(elevation_data %>% select(sample, zone), by = "sample")

otu_zones <- otu_binary %>%
  group_by(species, zone.x) %>%
  summarise(present = max(abundance), .groups = "drop") %>%
  pivot_wider(names_from = zone.x, values_from = present, values_fill = 0)

otu_zones$zones_present <- rowSums(otu_zones[, -1])

#now, create difference between stenotopic and eurytopic ones
stenotopic <- otu_zones %>% filter(zones_present == 1)
View(stenotopic)
eurytopic <- otu_zones %>% filter(zones_present > 1)

#create a unique sets for each zone
unique_forest <- otu_zones %>% filter(forest == 1 & zones_present == 1)
unique_transient<- otu_zones %>% filter(transient == 1 & zones_present == 1)
unique_alpine <- otu_zones %>% filter(alpine == 1 & zones_present == 1)

#create a sets concerning the shared OTUs
shared_forest_transient <- otu_zones %>% filter(forest == 1 & transient == 1 & zones_present == 2)
shared_forest_alpine <- otu_zones %>% filter(forest == 1 & alpine == 1 & zones_present == 2)
shared_transient_alpine <- otu_zones %>% filter(transient == 1 & alpine == 1 & zones_present == 2)
shared_all<-otu_zones %>% filter(forest == 1 & transient == 1 & alpine ==1 & zones_present == 3)

#making a Venn diagram
library(VennDiagram)
venn.plot <- draw.triple.venn(
  area1 = nrow(unique_forest) + nrow(shared_forest_transient) + nrow(shared_forest_alpine) + nrow(shared_all),
  area2 = nrow(unikátní_transient) + nrow(shared_forest_transient) + nrow(shared_transient_alpine) + nrow(shared_all),
  area3 = nrow(unique_alpine) + nrow(shared_forest_alpine) + nrow(shared_transient_alpine) + nrow(shared_all),
  n12 = nrow(shared_forest_transient) + nrow(shared_all),
  n23 = nrow(shared_transient_alpine) + nrow(shared_all),
  n13 = nrow(shared_forest_alpine) + nrow(shared_all),
  n123 = nrow(shared_all),
  category = c("Forest", "Transient", "Alpine"),
  fill = c("green", "blue", "grey"),
  alpha = 0.5,
  cex = 1.5
)

#_______________________________________________________________________

###PCoA###
#for PCoA we need the ape package
library(ape)
library(vegan)

otu_data_pcoa<-otu_data[,-1]
#Bray-Curtis dissimilarity counting
otu_dist<-vegdist(otu_data_pcoa, method="bray")

#PCoA analysis
pcoa_result<-pcoa(otu_dist)

#write down of the main result
head(pcoa_result$vectors)
summary(pcoa_result)

#upload datataset "vse2" containing metadata
vse2<-read_excel("vse2.xlsx")

vse2$sample <- sub("SAMPLE", "", vse2$sample)


#creating dataframe for the graph
pcoa_scores <- as.data.frame(pcoa_result$vectors)
pcoa_scores$zone <- vse2$zone
pcoa_scores$sample <- vse2$sample

#we use ggplot2 for visualization
library(ggplot2)

ggplot(pcoa_scores, aes(x = Axis.1, y = Axis.2, color = zone)) +
  geom_point(size = 3) +
  labs(
    title = "PCoA – Bray-Curtis",
    x = paste0("Axis 1 (", round(pcoa_result$values$Relative_eig[1] * 100, 1), "%)"),
    y = paste0("Axis 2 (", round(pcoa_result$values$Relative_eig[2] * 100, 1), "%)")
  ) +
  theme_minimal()
pcoa_result$values$Relative_eig

ggplot(pcoa_scores, aes(x = Axis.1, y = Axis.2, color = zone)) +
  geom_point(size = 3) +
  labs(
    title = "PCoA – Bray-Curtis",
    x = paste0("Axis 1 (", round(pcoa_result$values$Relative_eig[1] * 100, 1), "%)"),
    y = paste0("Axis 2 (", round(pcoa_result$values$Relative_eig[2] * 100, 1), "%)")
  ) +
  scale_color_manual(values = c(
    "forest" = "green3",
    "transient" = "coral",
    "alpine" = "steelblue1"
  )) +
  theme_minimal()

#_____________________________________________________________________________________

###NMDS###
#loading the otu table, it is the same one we already used in the PCOA
nmds_result<-metaMDS(otu_data_pcoa, distance="bray", k=2, trymax= 100)

#visualization
plot(nmds_result, type="n")
points(nmds_result, col=as.factor(vse2$zone), pch=19)
ordiellipse(nmds_result, groups=vse2$zone, col=1:3)
legend("topright", legend = levels(as.factor(vse2$zone)), col=1:3, 19)

#_____________________________________________________________________________________

###Taxonomy###

##Number of OTU reads based on the class

#loading taxonomy data
tax_data<-read_excel("taxonomy_data.xlsx")

#count number of reads for each class
tax_data %>%
  count(class, sort = TRUE)

#delete the specific rows with incertae sedis OTUs
tax_data_filtered<- tax_data %>%
  filter(!is.na(class),
         class != "Chlorophyta_cls_Incertae_sedis")

#visualization
tax_data_filtered %>%
  count(class) %>%
  ggplot(aes(x = reorder(class, n), y = n)) +
  geom_col(fill = "forestgreen") +
  coord_flip() +
  labs(title = "Number of OTUs basd on class",
       x = "Class",
       y = "Number of occurences") +
  theme_minimal()


###Heatmap of the 60 most abundant OTUs###
#Rebranding the OTUs to characters
otu_sums <- otu_top %>%
  column_to_rownames("cluster") %>%
  rowSums(na.rm = TRUE)

top60_ids <- names(sort(otu_sums, decreasing = TRUE))[1:60]

#Picking out the relevant OTUs
otu_top <- OTUtabu %>%
  filter(cluster %in% top60_ids)

top60_ids<-
  
  #Connecting names of species to OTUs  
  species_labels <- taxonomy_top60 %>%
  filter(cluster %in% top60_ids) %>%
  mutate(label = paste0(genus, " (", cluster, ")")) %>%
  arrange(match(cluster, otu_top$cluster)) %>%
  pull(label)

#Transfering the table to matrix for heatmap
otu_top_matrix <- otu_top %>%
  column_to_rownames("cluster") %>%
  as.matrix()

mode(otu_top_matrix) <- "numeric"
rownames(otu_top_matrix) <- species_labels

#Renaming the "Sample" variable
vse2$sample <- paste0("SAMPLE", vse2$sample)


#Altitude zones annotation colors
annotation_colors <- list(
  zone = c(
    "forest" = "forestgreen",
    "transient" = "deeppink",
    "alpine" = "dodgerblue"
  )
)

ann_col <- vse2 %>%
  filter(sample %in% colnames(otu_top_matrix)) %>%
  select(sample, zone) %>%
  column_to_rownames(var = "sample")


#Visualization
pheatmap(
  otu_top_matrix,
  scale = "row",
  annotation_col = ann_col,
  annotation_colors = annotation_colors,
  color = colorRampPalette(c("white", "orange", "red"))(100),
  main = "Heatmapa 60 nejhojnějších ",
  fontsize_row = 6
)



###Venn Diagram###
library("purrr")
library("readxl")
library("ggplot2")
library(tidyr)
library(dplyr)
library(VennDiagram)


#Connecting sample names with altitude#
elevation_data <- data.frame(
  sample = c("030", "029", "028", "009","010", "027", "008", "026", "016", "007","024","003","004","006","001","023","014","002","013","012","022","011","021","020","005","019","018","017"),
  elevation = c(850,	900,	950,	980,	1100,	1190,	1460,	1630, 1670,	1920,	1920,	1960,	1970,	1980,	2040,	2050,	2150,	2190,	2230,	2250,	2360,	2400,	2460,	2490,	2500,	2590,	2620,	2655)
)
elevation_data$zone <- cut(elevation_data$elevation, 
                           breaks = c(800, 1700, 2200, 2700), 
                           labels = c("forest", "transient", "alpine"), 
                           include.lowest = TRUE)
elevation_data1<-data_frame(metadata$elevation, metadata$sample)
elevation_datat1<-metadata[,-4]
elevation_datat1<-elevation_datat1[,-3]

#adding OTU table to the created zones

#transfering the OTU table to the long formate
otu_long<-pivot_longer(OTU.tab1, cols=starts_with("0"),names_to="sample", values_to="count")
otu_long<-merge(otu_long, elevation_data, by="sample")



#transfering the OTU table to the long formate
otu_long<-pivot_longer(OTU.tab1, cols=starts_with("0"),names_to="sample", values_to="count")
otu_long<-merge(otu_long, elevation_data, by="sample")

#now we need to filter if the OTU is present
presence_absence<-otu_long %>%
  group_by(No., zone) %>%
  summarise(present = any(count > 0), .groups = "drop")


#now count in how many zones is each OTU present
otu_per_zone <- presence_absence %>% 
  group_by(zone) %>% 
  summarise(unique_otu = sum(present))

print(otu_per_zone)

#OTUs present in all zones
otu_all_zones<-presence_absence %>%
  group_by(No.) %>%
  summarise(zones_present = sum(present)) %>%
  filter(zones_present ==3)

print(otu_all_zones)

#visualization of unique OTUs
ggplot(otu_per_zone, aes(x = zone, y = unique_otu, fill = zone)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Number of unique OTUs per altitudinal zone",
       x = "Altitudina zone", y = "Number of unique OTUs")


#transfering the presence/absence of OTU to binary (1=present, 0=absent)
otu_binary<-OTU.tab1
otu_binary<-otu_binary[,-1]
otu_binary<-map_df(otu_binary, as.numeric)

otu_binary[otu_binary > 0] <- 1
sample_names<-read_excel("sample_names.xlsx")
otu_binary <- cbind(Sample_ID = sample_names, otu_binary)
otu_binary <- otu_binary %>%
  pivot_longer(cols = -species, names_to = "sample", values_to = "abundance")

otu_binary<-merge(otu_binary, elevation_data, by="sample")

#merge by zones, so each OTU has the information about presence in each zone
otu_binary <- otu_binary %>%
  left_join(elevation_data %>% select(sample, zone), by = "sample")

otu_zones <- otu_binary %>%
  group_by(species, zone.x) %>%
  summarise(present = max(abundance), .groups = "drop") %>%
  pivot_wider(names_from = zone.x, values_from = present, values_fill = 0)

otu_zones$zones_present <- rowSums(otu_zones[, -1])

#now, create difference between stenotopic and eurytopic ones
stenotopic <- otu_zones %>% filter(zones_present == 1)
View(stenotopic)
eurytopic <- otu_zones %>% filter(zones_present > 1)

#create a unique sets for each zone
unique_forest <- otu_zones %>% filter(forest == 1 & zones_present == 1)
unique_transient<- otu_zones %>% filter(transient == 1 & zones_present == 1)
unique_alpine <- otu_zones %>% filter(alpine == 1 & zones_present == 1)

#create a sets concerning the shared OTUs
shared_forest_transient <- otu_zones %>% filter(forest == 1 & transient == 1 & zones_present == 2)
shared_forest_alpine <- otu_zones %>% filter(forest == 1 & alpine == 1 & zones_present == 2)
shared_transient_alpine <- otu_zones %>% filter(transient == 1 & alpine == 1 & zones_present == 2)
shared_all<-otu_zones %>% filter(forest == 1 & transient == 1 & alpine ==1 & zones_present == 3)

#making a Venn diagram
venn.plot <- draw.triple.venn(
  area1 = nrow(unique_forest) + nrow(shared_forest_transient) + nrow(shared_forest_alpine) + nrow(shared_all),
  area2 = nrow(unique_transient) + nrow(shared_forest_transient) + nrow(shared_transient_alpine) + nrow(shared_all),
  area3 = nrow(unique_alpine) + nrow(shared_forest_alpine) + nrow(shared_transient_alpine) + nrow(shared_all),
  n12 = nrow(shared_forest_transient) + nrow(shared_all),
  n23 = nrow(shared_transient_alpine) + nrow(shared_all),
  n13 = nrow(shared_forest_alpine) + nrow(shared_all),
  n123 = nrow(shared_all),
  category = c("Forest", "Transient", "Alpine"),
  fill = c("green", "blue", "grey"),
  alpha = 0.5,
  cex = 1.5
)