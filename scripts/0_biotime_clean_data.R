library(reshape2)
library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table) 
library(dplyr)
library(terra)
library(here)
library(glmmTMB)

# Add your folder
project_folder="path/to/downloaded_biotime_data.csv"
output_folder="path/to/outputfiles"

# Read data
rawdata = fread(paste0(project_folder,"biotime_v2_rawdata_2025.csv")) 
metadata = fread(paste0(project_folder,"biotime_v2_metadata_2025.csv"))

# PLOTS #

# First filter by Taxa and by REALM
taxon_selection = c("Birds", "Mammals", "Invertebrates", "Reptiles", "Amphibians")

meta_inv = metadata %>%
  filter(TAXA %in% taxon_selection & REALM == "Terrestrial")

studyid = meta_inv$STUDY_ID
fil_data = rawdata %>% filter(STUDY_ID %in% studyid)
taxon_df = fread (paste0(project_folder, "data/taxon_df_all.csv"))

fil_data_tax <- fil_data %>%
  left_join(taxon_df, by = c("valid_name" = "canonicalName"))

fil_data_tax = fil_data_tax %>% filter(!kingdom=="Plantae", rank=="SPECIES") %>%
  select(-c(usageKey, status, confidence, matchType, kingdomKey, speciesKey, phylumKey,
            classKey, orderKey, familyKey, genusKey, speciesKey, acceptedUsageKey))


## add latitude and longitude

add_plot_id <- function(data, class_type, cluster_size) {
  
  data_class = data %>% filter (class==class_type)
  spat_plots = vect (unique(data_class[,c("STUDY_ID","LATITUDE","LONGITUDE")]), geom = c("LONGITUDE","LATITUDE"),crs="+proj=longlat +datum=WGS84")
  spat_plots$PLOT_ID = NA
  studies = sort (unique(spat_plots$STUDY_ID))
  
  time1=Sys.time()
  for(i in studies){
    bidon=spat_plots[spat_plots$STUDY_ID==i,]
    if(nrow(bidon)>1){
      mati=distance(bidon)
      hc=hclust(mati)
      clus=cutree(hc, h = cluster_size)
      spat_plots$PLOT_ID[spat_plots$STUDY_ID==i]=clus
    }
  }
  time2=Sys.time()
  time2-time1
  
  coords=as.data.frame(geom(spat_plots))
  tab_plots=cbind(as.data.frame(spat_plots),data.frame(LATITUDE=coords$y,LONGITUDE=coords$x))
  data_class=merge(data_class,tab_plots,by=c("STUDY_ID","LATITUDE","LONGITUDE"),all=TRUE)
  #add a column to identify study and plots
  data_class$SP_ID = paste0(data_class$STUDY_ID,"-",data_class$PLOT_ID)
  return(data_class)
}

squamata_data   <- add_plot_id (data = fil_data_tax, class_type="Squamata", cluster_size = 2e3)
aves_data   <- add_plot_id (data = fil_data_tax, class_type="Aves", cluster_size = 5e3)
mammal_data <- add_plot_id (data = fil_data_tax, class_type="Mammalia", cluster_size = 5e3)
insect_data <- add_plot_id (data = fil_data_tax, class_type="Insecta", cluster_size = 1e3)
arac_data   <- add_plot_id (data = fil_data_tax, class_type="Arachnida", cluster_size = 1e3)
amph_data   <- add_plot_id (data = fil_data_tax, class_type="Amphibia", cluster_size = 2e3)

data_plot <- rbind(squamata_data, aves_data, mammal_data, 
                            insect_data, arac_data, amph_data)

data_plot$MODEL_ID = paste(data_plot$SP_ID, data_plot$valid_name, sep = "-")

# FILTERS #####

filter_data <- function(data) {
  
  # Filter consecutive years and duration
  meta_ts <- data %>%
    group_by(STUDY_ID, PLOT_ID, SP_ID, class) %>%
    summarise(
      Duration = n_distinct(YEAR),
      Min_year = min(YEAR),
      Max_year = max(YEAR),
      Consecutive = all(diff(sort(unique(
        YEAR
      ))) == 1)
    ) %>%
    filter (Duration >= 3 & Consecutive == TRUE)
  
  
  data = subset (data, SP_ID %in% meta_ts$SP_ID)
  data = subset (data, resolution == "species")
  
  # Infer the 0s
  mat <- dcast(
    data,
    STUDY_ID + PLOT_ID + SP_ID + YEAR + class ~ valid_name,
    value.var = "ABUNDANCE",
    fun.aggregate = sum,
    fill = 0
  )
  
  df_with_zeros = melt(
    mat,
    id.vars = c("STUDY_ID", "PLOT_ID", "SP_ID", "YEAR", "class"),
    variable.name = "valid_name",
    value.name = "ABUNDANCE"
  )
  # remove species all zeros across all years in a study
  df_with_zeros <- df_with_zeros %>%
    group_by(valid_name, STUDY_ID, SP_ID) %>%
    mutate(abund_sum = sum(ABUNDANCE)) %>%
    filter(abund_sum > 0) %>%
    ungroup()
  
  #filter species that have been seen in at least 2 times & with at least 50% of non zeros values
  df_with_zeros <- df_with_zeros %>%
    group_by(valid_name, STUDY_ID, SP_ID) %>%
    mutate(
      nb_year_detect = length(unique(YEAR[ABUNDANCE > 0])),
      perc_zeros = length(unique(YEAR[ABUNDANCE > 0])) / length(unique(YEAR))
    ) %>%
    filter(nb_year_detect > 1 & perc_zeros >= 0.5)
  
  df_with_zeros$MODEL_ID = paste(df_with_zeros$SP_ID,
                                 df_with_zeros$valid_name,
                                 sep = "-")
  
  length(unique(df_with_zeros$MODEL_ID))
  return(df_with_zeros)
  
}

squamat_data_filtered = filter_data(data=squamata_data)
aves_data_filtered = filter_data(data=aves_data)
mamm_data_filtered = filter_data(data=mammal_data)
inse_data_filtered = filter_data(data=insect_data)
arac_data_filtered = filter_data(data=arac_data)
amph_data_filtered = filter_data(data=amph_data)

data_plot_filtered <- rbind(squamat_data_filtered, aves_data_filtered, mamm_data_filtered, 
                       inse_data_filtered, arac_data_filtered, amph_data_filtered)

# centralize years in YEAR2  
data_plot_filtered = data_plot_filtered %>% 
  group_by(valid_name, MODEL_ID, class, STUDY_ID) %>%
  summarise (YEAR2 = YEAR-min(YEAR), 
             YEAR = YEAR,
             ABUNDANCE = ABUNDANCE)

data_plot_filtered$group <- factor(1)
data_plot_filtered$MODEL_ID


# Remove 17828 birds' time series to reduce computational time

set.seed(8)

calculation = data_plot_filtered %>% group_by(class) %>% 
  summarise(n_distinct(MODEL_ID))
n = calculation[3,2]-calculation[1,2]
n=n$`n_distinct(MODEL_ID)`
modelsaves = data_plot_filtered %>% filter (class=="Aves") %>% distinct (MODEL_ID)
models_aves_out = sample(modelsaves$MODEL_ID, n, replace=F)
data_biotime_models = data_plot_filtered %>% filter (!MODEL_ID %in% models_aves_out)


# Save data
fwrite(data_biotime_models, paste0(output_folder, "/data_filtered.csv"))

conversor= data_plot %>% select(LATITUDE, LONGITUDE, MODEL_ID) %>% distinct()
fwrite(conversor, paste0(output_folder, "/conversor.csv"))

