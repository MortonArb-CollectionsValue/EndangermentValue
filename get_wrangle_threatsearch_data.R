################################################################################

## get_wrangle_threatsearch_data.R
### Authors: Emily Beckman ### Date: 09/02/2021

### DESCRIPTION:
  # This script downloads, reads in, and wrangles BGCI ThreatSearch data
  #   <https://tools.bgci.org/threat_search.php>

################################################################################
# Download BGCI's ThreatSearch data by genus
################################################################################

# load packages
library("tidyverse")

# create list of target genera
gen <- c("Malus","Quercus","Tilia","Ulmus")

# set location you want the data saved to
main_dir <- "/Volumes/GoogleDrive/Shared drives/IMLS MFA/Endangerment Value"

# download ThreatSearch data
for(i in 1:length(gen)){
  download.file(paste0("https://tools.bgci.org/export_threatsearch.php?ftrFamily=&ftrGenus=",
    gen[[i]],"&ftrSpecies=&ftrInfraSpecName=&ftrBGCI_Scope=&ftrAssessmentYear=&ftrThreatened=&ftrPagerLimit=100000&action=Find&export=1"),
    destfile=file.path(main_dir,
      paste0(gen[[i]],"_BGCI-ThreatSearch-Results.csv")))
}

################################################################################
# Read in and wrangle ThreatSearch data
################################################################################

# read in
file_list <- list.files(main_dir, pattern = "BGCI-ThreatSearch-Results.csv",
  full.names = T)
file_dfs <- lapply(file_list, read.csv, header = T, na.strings = c("","NA"),
  colClasses = "character")
length(file_dfs) #4

# compile into dataframe
threat_df <- Reduce(rbind, file_dfs)
nrow(threat_df) #1577

# wrangle ThreatSearch data, depending on your purpose
  # separate out taxon name
threat_df <- threat_df %>%
  separate("Plant.Name",c("genus","species","infra_rank","infra_name"),remove=F)
  # create genus_species column
threat_df$genus_species <- paste(threat_df$genus,threat_df$species)
  # standardize Interpreted.Conservation.Status column
threat_df$Interpreted.Conservation.Status <-
  str_squish(threat_df$Interpreted.Conservation.Status)
table(threat_df$Interpreted.Conservation.Status)
  ## DEPENDING ON WHAT YOU'RE LOOKING FOR!...filter the rows
threat_df_sp <- threat_df %>%
  # remove infrataxa and hybrids
  filter(!grepl(" var. | subsp. | ssp. | x | spp. | sp. ",Plant.Name)) %>%
  # remove non-global assessments
  filter(!grepl("Not Global|Not global|Unknown|\\?",Scope))
nrow(threat_df_sp) #1053
  # remove duplicates (keep most recent assessment for each species)
threat_df_sp_unq <- threat_df_sp %>%
  arrange(desc(Assessment.Year),Interpreted.Conservation.Status,Source) %>%
  distinct(genus_species,.keep_all=T) %>%
  select(genus,genus_species,Interpreted.Conservation.Status,
    Published.Conservation.Status,Assessment.Year,Source)
nrow(threat_df_sp_unq) #580
# remove genera not in original list (not sure how they appeared, but they did
#   somehow during ThreatSearch query.. perhaps through ThreatSearch's synonymy?)
threat_df_sp_unq <- threat_df_sp_unq %>%
  filter(genus %in% gen)
nrow(threat_df_sp_unq) #575

# write file
write.csv(threat_df_sp_unq, file.path(main_dir,
  "Wrangled_BGCI_ThreatSearch.csv"), row.names = F)
