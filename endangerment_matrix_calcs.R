################################################################################

## endangerment_matrix_calcs.R
### Authors: Emily Beckman ### Date: 08/03/2021

### DESCRIPTION:
  # This script

### DATA IN:
  #

### DATA OUT:
  #

################################################################################
# Load libraries
################################################################################

# this little code chunk is from Shannon M Still !
  rm(list=ls())
  my.packages <- c('tidyverse','PerformanceAnalytics','ggplot2','ggrepel',
    'data.table')
  select <- dplyr::select
  # install.packages (my.packages) #Turn on to install current versions
  lapply(my.packages, require, character.only=TRUE)
    rm(my.packages)

################################################################################
# Set working directory
################################################################################

main_dir <- "/Volumes/GoogleDrive/Shared drives/IMLS MFA/Endangerment Value"

################################################################################
# Functions
################################################################################

##### FUNCTION FROM SEAN HOBAN ####

#Look at correlations among columns
#This function calculates correlations, outputs top correlated metrics
col_corr<-function(df_matrix){
	print("Correlation among <<measures>> themselves")
	print(sort(rowSums(cor(df_matrix[,cols_data],use="complete.obs")>.80),decreasing=T)[1:4])
	print(sort(rowMeans(cor(df_matrix[,cols_data],use="complete.obs")),decreasing=T)[1:4])
	print("Correlation among <<ranks>>- when put in order from measures")
	print(sort(rowSums(cor(apply(df_matrix[,cols_data],2,rank))>.80),decreasing=T)[1:4])
	print(sort(rowMeans(cor(apply(df_matrix[,cols_data],2,rank))),decreasing=T)[1:4])
	chart.Correlation(df_matrix[,cols_data], histogram = TRUE, method = "pearson")
	#To also add the one on ranking?
}

################################################################################
# MANUAL CHANGES REQUIRED: Assign values for scoring
################################################################################

## create dataframe for assigning scores to categorical columns
categories <- c("AB","PR",
                "EW","CR","EN","DD","VU","NT","LC","NE?",
                "Y","N",
                "C-A","C-B","C-C","C-D","C-E1","C-E2","C-E3","C-E4",
                "P-A1","P-A2","P-A3","P-B","P-A4","P-C","P-D","P-E",
                "N/A")
category_scores <- c(1, 0,
                    1, 0.835, 0.668, 0.501, 0.334, 0.167, 0, 0,
                    1, 0,
                    1, 0.5, 0.5, 0.5, 0.1, 0.1, 0.1, 0,
                    1, 0.5, 0.5, 0.5, 0.1, 0.1, 0.1, 0,
                    0)
vals <- data.frame(categories,category_scores)
vals

## create vector for assigning scores to quartiles 1-4
quartile_scores <- c(0, 0, 0.5, 1)

## select columns that need reverse log transformation and quantile scoring
log_cols <- c(5,6,7,8)
quantile_cols <- c(9,11)

## assign weight to each column (must all add up to 1)
#colnames(df)
col_wt <- c(0, 0.05, 0.2, 0.05, 0.1, 0.1, 0.2, 0.1, 0.1, 0, 0.1, 0)
if(sum(col_wt)!=1){ print("ERROR: !!THE COLUMN WEIGHTS MUST ADD UP TO ONE!!")}
print(sum(col_wt))

################################################################################
# Set up endangerment matrix scoring
################################################################################

# read in endangerment matrix
df_raw <- read.csv(file.path(main_dir,"endangerment_matrix_forR.csv"),
  header = T, na.strings = c("","NA"))
df <- df_raw

  ########## look at adding TRY Database (trait) data
    # read in TRY data
  trydb <- fread(file.path(main_dir,"TRY_Database_download",
    "16564_01092021010239","16564.txt"),
    header = T, sep = "\t", dec = ".", quote = "", data.table = T)
    nrow(trydb) #2636596
    # keep just records for target species and traits that have numbers (not
    #   sure where the other data came from?)
  spp <- unique(df$species)
  try_target <- trydb %>%
    filter(SpeciesName %in% spp) %>%
    filter(!is.na(TraitID)) %>%
    arrange(SpeciesName)
  nrow(try_target) #2764
  head(as.data.frame(try_target))
    # summarize available data
  try_summary <- data.frame()
  for(i in 1:length(unique(try_target$DataName))){
    trait_cols <- unique(try_target[,c(11,13)])[i]
    unq_vals <- unique(sort(try_target[
      which(try_target$DataName==unique(try_target$DataName)[i]),]$OrigValueStr))
    num_spp <- length(unique(try_target[
      which(try_target$DataName==unique(try_target$DataName)[i]),]$SpeciesName))
    which_spp <- unique(sort(try_target[
      which(try_target$DataName==unique(try_target$DataName)[i]),]$SpeciesName))
    add <- data.frame(
      TraitName = trait_cols[,1],
      DataName = trait_cols[,2],
      percent_target_sp = round(num_spp/length(spp)*100,2),
      which_target_sp = paste(which_spp,sep='',collapse='; '),
      unique_values = paste(unq_vals,sep='',collapse='; '))
    try_summary <- rbind(try_summary,add)
  }
  try_summary <- try_summary %>% arrange(TraitName,desc(percent_target_sp))
  try_summary[,1:3]
    # write file
  write.csv(try_summary, file.path(main_dir,"TRY_data_summary.csv"),
    row.names = F)
  ##########

# calculate correlations among exsitu data columns (raw values)
cols_data <- 5:8
col_corr(df)

# convert qualitative values to scores using ref vals created above
for(i in 1:nrow(vals)){
  df[df == vals[i,1]] <- vals[i,2]
}
# make everything numeric except species name
df[,2:ncol(df)] <- df[,2:ncol(df)] %>% mutate_if(is.character,as.numeric)
str(df)

# calculate correlations among all columns once scores are filled in
cols_data <- 2:12
col_corr(df)

## convert quantitative values to scores using equations
  # Collections-based columns:
  #   Log-transformed then scaled in reverse from 1 (no collections) to 0
  #   (max num collections across all species)
for(i in 1:length(log_cols)){
  ln_max <- max(log(df[,log_cols[i]]+1))
  ln_min <- min(log(df[,log_cols[i]]+1))
  for(j in 1:nrow(df)){
    df[j,log_cols[i]] <- 1-((log(df[j,log_cols[i]]+1)-ln_min)/(ln_max-ln_min))
  }
}
  # Climate change and pest/disease score columns:
  #   use scores entered at the beginning of the script
for(i in 1:length(quantile_cols)){
  has_val <- df[,quantile_cols[i]]>0
  quartiles <- quantile(df[has_val,quantile_cols[i]])
  print(quartiles)
  print(mean(df[has_val,quantile_cols[i]]))
  for(j in 1:nrow(df)){
    if (df[j,quantile_cols[i]]>=quartiles[4]){
      df[j,quantile_cols[i]] <- quartile_scores[4]
    } else if (df[j,quantile_cols[i]]>=quartiles[3] & df[j,quantile_cols[i]]<quartiles[4]){
      df[j,quantile_cols[i]] <- quartile_scores[3]
    } else if (df[j,quantile_cols[i]]>=quartiles[2] & df[j,quantile_cols[i]]<quartiles[3]){
      df[j,quantile_cols[i]] <- quartile_scores[2]
    } else if (df[j,quantile_cols[i]]<quartiles[2]){
      df[j,quantile_cols[i]] <- quartile_scores[1]
    }
  }
}
str(df)

################################################################################
# Calculate total score using all columns
################################################################################

# select only columns of interest
df <- df %>% select(species,rl_category,nativity_to_us,exsitu_sites_plantsearch,
  exsitu_wz_sites,exsitu_wz_accessions,climate_change_vul_class,
  pest_disease_vul_class)
colnames(df)
# look at correlations for just these selected columns
cols_data <- 2:8
col_corr(df)
  # for each genus
  gen_df <- df %>% filter(grepl("Malus ",species))
  col_corr(gen_df)
  gen_df <- df %>% filter(grepl("Quercus ",species))
  col_corr(gen_df)
  gen_df <- df %>% filter(grepl("Tilia ",species))
  col_corr(gen_df)
  gen_df <- df %>% filter(grepl("Ulmus ",species))
  col_corr(gen_df)

# create vector of column weights if all are weighted evenly
col_wt_even <- rep(x = 1/(ncol(df)-1), times = ncol(df))

# calculate total score for each species
  ## FUNCTION
calc_total_score <- function(selected_cols,df,col_wt){
  total_score <- vector(mode="numeric", length=nrow(df))
  for(i in 1:nrow(df)){
    temp <- 0
    for(j in 1:length(selected_cols)){
      temp <- sum(temp, df[i,selected_cols[j]]*col_wt[selected_cols[j]])
    }
    total_score[i] <- temp*100
  }
  return(total_score)
}
  ## calculate using all columns
col_use <- c(2:12)
total_WeightsAsIs <- calc_total_score(col_use,df,col_wt)
df <- cbind(df,total_WeightsAsIs)

################################################################################
# Sensitivity analysis
################################################################################

# calculate score after individually dropping columns we're unsure about
  # all columns, even weights
total_EvenWeights <- calc_total_score(col_use,df,col_wt_even)
df <- cbind(df,total_EvenWeights)

  # drop nativity
col_use <- c(2,4:8)
total_NoNativity <- calc_total_score(col_use,df,col_wt_even)
df <- cbind(df,total_NoNativity)
  # drop Potter
col_use <- c(2:6)
total_NoPotter <- calc_total_score(col_use,df,col_wt_even)
df <- cbind(df,total_NoPotter)
  # drop Potter overall scores only
#col_use <- c(2:8,10,12)
#total_NoOverallPotter <- calc_total_score(col_use,df,col_wt_even)
#df <- cbind(df,total_NoOverallPotter)
  # drop Potter vulnerability categories only
#col_use <- c(2:9,11)
#total_NoVulernPotter <- calc_total_score(col_use,df,col_wt_even)
#df <- cbind(df,total_NoVulernPotter)
  # drop wz accessions
col_use <- c(2:5,7:8)
total_NoWZAcc <- calc_total_score(col_use,df,col_wt_even)
df <- cbind(df,total_NoWZAcc)

  # use mean for species with no Potter data
    # select potter columns from raw data
df2 <- df_raw[,c(1,9:12)]
    # create categories again but with N/A = 999
categories <- c("C-A","C-B","C-C","C-D","C-E1","C-E2","C-E3","C-E4",
                "P-A1","P-A2","P-A3","P-B","P-A4","P-C","P-D","P-E",
                "N/A")
category_scores <- c(1, 0.5, 0.5, 0.5, 0.1, 0.1, 0.1, 0,
                    1, 0.5, 0.5, 0.5, 0.1, 0.1, 0.1, 0,
                    999)
vals <- data.frame(categories,category_scores)
    # replace categories with scores
for(i in 1:nrow(vals)){
  df2[df2 == vals[i,1]] <- vals[i,2]
}

  # calculate correlations among Potter columns that have scores
#have_Potter_val <- df2[which(df2$climate_change_score!="999"),]
#have_Potter_val <- have_Potter_val %>% mutate_if(is.character,as.numeric)
#cols_data <- 2:5
#col_corr(have_Potter_val)

    ### for categorical cols, replace N/A with mean score for the col
df2[df2[,3]=="999",3] <- round(mean(as.numeric(df2[df2[,3]!="999",3])),3)
df2[df2[,5]=="999",5] <- round(mean(as.numeric(df2[df2[,5]!="999",5])),5)
    # make everything numeric except species name
df2[,2:ncol(df2)] <- df2[,2:ncol(df2)] %>% mutate_if(is.character,as.numeric); str(df2)
    # convert quantitative cols to scores (0-1) using quartiles
quantile_cols <- c(2,4)
for(i in 1:length(quantile_cols)){
  has_val <- df2[,quantile_cols[i]]<100
  quartiles <- quantile(df2[has_val,quantile_cols[i]])
    ### for quantitative cols, replace N/A with mean value for the col
  df2[!has_val,quantile_cols[i]] <- mean(df2[has_val,quantile_cols[i]])
  for(j in 1:nrow(df2)){
    if (df2[j,quantile_cols[i]]>=quartiles[4]){
      df2[j,quantile_cols[i]] <- quartile_scores[4]
    } else if (df2[j,quantile_cols[i]]>=quartiles[3] & df2[j,quantile_cols[i]]<quartiles[4]){
      df2[j,quantile_cols[i]] <- quartile_scores[3]
    } else if (df2[j,quantile_cols[i]]>=quartiles[2] & df2[j,quantile_cols[i]]<quartiles[3]){
      df2[j,quantile_cols[i]] <- quartile_scores[2]
    } else if (df2[j,quantile_cols[i]]<quartiles[2]){
      df2[j,quantile_cols[i]] <- quartile_scores[1]
    }
  }
}
    # add new potter scores to other data
df2 <- cbind(df[,1:8],df2[,2:5])
df2 <- cbind(df2,df[,13:19])
str(df2)
    # calculate total scores
col_use <- c(2:12)
total_PotterMeans <- calc_total_score(col_use,df2,col_wt_even)
df <- cbind(df,total_PotterMeans)
  # drop Potter overall scores only
col_use <- c(2:8,10,12)
total_NoOverallPotterMeans <- calc_total_score(col_use,df2,col_wt_even)
df <- cbind(df,total_NoOverallPotterMeans)
  # drop Potter vulnerability categories only
col_use <- c(2:9,11)
total_NoVulernPotterMeans <- calc_total_score(col_use,df2,col_wt_even)
df <- cbind(df,total_NoVulernPotterMeans)

  # switch DD score with VU score (so DD is below VU not above)
df$rl_category[which(df$rl_category==0.501)] <- 999
df$rl_category[which(df$rl_category==0.334)] <- 0.501
df$rl_category[which(df$rl_category==999)] <- 0.334
    # calculate total scores
col_use <- c(2:12)
total_SwapDDVU <- calc_total_score(col_use,df,col_wt_even)
df <- cbind(df,total_SwapDDVU)

## convert total scores to ranks
col_to_rank <- df[,13:23]
for(i in 1:length(col_to_rank)){
  col_ranked <- as.vector(rank(col_to_rank[i],ties.method="average"))
  df <- cbind(df,col_ranked)
}
str(df)

# write file
write.csv(df, file.path(main_dir,"EndangermentMatrix_SensitivityAnalysis.csv"),
	row.names = F)








# set up dataframe for graphing
tests <- c("A_WeightsAsIs","B_EvenWeights","C_NoPresAbs","D_NoNativity",
  "E_NoPotter","F_NoOverallPotter","G_NoVulernPotter","H_PotterMeans",
  "I_NoOverallPotterMeans","J_NoVulernPotterMeans","K_SwapDDVU")
view_ranks <- data.frame()
for(i in 1:length(tests)){
  add <- data.frame(
    Test = tests[i],
    Species = df[,1],
    Rank = df[,(length(tests)+12+i)]) #the number here is the number of columns before totals
  view_ranks <- rbind(view_ranks,add)
}
view_ranks$Abbr <- paste0(substr(view_ranks$Species, 0, 1),".",substr(sub(".* ","",view_ranks$Species), 0, 3))
head(view_ranks)

# FROM: https://ibecav.github.io/slopegraph/
MySpecial <- list(
  # move the x axis labels up top
  scale_x_discrete(position = "top"),
  theme_bw(),
  # Format tweaks
  # Remove the legend
  theme(legend.position = "none"),
  # Remove the panel border
  theme(panel.border     = element_blank()),
  # Remove just about everything from the y axis
  theme(axis.title.y     = element_blank()),
  theme(axis.text.y      = element_blank()),
  theme(panel.grid.major.y = element_blank()),
  theme(panel.grid.minor.y = element_blank()),
  # Remove a few things from the x axis and increase font size
  theme(axis.title.x     = element_blank()),
  theme(panel.grid.major.x = element_blank()),
  theme(axis.text.x.top      = element_text(size=12)),
  # Remove x & y tick marks
  theme(axis.ticks       = element_blank()),
  # Format title & subtitle
  theme(plot.title       = element_text(size=14, face = "bold", hjust = 0.5)),
  theme(plot.subtitle    = element_text(hjust = 0.5))
)

ggplot(data = view_ranks, aes(x = Test, y = Rank, group = Species)) +
  geom_line(aes(color = Species, alpha = 1), size = 1) +
#  geom_point(aes(color = Type, alpha = .1), size = 4) +
  #geom_text_repel(data = view_ranks %>% filter(Test == "SwapDDVU"),
  #                aes(label = Abbr) ,
  #                hjust = "left",
  #                fontface = "bold",
  #                size = 3,
  #                nudge_x = -.45,
  #                direction = "y") +
  #geom_text_repel(data = view_ranks %>% filter(Test == "NoVulernPotterMeans"),
  #                aes(label = Abbr) ,
  #                hjust = "right",
  #                fontface = "bold",
  #                size = 3,
  #                nudge_x = .5,
  #                direction = "y") +
  geom_label(aes(label = Abbr),
             size = 2,
             label.padding = unit(0.02, "lines"),
             label.size = 0.0) +
  MySpecial +
  labs(
    title = "Sensitivity Analysis: Visualization of Rank Changes in the Endangerment Matrix",
    #subtitle = "Subtitle",
    caption = "Code from https://ibecav.github.io/slopegraph/"
  )
