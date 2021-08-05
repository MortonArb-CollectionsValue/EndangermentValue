################################################################################

## endangerment_sensitivity_and_correlations.R
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

rm(list=ls())
  my.packages <- c('tidyverse')
  select <- dplyr::select
# install.packages (my.packages) #Turn on to install current versions
lapply(my.packages, require, character.only=TRUE)
  rm(my.packages)

################################################################################
# Set working directory
################################################################################

main_dir <- "/Volumes/GoogleDrive/Shared drives/IMLS MFA/Endangerment Value"

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
colnames(df)
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

# convert qualitative values to scores using ref vals created above
for(i in 1:nrow(vals)){
  df[df == vals[i,1]] <- vals[i,2]
}
# make everything numeric except species name
df[,2:ncol(df)] <- df[,2:ncol(df)] %>% mutate_if(is.character,as.numeric)
str(df)

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
    } else if (df[j,quantile_cols[i]]>=quartiles[3]){
      df[j,quantile_cols[i]] <- quartile_scores[3]
    } else if (df[j,quantile_cols[i]]>=quartiles[2]){
      df[j,quantile_cols[i]] <- quartile_scores[2]
    } else {
      df[j,quantile_cols[i]] <- quartile_scores[1]
    }
  }
}
str(df)

################################################################################
# Calculate total score using all columns
################################################################################

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

  # drop presence/absence
col_use <- c(3:12)
total_noPresAbs <- calc_total_score(col_use,df,col_wt_even)
df <- cbind(df,total_noPresAbs)
  # drop nativity
col_use <- c(2:3,5:12)
total_noNativity <- calc_total_score(col_use,df,col_wt_even)
df <- cbind(df,total_noNativity)
  # drop Potter
col_use <- c(2:8)
total_noPotter <- calc_total_score(col_use,df,col_wt_even)
df <- cbind(df,total_noPotter)
  # drop Potter overall scores only
col_use <- c(2:8,10,12)
total_noOverallPotter <- calc_total_score(col_use,df,col_wt_even)
df <- cbind(df,total_noOverallPotter)
  # drop Potter vulnerability categories only
col_use <- c(2:9,11)
total_noVulernPotter <- calc_total_score(col_use,df,col_wt_even)
df <- cbind(df,total_noVulernPotter)






  # use mean for species with no Potter data

df2 <- df_raw[,c(1,9:12)]

for(i in 1:nrow(vals)){
  df2[df2 == vals[i,1]] <- vals[i,2]
}

df2 <- df2

df2[df2[,2]=="N/A",2] <- round(mean(as.numeric(df2[df2[,2]!="N/A",2])),2)
df2[df2[,4]=="N/A",4] <- round(mean(as.numeric(df2[df2[,4]!="N/A",4])),2)

# make everything numeric except species name
df2[,2:ncol(df2)] <- df2[,2:ncol(df2)] %>% mutate_if(is.character,as.numeric)
str(df2)
  # Climate change and pest/disease score columns:
  #   use scores entered at the beginning of the script
for(i in 1:length(quantile_cols)){
  has_val <- df2[,quantile_cols[i]]>0
  quartiles <- quantile(df2[has_val,quantile_cols[i]])
  df2[!has_val,quantile_cols[i]] <- mean(df2[has_val,quantile_cols[i]])
  for(j in 1:nrow(df2)){
    if (df2[j,quantile_cols[i]]>=quartiles[4]){
      df2[j,quantile_cols[i]] <- quartile_scores[4]
    } else if (df2[j,quantile_cols[i]]>=quartiles[3]){
      df2[j,quantile_cols[i]] <- quartile_scores[3]
    } else if (df2[j,quantile_cols[i]]>=quartiles[2]){
      df2[j,quantile_cols[i]] <- quartile_scores[2]
    } else {
      df2[j,quantile_cols[i]] <- quartile_scores[1]
    }
  }
}
str(df2)





df$climate_change_score[which(df$climate_change_score==0)] <- 0.5
df$climate_change_vul_class[which(df$climate_change_vul_class==0)] <- 0.5

df$pest_disease_score[which(df$pest_disease_score==0)] <- 0

col_use <- c(2:9,11)
total_PotterMeans <- calc_total_score(col_use,df,col_wt_even)
df <- cbind(df,total_PotterMeans)
