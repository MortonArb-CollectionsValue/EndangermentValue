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
total_NoPresAbs <- calc_total_score(col_use,df,col_wt_even)
df <- cbind(df,total_NoPresAbs)
  # drop nativity
col_use <- c(2:3,5:12)
total_NoNativity <- calc_total_score(col_use,df,col_wt_even)
df <- cbind(df,total_NoNativity)
  # drop Potter
col_use <- c(2:8)
total_NoPotter <- calc_total_score(col_use,df,col_wt_even)
df <- cbind(df,total_NoPotter)
  # drop Potter overall scores only
col_use <- c(2:8,10,12)
total_NoOverallPotter <- calc_total_score(col_use,df,col_wt_even)
df <- cbind(df,total_NoOverallPotter)
  # drop Potter vulnerability categories only
col_use <- c(2:9,11)
total_NoVulernPotter <- calc_total_score(col_use,df,col_wt_even)
df <- cbind(df,total_NoVulernPotter)

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






##### FROM SEAN HOBAN ####

#	RANKING

#One way to actually rank species is to identify those that most frequently are ranked in a given bunch, say in the top 10
species_ranked1<-data.frame(sp_names_wcoll,rowSums(apply(eco_geo_results[,2:(ncol(eco_geo_results)-2)],2,rank)<10))
#Another way to do actually rank species is to take the mean across rows in the rank order
species_ranked2<-data.frame(sp_names_wcoll,rowMeans(apply(eco_geo_results[,2:(ncol(eco_geo_results)-2)],2,rank)))
#species_ranked<-data.frame(sp_names_wcoll,rank(rowMeans(apply(eco_geo_results[,2:(ncol(eco_geo_results)-2)],2,rank))))
colnames(species_ranked1)<-c("sp","rank"); colnames(species_ranked2)<-c("sp","rank")
#It really doesn't matter which approach :)
cbind(species_ranked1[order(species_ranked1$rank),],species_ranked2[order(species_ranked2$rank,decreasing=T),])

#Geo 10
species_ranked_geo10<-data.frame(sp_names_wcoll,apply(eco_geo_results[,2:(ncol(eco_geo_results)-2)],2,rank)[,1])
#Geo 100
species_ranked_geo100<-data.frame(sp_names_wcoll,apply(eco_geo_results[,2:(ncol(eco_geo_results)-2)],2,rank)[,4])
#Eco 10
species_ranked_eco10<-data.frame(sp_names_wcoll,apply(eco_geo_results[,2:(ncol(eco_geo_results)-2)],2,rank)[,5])
#Eco USA 50
species_ranked_ecous50<-data.frame(sp_names_wcoll,apply(eco_geo_results[,2:(ncol(eco_geo_results)-2)],2,rank)[,10])
colnames(species_ranked_geo10)<-c("sp","rank-geo10"); colnames(species_ranked_geo100)<-c("sp","rank-geo100")
colnames(species_ranked_eco10)<-c("sp","rank-eco10"); colnames(species_ranked_ecous50)<-c("sp","rank-ecous50")
cbind(species_ranked_geo10[order(species_ranked_geo10$rank),],species_ranked_geo100[order(species_ranked_geo100$rank),],
	species_ranked_eco10[order(species_ranked_eco10$rank),],species_ranked_ecous50[order(species_ranked_ecous50$rank),])
#Maybe 10 is just too fine scale!!


#Trying to figure out if any of the measures are especially bad/ different
#Look for correlations among columns
#This function calculates correlations, outputs top correlated metrics
eco_geo_corr<-function(eg_matrix){
	print("Correlation among <<measures>> themselves")
	print(sort(rowSums(cor(eg_matrix[,cols_eco_geo],use="complete.obs")>.80),decreasing=T)[1:4])
	print(sort(rowMeans(cor(eg_matrix[,cols_eco_geo],use="complete.obs")),decreasing=T)[1:4])
	print("Correlation among <<ranks>>- when put in order from measures")
	print(sort(rowSums(cor(apply(eg_matrix[,cols_eco_geo],2,rank))>.80),decreasing=T)[1:4])
	print(sort(rowMeans(cor(apply(eg_matrix[,cols_eco_geo],2,rank))),decreasing=T)[1:4])
	chart.Correlation(eg_matrix[,cols_eco_geo], histogram = TRUE, method = "pearson")
	#To also add the one on ranking?
}
eco_geo_corr(eco_geo_results)
	#Note: Geo10 and Eco10 are most different- not well correlated
#Correlation among USA species
eco_geo_results_usa<-eco_geo_results[-which(is.na(eco_geo_results[,9])),]
eco_geo_corr(eco_geo_results_usa)
#Correlation just for within rare and not rare sets
eco_geo_results_LC<-eco_geo_results[eco_geo_results[,ncol(eco_geo_results)]=="LC",]
eco_geo_corr(eco_geo_results_LC)
	#Note: For geo10 and ecousa10, LC low correlation (probably too fine scale  for LC
eco_geo_results_TH<-eco_geo_results[eco_geo_results[,ncol(eco_geo_results)]!="LC",]
eco_geo_corr(eco_geo_results_TH)
	#Note: here it is eco10 and eco50 low correlation (probably too fine scale  for LC
#OVERALL: geo50, geo100, eco50, ecoUS50 have highest overall

#But should probably find highest pairs and drop one of each pair... iteratively
#This applies the function to subsets of: all, US, LC, Threatened
lapply(list(eco_geo_results,eco_geo_results_usa,eco_geo_results_LC,eco_geo_results_TH),eco_geo_corr)
