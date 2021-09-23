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
chart_folder <- "corr_charts"

################################################################################
# Functions
################################################################################

##### FUNCTION FROM SEAN HOBAN ####

#Look at correlations among columns
#This function calculates correlations, outputs top correlated metrics
col_corr<-function(df_matrix,cols_data){
	print("Correlation among <<measures>> themselves")
	print(sort(rowSums(cor(df_matrix[,cols_data],use="complete.obs")>.80),decreasing=T)[1:4])
	print(sort(rowMeans(cor(df_matrix[,cols_data],use="complete.obs")),decreasing=T)[1:4])
	print("Correlation among <<ranks>>- when put in order from measures")
	print(sort(rowSums(cor(apply(df_matrix[,cols_data],2,rank))>.80),decreasing=T)[1:4])
	print(sort(rowMeans(cor(apply(df_matrix[,cols_data],2,rank))),decreasing=T)[1:4])
	chart.Correlation(df_matrix[,cols_data], histogram = TRUE, method = "pearson")
	#To also add the one on ranking?
}

####### VISUALIZATION FROM https://ibecav.github.io/slopegraph/ ######

# fomatting for the slope graph
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

################################################################################
# MANUAL CHANGES REQUIRED: Assign values for scoring
################################################################################

## read in endangerment matrix
df_raw <- read.csv(file.path(main_dir,"endangerment_matrix_forR.csv"),
  header = T, na.strings = c("","NA"))
df_raw <- df_raw %>% arrange(species)
df <- df_raw
colnames(df)

## create dataframe for assigning scores to categorical columns
categories <- c("AB","PR",
                "EW","CR","EN","DD","VU","NT","DD*","LC","NE",
                "Y","N",
                "C-A","C-B","C-C","C-D","C-E1","C-E2","C-E3","C-E4",
                "P-A1","P-A2","P-A3","P-B","P-A4","P-C","P-D","P-E",
                "N/A")
category_scores <- c(1, 0,
                    1, 0.835, 0.668, 0.501, 0.334, 0.167, 0.167, 0, 0,
                    1, 0,
                    1, 0.5, 0.5, 0.5, 0.1, 0.1, 0.1, 0,
                    1, 0.5, 0.5, 0.5, 0.1, 0.1, 0.1, 0,
                    NA)
vals <- data.frame(categories,category_scores)
vals

## create vector for assigning scores to quartiles 1-4
quartile_scores <- c(0, 0, 0.5, 1)

## select columns that need reverse log transformation and quantile scoring
log_cols <- c(5,6,7,8) # ex situ collections data
quantile_cols <- c(9,11) # pest/disease and climate change overall scores

## assign weight to each column (must all add up to 1)
#colnames(df)
  # "species"                  XX"presence_absence"XX        "rl_category"
  # "nativity_to_us"           "exsitu_sites_plantsearch" XX"exsitu_sites_survey"XX
  # "exsitu_wz_sites"          "exsitu_wz_accessions"     XX"climate_change_score"XX
  # "climate_change_vul_class" XX"pest_disease_score"XX       "pest_disease_vul_class"
col_wt <- c(0, 0, 0.3, 0.05, 0.1, 0, 0.05, 0.25, 0, 0.125, 0, 0.125, 0.125, 0.125)
#if(sum(col_wt)!=1){ print("ERROR: !!THE COLUMN WEIGHTS MUST ADD UP TO ONE!!")}
#print(sum(col_wt))

## create vector of column weights if all are weighted evenly
col_using <- 7 # number of columns you're using (non-zero weights)
1/col_using # use this in vector below, for all non-zero weighted columns
col_wt_even <- c(0, 0, 0.1429, 0.1429, 0.1429, 0, 0.1429, 0.1429, 0, 0.1429, 0, 0.1429, 0.1429, 0.1429)

## correlations to run
  # all columns
all_col <- 2:12
  # ex situ data
exsitu_col <- 5:8
  # selected columns
sel_col <- c(3:5,7:8,10,12)
sel_col_means <- c(3:5,7:8,13:14)

## sensitivity tests to run
  # drop nativity
no_nativity <- c(3,5,7:8,13:14)
  # drop Potter
no_potter <- c(3:5,7:8)
  # drop wz sites
no_wzsite <- c(3:5,8,13:14)
  # drop wz accessions
no_wzacc <- c(3:5,7,13:14)

##

################################################################################
# Set up endangerment matrix scoring
################################################################################

# calculate correlations among exsitu data columns (raw values)
chart_path = file.path(main_dir,chart_folder,"Exsitu_cols-correlation_matrix.png")
png(height = 1000, width = 1000, file = chart_path, type = "cairo")
col_corr(df,exsitu_col)
dev.off()

# convert qualitative values to scores using ref vals created above
for(i in 1:nrow(vals)){
  df[df == vals[i,1]] <- vals[i,2]
}
# make everything numeric except species name
df[,2:ncol(df)] <- df[,2:ncol(df)] %>% mutate_if(is.character,as.numeric)
str(df)

## add columns with Potter mean filled in for NA
df$climate_change_vul_class_means <- df$climate_change_vul_class
df$climate_change_vul_class_means[is.na(df$climate_change_vul_class_means)] <-
  mean(df$climate_change_vul_class[!is.na(df$climate_change_vul_class)])
df$pest_disease_vul_class_means <- df$pest_disease_vul_class
df$pest_disease_vul_class_means[is.na(df$pest_disease_vul_class_means)] <-
  mean(df$pest_disease_vul_class[!is.na(df$pest_disease_vul_class)])

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
  has_val[which(is.na(has_val))] <- FALSE
  quartiles <- quantile(df[has_val,quantile_cols[i]],na.rm=T)
  print(quartiles)
  print(mean(df[has_val,quantile_cols[i]]))
  for(j in 1:nrow(df)){
    if (df[j,quantile_cols[i]]>=quartiles[4] & !is.na(df[j,quantile_cols[i]])){
      df[j,quantile_cols[i]] <- quartile_scores[4]
    } else if (df[j,quantile_cols[i]]>=quartiles[3] & df[j,quantile_cols[i]]<quartiles[4] & !is.na(df[j,quantile_cols[i]])){
      df[j,quantile_cols[i]] <- quartile_scores[3]
    } else if (df[j,quantile_cols[i]]>=quartiles[2] & df[j,quantile_cols[i]]<quartiles[3] & !is.na(df[j,quantile_cols[i]])){
      df[j,quantile_cols[i]] <- quartile_scores[2]
    } else if (df[j,quantile_cols[i]]<quartiles[2] & !is.na(df[j,quantile_cols[i]])){
      df[j,quantile_cols[i]] <- quartile_scores[1]
    }
  }
}
str(df)

# calculate correlations among all columns once scores are filled in
chart_path = file.path(main_dir,chart_folder,"All_cols_scored-correlation_matrix.png")
png(height = 1500, width = 1500, file = chart_path, type = "cairo")
col_corr(df,all_col)
dev.off()

################################################################################
# Calculate total score using all columns
################################################################################

# look at correlations for selected columns
chart_path = file.path(main_dir,chart_folder,"Selected_cols_scored-correlation_matrix.png")
png(height = 1000, width = 1000, file = chart_path, type = "cairo")
col_corr(df,sel_col_means)
dev.off()
  # for each genus
  gen_df <- df %>% filter(grepl("Malus ",species))
    chart_path = file.path(main_dir,chart_folder,"Malus-correlation_matrix.png")
    png(height = 1000, width = 1000, file = chart_path, type = "cairo")
    col_corr(gen_df,sel_col_means)
    dev.off()
  gen_df <- df %>% filter(grepl("Quercus ",species))
    chart_path = file.path(main_dir,chart_folder,"Quercus-correlation_matrix.png")
    png(height = 1000, width = 1000, file = chart_path, type = "cairo")
    col_corr(gen_df,sel_col_means)
    dev.off()
  gen_df <- df %>% filter(grepl("Tilia ",species))
    chart_path = file.path(main_dir,chart_folder,"Tilia-correlation_matrix.png")
    png(height = 1000, width = 1000, file = chart_path, type = "cairo")
    col_corr(gen_df,sel_col_means)
    dev.off()
  gen_df <- df %>% filter(grepl("Ulmus ",species))
    chart_path = file.path(main_dir,chart_folder,"Ulmus-correlation_matrix.png")
    png(height = 1000, width = 1000, file = chart_path, type = "cairo")
    col_corr(gen_df,sel_col_means)
    dev.off()

# calculate total score for each species (FUNCTION)
  ## first way, using mean for NA in pest/climate columns and scaling 1-100
calc_total_score_Means <- function(selected_cols,df,col_wt){
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
  ## second way, using division to remove pest/climate columns if NA; no scale
calc_total_score_NoNA <- function(selected_cols,df,col_wt){
  total_score <- vector(mode="numeric", length=nrow(df))
  for(i in 1:nrow(df)){
    temp <- 0
    na_count <- 0
    for(j in 1:length(selected_cols)){
      if(is.na(df[i,selected_cols[j]])){
        na_count <- na_count + 1
      } else {
        temp <- sum(temp, df[i,selected_cols[j]]*col_wt[selected_cols[j]])
      }
    }
      # divide score by number of values considers (removes NA columns)
    temp <- temp/(length(selected_cols)-na_count)
    total_score[i] <- temp#*100
  }
  return(total_score)
}

  ## calculate using desired columns
    # Potter means for NA
total_WeightsAsIs_MeanNA <- calc_total_score_Means(sel_col_means,df,col_wt)
df <- cbind(df,total_WeightsAsIs_MeanNA)
    # remove NA columns (divide sum by number of non-NA columns)
total_WeightsAsIs_NoNA <- calc_total_score_NoNA(sel_col,df,col_wt)
df <- cbind(df,total_WeightsAsIs_NoNA)
    # Potter with zeros for NA
df$pest_disease_vul_class[is.na(df$pest_disease_vul_class)] <- 0
df$climate_change_vul_class[is.na(df$climate_change_vul_class)] <- 0
total_WeightsAsIs_ZeroNA <- calc_total_score_Means(sel_col,df,col_wt)
df <- cbind(df,total_WeightsAsIs_ZeroNA)

########
## can look at difference between these ranks
df_test <- df
col_to_rank <- df_test[,15:17]
  # assign ranks
for(i in 1:length(col_to_rank)){
  col_ranked <- frankv(col_to_rank[i],ties.method="min",order=-1)
  df_test <- cbind(df_test,col_ranked)
}
str(df_test)
  # view slope graph
    # set up dataframe for graphing
df_graph1 <- data.frame(
    Test = "1_MeanNA",
    Species = df_test[,1],
    Rank = df_test[,18])
df_graph2 <- data.frame(
    Test = "2_RemoveNA",
    Species = df_test[,1],
    Rank = df_test[,19])
df_graph3 <- data.frame(
    Test = "3_MeanNA",
    Species = df_test[,1],
    Rank = df_test[,18])
df_graph4 <- data.frame(
    Test = "4_ZeroNA",
    Species = df_test[,1],
    Rank = df_test[,20])
view_ranks <- Reduce(rbind,list(df_graph1,df_graph2,df_graph3,df_graph4))
view_ranks$Abbr <- paste0(substr(view_ranks$Species, 0, 1),".",substr(sub(".* ","",view_ranks$Species), 0, 3))
head(view_ranks)
    # create chart
ggplot(data = view_ranks, aes(x = Test, y = Rank, group = Species)) +
  geom_line(aes(color = Species, alpha = 1), size = 1) +
  geom_label(aes(label = Abbr),
             size = 2,
             label.padding = unit(0.02, "lines"),
             label.size = 0.0) +
  MySpecial +
  labs(title = "Sensitivity Analysis: Visualization of Rank Changes in the Endangerment Matrix")
########


################################################################################
# Sensitivity analysis
################################################################################

# see how much the weighting affects the scores

  # all columns, even weights
total_EvenWeights <- calc_total_score_Means(sel_col_means,df,col_wt_even)
df <- cbind(df,total_EvenWeights)

# calculate score after individually dropping columns we're unsure about

  # drop nativity
total_NoNativity <- calc_total_score_Means(no_nativity,df,col_wt)
df <- cbind(df,total_NoNativity)
  # drop Potter
total_NoPotter <- calc_total_score_Means(no_potter,df,col_wt)
df <- cbind(df,total_NoPotter)
  # drop Potter overall scores only
#total_NoOverallPotter <- calc_total_score(col_use,df,col_wt)
#df <- cbind(df,total_NoOverallPotter)
  # drop Potter vulnerability categories only
#total_NoVulernPotter <- calc_total_score(col_use,df,col_wt)
#df <- cbind(df,total_NoVulernPotter)
  # drop wz sites
total_NoWZSite <- calc_total_score_Means(no_wzsite,df,col_wt)
df <- cbind(df,total_NoWZSite)
  # drop wz accessions
total_NoWZAcc <- calc_total_score_Means(no_wzacc,df,col_wt)
df <- cbind(df,total_NoWZAcc)

## convert total scores to ranks
col_to_rank <- df[,(ncol(df)-(ncol(df)-ncol(df_raw))+3):ncol(df)]
for(i in 1:length(col_to_rank)){
  col_ranked <- frankv(col_to_rank[i],ties.method="min",order=-1)
  df <- cbind(df,col_ranked)
}
str(df)

# write file
write.csv(df, file.path(main_dir,
  "EndangermentMatrix_SensitivityAnalysis_10-22-21.csv"), row.names = F)


## FROM SEAN HOBAN:
#####################
#	RANKING			#
#####################

all_ranks<-df[,23:30]
colnames(all_ranks)<-colnames(df[,15:22])
sp_names_wcoll<-df[,1]
#Could take the mean or majority decision...
#There are two ways to get agreement across all of them
#One way to actually rank species is to identify those that most frequently are ranked in a given bunch, say in the top 10
species_ranked1<-data.frame(sp_names_wcoll,rowSums(all_ranks<10))
#Another way to do actually rank species is to take the mean across rows in the rank order
species_ranked2<-data.frame(sp_names_wcoll,rowMeans(all_ranks))
#species_ranked<-data.frame(sp_names_wcoll,rank(rowMeans(apply(eco_geo_results[,2:(ncol(eco_geo_results)-2)],2,rank))))
colnames(species_ranked1)<-c("sp","rank"); colnames(species_ranked2)<-c("sp","rank")
#It really doesn't matter which approach :)
cbind(species_ranked1[order(species_ranked1$rank),],species_ranked2[order(species_ranked2$rank,decreasing=T),])

#But let's examine more closely the individual columns and how they differ
#Could look at those that might be most different and figure out why
#Examine them by eye
all_ranks <- cbind(sp_names_wcoll,all_ranks)
cbind(all_ranks[order(all_ranks[,2]),1],
	    all_ranks[order(all_ranks[,3]),1],
	    all_ranks[order(all_ranks[,4]),1],
	    all_ranks[order(all_ranks[,5]),1],
	    all_ranks[order(all_ranks[,6]),1],
	    all_ranks[order(all_ranks[,7]),1],
	    all_ranks[order(all_ranks[,8]),1],
	    all_ranks[order(all_ranks[,9]),1]
    )

#TO DO ADD IN EMILY CODE FOR LINES

#How to get the difference in new rank from old rank
#match(species_ranked_geo50[order(species_ranked_geo50$rank),1],species_ranked_geo100[order(species_ranked_geo100$rank),1])-1:41
count_changes<-function(list_ranks,base=1){
	base_order<-list_ranks[[base]][order(list_ranks[[base]]$rank),1]
	examine_changes<-data.frame(base_order)
	for (i in 1:length(list_ranks)){
		ranks_diff<-match(base_order,list_ranks[[i]][order(list_ranks[[i]]$rank),1])
		examine_changes<-cbind(examine_changes,match(base_order,list_ranks[[i]][order(list_ranks[[i]]$rank),1])-1:length(base_order))
		colnames(examine_changes)[i+1]<-names(list_ranks[[i]][2])
	}
	#examine_changes<-examine_changes[,-2]
	examine_changes
}
#compare to geo50
examine_changes<-count_changes(list(data.frame("sp"=all_ranks[,1],"rank_WeightsAsIs_MeanNA"=as.numeric(all_ranks[,2])),
                                    data.frame("sp"=all_ranks[,1],"rank_WeightsAsIs_NoNA"=as.numeric(all_ranks[,3])),
                                    data.frame("sp"=all_ranks[,1],"rank_WeightsAsIs_ZeroNA"=as.numeric(all_ranks[,4])),
                                    data.frame("sp"=all_ranks[,1],"rank_EvenWeights"=as.numeric(all_ranks[,5])),
                                    data.frame("sp"=all_ranks[,1],"rank_NoNativity"=as.numeric(all_ranks[,6])),
                                    data.frame("sp"=all_ranks[,1],"rank_NoPotter"=as.numeric(all_ranks[,7])),
                                    data.frame("sp"=all_ranks[,1],"rank_NoWZSite"=as.numeric(all_ranks[,8])),
                                    data.frame("sp"=all_ranks[,1],"rank_NoWZAcc"=as.numeric(all_ranks[,9]))
                                  ))
colSums(abs(examine_changes[,-1])>5)
colSums(abs(examine_changes[,-1])>10)
colSums(abs(examine_changes[,-1])>20)








## OLD CHUNKS

#  # use mean for species with no Potter data
#    # select potter columns from raw data
#df2 <- df_raw[,c(1,9:12)]
#    # create categories again but with N/A = 999
#categories <- c("C-A","C-B","C-C","C-D","C-E1","C-E2","C-E3","C-E4",
#                "P-A1","P-A2","P-A3","P-B","P-A4","P-C","P-D","P-E",
#                "N/A")
#category_scores <- c(1, 0.5, 0.5, 0.5, 0.1, 0.1, 0.1, 0,
#                    1, 0.5, 0.5, 0.5, 0.1, 0.1, 0.1, 0,
#                    999)
#vals <- data.frame(categories,category_scores)
#    # replace categories with scores
#for(i in 1:nrow(vals)){
#  df2[df2 == vals[i,1]] <- vals[i,2]
#}
#  # calculate correlations among Potter columns that have scores
#have_Potter_val <- df2[which(df2$climate_change_score!="999"),]
#have_Potter_val <- have_Potter_val %>% mutate_if(is.character,as.numeric)
#cols_data <- 2:5
#col_corr(have_Potter_val)
#    ### for categorical cols, replace N/A with mean score for the col
#df2[df2[,3]=="999",3] <- round(mean(as.numeric(df2[df2[,3]!="999",3])),3)
#df2[df2[,5]=="999",5] <- round(mean(as.numeric(df2[df2[,5]!="999",5])),5)
#    # make everything numeric except species name
#df2[,2:ncol(df2)] <- df2[,2:ncol(df2)] %>% mutate_if(is.character,as.numeric); str(df2)
    # convert quantitative cols to scores (0-1) using quartiles
#quantile_cols <- c(2,4)
#for(i in 1:length(quantile_cols)){
#  has_val <- df2[,quantile_cols[i]]<100
#  quartiles <- quantile(df2[has_val,quantile_cols[i]])
#    ### for quantitative cols, replace N/A with mean value for the col
#  df2[!has_val,quantile_cols[i]] <- mean(df2[has_val,quantile_cols[i]])
#  for(j in 1:nrow(df2)){
#    if (df2[j,quantile_cols[i]]>=quartiles[4]){
#      df2[j,quantile_cols[i]] <- quartile_scores[4]
#    } else if (df2[j,quantile_cols[i]]>=quartiles[3] & df2[j,quantile_cols[i]]<quartiles[4]){
#      df2[j,quantile_cols[i]] <- quartile_scores[3]
#    } else if (df2[j,quantile_cols[i]]>=quartiles[2] & df2[j,quantile_cols[i]]<quartiles[3]){
#      df2[j,quantile_cols[i]] <- quartile_scores[2]
#    } else if (df2[j,quantile_cols[i]]<quartiles[2]){
#      df2[j,quantile_cols[i]] <- quartile_scores[1]
#    }
#  }
#}
#    # add new potter scores to other data
#df_new <- full_join(df[,1:8],df2[,c(1,3,5)])
#df <- cbind(df_new,df[13:18])
#str(df)
#colnames(df)
#    # calculate total scores
#total_PotterMeans <- calc_total_score(col_use,df,col_wt)
#df <- cbind(df,total_PotterMeans)
  # drop Potter overall scores only
#col_use <- c(2:8,10,12)
#total_NoOverallPotterMeans <- calc_total_score(col_use,df2,col_wt)
#df <- cbind(df,total_NoOverallPotterMeans)
  # drop Potter vulnerability categories only
#col_use <- c(2:9,11)
#total_NoPotterMeans <- calc_total_score(col_use,df2,col_wt)
#df <- cbind(df,total_NoPotterMeans)

#  # switch DD score with VU score (so DD is below VU not above)
#df$rl_category[which(df$rl_category==0.501)] <- 999
#df$rl_category[which(df$rl_category==0.334)] <- 0.501
#df$rl_category[which(df$rl_category==999)] <- 0.334
#    # calculate total scores
#total_SwapDDVU <- calc_total_score(col_use,df,col_wt)
#df <- cbind(df,total_SwapDDVU)

##### OPTIONAL GRAPHING #####

# set up dataframe for graphing
#tests <- c(
#  "A_STANDARD",
#  "B_EvenWeights",
#  "C_STANDARD",
#  "D_NoNativity",
#  "E_STANDARD",
#  "F_NoPotter",
#  "G_STANDARD",
#  "H_NoWZSite",
#  "I_STANDARD",
#  "J_NoWZAcc",
#  "K_STANDARD",
#  "L_PotterMeans",
#  "M_STANDARD",
#  "N_SwapDDVU"
#)
#colnames(df)
#standard <- data.frame(
#    Test = "STANDARD",
#    Species = df[,1],
#    Rank = df[,19]) # first column with ranks
#view_ranks <- data.frame()
#num_tests <- 7 # number of non-"STANDARD" tests
#for(i in 1:num_tests){
#  print(tests[i+i])
#  add <- data.frame(
#    Test = tests[i+i],
#    Species = df[,1],
#    Rank = df[,19+i]) # first column with ranks
#  standard$Test <- tests[i+i-1]
#  view_ranks <- Reduce(rbind,list(view_ranks,standard,add))
#}
#view_ranks$Abbr <- paste0(substr(view_ranks$Species, 0, 1),".",substr(sub(".* ","",view_ranks$Species), 0, 3))
#head(view_ranks)

# create and save chart
#chart_path = file.path(main_dir,chart_folder,"Slope_graph.png")
#  png(height = 1800, width = 1800, file = chart_path, type = "cairo")

#ggplot(data = view_ranks, aes(x = Test, y = Rank, group = Species)) +
#  geom_line(aes(color = Species, alpha = 1), size = 1) +
#  geom_point(aes(color = Type, alpha = .1), size = 4) +
#  #geom_text_repel(data = view_ranks %>% filter(Test == "SwapDDVU"),
#  #                aes(label = Abbr) ,
#  #                hjust = "left",
#  #                fontface = "bold",
#  #                size = 3,
#  #                nudge_x = -.45,
#  #                direction = "y") +
#  #geom_text_repel(data = view_ranks %>% filter(Test == "NoVulernPotterMeans"),
#  #                aes(label = Abbr) ,
#  #                hjust = "right",
#  #                fontface = "bold",
#  #                size = 3,
#  #                nudge_x = .5,
#  #                direction = "y") +
#  geom_label(aes(label = Abbr),
#             size = 2,
#             label.padding = unit(0.02, "lines"),
#             label.size = 0.0) +
#  MySpecial +
#  labs(
#    title = "Sensitivity Analysis: Visualization of Rank Changes in the Endangerment Matrix"#,
#    #subtitle = "Subtitle",
#    #caption = "Code from https://ibecav.github.io/slopegraph/"
#  )
##dev.off()
