
## endangerment_matrix_calcs-ByGenus.R
### Authors: Emily Beckman ### Date: 09/24/2021

### DESCRIPTION:
  # This script

### DATA IN:
  #

### DATA OUT:
  #

################################################################################
# Load libraries
################################################################################

# this code chunk is from Shannon M Still !
  rm(list=ls())
  my.packages <- c('tidyverse','PerformanceAnalytics','ggplot2','ggrepel',
    'data.table','tidytext')
  select <- dplyr::select
  # install.packages (my.packages) #Turn on to install current versions
  lapply(my.packages, require, character.only=TRUE)
    rm(my.packages)

################################################################################
# Set working directory
################################################################################

main_dir <- "/Volumes/GoogleDrive/Shared drives/IMLS MFA/Endangerment Value/script_outputs"
chart_folder <- "corr_charts"

################################################################################
# Functions
################################################################################

### FROM SEAN HOBAN ###
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

### Functions for scoring endangerment matrix ###
# calculate total score for each species
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
  ## -- not using right now --
#calc_total_score_NoNA <- function(selected_cols,df,col_wt){
#  total_score <- vector(mode="numeric", length=nrow(df))
#  for(i in 1:nrow(df)){
#    temp <- 0
#    na_count <- 0
#    for(j in 1:length(selected_cols)){
#      if(is.na(df[i,selected_cols[j]])){
#        na_count <- na_count + 1
#      } else {
#        temp <- sum(temp, df[i,selected_cols[j]]*col_wt[selected_cols[j]])
#      }
#    }
#      # divide score by number of values considers (removes NA columns)
#    temp <- temp/(length(selected_cols)-na_count)
#    total_score[i] <- temp#*100
#  }
#  return(total_score)
#}

# list of target species for genetic dimension
seans_spp <- c("Ulmus rubra","Ulmus americana","Quercus laurifolia",
  "Tilia platyphyllos","Quercus faginea","Tilia americana","Malus angustifolia",
  "Quercus phellos","Quercus laevis","Quercus emoryi","Quercus falcata",
  "Quercus virginiana","Quercus texana","Quercus minima",
  "Quercus muehlenbergii","Quercus turbinella","Quercus vacciniifolia",
  "Malus ioensis","Quercus chrysolepis","Quercus chenii","Tilia japonica",
  "Quercus multinervis","Quercus palmeri","Quercus austrina","Quercus lobata",
  "Quercus castaneifolia","Tilia paucicostata","Quercus sadleriana",
  "Quercus havardii","Quercus arkansana","Quercus oglethorpensis",
  "Malus trilobata","Ulmus chenmoui","Ulmus gaussenii","Ulmus elongata",
  "Quercus pontica","Quercus georgiana","Quercus acerifolia","Malus sieversii",
  "Quercus boyntonii","Malus spontanea")

################################################################################
# MANUAL CHANGES REQUIRED: Assign values for scoring
################################################################################

## read in endangerment matrix
en_df_raw <- read.csv(file.path(main_dir,"endangerment_matrix_forR.csv"),
  header = T, na.strings = c("","NA"))
en_df_raw <- en_df_raw %>% arrange(species)
en_df <- en_df_raw
colnames(en_df)

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
log_cols <- c(6:9) # ex situ collections data
quantile_cols <- c(10,12) # pest/disease and climate change overall scores

## assign weight to each column (must all add up to 1)
#colnames(en_df)
  # "species"     "group"      XX"presence_absence"XX        "rl_category"
  # "nativity_to_us"           "exsitu_sites_plantsearch" XX"exsitu_sites_survey"XX
  # "exsitu_wz_sites"          "exsitu_wz_accessions"     XX"climate_change_score"XX
  # "climate_change_vul_class" XX"pest_disease_score"XX       "pest_disease_vul_class"
col_wt <- c(0, 0, 0, 0.3, 0.05, 0.1, 0, 0.05, 0.25, 0, 0.125, 0, 0.125, 0.125, 0.125)
#if(sum(col_wt)!=1){ print("ERROR: !!THE COLUMN WEIGHTS MUST ADD UP TO ONE!!")}
#print(sum(col_wt))

## create vector of column weights if all are weighted evenly
col_using <- 7 # number of columns you're using (non-zero weights)
1/col_using # use this in vector below, for all non-zero weighted columns
col_wt_even <- c(0, 0, 0, 0.1429, 0.1429, 0.1429, 0, 0.1429, 0.1429, 0, 0.1429, 0, 0.1429, 0.1429, 0.1429)

## correlations to run
  # all columns
all_col <- 3:13
  # ex situ data
exsitu_col <- 6:9
  # selected columns
sel_col <- c(4:6,8:9,14:15) # using Potter vulnerability class with means for NA
  # nicknames for easier viewing of results when we drop each column one-by-one
sel_col_nicknames <- c("no_rl","no_nativity","no_plantsearch","no_wzsites",
  "no_wzacc","no_climchange","no_pestdisease")
#sel_col_PotterZero <- c(4:6,8:9,11,13) # using Potter vul class with 0 for NA

## threashold for using Potter columns (minimum percent filled to use col)
potter_thresh <- 0.15

################################################################################
# START LOOPING THROUGH GENERA/GROUPS
################################################################################

# create list of dataframes, one for each genus or other groupping, for
#   looping through
all <- en_df %>% mutate(group="AllSpp") %>% select(species,group,everything())
malus <- en_df %>% filter(grepl("Malus ",species)) %>% mutate(group="Malus") %>% select(species,group,everything())
quercus <- en_df %>% filter(grepl("Quercus ",species)) %>% mutate(group="Quercus") %>% select(species,group,everything())
tilia <- en_df %>% filter(grepl("Tilia ",species)) %>% mutate(group="Tilia") %>% select(species,group,everything())
ulmus <- en_df %>% filter(grepl("Ulmus ",species)) %>% mutate(group="Ulmus") %>% select(species,group,everything())
seans <- en_df %>% filter(species %in% seans_spp) %>% mutate(group="GeneticDimSpp") %>% select(species,group,everything())
  # list of groups to loop through
genera <- list(all,malus,quercus,tilia,ulmus,seans)

##
#### start looping...
##
for(g in 1:length(genera)){
#g <- 1

  ################################################################################
  # Set up endangerment matrix scoring
  ################################################################################

  # calculate correlations among exsitu data columns (raw values)
  chart_path = file.path(main_dir,chart_folder,
    paste0(gsub(" .*$","",genera[[g]]$group[1]),"-Exsitu_cols-correlation_matrix.png"))
  png(height = 1000, width = 1000, file = chart_path, type = "cairo")
  col_corr(genera[[g]],exsitu_col)
  dev.off()

  # convert qualitative values to scores using ref vals created above
  for(i in 1:nrow(vals)){
    genera[[g]][genera[[g]] == vals[i,1]] <- vals[i,2]
  }
  # make everything numeric except species name and group name
  genera[[g]][,3:ncol(genera[[g]])] <- genera[[g]][,3:ncol(genera[[g]])] %>% mutate_if(is.character,as.numeric)
  str(genera[[g]])

  ## convert quantitative values to scores using equations
    # Collections-based columns:
    #   Log-transformed then scaled in reverse from 1 (no collections) to 0
    #   (max num collections across all species)
  for(i in 1:length(log_cols)){
    ln_max <- max(log(genera[[g]][,log_cols[i]]+1))
    ln_min <- min(log(genera[[g]][,log_cols[i]]+1))
    for(j in 1:nrow(genera[[g]])){
      genera[[g]][j,log_cols[i]] <- 1-((log(genera[[g]][j,log_cols[i]]+1)-ln_min)/(ln_max-ln_min))
    }
  }

  # score Potter columns
    #first check percent of species with data, and compare to threshold you set
    potter_thresh_g <- nrow(genera[[g]][which(!is.na(genera[[g]]$climate_change_vul_class)),]) / nrow(genera[[g]])
  if(potter_thresh < potter_thresh_g){
    # add columns with Potter mean filled in for NA
    genera[[g]]$climate_change_vul_class_means <- genera[[g]]$climate_change_vul_class
    genera[[g]]$climate_change_vul_class_means[is.na(genera[[g]]$climate_change_vul_class_means)] <-
      mean(genera[[g]]$climate_change_vul_class[!is.na(genera[[g]]$climate_change_vul_class)])
    genera[[g]]$pest_disease_vul_class_means <- genera[[g]]$pest_disease_vul_class
    genera[[g]]$pest_disease_vul_class_means[is.na(genera[[g]]$pest_disease_vul_class_means)] <-
      mean(genera[[g]]$pest_disease_vul_class[!is.na(genera[[g]]$pest_disease_vul_class)])
    # Climate change and pest/disease score columns:
    #   use scores entered at the beginning of the script
    for(i in 1:length(quantile_cols)){
      has_val <- genera[[g]][,quantile_cols[i]]>0
      has_val[which(is.na(has_val))] <- FALSE
      quartiles <- quantile(genera[[g]][has_val,quantile_cols[i]],na.rm=T)
      print(quartiles)
      print(mean(genera[[g]][has_val,quantile_cols[i]]))
      for(j in 1:nrow(genera[[g]])){
        if (genera[[g]][j,quantile_cols[i]]>=quartiles[4] & !is.na(genera[[g]][j,quantile_cols[i]])){
          genera[[g]][j,quantile_cols[i]] <- quartile_scores[4]
        } else if (genera[[g]][j,quantile_cols[i]]>=quartiles[3] & genera[[g]][j,quantile_cols[i]]<quartiles[4] & !is.na(genera[[g]][j,quantile_cols[i]])){
          genera[[g]][j,quantile_cols[i]] <- quartile_scores[3]
        } else if (genera[[g]][j,quantile_cols[i]]>=quartiles[2] & genera[[g]][j,quantile_cols[i]]<quartiles[3] & !is.na(genera[[g]][j,quantile_cols[i]])){
          genera[[g]][j,quantile_cols[i]] <- quartile_scores[2]
        } else if (genera[[g]][j,quantile_cols[i]]<quartiles[2] & !is.na(genera[[g]][j,quantile_cols[i]])){
          genera[[g]][j,quantile_cols[i]] <- quartile_scores[1]
        }
      }
    }
  } else { # remove Potter columns if % filled is below threshold
    genera[[g]]$climate_change_vul_class <- 0
    genera[[g]]$pest_disease_vul_class <- 0
    genera[[g]]$climate_change_score <- 0
    genera[[g]]$pest_disease_score <- 0
    genera[[g]]$climate_change_vul_class_means <- 0
    genera[[g]]$pest_disease_vul_class_means <- 0
  }

  str(genera[[g]])

  # calculate correlations among all columns once scores are filled in
  chart_path = file.path(main_dir,chart_folder,
    paste0(gsub(" .*$","",genera[[g]]$group[1]),"-All_cols_scored-correlation_matrix.png"))
  png(height = 1500, width = 1500, file = chart_path, type = "cairo")
  col_corr(genera[[g]],all_col)
  dev.off()

  # look at correlations for selected columns
  chart_path = file.path(main_dir,chart_folder,
        paste0(gsub(" .*$","",genera[[g]]$group[1]),"-Selected_cols_scored-correlation_matrix.png"))
  png(height = 1000, width = 1000, file = chart_path, type = "cairo")
  col_corr(genera[[g]],sel_col)
  dev.off()

  ################################################################################
  # Sensitivity analysis
  ################################################################################

  # see how much the weighting and inclusion/removal of columns affects the scores

  ### CALCULATE TOTAL SCORES FOR DIFFERENT TESTS ###

    # Potter with zeros for NA
  #genera[[g]]$pest_disease_vul_class[is.na(genera[[g]]$pest_disease_vul_class)] <- 0
  #genera[[g]]$climate_change_vul_class[is.na(genera[[g]]$climate_change_vul_class)] <- 0
  #PotterZeroNA <- calc_total_score_Means(sel_col_PotterZero,genera[[g]],col_wt)
  #genera[[g]] <- cbind(genera[[g]],PotterZeroNA)

  ## BASE TEST: Weighted cols and Potter means for NAs
  wt_all <- calc_total_score_Means(sel_col,genera[[g]],col_wt)
  genera[[g]] <- cbind(genera[[g]],wt_all)

  # recalculate total score for different tests, using weighted score
  for(i in 1:length(sel_col)){
    test_col <- sel_col[-i]
    test <- calc_total_score_Means(test_col,genera[[g]],col_wt)
    genera[[g]] <- cbind(genera[[g]],test)
    colnames(genera[[g]])[colnames(genera[[g]])=="test"] <-
      paste0("wt_",sel_col_nicknames[i])
  }
  head(genera[[g]])

  ## BASE TEST 2: Even weight for each col and Potter means for NAs
  ev_all <- calc_total_score_Means(sel_col,genera[[g]],col_wt_even)
  genera[[g]] <- cbind(genera[[g]],ev_all)

  # recalculate total score for different tests, using even weights
  for(i in 1:length(sel_col)){
    test_col <- sel_col[-i]
    test <- calc_total_score_Means(test_col,genera[[g]],col_wt_even)
    genera[[g]] <- cbind(genera[[g]],test)
    colnames(genera[[g]])[colnames(genera[[g]])=="test"] <-
      paste0("ev_",sel_col_nicknames[i])
  }
  head(genera[[g]])

  # visualize tests

    # format dataframe for visualizing
  vis <- data.frame(species="start",test="start",values=0)
  for(i in 1:(length(genera[[g]])-15)){ #subtract num rows before total scores
    add_vis <- data.frame(
                species=genera[[g]][,1],
                test=rep(colnames(genera[[g]][i+15]),n=nrow(genera[[g]])),
                values=genera[[g]][,i+15] #num columns before scores
              )
    vis <- rbind(vis,add_vis)
  }; vis <- vis[-1,]

    ## multi-facet scatter plot
  vis %>%
    #group_by(values) %>%
    #top_n(15) %>%
    #ungroup %>%
    mutate(test = as.factor(test),
           species = reorder_within(species, values, test)) %>%
    ggplot(aes(species, values, fill = test)) +
    geom_point(show.legend = FALSE) +
    facet_wrap(~test, scales = "free_y") +
    coord_flip() +
    scale_x_reordered() +
    scale_y_continuous(expand = c(0,0)) +
    labs(y = "Endangerment Score",
         x = "Species",
         title = "Endangerment Matrix Sensitivity Analysis")
  ggsave(file.path(main_dir,chart_folder,
    paste0(gsub(" .*$","",genera[[g]]$group[1]),"-Selected_cols_scored-sensitivity_analysis_scatterplots.png")),
    width=20,height=10)

    ## jitter plot visualization
  vis %>%
    mutate(test = as.factor(test),
           species = reorder_within(species, values, test)) %>%
    ggplot(aes(x = test, y = values)) +
    geom_jitter(alpha = 0.3, color = "tomato") +
    geom_boxplot(alpha = 0) +
    labs(y = "Endangerment Score",
         x = "Sensitivity Test",
         title = "Endangerment Matrix Sensitivity Analysis")
  ggsave(file.path(main_dir,chart_folder,
    paste0(gsub(" .*$","",genera[[g]]$group[1]),"-Selected_cols_scored-sensitivity_analysis_boxplots.png")),
    width=20,height=10)

  ### ADD RANKS FOR DIFFERENT TESTS ###

  # arrange by 'base' test before adding ranks
  genera[[g]] <- genera[[g]] %>% arrange(desc(wt_all))

  ## convert total scores to ranks
  col_to_rank <- genera[[g]][,(ncol(genera[[g]])-(ncol(genera[[g]])-ncol(en_df_raw))+4):ncol(genera[[g]])]
  for(i in 1:length(col_to_rank)){
    col_ranked <- frankv(col_to_rank[i],ties.method="min",order=-1)
    genera[[g]] <- cbind(genera[[g]],col_ranked)
    colnames(genera[[g]])[colnames(genera[[g]])=="col_ranked"] <-
      paste0(names(col_to_rank[i]),"_r")
  }
  str(genera[[g]])


  ## BASE CODE FROM SEAN HOBAN; edited for this purpose:
  #####################
  #	QUANTIFYING SENSITIVITY	#
  #####################

  all_ranks<-genera[[g]][,32:47]
  colnames(all_ranks)<-colnames(genera[[g]][,16:31])
  species<-genera[[g]][,1]
  #Could take the mean or majority decision...
  #There are two ways to get agreement across all of them
  #One way to actually rank species is to identify those that most frequently are ranked in a given bunch, say in the top 10
  species_ranked1<-data.frame(species,rowSums(all_ranks<10))
  #Another way to do actually rank species is to take the mean across rows in the rank order
  species_ranked2<-data.frame(species,rowMeans(all_ranks))
  #species_ranked<-data.frame(species,rank(rowMeans(apply(eco_geo_results[,2:(ncol(eco_geo_results)-2)],2,rank))))
  colnames(species_ranked1)<-c("sp","rank"); colnames(species_ranked2)<-c("sp","rank")
  #It really doesn't matter which approach :)
  cbind(species_ranked1[order(species_ranked1$rank),],species_ranked2[order(species_ranked2$rank,decreasing=T),])

  #But let's examine more closely the individual columns and how they differ
  #Could look at those that might be most different and figure out why
  #Examine them by eye
  #all_ranks <- cbind(species,all_ranks)
    #cbind(all_ranks[order(all_ranks[,2]),1],
    #	    all_ranks[order(all_ranks[,3]),1],
    #	    all_ranks[order(all_ranks[,4]),1],
    #	    all_ranks[order(all_ranks[,5]),1],
    #	    all_ranks[order(all_ranks[,6]),1],
    #	    all_ranks[order(all_ranks[,7]),1],
    #	    all_ranks[order(all_ranks[,8]),1],
    #	    all_ranks[order(all_ranks[,9]),1]
    #    )
  # calculate number of changes based on different thresholds

  ### count changes compared to base test

    # list of dataframes, one for each test; input for analysis
  list_for_analysis <- list()
  for(i in 1:(length(all_ranks))){
    t <- data.frame("sp"=species,"temp"=as.numeric(all_ranks[,i]))
    colnames(t)[colnames(t)=="temp"] <- paste0("rank_",names(all_ranks[i]))
    list_for_analysis[[i]] <- t
  }
  examine_changes<-count_changes(list_for_analysis)
  ## use number of species to calculate cutoffs (percents) for ranks shifting up and down
  num_spp <- nrow(genera[[g]])
  cutoff1 <- round((num_spp * 0.05),0)
  cutoff2 <- round((num_spp * 0.10),0)
  cutoff3 <- round((num_spp * 0.15),0)
  cutoff4 <- round((num_spp * 0.20),0)

  rank_changes <- cbind(
    greater.than.5 = colSums(abs(examine_changes[,-1])>cutoff1),
    greater.than.10 = colSums(abs(examine_changes[,-1])>cutoff2),
    greater.than.15 = colSums(abs(examine_changes[,-1])>cutoff3),
    greater.than.20 = colSums(abs(examine_changes[,-1])>cutoff4)
  ); rank_changes

  write.csv(rank_changes, file.path(main_dir,
    paste0(gsub(" .*$","",genera[[g]]$group[1]),
    "-SensitivityAnalysis_rank_change_calcs.csv")), row.names = T)

  ### count number of times a species shows up in the top 'bunch'
  ncol(all_ranks)

  # top 5 spots
  species_ranked1<-data.frame(species,rowSums(all_ranks<5))
  head(species_ranked1[order(species_ranked1[,2],decreasing=T),],n=20)

  # top 10 spots
  species_ranked2<-data.frame(species,rowSums(all_ranks<10))
  head(species_ranked2[order(species_ranked2[,2],decreasing=T),],n=30)

  # top 20 spots
  species_ranked3<-data.frame(species,rowSums(all_ranks<20))
  head(species_ranked3[order(species_ranked3[,2],decreasing=T),],n=50)

  # join all together for saving
  genera[[g]] <- Reduce(full_join,list(genera[[g]],species_ranked1,species_ranked2,species_ranked3))
  head(genera[[g]])

  # write file
  write.csv(genera[[g]], file.path(main_dir,
    paste0(gsub(" .*$","",genera[[g]]$group[1]),
    "-EndangermentMatrix_SensitivityAnalysis.csv")), row.names = F)

}








## OLD BIT FOR CREATING SLOPE GRAPH VISUALIZATION

### VISUALIZATION HELPER FUNCTION FROM https://ibecav.github.io/slopegraph/ ###
# fomatting for the slope graph
#MySpecial <- list(
#  # move the x axis labels up top
#  scale_x_discrete(position = "top"),
#  theme_bw(),
#  # Format tweaks
#  # Remove the legend
#  theme(legend.position = "none"),
#  # Remove the panel border
#  theme(panel.border     = element_blank()),
#  # Remove just about everything from the y axis
#  theme(axis.title.y     = element_blank()),
#  theme(axis.text.y      = element_blank()),
#  theme(panel.grid.major.y = element_blank()),
#  theme(panel.grid.minor.y = element_blank()),
#  # Remove a few things from the x axis and increase font size
#  theme(axis.title.x     = element_blank()),
#  theme(panel.grid.major.x = element_blank()),
#  theme(axis.text.x.top      = element_text(size=12)),
#  # Remove x & y tick marks
#  theme(axis.ticks       = element_blank()),
#  # Format title & subtitle
#  theme(plot.title       = element_text(size=14, face = "bold", hjust = 0.5)),
#  theme(plot.subtitle    = element_text(hjust = 0.5))
#)

#  ########
#  ## visualize sensitivity
#  en_df_test <- genera[[g]]
#    # view slope graph
#      # set up dataframe for graphing
#  graph1 <- data.frame(
#      Test = "1_WeightsAsIs",
#      Species = en_df_test[,1],
#      Rank = en_df_test[,23])
#  graph2 <- data.frame(
#      Test = "2_EvenWeights",
#      Species = en_df_test[,1],
#      Rank = en_df_test[,26])
#  graph3 <- data.frame(
#      Test = "1_WeightsAsIs",
#      Species = en_df_test[,1],
#      Rank = en_df_test[,23])
#  graph4 <- data.frame(
#      Test = "2_NoNativity",
#      Species = en_df_test[,1],
#      Rank = en_df_test[,27])
#view_ranks <- rbind(graph1,graph2)
#view_ranks <- rbind(graph3,graph4)
#  view_ranks <- Reduce(rbind,list(graph1,graph2,graph3,graph4))
#  view_ranks$Abbr <- paste0(substr(view_ranks$Species, 0, 1),".",substr(sub(".* ","",view_ranks$Species), 0, 3))
#  head(view_ranks)
#      # create chart
#        # I have not bothered to save these because the labels don't show
#        #   up in the saved version :(
#  ggplot(data = view_ranks, aes(x = Test, y = Rank, group = Species)) +
#    geom_line(aes(color = Species, alpha = 1), size = 1) +
#    geom_label(aes(label = Abbr),
#               size = 3,
#               label.padding = unit(0.02, "lines"),
#               label.size = 0.0) +
#    MySpecial +
#    labs(title = "Visualization of Rank Changes")
#  ########
