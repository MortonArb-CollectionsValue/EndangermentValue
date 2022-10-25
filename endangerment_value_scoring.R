################################################################################

## endangerment_matrix_scoring.R

### Authors: Emily Beckman Bruns & Sean Hoban
### Funding: Institude of Museum and Library Services (IMLS MFA program grant
###          MA-30-18-0273-18 to The Morton Arboretum).
### Last Updated: 24 October 2022

### R version 4.1.3

### DESCRIPTION:
  # This script calculates the total endangerment value score for each target
  # species. It also calculates correlation among columns, performs a
  # sensitivity analysis, and calculates summary results. All analyses can be
  # performed for subsets of species individually (e.g., by genus) and/or for
  # all species together.

### FOLDER STRUCTURE:
  # endangerment_value_scoring
    # correlation_analysis
    # sensitivity_analysis
      # visualizations
      # summary_by_metric

### DATA INPUTS:
  # endangerment_matrix_forR.csv
  # str()...
    # species                 : chr  "Malus angustifolia" "Malus coronaria" "Malus fusca" "Malus ioensis" ...
    # extinction_risk         : chr  "LC" "LC" "LC" "LC" ...
    # in_natural_country      : chr  "Yes" "Yes" "Yes" "Yes" ...
    # exsitu_sites_plantsearch: int  24 59 59 51 31 186 16 40 17 176 ...
    # exsitu_sites_with_wild  : int  4 11 15 8 19 39 12 14 5 27 ...
    # num_wild_accessions     : int  8 19 22 11 86 321 18 55 9 78 ...
    # climate_change_vul_class: chr  "C" "A" "E2" "C" ...
    # pest_disease_vul_class  : chr  "B" "B" "E" "C" ...

### DATA OUTPUTS:
  # script_outputs folder...
  # For each group (all spp., each genus, or other grouping you've defined):
    # [GROUP]-EndangermentMatrix_SensitivityAnalysis.csv
    #   ^ Full sensitivity analysis with...
    #       1) raw scores for each column,
    #       2) total score calculated in each test
    #       3) species rank in each test
    #       4) final summary that counts the number of times each species
    #          is in the top 'bunch' (10th, 3rd, 4th) through all tests
    # [GROUP]-SensitivityAnalysis_rank_change_calcs.csv
    #   ^ Number of species moving more than X spots (5%, 10%, 15%, 20%)
    #     up or down in rank for each test
  # corr_charts folder...
  # For each group (all spp., each genus, or other grouping you've defined):
    # [GROUP]-All_cols_scored-correlation_matrix.png
    #   ^ Correlation matrix for all columns
    # [GROUP]-Exsitu_cols_scored-correlation_matrix.png
    #   ^ Correlation matrix for just the ex situ columns
  # sensitivity_charts folder...
  # For each group (all spp., each genus, or other grouping you've defined):
    # [GROUP]-Selected_cols_scored-sensitivity_analysis_line.png
    #   ^ One scatter plot that shows the 'base' total score for each species,
    #     color-coded by genus when multiple groups are included
    # [GROUP]-Selected_cols_scored-sensitivity_analysis_boxplots.png
    #   ^ Jitter box plot for each sensitivity test, showing spp & total scores
    # [GROUP]-Selected_cols_scored-sensitivity_analysis_scatterplots.png
    #   ^ Individual scatter plot for each sensitivity test, showing spp & total scores

################################################################################
# Load libraries
################################################################################

# this code chunk is from Shannon M Still !
  rm(list=ls())
  my.packages <- c('tidyverse','PerformanceAnalytics','ggplot2','ggrepel',
    'data.table','tidytext')
  # install.packages (my.packages) #Turn on to install current versions
  lapply(my.packages, require, character.only=TRUE)
    rm(my.packages)
  select <- dplyr::select

################################################################################
# Set working directory
################################################################################

main_dir <- "/Volumes/GoogleDrive-103729429307302508433/Shared drives/IMLS MFA/Endangerment Value/endangerment_value_scoring"

################################################################################
# Functions
################################################################################

# look at correlations among columns
col_corr <- function(df_matrix,cols_data){
	print("Correlation among <<measures>> themselves")
	print(sort(rowSums(cor(df_matrix[,cols_data],use="complete.obs")>.80),decreasing=T)[1:length(sel_col)])
	print(sort(rowMeans(cor(df_matrix[,cols_data],use="complete.obs")),decreasing=T)[1:length(sel_col)])
	print("Correlation among <<ranks>>- when put in order from measures")
	print(sort(rowSums(cor(apply(df_matrix[,cols_data],2,rank))>.80),decreasing=T)[1:length(sel_col)])
	print(sort(rowMeans(cor(apply(df_matrix[,cols_data],2,rank))),decreasing=T)[1:length(sel_col)])
	chart.Correlation(df_matrix[,cols_data], histogram = TRUE, method = "pearson")
}

# get the difference in new rank versus old rank
count_changes <- function(list_ranks,base=1){
	base_order<-list_ranks[[base]][order(list_ranks[[base]]$rank),1]
	examine_changes<-data.frame(base_order)
	for (i in 1:length(list_ranks)){
		ranks_diff <- match(base_order,list_ranks[[i]][order(list_ranks[[i]]$rank),1])
		examine_changes <- cbind(examine_changes,match(base_order,list_ranks[[i]][order(list_ranks[[i]]$rank),1])-1:length(base_order))
		colnames(examine_changes)[i+1]<-names(list_ranks[[i]][2])
	}
	examine_changes
}

# calculate total score for each species in the endangerment value matrix
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

# fomatting for the slope graph, via https://ibecav.github.io/slopegraph/
special_format <- list(
  # move the x axis labels up top
  scale_x_discrete(position = "top"),
  theme_bw(),
  # Format tweaks
  # Remove the legend
  theme(legend.position = "none"),
  # Remove the panel border
  theme(panel.border = element_blank()),
  # Remove just about everything from the y axis
  theme(axis.title.y = element_blank()),
  theme(axis.text.y = element_blank()),
  theme(panel.grid.major.y = element_blank()),
  theme(panel.grid.minor.y = element_blank()),
  # Remove a few things from the x axis and increase font size
  theme(axis.title.x = element_blank()),
  theme(panel.grid.major.x = element_blank()),
  theme(axis.text.x.top = element_text(size=12)),
  # Remove x & y tick marks
  theme(axis.ticks = element_blank()),
  # Format title & subtitle
  theme(plot.title = element_text(size=14, face = "bold", hjust = 0.5)),
  theme(plot.subtitle = element_text(hjust = 0.5))
)

# create slope graph that visualizes species rank switching between two tests
slope_graph_vis <- function(df,base_col,base_label,test_col,test_label){
  # left panel in visualization ('base' test)
  graph_base <- data.frame(
    Test = base_label,
    Species = en_df_test[,1],
    Rank = en_df_test[,base_col] )
  # right panel in visualization (comparison test)
  graph_compare <- data.frame(
    Test = test_label,
    Species = df[,1],
    Rank = df[,test_col] )
  # bind the two tests together
  view_ranks <- rbind(graph_base,graph_compare)
  # create species name abbreviation for viewing
  view_ranks$Abbr <- paste0(substr(view_ranks$Species, 0, 1),".",
    substr(sub(".* ","",view_ranks$Species), 0, 15))
  # create chart
  ggplot(data = view_ranks, aes(x = Test, y = Rank, group = Species)) +
    geom_line(aes(color = Species, alpha = 1), size = 2) +
    geom_label(aes(label = Abbr),
                   size = 4,
                   label.padding = unit(0.02, "lines"),
                   label.size = 0.0) +
    special_format +
    labs(title = "Visualization of Rank Changes")
}

################################################################################
# MANUAL CHANGES REQUIRED: Assign values for scoring and create species subsets
################################################################################


## read in endangerment value matrix formatted for R
en_df_raw <- read.csv(file.path(main_dir,"endangerment_matrix_forR.csv"),
  header = T, na.strings = c("","NA"))
en_df_raw <- en_df_raw %>% arrange(species)
en_df <- en_df_raw


## IF USING POTTER ET AL 2017/2019:
# add "C-" in front of climate change vals & "P-" in front of pest/disease vals
en_df$climate_change_vul_class[which(en_df$climate_change_vul_class != "No data")] <-
  paste0("C-",en_df$climate_change_vul_class[which(en_df$climate_change_vul_class != "No data")])
en_df$pest_disease_vul_class[which(en_df$pest_disease_vul_class != "No data")] <-
  paste0("P-",en_df$pest_disease_vul_class[which(en_df$pest_disease_vul_class != "No data")])


## IF DESIRED, CREATE SUBSETS OF TARGET SPECIES FOR SCORING SEPARATELY
# OPTION 1) group by genus
  # add column with genus name, for groupping ("Additional pieces" warning is ok)
  en_df <- en_df %>% separate("species","genus",sep=" ",remove=F)
# OPTION 2) manually select target species to group
  # e.g., list of target species we included in a separate genetic analysis
  special_group <- c("Ulmus rubra","Ulmus americana","Quercus laurifolia",
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
## create subsets
all <- en_df %>% mutate(group="AllSpp") %>% select(species,group,everything())
malus <- en_df %>% filter(grepl("Malus ",species)) %>% mutate(group="Malus") %>% select(species,group,everything())
quercus <- en_df %>% filter(grepl("Quercus ",species)) %>% mutate(group="Quercus") %>% select(species,group,everything())
tilia <- en_df %>% filter(grepl("Tilia ",species)) %>% mutate(group="Tilia") %>% select(species,group,everything())
ulmus <- en_df %>% filter(grepl("Ulmus ",species)) %>% mutate(group="Ulmus") %>% select(species,group,everything())
special <- en_df %>% filter(species %in% special_group) %>% mutate(group="GeneticDimSpp") %>% select(species,group,everything())
## create list of dataframes, one for each genus or other groupping, for
##  looping through to calculate scores, correlations, etc.
groups <- list(all,malus,quercus,tilia,ulmus,special)


## create dataframe for assigning scores to categorical columns
categories <- c(
    # IUCN Red List categories
    "EW","CR","EN","VU","NT","LC","DD","NE",
      # p is used to mark assessments that haven't been published yet;
      # * is used to mark species that have taxonomic notes
      # these are not scored differently
    "EWₚ","CRₚ","ENₚ","VUₚ","NTₚ","LCₚ","DDₚ","NE*",
    # nativity to country of interest
    "Yes","No",
    # climate change vulnerability classes (Potter et al. 2017)
    "C-A","C-B","C-C","C-D","C-E1","C-E2","C-E3","C-E4",
    # pest/disease vulnerability classes (Potter et al. 2019)
    "P-A1","P-A2","P-A3","P-B","P-A4","P-C","P-D","P-E",
    # when no data are available
    "No data")
## now provide a score for each category listed above
category_scores <- c(1, 0.8, 0.6, 0.4, 0.2, 0, 0.4, 0.2,
                     1, 0.8, 0.6, 0.4, 0.2, 0, 0.4, 0.2,
                     1, 0,
                     1, 0.5, 0.5, 0.5, 0.1, 0.1, 0.1, 0,
                     1, 0.5, 0.5, 0.5, 0.1, 0.1, 0.1, 0,
                     NA)
## create the df of categories and scores
vals <- data.frame(categories,category_scores)
vals
## look at df structure to help fill in next sections
str(all)
colnames(all)


## assign weight to each column (must add up to 1 if you want the final
##   scores to be scaled from 0 to 100)
## put zero if the column won't be scored (e.g., 'species' column)
col_wt <- c(0, 0, 0, 0.3, 0.05, 0.1, 0.05, 0.25, 0.125, 0.125)
if(sum(col_wt)!=1){ print("STOP !! THE COLUMN WEIGHTS SHOULD ADD UP TO ONE !!")}

## create vector of column weights when all are weighted evenly
col_using <- 7 # number of columns you're using (non-zero weights)
1/col_using # use this in the vector below, for all non-zero weighted columns
col_wt_even <- c(0, 0, 0, 0.1429, 0.1429, 0.1429, 0.1429, 0.1429, 0.1429, 0.1429)

## select columns that need reverse log transformation
# ex situ collections data columns
log_cols <- c(6:8)

## select columns you want to examine for correlations
# right now we are just using all columns
sel_col <- c(4:10)
# column nicknames for easier viewing of results when we drop each column
sel_col_nicknames <- c("no_extinction-risk","no_natural-country",
  "no_exsitu-sites","no_exsitu-with-wild","no_num-accessions",
  "no_climate-change","no_pest-disease")

## threashold for using Potter columns (minimum percent filled to use col)
# used for groups (e.g. genera) that have very little data; value is
# arbitrary and should be tested further
potter_thresh <- 0.15


## create dataframes to collect final sensitivity analysis rank stats
top_tenth <- data.frame(species="start",group="start",count_times_in_top10=0)
top_fourth <- data.frame(species="start",group="start",count_times_in_top25=0)
top_third <- data.frame(species="start",group="start",count_times_in_top33=0)


################################################################################
# START LOOPING THROUGH GROUPS
################################################################################

for(g in 1:length(groups)){


  ##############################################################################
  # Set up endangerment value matrix scoring
  ##############################################################################

  ## convert categories to scores using reference values created above
  for(i in 1:nrow(vals)){
    groups[[g]][groups[[g]] == vals[i,1]] <- vals[i,2]
  }

  ## make all scored columns numeric
  groups[[g]][,sel_col] <- groups[[g]][,sel_col] %>%
    mutate_if(is.character,as.numeric)
  str(groups[[g]])

  ## convert quantitative values to scores using equations
    # collections-based columns:
    #   Log-transformed then scaled in reverse from 1 (no collections) to
    #   0 (max num collections across all species)
  for(i in 1:length(log_cols)){
    ln_max <- max(log(groups[[g]][,log_cols[i]]+1))
    ln_min <- min(log(groups[[g]][,log_cols[i]]+1))
    for(j in 1:nrow(groups[[g]])){
      groups[[g]][j,log_cols[i]] <-
        1-((log(groups[[g]][j,log_cols[i]]+1)-ln_min)/(ln_max-ln_min))
    }
  }

  ## score Potter columns (climate change & pest/disease)
    # first check percent of species with data, and compare to threshold you set
    potter_thresh_g <- nrow(groups[[g]][which(!is.na(
      groups[[g]]$climate_change_vul_class)),]) / nrow(groups[[g]])
  if(potter_thresh < potter_thresh_g){
    # fill NA with mean
    cc_mean <- mean(groups[[g]]$climate_change_vul_class[!is.na(groups[[g]]$climate_change_vul_class)])
    pd_mean <- mean(groups[[g]]$pest_disease_vul_class[!is.na(groups[[g]]$pest_disease_vul_class)])
    groups[[g]]$climate_change_vul_class[is.na(groups[[g]]$climate_change_vul_class)] <- cc_mean
    groups[[g]]$pest_disease_vul_class[is.na(groups[[g]]$pest_disease_vul_class)] <- pd_mean
  } else { # remove Potter columns if % filled is below threshold
    groups[[g]]$climate_change_vul_class <- 0
    groups[[g]]$pest_disease_vul_class <- 0
  }
  str(groups[[g]])


  ##############################################################################
  # Correlation analysis
  #   calculate correlations among metrics and create correlation chart
  ##############################################################################

  # create folder for output charts
  if(!dir.exists(file.path(main_dir,"correlation_analysis")))
    dir.create(file.path(main_dir,"correlation_analysis"), recursive=T)

  # calculate correlations among all columns
  chart_path = file.path(main_dir,"correlation_analysis",
    paste0(gsub(" .*$","",groups[[g]]$group[1]),
    "_All-cols-scored_correlation-matrix.png"))
  png(height = 1500, width = 1500,
    file = chart_path, type = "cairo")
  col_corr(groups[[g]],sel_col)
  dev.off()


  ##############################################################################
  # Sensitivity analysis
  #   see how much the weighting and removal of columns affects the scores
  ##############################################################################

  # create folders for outputs
    # for all sensitivity analysisoutputs
  if(!dir.exists(file.path(main_dir,"sensitivity_analysis")))
    dir.create(file.path(main_dir,"sensitivity_analysis"), recursive=T)
    # subfolder for charts (visualizations)
  if(!dir.exists(file.path(main_dir,"sensitivity_analysis","visualizations")))
    dir.create(file.path(main_dir,"sensitivity_analysis","visualizations"), recursive=T)
    # subfolder for summary stats by metric (not by species)
  if(!dir.exists(file.path(main_dir,"sensitivity_analysis","summary_by_metric")))
    dir.create(file.path(main_dir,"sensitivity_analysis","summary_by_metric"), recursive=T)


  ### CALCULATE SCORES WITHIN DIFFERENT TESTS (e.g., removing metrics) ###


  ## IF USING POTTER ET AL 2017/2019:
  # create Potter cols for testing, where NA is replaced with 0 instead of mean
  groups[[g]]$climate_change_vul_class_zero <- as.numeric(gsub(
    cc_mean,0,groups[[g]]$climate_change_vul_class))
  groups[[g]]$pest_disease_vul_class_zero <- as.numeric(gsub(
    pd_mean,0,groups[[g]]$pest_disease_vul_class))


  ## 'base' test: relative weighting for metrics & Potter means for NAs
  wt_all <- calc_total_score(sel_col,groups[[g]],col_wt)
  groups[[g]] <- cbind(groups[[g]],wt_all)
  # recalculate total score for different tests, using weighted score
  for(i in 1:length(sel_col)){
    test_col <- sel_col[-i]
    test <- calc_total_score(test_col,groups[[g]],col_wt)
    groups[[g]] <- cbind(groups[[g]],test)
    colnames(groups[[g]])[colnames(groups[[g]])=="test"] <-
      paste0("wt_",sel_col_nicknames[i])
  }
  ## IF USING POTTER ET AL 2017/2019:
  # test using 0 instead of NA when no data
   # select columns
  sel_col_PotterZero <- c(4:8,11:12)
   # selected column weights
  col_wt_PotterZero <- c(0, 0, 0, 0.3, 0.05, 0.1, 0.05, 0.25, 0, 0, 0.125, 0.125)
  groups[[g]]$wt_zero_potter <-
    calc_total_score(sel_col_PotterZero,groups[[g]],col_wt_PotterZero)
  head(groups[[g]])


  ## now test using the same weight for every column
  ev_all <- calc_total_score(sel_col,groups[[g]],col_wt_even)
  groups[[g]] <- cbind(groups[[g]],ev_all)
  # recalculate total score for different tests, using even weights
  for(i in 1:length(sel_col)){
    test_col <- sel_col[-i]
    test <- calc_total_score(test_col,groups[[g]],col_wt_even)
    groups[[g]] <- cbind(groups[[g]],test)
    colnames(groups[[g]])[colnames(groups[[g]])=="test"] <-
      paste0("ev_",sel_col_nicknames[i])
  }
  ## IF USING POTTER ET AL 2017/2019:
  # test using 0 instead of NA when no data
   # selected column weights
  col_ev_PotterZero <- c(0, 0, 0, 0.1429, 0.1429, 0.1429, 0.1429, 0.1429, 0, 0, 0.1429, 0.1429)
  groups[[g]]$wt_zero_potter <-
    calc_total_score(sel_col_PotterZero,groups[[g]],col_ev_PotterZero)
  head(groups[[g]])


  ### VISUALIZE RESULTS FROM TESTS ###

  ###
  ### THE REST OF THIS SCRIPT NEEDS SOME MANUAL EDITS TO RUN
  ### places that need manual edits are marked with '!!'
  ###

  ## format dataframe for visualizing
  vis <- data.frame(species="start",test="start",values=0)
  for(i in 1:(length(groups[[g]])-12)){ #!!subtract num cols before total scores
    add_vis <- data.frame(
                species=groups[[g]][,1],
                test=rep(colnames(groups[[g]][i+12]),n=nrow(groups[[g]])),#!!
                values=groups[[g]][,i+12] #!!num columns before scores
              )
    vis <- rbind(vis,add_vis)
  }; vis <- vis[-1,]

  ## multi-facet scatter plot
  # not optimized so you can read all species names
  vis %>%
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
  ggsave(file.path(main_dir,"sensitivity_analysis","visualizations",
    paste0(gsub(" .*$","",groups[[g]]$group[1]),
    "_Sensitivity-analysis_scatterplots.png")),width=20,height=10)

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
  ggsave(file.path(main_dir,"sensitivity_analysis","visualizations",
    paste0(gsub(" .*$","",groups[[g]]$group[1]),
    "_Sensitivity-analysis_boxplots.png")),width=30,height=10)


  ##############################################################################
  # Synthesis: summary of sensitivity analysis
  #   see which species are ranked highly through all tests
  ##############################################################################


  ### ADD RANKS FOR DIFFERENT TESTS ###


  ## sort by 'base' test (high to low score) before adding ranks
  groups[[g]] <- groups[[g]] %>% arrange(desc(wt_all))

  ## add columns where total scores are converted to ranks
  col_to_rank <- groups[[g]][,(ncol(groups[[g]])-(ncol(groups[[g]])-ncol(en_df_raw))+5):ncol(groups[[g]])]
  for(i in 1:length(col_to_rank)){
    col_ranked <- frankv(col_to_rank[i],ties.method="min",order=-1)
    groups[[g]] <- cbind(groups[[g]],col_ranked)
    colnames(groups[[g]])[colnames(groups[[g]])=="col_ranked"] <-
      paste0(names(col_to_rank[i]),"_r")
  }
  str(groups[[g]])


  ### QUANTIFY SENSITIVITY ###


  ## first let's look at two different methods and compare them
  # select columns we'll use for the analyses
  colnames(groups[[g]])
  all_ranks <- groups[[g]][,30:46] #!!select rank columns (end with "_r")
  colnames(all_ranks) <- colnames(groups[[g]][,13:29]) #!!select total score cols
  species <- groups[[g]][,1]
  # one way to find highest priority species is to identify those that are most
  #   frequently ranked in a given bunch, say in the top 10
  species_ranked1 <- data.frame(species,rowSums(all_ranks<10))
  colnames(species_ranked1) <- c("sp","freq_in_top_10")
  # another way is to take the mean across rows in the rank order
  species_ranked2 <- data.frame(species,rowMeans(all_ranks))
  colnames(species_ranked2)<-c("sp","mean_of_ranks")
  # it really doesn't matter which approach :)
  cbind(species_ranked1[order(species_ranked1$freq_in_top_10,decreasing=T),],
        species_ranked2[order(species_ranked2$mean_of_ranks),])


  ## calculate changes during senstivity analysis, based on different thresholds
  # list of dataframes, one for each test; input for analysis
  list_for_analysis <- list()
  for(i in 1:(length(all_ranks))){
    t <- data.frame("sp"=species,"temp"=as.numeric(all_ranks[,i]))
    colnames(t)[colnames(t)=="temp"] <- paste0("rank_",names(all_ranks[i]))
    list_for_analysis[[i]] <- t
  }
  # use the function defned above to count the rank changes,
  #   i.e., compare the species 'base' test rank to its rank in each
  #        additional test
  examine_changes <- count_changes(list_for_analysis)
  # use number of target species to calculate thresholds (percents) for ranks
  #   shifting up and down
  num_spp <- nrow(groups[[g]])
  cutoff1 <- round((num_spp * 0.05),0) #!! can change percent
  cutoff2 <- round((num_spp * 0.10),0) #!! can change percent
  cutoff3 <- round((num_spp * 0.15),0) #!! can change percent
  cutoff4 <- round((num_spp * 0.20),0) #!! can change percent
  # calculate number of species that have a rank change greater
  #   than the threshold, within each test in the sensitivity analysis
  rank_changes <- cbind(
    greater.than.5 = colSums(abs(examine_changes[,-1])>cutoff1),
    greater.than.10 = colSums(abs(examine_changes[,-1])>cutoff2),
    greater.than.15 = colSums(abs(examine_changes[,-1])>cutoff3),
    greater.than.20 = colSums(abs(examine_changes[,-1])>cutoff4)
  ); rank_changes
  # write file with restuls
  write.csv(rank_changes, file.path(main_dir,"sensitivity_analysis",
    "summary_by_metric", paste0(gsub(" .*$","",groups[[g]]$group[1]),
    "_Sensitivity-analysis_num-species-rank-change-thresholds.csv")), row.names = T)


  ## calculate the number of tests where each species has a rank change greater
  ##   than the threshold
  ncol(all_ranks)
  # top 10% (tenth) of spots
  species_ranked1 <- data.frame(
    species, group=rep(groups[[g]]$group[1]),
    count_times_in_top10=rowSums(all_ranks<(round(nrow(groups[[g]])/10,0))))
  species_ranked1 <- species_ranked1 %>% arrange(desc(count_times_in_top10))
  top_tenth <- rbind(top_tenth,species_ranked1)
  head(species_ranked1[order(species_ranked1[,2],decreasing=T),],n=20)
  # top 25% (fourth) of spots
  species_ranked2 <- data.frame(
    species, group=rep(groups[[g]]$group[1]),
    count_times_in_top25=rowSums(all_ranks<(round(nrow(groups[[g]])/4,0))))
  species_ranked2 <- species_ranked2 %>% arrange(desc(count_times_in_top25))
  top_fourth <- rbind(top_fourth,species_ranked2)
  head(species_ranked2[order(species_ranked2[,2],decreasing=T),],n=30)
  # top 33% (third) of spots
  species_ranked3 <- data.frame(
    species, group=rep(groups[[g]]$group[1]),
    count_times_in_top33=rowSums(all_ranks<(round(nrow(groups[[g]])/3,0))))
  species_ranked3 <- species_ranked3 %>% arrange(desc(count_times_in_top33))
  top_third <- rbind(top_third,species_ranked3)
  head(species_ranked3[order(species_ranked3[,2],decreasing=T),],n=50)
  # join all together for saving
  groups[[g]] <- Reduce(full_join,list(
    groups[[g]],species_ranked1,species_ranked2,species_ranked3))
  head(groups[[g]])
  # write file
  write.csv(groups[[g]], file.path(main_dir,"sensitivity_analysis",
    paste0(gsub(" .*$","",groups[[g]]$group[1]),
    "_Endangerment-matrix_sensitivity-analysis-and-synthesis.csv")),
    row.names = F)


  ### VISUALIZE SOME RANK CHANGE RESULTS ###


  ## slope graph visualization
  colnames(groups[[g]])
  en_df_test <- groups[[g]] %>% arrange(wt_all_r)
  # create charts using function defined at the beginning of the script;
  #   inputs are: dataframe, column num for 'base' test, label for 'base' test,
  #               column num for comparison test, label for comparison test
  # select whichever tests you'd like to compare visually
  slope_graph_vis(en_df_test,30,"*All_Columns",39,"All_Same_Weight")
  slope_graph_vis(en_df_test,30,"*All_Columns",31,"No_Extinction-risk")
  slope_graph_vis(en_df_test,30,"*All_Columns",32,"No_Natural-dist-in-US")
  # I have not bothered to save these charts because the labels don't show
  #   up in the saved version :(

}
