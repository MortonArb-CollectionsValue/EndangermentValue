

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
