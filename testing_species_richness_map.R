### Author: Emily Beckman  ###  Date: 12/16/2020

### DESCRIPTION:
  #

### INPUTS:
	#

### OUTPUTS:
	#


################################################################################
# Load libraries
################################################################################

my.packages <- c("leaflet","raster","sp","rgeos","dplyr","rgdal",#"knitr",
	"RColorBrewer","tidyverse")

#install.packages (my.packages) # turn on to install current versions
lapply(my.packages, require, character.only=TRUE)
rm(my.packages)

################################################################################
# Working directory
################################################################################

local_dir <- "./Desktop/work"
data_dir <- "/Volumes/GoogleDrive/Shared drives/IMLS MFA/occurrence_points"
output_dir <- "/Volumes/GoogleDrive/Shared drives/IMLS MFA/Endangerment Value/script_outputs"

################################################################################
# Country-level richness maps
################################################################################

### READ IN COUNTRIES SHAPEFILE

# https://www.arcgis.com/home/item.html?id=2ca75003ef9d477fb22db19832c9554f
#world_country_shp <- readOGR("countries_shp/countries.shp")
# https://hub.arcgis.com/datasets/252471276c9941729543be8789e06e12_0?geometry=23.192%2C13.203%2C-13.370%2C79.425
world_country_shp <- readOGR(file.path(local_dir,
	"UIA_World_Countries_Boundaries-shp/World_Countries__Generalized_.shp"))

	# if needed, look at country names and codes
country_list <- paste(world_country_shp@data$COUNTRY,world_country_shp@data$ISO,sep="; ")
sort(country_list)

### READ IN SPECIES COUNTRY-LEVEL DISTRIBUTION SPREADSHEET

ctry_dist <- read.csv(file.path(data_dir,"inputs","known_distribution",
	"target_taxa_with_native_dist_12-16-21.csv"),
	as.is=T, na.strings=c("","NA"), colClasses="character")
head(ctry_dist)
ctry_dist$all_native_dist_iso2 <- gsub("; ",",",ctry_dist$all_native_dist_iso2)
ctry_dist <- ctry_dist %>% rename(countryCodes = all_native_dist_iso2)

### SPECIES RICHNESS CALCULATIONS

# function to create richness table and join to polygon data
richness.poly.countries <- function(df,polygons){
	# see max number of country codes for one species
	count_codes <- sapply(df$countryCodes,function(x) str_count(x, pattern = ","))
	max_codes <- str_count(df$countryCodes, pattern = ",")
	# create array of separated country codes
	ISO <- str_split_fixed(df$countryCodes, ",", n = max_codes+1)
	# sum to calculate richness
	richness <- as.data.frame(table(ISO))
	#richness <- richness[-1,]
	print(richness)
	# merge polygons with species richness data
	merged <- merge(polygons,richness)
	merged@data$Freq[which(is.na(merged@data$Freq))] <- 0
	merged <- merged[merged@data$Freq > 0,]
	return(merged)
}

# country-level richness for ALL species
map_countries <- richness.poly.countries(ctry_dist,world_country_shp)

# country-level richness for TOP 20 ENDANGERMENT DIMENSION species
top_endangerment <- c("Malus komarovii","Malus armeniacaefolia","Malus turkmenorum","Ulmus elongata",
"Quercus boyntonii","Malus honanensis","Ulmus pseudopropinqua","Ulmus microcarpa",
"Ulmus canescens","Ulmus prunifolia","Ulmus lamellosa","Quercus arkansana",
"Malus angustifolia","Malus praecox","Malus leiocalyca","Quercus austrina",
"Quercus sadleriana","Quercus acrodonta","Quercus cocciferoides",
"Ulmus elliptica")
ctry_dist_en <- ctry_dist %>%
	filter(ctry_dist$species_name_acc %in% top_endangerment)
map_countries_en <- richness.poly.countries(ctry_dist_en,world_country_shp)

# country-level richness for SEANS GENETIC DIMENSION TARGET species
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
ctry_dist_ge <- ctry_dist %>%
	filter(ctry_dist$species_name_acc %in% seans_spp)
map_countries_ge <- richness.poly.countries(ctry_dist_ge,world_country_shp)

### CREATE MAPS

# maping function
map.countries <- function(countries,pal,legend_text,legend_labels){
	map <- leaflet() %>%
		addProviderTiles("CartoDB.PositronNoLabels") %>%
		addPolygons(data = world_country_shp,
			color = "grey", weight = 0.6, opacity = 1,
			fillColor = "white",fillOpacity = 1) %>%
		addPolygons(data = countries,
			color = "grey", weight = 1, opacity = 1,
			fillColor = ~pal(countries@data$Freq),
			fillOpacity = 1) %>%
		addLegend(values = countries@data$Freq,
			pal = pal, opacity = 1,
			title = legend_text,
			labFormat = function(type, cuts, p) {paste0(legend_labels)},
			position = "bottomleft")
	return(map)
}

# ALL species
	# create color bins and labels
	hist(map_countries@data$Freq,breaks=90,xlim=c(0,80),ylim=c(0,15))
	bins <- c(0,1,3,5,10,15,20,25,60,Inf)
	labels <- c("0","1-2","3-4","5-9","10-14","15-19","20-39","25-30","60+")
	# create color palette
	#display.brewer.all()
	palette_country <- colorBin(palette = "PuRd", bins = bins,
		domain = map_countries@data$Freq, reverse = F, na.color = "white")
	# create map
	legend <- paste0("Number of native","<br/>","desiderata species")
	map_richness <- map.countries(map_countries,palette_country,legend,labels)
	map_richness
  # save map
  htmlwidgets::saveWidget(map_richness,
		file.path(output_dir,"CountryRichness_Desiderata_leafletmap.html"))

# TOP ENDANGERMENT species
	# create color bins and labels
	hist(map_countries_en@data$Freq,breaks=20,xlim=c(0,15),ylim=c(0,20))
	bins <- c(0,1,2,3,4,5,6,7,10,Inf)
	labels <- c("0","1","2","3","4","5","6","7","10+")
	# create color palette
	#display.brewer.all()
	palette_country <- colorBin(palette = "PuRd", bins = bins,
		domain = map_countries_en@data$Freq, reverse = F, na.color = "white")
	# create map
	legend <- paste0("Number of native,","<br/>","top endangerment species")
	map_richness_en <- map.countries(map_countries_en,palette_country,legend,labels)
	map_richness_en
  # save map
  htmlwidgets::saveWidget(map_richness_en, "CountryRichness_Top20Endangerment_leafletmap.html")

# SEANS GENETIC TARGET species
	# create color bins and labels
	hist(map_countries_ge@data$Freq,breaks=30,xlim=c(0,15),ylim=c(0,50))
	bins <- c(0,1,2,3,4,5,6,8,25,Inf)
	labels <- c("0","1","2","3","4","5","6-7","8","25+")
	# create color palette
	#display.brewer.all()
	palette_country <- colorBin(palette = "PuRd", bins = bins,
		domain = map_countries_ge@data$Freq, reverse = F, na.color = "white")
	# create map
	legend <- paste0("Number of native,","<br/>","target genetic species")
	map_richness_ge <- map.countries(map_countries_ge,palette_country,legend,labels)
	map_richness_ge
  # save map
  htmlwidgets::saveWidget(map_richness_en, "CountryRichness_TargetGenetic_leafletmap.html")


################################################################################
# Occurrence-level richness maps
################################################################################

















## NOT USING NOW -- FROM OTHER SCRIPT:



################################################################################
# MEXICO PROPAGATION MANUAL
################################################################################

# set working directory
setwd("./Desktop/work")

####
##### MEXICO BY STATE
####

### READ IN POLYGONS

# https://www.arcgis.com/home/item.html?id=ac9041c51b5c49c683fbfec61dc03ba8
mx_states_shp <- readOGR("mexstates/mexstates.shp")

### READ IN DISTRIBUTION AND GARDEN LOCATION SPREADSHEETS

# read in country/state distribution information
dist <- read.csv("Mexico2020_RevisadoSusanaValencia_andAllenCoombes.csv",
	as.is=T,na.strings=c("","NA"),colClasses="character",fileEncoding="UTF-8")
dist <- dist %>%
	rename(MX_states = Mexico.state.distribution,
				 RL_countries = Country.distribution,
			 	 RL_category = IUCN.Red.List.category)
	# remove accents
dist$MX_states <- stringi::stri_trans_general(dist$MX_states, "Latin-ASCII")
	# remove any double spaces
dist$MX_states <- str_squish(dist$MX_states)
# read in garden location points
mx_gardens <- read.csv("mexico_gardens.csv", as.is=T, na.strings=c("","NA"))

### CALCULATE STATE RICHNESS

# function to create richness table and join to polygon data
richness.poly.states <- function(df,polygons){
	l <- sapply(df$MX_states,function(x) str_count(x, pattern = ","))
	ADMIN_NAME <- str_split_fixed(df$MX_states, ", ", n = (max(l)+1))
	richness <- as.data.frame(table(ADMIN_NAME))
	richness <- richness[-1,]
  richness <- richness %>%
	  mutate(ADMIN_NAME = recode(ADMIN_NAME,
			"Mexico Distrito Federal" = "Distrito Federal",
			"Mexico State" = "Mexico"))
	print(richness)
	merged <- merge(polygons,richness)
	merged@data$Freq[which(is.na(merged@data$Freq))] <- 0
	return(merged)
}

# find state richness for ALL species
state_dist <- dist %>% filter(!is.na(MX_states))
map_states <- richness.poly.states(state_dist,mx_states_shp)
# find state richness for THREATENED species
state_dist_th <- dist %>% filter(!is.na(MX_states) &
	(RL_category=="CR" | RL_category=="EN" | RL_category=="VU"))
map_states_th <- richness.poly.states(state_dist_th,mx_states_shp)
# find state richness for ALL POSSIBLY THREATENED species(includes NT and DD)
#state_dist_pth <- dist %>% filter(!is.na(MX_states) &
#	(RL_category=="CR" | RL_category=="EN" | RL_category=="VU" |
#	 RL_category=="NT" | RL_category=="DD"))
#map_states_pth <- richness.poly.states(state_dist_pth,mx_states_shp)

# find number of gardens in each state
proj <- map_states@proj4string
latlong <- mx_gardens %>% select(long,lat)
spdf <- SpatialPointsDataFrame(latlong, mx_gardens, proj4string = proj)
proj_df <- spTransform(spdf,proj)
join <- raster::intersect(proj_df,mx_states_shp)
as.data.frame(join) %>% count(ADMIN_NAME)

### MAP

# map function
map.state.richness <- function(shapes,my_palette,labels,legend_txt1,
	legend_txt2){
	map <- leaflet() %>%
		addProviderTiles(
			"CartoDB.PositronNoLabels",
			options = providerTileOptions(maxZoom = 10)) %>%
		addPolygons(
			data = shapes,
			label = ~ADMIN_NAME,
			color = "grey",
			weight = 2,
			opacity = 0.5,
			fillColor = ~my_palette(shapes@data$Freq),
			fillOpacity = 0.8) %>%
		addCircleMarkers(
			data = mx_gardens,
			lat = ~lat, lng = ~long,
			label = ~garden,
			color = "black",
			fillColor = "#0fb909",
			radius = 3,
			weight = 2,
			stroke = T,
			fillOpacity = 1) %>%
		#addControl(
		#	title,
		#	position = "topright") %>%
  	addLegend(
    	labels = legend_txt2,
    	colors = "#0fb909; width: 10px; height: 10px; border-radius: 50%",
    	position = "bottomleft",
    	opacity = 1) %>%
		addLegend(
			pal = my_palette,
			values = shapes@data$Freq,
			opacity = 0.7,
			title = legend_txt1,
			labFormat = function(type, cuts, p) {paste0(labels)},
			position = "bottomleft") #%>%
		#leafem::addStaticLabels(
		#	.,
		#	data = shapes,
		#	label = shapes@data$Freq,
		#	style = list("font-weight"="bold","font-size"="14px"))
}

# ALL SPECIES
#legend_txt1 <- "Number of native oak species"
#legend_txt2 <- "Botanic gardens"
legend_txt1 <- "Número de especies <br/> de encinos nativos"
legend_txt2 <- "Jardines botanicos"
bins <- c(0,1,11,16,21,26,31,36,41,46,50,Inf)
label_state <- c("0","1-10","11-15","16-20","21-25","26-30","31-35","36-40",
	"41-45","46-49","50+")
palette_state <- colorBin(palette = "RdYlBu", bins = bins,
	domain = map_states@data$Freq, reverse = T, na.color = "transparent")
map <- map.state.richness(map_states,palette_state,label_state,legend_txt1,
	legend_txt2); map
#htmlwidgets::saveWidget(map, file = "Quercus_richness_mexico.html")

# THREATENED
#legend_txt1 <- "Number of native <br/> threatened oak <br/> species"
#legend_txt2 <- "Botanic gardens"
legend_txt1 <- "Número de especies <br/> de robles nativos <br/> amenazados"
legend_txt2 <- "Jardines Botanicos"
bins <- c(0,1,2,3,4,5,6,7,8)
label_state <- c("0","1","2","3","4","5","6","7")
palette_state <- colorBin(palette = "YlOrRd", bins = bins,
	domain = map_states_th@data$Freq, reverse = F, na.color = "transparent")
map <- map.state.richness(map_states_th,palette_state,label_state,legend_txt1,
	legend_txt2); map

# POSSIBLY THREATENED
#legend_txt1 <- "Number of native <br/> likely threatened <br/> oak species"
#bins <- c(0,1,3,6,9,12,15,19,21,Inf)
#label_state <- c("0","1-2","3-5","6-8","9-11","12-14","15-17","18-20","21+")
#palette_state <- colorBin(palette = "YlOrRd", bins = bins,
#	domain = map_states_th@data$Freq, reverse = F, na.color = "transparent")
#map <- map.state.richness(map_states_pth,palette_state,label_state,legend_txt1,
#	legend_txt2); map



####
##### MESOAMERICA BY COUNTRY
####

# set working directory
setwd("./Desktop/work")

### READ IN POLYGONS

	# https://www.arcgis.com/home/item.html?id=2ca75003ef9d477fb22db19832c9554f
#world_country_shp <- readOGR("countries_shp/countries.shp")
	# https://hub.arcgis.com/datasets/252471276c9941729543be8789e06e12_0?geometry=23.192%2C13.203%2C-13.370%2C79.425
world_country_shp <- readOGR("UIA_World_Countries_Boundaries-shp/World_Countries__Generalized_.shp")

# select only countries of interest
#select_countries <- c("Mexico","Guatemala","Belize","El Salvador","Honduras",
#	"Nicaragua","Costa Rica","Panama")
#target_countries <- world_country_shp[world_country_shp@data$COUNTRY %in%
#	select_countries,]
#head(target_countries)

### READ IN DISTRIBUTION DATA

# read in country/state distribution information
dist <- read.csv("Oaks of the World Master List - Oak spp. list w_ categories (for sharing).csv",
	as.is=T,na.strings=c("","NA"),colClasses="character",fileEncoding="UTF-8")
dist <- dist %>%
	rename(RL_countries = Country.distribution,
			 	 RL_category = IUCN.Red.List.category)
	# remove any extra spaces
dist$RL_countries <- str_squish(dist$RL_countries)
dist$RL_countries <- gsub(" ","",dist$RL_countries)

### CALCULATE STATE RICHNESS

# function to create richness table and join to polygon data
richness.poly.countries <- function(df,polygons){
	l <- sapply(df$RL_countries,function(x) str_count(x, pattern = ","))
	ISO <- str_split_fixed(df$RL_countries, ",", n = (max(l)+1))
	richness <- as.data.frame(table(ISO))
	richness <- richness[-1,]
	print(richness)
	merged <- merge(polygons,richness)
	merged@data$Freq[which(is.na(merged@data$Freq))] <- 0
	merged <- merged[merged@data$Freq > 0,]
	return(merged)
}
# find country richness
map_countries <- richness.poly.countries(dist,world_country_shp)
# find threatened species country richness
dist_th <- dist %>%
	filter(RL_category == "VU" | RL_category == "EN" | RL_category == "CR")
map_countries_th <- richness.poly.countries(dist_th,world_country_shp)

### MAP

map.countries <- function(countries,pal,legend_text,legend_labels){
	map <- leaflet() %>%
		addProviderTiles("CartoDB.PositronNoLabels") %>%
		addPolygons(data = world_country_shp,
			color = "grey", weight = 0.6, opacity = 1,
			fillColor = "white",fillOpacity = 1) %>%
		addPolygons(data = countries,
			color = "grey", weight = 1, opacity = 1,
			fillColor = ~pal(countries@data$Freq),
			fillOpacity = 1) %>%
		addLegend(values = countries@data$Freq,
			pal = pal, opacity = 1,
			title = legend_text,
			labFormat = function(type, cuts, p) {paste0(legend_labels)},
			position = "bottomleft") %>%
		leafem::addStaticLabels(
			.,
			data = countries,
			label = countries@data$Freq,
			style = list("font-weight"="bold","font-size"="14px"))
	return(map)
}

# ALL species
	# create color bins and labels
	#hist(map_countries@data$Freq,breaks=90,xlim=c(0,200),ylim=c(0,25))
	bins <- c(0,1,10,15,20,30,60,100,Inf)
	labels <- c("0","1-9","10-14","15-19","20-29","30-59","60-99","100+")
	# create color palette
	#display.brewer.all()
	palette_country <- colorBin(palette = "YlOrRd", bins = bins,
		domain = map_countries@data$Freq, reverse = F, na.color = "white")
	# create map
	legend <- paste0("Number of native","<br/>","oak species")
	map_richness <- map.countries(map_countries,palette_country,legend,labels)
	map_richness
  # save map
  #htmlwidgets::saveWidget(map_richness, "GlobalQuercusRichness_leaflet_map.html")

# THREATENED species
	# create color bins and labels
	#hist(map_countries_th@data$Freq,breaks=31,xlim=c(0,40),ylim=c(0,25))
	bins <- c(0,1,2,3,4,5,10,30,Inf)
	labels <- c("0","1","2","3","4","5-9","10-29","30+")
	# create color palette
	#display.brewer.all()
	palette_country <- colorBin(palette = "YlOrRd", bins = bins,
		domain = map_countries_th@data$Freq, reverse = F, na.color = "white")
	# create map
	legend <- paste0("Number of native,","<br/>","threatened oak species")
	map_richness_th <- map.countries(map_countries_th,palette_country,legend,labels)
	map_richness_th
  # save map
  #htmlwidgets::saveWidget(map_richness_th, "GlobalThreatenedQuercusRichness_leaflet_map.html")





################################################################################
# GAP ANALYSIS 2.0
################################################################################

### SET WORKING DIRECTORY

local_dir <- "./Desktop/work"

### READ IN BOUNDARIES SHAPEFILE

# define initial projection of points (usually WGS 84)
wgs.proj <- CRS("+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84
	+no_defs +towgs84=0,0,0")

# U.S. counties
	# https://www.census.gov/geographies/mapping-files/time-series/geo/carto-boundary-file.html
us_ctys <- readOGR(file.path(local_dir,"cb_2018_us_county_20m/cb_2018_us_county_20m.shp"))
#target_states <- c("AL","AZ","AR","CA","CO","CT","DE","FL","GA","ID",
#  "IL","IN","IA","KS","KY","LA","ME","MD","MA","MI","MN","MS","MO","MT","NE",
#  "NV","NH","NJ","NM","NY","NC","ND","OH","OK","OR","PA","RI","SC","SD","TN",
#  "TX","UT","VT","WV","WA","VA","WI","WY")
#us_states <- us_states[us_states@data$STUSPS %in% target_states,]
us_ctys.wgs <- spTransform(us_ctys,wgs.proj)
# create FIPS column to match to species distribution data
us_ctys@data$FIPS <- paste0("US",us_ctys@data$STATEFP,us_ctys@data$COUNTYFP)

# U.S. states
us_states <- readOGR(file.path(local_dir,"cb_2018_us_state_5m/cb_2018_us_state_5m.shp"))

# world countries

# https://hub.arcgis.com/datasets/252471276c9941729543be8789e06e12_0?geometry=23.192%2C13.203%2C-13.370%2C79.425
world_country_shp <- readOGR(file.path(local_dir,"UIA_World_Countries_Boundaries-shp/World_Countries__Generalized_.shp"))

### READ IN SPECIES DISTRIBUTION SPREADSHEET

dist <- read.csv(file.path(local_dir,"GA2_BONAP_USDA_countyDist.csv"),
	as.is=T, na.strings=c("","NA"), colClasses="character")
# create genus column and concatenate all county codes for each accepted species name
dist_collapse <- dist %>%
	separate("taxon_name_acc","genus",sep=" ",remove=F) %>%
	group_by(taxon_name_acc) %>%
  mutate(all_FIPS = paste(unique(FIPS), collapse=",")) %>%
	distinct(genus,taxon_name_acc,all_FIPS) %>%
	ungroup()
head(dist_collapse)

carya <- dist_collapse %>% filter(genus == "Carya")
fagus <- dist_collapse %>% filter(genus == "Fagus")
gymnocladus <- dist_collapse %>% filter(genus == "Gymnocladus")
juglans <- dist_collapse %>% filter(genus == "Juglans")
laurels <- dist_collapse %>% filter(genus == "Lindera" |
																		genus == "Persea" |
																		genus == "Sassafras")
pinus <- dist_collapse %>% filter(genus == "Pinus")
taxus <- dist_collapse %>% filter(genus == "Taxus")

### SPECIES RICHNESS CALCULATIONS

# function to create richness table and join to polygon data
richness.poly.counties <- function(df,polygons){
	# see max number of country codes for one species
	count_codes <- sapply(df$all_FIPS,function(x) str_count(x, pattern = ","))
	# create array of separated country codes
	FIPS <- str_split_fixed(df$all_FIPS, ",", n = (max(count_codes)+1))
	# sum to calculate richness
	richness <- as.data.frame(table(FIPS))
	#richness <- richness[-1,]
	print(richness)
	# merge polygons with species richness data
	merged <- merge(polygons,richness)
	merged@data$Freq[which(is.na(merged@data$Freq))] <- 0
	merged <- merged[merged@data$Freq > 0,]
	return(merged)
}

### CREATE MAPS

# maping function
map.counties <- function(counties,pal,legend_text,legend_labels){
	map <- leaflet() %>%
		addProviderTiles("CartoDB.PositronNoLabels") %>%
		addPolygons(data = counties,
			color = "grey", weight = 1, opacity = 0.9,
			fillColor = ~pal(counties@data$Freq),
			fillOpacity = 0.8) %>%
		#addPolygons(data = world_country_shp,
		#	color = "grey", weight = 2, opacity = 0.7,
		#	fillOpacity = 0) %>%
		addPolygons(data = us_states,
			color = "#636363", weight = 1.2, opacity = 1,
			fillOpacity = 0) %>%
		addLegend(values = counties@data$Freq,
			pal = pal, opacity = 0.8,
			title = legend_text,
			labFormat = function(type, cuts, p) {paste0(legend_labels)},
			position = "bottomright") %>%
		setView(-97, 38, zoom = 5)
	return(map)
}
### "COMMAND+" three times to make scale bar larger

legend <- "Number of species"

## Juglans
	# calculate counry-level richness
	map_counties <- richness.poly.counties(juglans,us_ctys)
	# create color bins and labels
	sort(unique(map_counties@data$Freq))
	bins <- c(1,2,3,Inf)
	labels <- c("1","2","3")
	# create color palette
	#display.brewer.all()
	palette_country <- colorBin(palette = "YlOrRd", bins = bins,
		domain = map_counties@data$Freq, reverse = F, na.color = "white")
	# create map
	map_richness <- map.counties(map_counties,palette_country,legend,labels)
	map_richness
## Carya
	map_counties <- richness.poly.counties(carya,us_ctys)
	sort(unique(map_counties@data$Freq))
	bins <- c(1,2,3,4,5,6,7,8,9,Inf)
	labels <- c("1","2","3","4","5","6","7","8","9")
	palette_country <- colorBin(palette = "YlOrRd", bins = bins,
		domain = map_counties@data$Freq, reverse = F, na.color = "white")
	map_richness <- map.counties(map_counties,palette_country,legend,labels)
	map_richness
## Laurels
	map_counties <- richness.poly.counties(laurels,us_ctys)
	sort(unique(map_counties@data$Freq))
	bins <- c(1,2,3,4,5,Inf)
	labels <- c("1","2","3","4","5")
	palette_country <- colorBin(palette = "YlOrRd", bins = bins,
		domain = map_counties@data$Freq, reverse = F, na.color = "white")
	map_richness <- map.counties(map_counties,palette_country,legend,labels)
	map_richness
## Pinus
	map_counties <- richness.poly.counties(pinus,us_ctys)
	#hist(map_counties@data$Freq,breaks=90,xlim=c(0,15),ylim=c(0,50))
	bins <- c(1,2,3,4,5,6,7,8,9,10,11,12,Inf)
	labels <- c("1","2","3","4","5","6","7","8","9","10","11","12")
	palette_country <- colorBin(palette = "YlOrRd", bins = bins,
		domain = map_counties@data$Freq, reverse = F, na.color = "white")
	map_richness <- map.counties(map_counties,palette_country,legend,labels)
	map_richness
## Taxus
	map_counties <- richness.poly.counties(taxus,us_ctys)
	sort(unique(map_counties@data$Freq))
	bins <- c(1,Inf)
	labels <- c("1")
	palette_country <- colorBin(palette = "YlOrRd", bins = bins,
		domain = map_counties@data$Freq, reverse = F, na.color = "white")
	map_richness <- map.counties(map_counties,palette_country,legend,labels)
	map_richness


## Fagus
	map_counties <- richness.poly.counties(fagus,us_ctys)
	sort(unique(map_counties@data$Freq))
	bins <- c(1,Inf)
	labels <- c("")
	palette_country <- colorBin(palette = "YlOrRd", bins = bins,
		domain = map_counties@data$Freq, reverse = F, na.color = "white")
	legend <- paste0("Native county of occurrence","<br/>"," for American beech")
	map_richness <- map.counties(map_counties,palette_country,legend,labels)
	map_richness
## Gymnocladus
	map_counties <- richness.poly.counties(gymnocladus,us_ctys)
	sort(unique(map_counties@data$Freq))
	bins <- c(1,Inf)
	labels <- c("")
	palette_country <- colorBin(palette = "YlOrRd", bins = bins,
		domain = map_counties@data$Freq, reverse = F, na.color = "white")
	legend <- paste0("Native county of occurrence","<br/>"," for Kentucky coffeetree")
	map_richness <- map.counties(map_counties,palette_country,legend,labels)
	map_richness

# version when map needs to be zoomed out
map.counties <- function(counties,pal,legend_text,legend_labels){
	map <- leaflet() %>%
		addProviderTiles("CartoDB.PositronNoLabels") %>%
		addPolygons(data = counties,
			color = "grey", weight = 0.5, opacity = 0.9,
			fillColor = ~pal(counties@data$Freq),
			fillOpacity = 0.8) %>%
		#addPolygons(data = world_country_shp,
		#	color = "grey", weight = 2, opacity = 0.7,
		#	fillOpacity = 0) %>%
		addPolygons(data = us_states,
			color = "#636363", weight = 1, opacity = 1,
			fillOpacity = 0) %>%
		addLegend(values = counties@data$Freq,
			pal = pal, opacity = 0.8,
			title = legend_text,
			labFormat = function(type, cuts, p) {paste0(legend_labels)},
			position = "bottomright") %>%
		setView(-97, 38, zoom = 5)
	return(map)
}
