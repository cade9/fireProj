# Purpose: Query the landsat USGS data for dates related to fire events 
# Select imagery 1 month before and 1 month after 
# Status: 
# Author: Christiana Ade
# Date: 8/8/2019
# Modified: 9/23/2019
# changed to filter by the correct year and seasons
# changed to filter by lat long points 
# **Requires** 
# 1)  
####################################################################################
## require packages
require(raster)
require(rgdal)
require(stringr)
require(lubridate)
require(tidyverse)
require(rgeos)
require(spdplyr)
require(spatstat)
require(getSpatialData)
require(spatialEco) # for point.in.poly tool
require(naniar)
require(birk)

#setwd("C:/Users/cade/Box/Ade_Xu")
setwd("Z:/Cade/CalFire")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####             START USER DEFINED VARIABLES       #### 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# The following variables should be defined by the user each time.
#### INPUT FILES ####
# 1) Spreadsheet detailing Fire events 
fireDat <- read_csv("./Data/CSV/2017-2018_california_wildfire_ca.csv") 
# 2) landsat 8 grid shapefile over california
l8grid <- readOGR("./Data/Vector/Grids/L8_WRS2_descending_intersectCali.shp")
# 3) 2017 shapefile 
#load("C:/Users/cade/Box/Ade_Xu/Data/Vector/shapefiles/USGS_2017_2018/2017shaplefile.rda")
# 4) 2018 shapefile
#load("C:/Users/cade/Box/Ade_Xu/Data/Vector/shapefiles/USGS_2017_2018/2018shapefile.rda")
# 5) Grouped shapefiles determined in step1 
#load("C:/Users/cade/Box/Ade_Xu/Data/Vector/shapefiles/USGS_2017_2018/groupedShp.rda")

#### OUTPUT FILES ####
outOrderList <- "./OrderList_Updated_Sep23.csv"
outUpdateFireInfo <- "./updated_fire_info_Sep23.csv"

#### USER DEFINED FUNCTIONS ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####             END USER DEFINED VARIABLES        ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####             START SCRIPT                      ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##################### PART 1: Create Point File from Fires #####################################
### filter fire csv for NA's
fireDat <- fireDat %>% 
  filter(!is.na(long)) %>%
  filter(!is.na(lat))

### Create spatial points data frame ####
firePoint <-  SpatialPointsDataFrame(fireDat[, c("long", "lat")], data = fireDat[1:ncol(fireDat)])


######################################
# QINGQING is this the correct datumn for the lat long points?
######################################

# change the coordinates to lat long 
proj4string(firePoint) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

##################### PART 2: Query USGS for the fire dates related to lat long #####################################
# create null order list
orderList <- NULL
updateFire <- NULL

# log into earth explorer 
login_USGS("c.ade92") 

#### Start Loop ####
for (i in 1:nrow(firePoint)){
  #### Step 1: Determine which path row the point for each the fire intersects ####
  
  # filter by fire entry
  fireE <- firePoint %>%
    slice(i) #i
  
  # find out which path row this overlaps with
  oLap <- point.in.poly(fireE,l8grid)
  
  ## create l8poly by merging any of the boxes together union
  l8Poly <- gUnaryUnion(l8grid %>%
    filter(PATH %in% oLap$PATH) %>%
    filter(ROW %in% oLap$ROW))
  
  #### Set the aoi for the landsat query ####
  # make sure the shapefile is in the projection? 
  set_aoi(l8Poly)
  
  #### Step 2: Determine which dates to query based on fire start and containment date ####
  # start and end of the fire
  fireStart<- dmy(first(oLap@data$Start.date))
  fireContain <- dmy(first(oLap$Containment.date))
  
  # anniversary start date of fire (fire date minus one year)
  aniDate <- fireStart - years(1)
  
  # determine the start date and end date for the query
  # subtract a year from the aniversary date
  start1 <- as.character(aniDate - year(1))
  
  ######################################
  # QINGQING if you want more than three months past the containment date go for it 
  # since it is 2017 and 2018 you could honestly just put this year?
  ######################################
  # add three months past the containment date 
  end1 <- as.character(fireContain + months(3))
  
  #### Enter in information into the API ####
  # create the time range for the search
  time_range <-  c(start1,end1)


  #### Step 3: Query the landsat data base ####
  
  # product names for landsat available through the getSpatialData package
  query <- getLandsat_query(time_range = time_range, name = "LANDSAT_8_C1", 
                            username = "c.ade92", password = "tuft9are25")
  
  
  #### Step 4: Define function for query database ####
  ######################################
  # QINGQING THIS IS WHERE YOU SHOULD INCLUDE ALL THE OTHER FILTER STATEMENTS 
  ######################################
  filtQuery <- function(df) {
    # determine which is the closest to the aniversary date
    minStart <- which.closest(df$acquisitionDate,aniDate)
    # which date is the closest to the containment date
    minEnd <- which.closest(df$acquisitionDate,fireContain)
    
    ## if it identifies the closest date to fire containment is before or equal to the day
    # the fire started change to the next available image
    x <- df$acquisitionDate[minEnd]
    if (x <= fireStart) {
      minEnd <- minEnd +1
    }
    ###################################################
    ##### SHOULD INCLUDE ALL OTHER FILTERS IN HERE #####
    ##################################################
    
    # filter the query
    queryFil <- df %>% slice(c(minStart,minEnd))
    return(queryFil)
  }
  
  #### Step 5: Filter the larger query list to the dates we want for pre and post fire ####

  query2 <- query %>%
    # filter for the path rows identified before
    filter(WRSPath %in% oLap$PATH) %>%
    filter(WRSRow %in% oLap$ROW) %>%
    # filter for a certain amount of landcover to be less than 50%
    filter(LandCloudCover <= 50) %>%
    # difference between start and containment date and the landsat 8 days
    mutate(acquisitionDate = ymd(acquisitionDate)) %>%
    # replace anything that is a numeric column with NA if there is a 0
    #replace_with_na_if(.predicate = is.numeric, condition = ~.x <0 ) %>%
    # group by path and row
    group_by(WRSPath,WRSRow) %>%
    # nest th data
    nest() %>%
    # apply the querying function to get the correct imagery for the different path rows
    # available related to the fire
    mutate(data = map(data,~filtQuery(.))) %>%
    # unest and ungroup data
    unnest(data) %>%
    ungroup(WRSPath,WRSRow) %>%
 
    #### Step 6: Update the orderList which will be entered into getLandsat_data() at the end of this loop #### 
    # rbind order to list to correct information
    orderList <- rbind(orderList,query2)
  
    #### Step 7: Update fire information #### 
    # add update the fireList
    fireInfo <- query2 %>%
    # add a fire name to the query
      mutate(FireName = fireE@data$Name) %>%
    # add the information contained within the shapefile to the orderlist
      left_join(fireE@data, by = c("FireName" = "Name"))
    updateFireInfo <- rbind(updateFireInfo,fireInfo)

}
##################### PART 3: write out csv files and ORDER DATA #####################################
##### Write out csv information ####
# order list
write_csv(orderList,outOrderList)
# updated fire information that includes the fire info in the orignial csv and the scenes that were ordered
# for that date 
write_csv(updateFire,outUpdateFireInfo)

###### ORDER DATA #######
files <- getLandsat_data(orderList[c(1:nrow(orderList)),], level = "sr", source = "auto",
                         dir_out = "C:\\Users\\cade\\Documents")  



