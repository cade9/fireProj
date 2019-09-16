# Purpose: Query the landsat USGS data for dates related to fire events 
# Select imagery 1 month before and 1 month after 
# Status: 
# Author: Christiana Ade
# Date: 8/8/2019
# Modified: 9/2/2019
# changed to filter by the correct year and seasons
# **Requires** 
# 1)  
####################################################################################
## require packages
require(raster)
require(rgdal)
require(stringr)
require(lubridate)
require(espa.tools)
require(tidyverse)
require(rgeos)
require(spdplyr)
require(spatstat)
require(getSpatialData)
# webshot::install_phantomjs()
library(htmlwidgets)
library(webshot)
require(naniar)
require(birk)
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
load("C:/Users/cade/Box/Ade_Xu/Data/Vector/shapefiles/USGS_2017_2018/2017shaplefile.rda")
# 4) 2018 shapefile
load("C:/Users/cade/Box/Ade_Xu/Data/Vector/shapefiles/USGS_2017_2018/2018shapefile.rda")
# 5) Grouped shapefiles determined in step1 
load("C:/Users/cade/Box/Ade_Xu/Data/Vector/shapefiles/USGS_2017_2018/groupedShp.rda")

#### OUTPUT FILES ####
#### USER DEFINED FUNCTIONS ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####             END USER DEFINED VARIABLES        ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####             START SCRIPT                      ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


orderList <- NULL
updateFire <- NULL

# log into earth explorer 
login_USGS("c.ade92") 

# need to select the ones that just have one entry 
groups17Fil <- groups17[lengths(groups17) == 1]

# 35 -- seems to have an issue

for (fl in 36:length(groups17Fil)){

  # identifing the list element that matches the name of that fire (fl)
  fir <- groups17Fil[fl]
  # name of the fire
  n17 <- names(fir)
  # identifing the names of the shapefiles associated with that fire
  shpFile <- fir[[n17]]
  # selecting the shapefiles
  f <- shape2017[shpFile]

  ## need to merge the shapefiles if length(f) is greater than 1

  # this needs to be t
  f2 <- f[[1]]

  # make sure the l8 grid matches that of the shapefile of the fire perimeter
  l8grid_trans <- spTransform(l8grid,crs(f2))
    # determine where they overlap
  # ~! there def has to be a way to do this in one step
  # oLap returns the info from the l8grid file which allows us to select the path and row
  # ~! what if it is more than one
  oLap <- over(f2,l8grid_trans)
  # need the shape file to determine the aoi later
  l8poly <- l8grid_trans %>% filter(PATH == oLap$PATH & ROW == oLap$ROW) # should actually be in this position
  # spTransform(f2, "+proj=longlat +datum=WGS84 +no_defs")

  # filter the dataset
  filDat <- fireDat %>% filter(Name == n17)

  # start and end of the fire
  fireStart<- dmy(filDat$`Start date`)
  fireContain <- dmy(filDat$`Containment date`)
  # anniversary start date of fire (fire date minus one year)
  aniDate <- fireStart - years(1)

  # determine the start date and end date and subtract and add 2 months
  start1 <- as.character( fireStart - months(13))
  end1 <- as.character( fireContain + months(2))

  #### Enter in information into the API ####
  # create the time range for the search
  time_range <-  c(start1,end1)

  # make sure the shapefile is in the projection? Does the function do that automatically
  set_aoi(l8poly)
  # product names for landsat available through the getSpatialData package
  product_names <- getLandsat_names() # you need to be logged on for this part to work
  # we are selecting product names 4 which is "LANDSAT_8_C1", but we might need to hard code this
  query <- getLandsat_query(time_range = time_range, name = "LANDSAT_8_C1", 
                            username = "c.ade92", password = "tuft9are25") %>%
    # filter by the correct path
    filter(WRSPath == as.vector(oLap$PATH) & WRSRow == as.vector(oLap$ROW))

  # filter for landcloud cover to be less than 50
  query2 <- query %>%
    filter(LandCloudCover <= 50) %>%
    # difference between start and containment date and the landsat 8 days
    #~! what if several of them are cloudy and we need more months??
    mutate(acquisitionDate = ymd(acquisitionDate),
          # startDif = as.numeric(acquisitionDate - fireStart, units ="days"),
           # containDif = as.numeric(acquisitionDate - fireContain, units = "days"),
           yr = year(acquisitionDate),
           mon = format(acquisitionDate,"%m-%d")) %>%
    replace_with_na_if(.predicate = is.numeric, condition = ~.x <0 )
  
  # determine which is the closest to the start date
  minStart <- which.closest(query2$acquisitionDate,aniDate)
  # which date is the closest to the containment date
  minEnd <- which.closest(query2$acquisitionDate,fireContain)
  # want to know what the actual date is though
  x <- query2$acquisitionDate[minEnd]
  
    if (x <= fireStart) {
      minEnd <- minEnd +1
    } 
  
  # filter the query
    m <- query2 %>% slice(c(minStart,minEnd))
    
    # rbind order to list to correct information
   orderList <- rbind(orderList,m)
   # add update the fireList
   filDat$preImg <- m$entityId[1]
   filDat$postImg <- m$entityId[2]
   updateFire <- rbind(updateFire,filDat)
}
  ## order the file list and make a file with the order number?  ## or you could just order the things individually

files <- getLandsat_data(orderList[c(61:65),], level = "sr", source = "auto",
                           dir_out = "C:\\Users\\cade\\Documents")  


write_csv(orderList,"./OrderList_Updated_Sep2.csv")
write_csv(updateFire,"./updated_fire_info_Sep2.csv")

## might want to include what happens if the closest image is the one of the fire
# if (m$startDif[1] == 0){
#   minStart2 <- minStart + 1
#   minEnd2 <- minEnd +1
#   m <- query2 %>% slice(c(minStart,minStart2,minEnd, minEnd2))
# } else {
#   print("poop")
# }


# espa-C.ade92@gmail.com-08192019-193031-603'
# 
# 
# 'espa-C.ade92@gmail.com-08152019-005724-662'
#   #   > saveWidget(z, "./temp.html")
#   #   > webshot("./temp.html","temp.png")


