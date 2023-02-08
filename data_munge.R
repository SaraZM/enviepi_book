load("data/uk_o3.rda")
load("data/uk_temp.rda")
load("data/uk_wind.rda")
load("data/uk_pollutant_coords.rda")
load("data/uk_pollutant_date.rda")

library(mvtnorm)

# uk_o3_uni <- data.frame(longitude = coords$longitude,
#                         latitude = coords$latitude,
#                         ozone = o3[, 1],
#                         temp = temp[,1],
#                         wind = wind[,1])
# 
# write.csv(uk_o3_uni, "data/uk_ozone_uni.csv")

uk_o3_one_site <- data.frame(ozone = o3[1,],
                        temp = temp[1, ],
                        wind = wind[1, ],
                        date = date)

write.csv(uk_o3_one_site, "data/uk_ozone_one_site.csv", row.names = FALSE)

# Chapter 10 ----
library(dplyr)
library(readxl)
library(spdep)
library(data.table)

## December benzene -----
VOCs <- read_excel("data/VOCDec2005.xlsx") 
VOCS_IDs <- read_excel("data/VOClogsheetDec2005.xlsx") 

VOCs_coords <- VOCS_IDs[,c("ID", "X", "Y")]
benzene <- VOCs[,c("ID", "Benzene")]
# avoid dealing with replicates
benzene <- benzene[unique(benzene$ID),]

benzene_coords <- merge(benzene, VOCs_coords, by = "ID") %>% 
  distinct(X, Y, .keep_all = TRUE)

# Note: the coordinates for site 26 are wrong, they were fixed for the article 
# but I just removed them
benzene_coords <- benzene_coords[-21,]
# Also, in order to make the code cleaner in the chapter, I will transform the 
# coordinates to lon lat (for the map)

benzene_utm <-
  SpatialPointsDataFrame(
    coords = benzene_coords[, c("X", "Y")],
    data = benzene_coords[, c("ID", "Benzene")],
    proj4string = CRS("+proj=utm +zone=18 +ellps=WGS72")
  ) 

benzene_geo <- spTransform(benzene_utm, CRS("+proj=longlat +datum=WGS84"))
benzene_geo_df <- as.data.frame(benzene_geo)

colnames(benzene_geo_df) <- c("ID", "Benzene", "lon", "lat")

write.csv(benzene_geo_df, "data/montreal_benzene.csv", row.names = FALSE)


## California temperature ----
cal_sites <- read.csv("data/metadataCA.txt", sep = "\t")
cal_temp <- read.csv("data/MaxCaliforniaTemp.csv")

cal_temp_oneday <- cal_temp[cal_temp$X == "20120401",]
cal_temp_oneday <- as.data.frame(t(cal_temp_oneday[c(2:14,16:19)]))
cal_temp_oneday <- tibble::rownames_to_column(cal_temp_oneday, "Location")

colnames(cal_temp_oneday) <- c("Location", 'MaxTemp')

cal_sites$MaxTemp <- cal_temp_oneday[,2]

write.csv(cal_sites, "data/CalTempData.csv", row.names = FALSE)

## April benzene
VOCs <- read_excel("data/VOCApr2006.xlsx") 
VOCS_IDs <- read_excel("data/VOClogsheetApr2006.xlsx") 

VOCs_coords <- VOCS_IDs[,c("ID", "X", "Y")]
benzene <- VOCs[,c("ID", "Benzene")]
# avoid dealing with replicates
benzene <- benzene[unique(benzene$ID),]

benzene_coords <- merge(benzene, VOCs_coords, by = "ID") %>% 
  distinct(X, Y, .keep_all = TRUE)

# Note: the coordinates for site 26 are wrong, they were fixed for the article 
# but I just removed them
benzene_coords <- benzene_coords[-c(4,5,23,92,104),]
# Also, in order to make the code cleaner in the chapter, I will transform the 
# coordinates to lon lat (for the map)

benzene_utm <-
  SpatialPointsDataFrame(
    coords = benzene_coords[, c("X", "Y")],
    data = benzene_coords[, c("ID", "Benzene")],
    proj4string = CRS("+proj=utm +zone=18 +ellps=WGS72")
  ) 

benzene_geo <- spTransform(benzene_utm, CRS("+proj=longlat +datum=WGS84"))
benzene_geo_df <- as.data.frame(benzene_geo)

colnames(benzene_geo_df) <- c("ID", "Benzene", "lon", "lat")

write.csv(benzene_geo_df[,-c( 1)], 
          "data/montreal_benzene_apr.csv", row.names = FALSE)

# Chapter 11 -----

## Daily data ----

ozone_day <- read.csv("data/O3_2013_LA_sites.csv")

site_day <- ozone_day[ozone_day$AQS_SITE_ID == "06-037-9033",] %>% 
  select(Date, Daily.Max.8.hour.Ozone.Concentration)


colnames(site_day) <- c("date", "max.ozone")

#write.csv(site_day, "data/LA_daily.csv", row.names = FALSE)

## Hourly data -----

ozone_hour <- read.csv("data/O3_2013_Hourly_LA_Site840060370113.csv")

site_hour <- select(ozone_hour, datetime, value)

colnames(site_hour) <- c("datetime", "ozone")

site_hour$time <- 1:nrow(site_hour)
# site_hour <- tidyr::separate(site_hour, datetime, c("date", "time"), sep = "T")


write.csv(site_hour, "data/LA_ozone_hourly.csv", row.names = FALSE)
