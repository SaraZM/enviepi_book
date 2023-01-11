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
