load("data/uk_o3.rda")
load("data/uk_temp.rda")
load("data/uk_wind.rda")
load("data/uk_pollutant_coords.rda")
load("data/uk_pollutant_date.rda")

library(mvtnorm)

uk_o3_uni <- data.frame(longitude = coords$longitude,
                        latitude = coords$latitude,
                        ozone = o3[, 1],
                        temp = temp[,1],
                        wind = wind[,1])

write.csv(uk_o3_uni, "data/uk_ozone_uni.csv")


