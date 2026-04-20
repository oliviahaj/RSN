library(daymetr)
library(purrr)
library(dplyr)
library(tibble)
library(lubridate)
library(dplyr)


## downloaded data
site_lat_long <- read.csv('/Users/olhajek/Desktop/RSN/RSN_proj/2026 Data/RSN_Site/RSN_latlong.csv')


#--- download daymet data ---#
output <- c()
for(i in 1:nrow(site_lat_long)){
  temp_daymet <- download_daymet(
    lat = site_lat_long$lat[i],
    lon = site_lat_long$long[i],
    start = 1980,
    end = 2025
  ) 
  temp_daymet_data <- temp_daymet$data %>% add_column(SITE = site_lat_long$SITE[i],
                                                      .before = 1)
  output <- rbind(output, temp_daymet_data)
  rm(temp_daymet_data)
  cat(paste(site_lat_long$SITE[i],"check\n"))
}

## fix names. 
namefix <- c("daylight_secday"= "dayl..s."  ,
             "precip_mmday" =  "prcp..mm.day.",
             "shortwave_rad_Wm2" = "srad..W.m.2."  ,
             "snow_water_eq_kgm2" = "swe..kg.m.2."  ,
             "tmax_degC" = "tmax..deg.c." ,
             "tmin_degC" = "tmin..deg.c." ,
             "vp_Pa" = "vp..Pa."  )


output2 <- output %>% rename(all_of(namefix))
# add mean temp. 
output2$tmean_degC <- (output2$tmax_degC + output2$tmin_degC)/2

## add date and month.
startDate <- as.Date("1980-01-01")
endDate <- as.Date("2025-12-31")
endDate - startDate
datedf <- data.frame(date = startDate + 0:16070)
datedf$year <- format(datedf$date, format = "%Y") %>% as.integer()
datedf$month <- format(datedf$date, format = "%m") %>% as.integer()
datedf$yday <- lubridate::yday(datedf$date)

# Daily weather by site from 1980-2025
weather <- left_join(output2, datedf, by = c("year", "yday" ))

# Summarize to monthly, annual, and mean annual--------
# Replace BFY11 with BFY6 because it is close
bad_site  <- "BFY11"
good_site <- "BFY6"

replacement <- weather %>%
  filter(SITE == good_site) %>%
  mutate(SITE = bad_site)

weather.2 <- weather %>%
  filter(SITE != bad_site) %>%
  rbind(replacement)

monthly <- weather.2 %>% group_by(SITE, year, month) %>%
  summarize(mean_tmax_degC = mean(tmax_degC),
            mean_tmin_degC = mean(tmin_degC),
            mean_tmean_degC = mean(tmean_degC),
            precip_mmmonth = sum(precip_mmday), 
            mean_vp_pa = mean(vp_Pa)) 
annual <- monthly %>% group_by(SITE, year) %>%
  summarize(mean_tmax_degC = mean(mean_tmax_degC),
            mean_tmin_degC = mean(mean_tmin_degC),
            mean_tmean_degC = mean(mean_tmean_degC),
            precip_mmyear = sum(precip_mmmonth), 
            ann_vp = mean(mean_vp_pa)
  )
meanannual <- annual %>% group_by(SITE) %>%
  summarize(mean_tmax_degC = mean(mean_tmax_degC),
            mean_tmin_degC = mean(mean_tmin_degC),
            mean_tmean_degC = mean(mean_tmean_degC),
            mean_precip_mmyear = mean(precip_mmyear)
  )
write.csv(monthly, file = "/Users/olhajek/Desktop/RSN/RSN_proj/2026 Data/RSN_Site/daymet_monthly_weather.csv" , row.names = FALSE)
write.csv(annual, file = "/Users/olhajek/Desktop/RSN/RSN_proj/2026 Data/RSN_Site/daymet_annual_weather.csv" , row.names = FALSE)
write.csv(meanannual, file = "/Users/olhajek/Desktop/RSN/RSN_proj/2026 Data/RSN_Site/daymet_meanannual_weather.csv" , row.names = FALSE)

## BFY11 is weird! shoudl replace this with nearest other weather...likely BFY6 seems pretty close!! 

## Make a few quick visualizations
library(tidyverse)
meanannual <- left_join(meanannual, site_lat_long)
ggplot(meanannual, aes(mean_tmean_degC, mean_precip_mmyear, color = lat))+
  geom_point()

ggplot(annual, aes(year, mean_tmean_degC, group = SITE)) +
  geom_line()

ggplot(annual, aes(year, precip_mmyear, group = SITE)) +
  geom_line()


ggplot(annual, aes(year, mean_tmax_degC, group = SITE)) +
  geom_line()

ggplot(annual, aes(year, mean_tmin_degC, group = SITE)) +
  geom_line()

ggplot(annual, aes(year, ann_vp, group = SITE)) +
  geom_line()
