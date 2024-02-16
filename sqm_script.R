#####Sensor Data Processor#######

#Create Database for sqm data####
#Run sensor data processor python script in folder
# Right click on folder SensorDataProcessor (folder) open new terminal > pipenv shell > python main.py . > updates existing database#

#setworking directory 
  setwd("~/Documents/OU/workspace/SensorDataProcessor")

#Pull data from SQL database####
  library(DBI)
  mydb <- dbConnect(RSQLite::SQLite(), 
                  "/Users/Documents/OU/workspace/SensorDataProcessor/Data.db")
  dbListTables(mydb)

#Select hours 23:00-04:00 
  #Split data frames into evening (2300-0000) and morning hours (0000-0400)
    # This is necessary if you want to analyze by night as opposed to calendar day. 

    Midnightsp22 <- dbGetQuery(mydb, "select *
      from sensor_data 
      where (local_time between '2022-04-22' and '2022-05-13')
      and (strftime('%H:%M', local_time) between '23:00' and '24:00')")
    
    Morningsp22 <- dbGetQuery(mydb, "select *
      from sensor_data 
      where (local_time between '2022-04-23' and '2022-05-14')
      and (strftime('%H:%M', local_time) between '00:00:00' and '04:00')")
    
    MidnightLights <- dbGetQuery(mydb, "select * 
      from sensor_data 
      where (local_time NOT between '2022-04-22' and '2022-05-13')
      and (strftime('%H:%M', local_time) between '23:00' and '24:00')")
    
    MorningLights <- dbGetQuery(mydb, "select * 
      from sensor_data 
      where (local_time NOT between '2022-04-23' and '2022-05-14')
      and (strftime('%H:%M', local_time) between '00:00:00' and '04:00')")

dbDisconnect(mydb) 

#Create dataframe with all obs####
  #Assign light treatment for Midnight hours
MidnightLights$Treatment <- "Lights On"

MidnightLights$Treatment[which(MidnightLights$Treatment<="Lights On")] <- "Lights On"
Midnightsp22$Treatment <- "Lights Out"
Midnightsp22$Treatment[which(Midnightsp22$Treatment<="Lights Out")] <- "Lights Out"


df_evening <- (rbind(MidnightLights,Midnightsp22))

#light treatment for Morning hours
MorningLights$Treatment <- "Lights On"
MorningLights$Treatment[which(MorningLights$Treatment<="Lights On")] <- "Lights On"

Morningsp22$Treatment <- "Lights Out"
Morningsp22$Treatment[which(Morningsp22$Treatment<="Lights Out")] <- "Lights Out"

df_morning <- (rbind(MorningLights,Morningsp22))


#seperate local date for each df
library(fasttime)
library(data.table)
is.data.frame(df_evening) == TRUE

# Convert data frames to data tables and create a new Date column
dt_evening <- as.data.table(df_evening, keep.rownames=FALSE, key=NULL, check.names=FALSE)
dt_evening[, Date := as.Date(substr(local_time, 0, 10))]

dt_morning <- as.data.table(df_morning, keep.rownames=FALSE, key=NULL, check.names=FALSE)
dt_morning[, Date := as.Date(substr(local_time, 0, 10))]

# create list of all the location in our data
locations <- unique(dt_morning$location)

#ONLY run this ONCE for morning datatable
  #Assigns 0000-0400 the previous calendar data to keep continuity of nights
dt_morning[, Date := Date - 1]

#combine the two
dt_together <- rbind(dt_morning,dt_evening)

#convert table with new column back to a to data frame
df_all <- as.data.frame(dt_together)
write.csv(df_all, "df_all.csv", row.names = F)

#Run df_all in python for cloud script then upload
# insure tcdc.2022.nc and tcdc.2023.nc are downloaded from NARR website and in SenosrDataProcessor file
# Right click on folder SensorDataProcessor open new terminal >
# Check working directory if needed with cwd
# Pip install netCDF4 >
# Pip install pandas >
# python3 extract_NARR_cloud-NEW.py

#Run steps at begining of file to get df_all_cloud
df_all_cloud <- read.csv("df_all_cloud.csv")

df_all_cloud$utc_time <- as.POSIXct(df_all_cloud$utc_time)


#Append weather -- closest index #####   
setwd("~/Documents/OU/workspace/SensorDataProcessor")
dal_weather_full <- read.csv("dal_weather_full.csv")
hou_weather_full <- read.csv("hou_weather_full.csv")
dal_precip <- read.csv("dal_precip.csv")
hou_precip <- read.csv("hou_precip.csv")

dal_weather <- dal_weather_full[dal_weather_full$Latitude == "32.85416",]

hou_weather <- hou_weather_full[hou_weather_full$Latitude == "29.9844",]

colnames(dal_weather)[6] = "utc_time"
colnames(hou_weather)[6] = "utc_time"

dal_weather$utc_time <- as.POSIXct(dal_weather$utc_time, format = "%Y-%m-%dT%H:%M:%S")
hou_weather$utc_time <- as.POSIXct(hou_weather$utc_time, format = "%Y-%m-%dT%H:%M:%OS")


#Rename weather variables
colnames(dal_weather)[9] = "pressure_Pa"
colnames(dal_weather)[10] = "rel_humidity"
colnames(dal_weather)[11] = "temp_F"
colnames(dal_weather)[12] = "wind_direction_deg"
colnames(dal_weather)[14] = "wind_speed_m_s"

colnames(hou_weather)[9] = "pressure_Pa"
colnames(hou_weather)[10] = "rel_humidity"
colnames(hou_weather)[11] = "temp_F"
colnames(hou_weather)[12] = "wind_direction_deg"
colnames(hou_weather)[14] = "wind_speed_m_s"

#Closest Index weather
#DALLAS
# Find closest timestamps from 'dallas_cloud' in 'dal_weather'
closest_indices <- sapply(dallas_cloud$utc_time, function(x) {
  closest_index <- which.min(abs(x - dal_weather$utc_time))
  return(closest_index)
})

# Merge
dal_result <- cbind(dallas_cloud, dal_weather[closest_indices, ])

#HOUSTON
# Find closest timestamps from 'dallas_cloud' in 'dal_weather'
closest_indices <- sapply(houston_cloud$utc_time, function(x) {
  closest_index <- which.min(abs(x - hou_weather$utc_time))
  return(closest_index)
})

# Merge
hou_result <- cbind(houston_cloud, hou_weather[closest_indices, ])

#write csvs
write.csv(dal_result, "dal_result.csv", row.names = F)
write.csv(hou_result, "hou_result.csv", row.names = F)
#Clean up and merge result
setwd("~/Documents/OU/workspace/SensorDataProcessor")
dal_result <- read.csv("dal_result.csv")
hou_reslut <- read.csv("hou_result.csv")

columns_to_remove <- c("X", "Solar.Radiation..W.m2.", "Soil.Volumetric.Water.Content.5.cm....", "X.1")
dal_result_clean <- dal_result[, -which(names(dal_result) %in% columns_to_remove)]

columns_to_remove <- c("X", "Precipitation.5min..in.", "X.1")
hou_result_clean <- hou_result[,-which(names(hou_result) %in% columns_to_remove)]

  #Add precipitation#####
#Dallas
#   setwd("~/Documents/OU/workspace/SensorDataProcessor/dallas_weather")
#   dal_22_precip <- read.csv("dal_apr22_apr23_precip.csv")
#   dal_23_precip <- read.csv("dal_apr23_jan23_precip.csv")
#   dal_precip <- rbind(dal_22_precip, dal_23_precip)
# #Houston
#   setwd("~/Documents/OU/workspace/SensorDataProcessor/houston_weather")
#   hou_22_precip <- read.csv("hou_apr22_apr23_precip.csv")
#   hou_23_precip <- read.csv("hou_apr23_jan23_precip.csv")
#   hou_precip <- rbind(hou_22_precip, hou_23_precip)  
#   
#   
#   setwd("~/Documents/OU/workspace/SensorDataProcessor")
#   write.csv(hou_precip, "hou_precip.csv", row.names=F)
#   write.csv(dal_precip, "dal_precip.csv", row.names = F)

setwd("~/Documents/OU/workspace/SensorDataProcessor")
dal_precip <- read.csv("dal_precip.csv")
hou_precip <- read.csv("hou_precip.csv")

#Rename columns
names(dal_precip)[7]<-"utc_time"
dal_precip$utc_time <- as.POSIXct(dal_precip$utc_time, format = "%Y-%m-%dT%H:%M:%OS")
dal_result_clean$utc_time <- as.POSIXct(dal_result_clean$utc_time, format = "%Y-%m-%d %H:%M:%S")

names(hou_precip)[7]<-"utc_time"
hou_precip$utc_time <- as.POSIXct(hou_precip$utc_time, format = "%Y-%m-%dT%H:%M:%OS")
hou_result_clean$utc_time <- as.POSIXct(hou_result_clean$utc_time, format = "%Y-%m-%dT%H:%M:%S")

#DALLAS
#Closests index
closest_indices <- sapply(dal_result_clean$utc_time, function(x) {
  closest_index <- which.min(abs(x - dal_precip$utc_time))
  return(closest_index)
})

# Merge 
dal_result_p <- cbind(dal_result_clean, dal_precip[closest_indices, ])

#HOUSTON
#closest index
closest_indices <- sapply(hou_result_clean$utc_time, function(x) {
  closest_index <- which.min(abs(x - hou_precip$utc_time))
  return(closest_index)
})
# Merge 
hou_result_p <- cbind(hou_result_clean, hou_precip[closest_indices, ])

#rename precipitation columns
names(dal_result_p)[names(dal_result_p) == "PRECIP..IN."] <- "precip"
names(hou_result_p)[names(hou_result_p) == "PRECIP..IN."] <- "precip"

#Delete duplicate utc_time columns
columns_to_remove <- c("utc_time.1", "utc_time.2", "Peak_Wind.Direction..Degrees.",
                       "Peak_Wind.Speed..mph.","wind_direction_deg","Wind.Gust.10m..kn.",
                       "wind_speed_m_s")
dal_result_p <- dal_result_p[, -which(names(dal_result_p) %in% columns_to_remove)]

hou_result_p <- hou_result_p[, -which(names(hou_result_p) %in% columns_to_remove)]

#save as .csv
write.csv(dal_result_p, "dal_result_p.csv", row.names = F)
write.csv(hou_result_p, "hou_result_p.csv", row.names = F)
  #add pm-- ljoin#####
#load output from last step
setwd("~/Documents/OU/workspace/SensorDataProcessor")
dal_result_p <- read.csv("dal_result_p.csv")
hou_result_p <- read.csv("hou_result_p.csv")

# #load original pm data
# setwd("~/Documents/OU/workspace/SensorDataProcessor/dallas_weather")
# 
# dallas_pm_2022 <- read.csv ("dallas_pm_2022.csv")
# dallas_pm_2023 <- read.csv ("dallas_pm_2023.csv")
# dallas_pm <- rbind(dallas_pm_2022, dallas_pm_2023)
# 
# setwd("~/Documents/OU/workspace/SensorDataProcessor/houston_weather")
# 
# houston_pm_2022 <- read.csv ("houston_pm_2022.csv")
# houston_pm_2023 <- read.csv ("houston_pm_2023.csv")
# houston_pm <- rbind(houston_pm_2022, houston_pm_2023)
# 
# setwd("~/Documents/OU/workspace/SensorDataProcessor")
# write.csv(dallas_pm, "dallas_pm.csv", row.names = F)
# write.csv(houston_pm, "houston_pm.csv", row.names = F)

#load simplified pm data
setwd("~/Documents/OU/workspace/SensorDataProcessor")
dallas_pm <- read.csv("dallas_pm.csv")
houston_pm <- read.csv("houston_pm.csv")

colnames(dallas_pm)[5] = "pm"
colnames(houston_pm)[5] = "pm"

#units of pm are daily mean in ug / m3

#convert DATE objects
library(lubridate)
dal_result_p$Date <- ymd(dal_result_p$Date)
dallas_pm$Date <- ymd(dallas_pm$Date)

average_dallas_pm <- dallas_pm %>%
  group_by(Date) %>%
  summarise_at(vars(pm), list(pm=mean))

houston_pm$Date <- ymd(houston_pm$Date) 
hou_result_p$Date <- ymd(hou_result_p$Date)

average_houston_pm <- houston_pm %>%
  group_by(Date) %>%
  summarise_at(vars(pm), list(pm=mean))

#left join

dal_result_pm <- left_join(dal_result_p, average_dallas_pm, by = "Date")
hou_result_pm <- left_join(hou_result_p,average_houston_pm, by = "Date")

#Add city and save file
dal_result_pm$City <- "Dallas"
hou_result_pm$City <- "Houston"

setwd("~/Documents/OU/workspace/SensorDataProcessor")
write.csv(dal_result_pm, "dal_result_pm.csv", row.names=F)
write.csv(hou_result_pm, "hou_result_pm.csv", row.names=F)

  #add wind -- closest index#####
setwd("~/Documents/OU/workspace/SensorDataProcessor/dallas_weather")
dal_wind <- read.csv("dal_wind.csv")

setwd("~/Documents/OU/workspace/SensorDataProcessor/houston_weather")
hou_wind <- read.csv("hou_wind.csv")

setwd("~/Documents/OU/workspace/SensorDataProcessor")
dal_result_pm <- read.csv("dal_result_pm.csv")
hou_result_pm <- read.csv("hou_result_pm.csv")


colnames(dal_wind)[2] = "utc_time"
colnames(hou_wind)[2] = "utc_time"


#select only these two because on the correct time interval
dal_wind <- dal_wind[, c("utc_time", "wind_speed_m_s", "wind_direction_deg")]
hou_wind <- hou_wind[,c("utc_time", "wind_speed_m_s", "wind_direction_deg")]

dal_wind$utc_time <- as.POSIXct(dal_wind$utc_time, format = "%Y-%m-%dT%H:%M:%S")
hou_wind$utc_time <- as.POSIXct(hou_wind$utc_time, format = "%Y-%m-%dT%H:%M:%OS")

dal_result_pm$utc_time <- as.POSIXct(dal_result_pm$utc_time, format = "%Y-%m-%d %H:%M:%S")
hou_result_pm$utc_time <- as.POSIXct(hou_result_pm$utc_time, format = "%Y-%m-%d %H:%M:%OS")

#matches closest time interval wind speed and direction only work  
closest_indices <- sapply(dal_result_pm$utc_time, function(x) {
  closest_index <- which.min(abs(x - dal_wind$utc_time))
  return(closest_index)
})
# Merge
dal_final_result <- cbind(dal_result_pm, dal_wind[closest_indices, ])

#HOUSTON
# Find closest timestamps from 'hou_result_pm' in 'dal_weather'
closest_indices <- sapply(hou_result_pm$utc_time, function(x) {
  closest_index <- which.min(abs(x - hou_wind$utc_time))
  return(closest_index)
})

# Merge
hou_final_result <- cbind(hou_result_pm, hou_wind[closest_indices, ])




##Filter out missing entries + trim 95%######
#Dallas time trim
dal_final_result <- dal_final_result[dal_final_result$local_time < "2024-01-11T23:00:05.000",]
#Houston Time Trim
hou_final_result <- hou_final_result[hou_final_result$local_time < "2024-01-11T23:00:05.000",]

#combine
dal_hou_msas <- rbind(dal_final_result, hou_final_result)
dal_hou_msas <- dal_hou_msas[, !colnames(dal_hou_msas) %in% "utc_time.1"]

#complete cases   
# columns_to_check <- c("pressure_Pa", "rel_humidity", "precip", "pm")
# 
# dal_hou_msas_filtered <- dal_hou_msas[complete.cases(dal_hou_msas[, columns_to_check]), ]
#Don't have to do this if not doing models...
#eliminates 360 results

#trim upper lower 95%
library(dplyr)

# Calculate quantiles by location and date
quantiles_by_date <- dal_hou_msas %>%
  group_by(location, Date) %>%
  summarise(p5_msas = quantile(msas, probs = 0.05, na.rm = TRUE),
            p95_msas = quantile(msas, probs = 0.95, na.rm = TRUE))

# Join the quantiles back to the original data and filter
dal_hou_msas_filtered <- dal_hou_msas %>%
  inner_join(quantiles_by_date, by = c("location", "Date")) %>%
  filter(msas >= p5_msas & msas <= p95_msas)

#erase obs. over 18
dal_hou_msas_filtered <- dal_hou_msas_filtered[dal_hou_msas_filtered$msas <18,] 
#eliminates another 278 results

hist(dal_hou_msas_filtered$msas, main= "95th Filtered Histogram")



#Summarize mean by night####
library(dplyr)

summary_data <- dal_hou_msas_filtered %>%
  group_by(location, Date) %>%
  summarise(msas_mean = mean(msas),
            msas_sd = sd(msas))
summary_data$Date <- as.Date(summary_data$Date)

#Data viz#####
#examples using dal_hou_msas_filtered and summary_data

#boxplot msas by location
boxplot(msas ~ location, data = dal_hou_msas_filtered, col = "lightblue", main = "Sky Brightness by Location")

library(ggplot2)
# Plot summary_data msas over time with facet_wrap and error bars
ggplot(summary_data, aes(x = Date, y = msas_mean, color = location)) +
  geom_point() +
  geom_errorbar(aes(ymin = msas_mean - msas_sd, ymax = msas_mean + msas_sd), width = 0.2) +
  facet_wrap(~ location, scales = "free") +
  coord_cartesian(ylim = c(12, 18)) +
  labs(x = "Date", y = "Sky Brightness (msas)", title = "75th Quantile Mean Sky Brightness Over Time by Location with SD")+
  scale_x_date(date_labels = "%b", date_breaks = "1 month") +
  theme_minimal()

#No SD bars
ggplot(summary_data, aes(x = Date, y = msas_mean, color = location)) +
  geom_point() +
  facet_wrap(~ location, scales = "free") +
  coord_cartesian(ylim = c(12, 18)) +
  labs(x = "Date", y = "Sky Brightness (msas)", title = "75th Quantile Mean Sky Brightness Over Time by Location")+
  scale_x_date(date_labels = "%b", date_breaks = "1 month") +
  theme_minimal() 


#two locations for comparison
compare_locals <- dal_hou_msas_filtered[dal_hou_msas_filtered$location == "City"|
                                          dal_hou_msas_filtered$location == "Audubon",]
compare_locals$Date <- as.Date(compare_locals$Date)

ggplot(data = compare_locals, aes(x = Date, y = msas)) +
  geom_point(aes(color=location)) +
  scale_color_manual(values = c("#FFCC80", "#6A0DAD")) +
  theme_classic() +
  theme(legend.position = "none") +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1, vjust = 0.5, size = 10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "City Hall vs. Audubon Center Brightness Over Time",
    y = "Brightness (MSAS)",
    x = " "
  ) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_date(date_breaks = "1 month", date_labels = "%B")+ 
  geom_vline(xintercept = as.Date(c("2022-09-05", "2022-10-29")), 
             color = "red", linetype = "dashed", linewidth = 1) #add lines with dates of interest

