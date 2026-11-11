library(dplyr)
library(zoo)

# Raw data not on Github because too large. See https://zenodo.org/records/3768080

## data preprocessing

# loading
data <- read.delim("./data/muskox_tracking_data_winter_20200417.txt", header = TRUE)

# renaming
data$ID = data$burst_id

# selecting relevant variables only
data <- data %>% select(
  ID, step, angle, x, y, datetime, tday, month, julian
)

# time formatting
data$datetime <- as.POSIXct(data$datetime, format = "%Y-%m-%d %H:%M:%S")

# linearly interpolate missing locations x,y (missing covariate values not possible)
data[c("x", "y")] <- data.frame(lapply(data[c("x", "y")], na.approx, na.rm = FALSE))

## save
saveRDS(data, "./data/muskox_winter.rds")
