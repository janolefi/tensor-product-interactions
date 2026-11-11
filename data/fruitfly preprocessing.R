library(dplyr)

# data = readRDS("./data/fruitflies_genotype1_raw.rds")
data = readRDS("./data/fruitflies_genotype0_raw.rds")


## remove columns with sine and cosine
data = data %>% select(-contains("sin"), -contains("cos"))

## set activity while startling to NA
data$activity[data$startling == 1 | data$startling_evening == 1] = NA

## on 24 hour scale
data$tod = data$tod / 2

# reorder and deselect columns
data = data %>%
  select(
    ID, ID_animal, time, genotype, activity, condition, tod, light,
  )

# dropping unused levels
data$ID = droplevels(data$ID)
data$ID_animal = droplevels(data$ID_animal)

data = data %>%
  select(
    aniID = ID_animal, trackID = ID, time, genotype, activity, condition, tod, light,
  )

aniIDs = unique(data$aniID)

# using a subset only for genotype 1
set.seed(123)
subIDs = sample(aniIDs, 35, replace = FALSE)

data = data %>% filter(aniID %in% subIDs)

data$aniID = droplevels(data$aniID)
data$trackID = droplevels(data$trackID)

# saveRDS(data, file = "./data/fruitflies_genotype1.rds")
saveRDS(data, file = "./data/fruitflies_genotype0.rds")
