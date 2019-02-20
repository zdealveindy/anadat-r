# Importing data into R and preparing them for analysis

# working directory
getwd ()  # get working directory:  "C:/Users/Bob/Documents"
setwd (dir = "c:/Users/Bob/Documents/quickRintro/data/")

# traits <- read.table (file = 'vltava.traits.txt', sep = '\t', header = TRUE, row.names = 1)

traits <- read.delim (file = 'vltava.traits.txt', row.names = 1)  # reading data from working directory (need to be set by setwd function)
# traits <- read.delim (file = 'c:/Users/Bob/Documents/quickRintro/data/vltava.traits.txt', row.names = 1) # specifying the whole path to the file (no need to set working directory by setwd function)

summary (traits)  # summarise the dataframe (also see how many missing values are there)

hist (traits$plant.height)  # histogram of raw plant.height data
hist (log (traits$plant.height)) # histogram of log-transformed plant.height data

traits$plant.height_t <- traits$plant.height  # creating new variable which will store the values after removing outliers and transformation

traits$plant.height_t [traits$plant.height_t > 6] <- NA  # replacing outliers by NA values
hist (traits$plant.height_t)
traits$plant.height_t <- log (traits$plant.height_t)  # log-transforming the variable (already without outliers)
hist (traits$plant.height_t)  # histogram of final distribution
