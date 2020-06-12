##Code Snippets Usefull
readFun <- function( filename ) {
  
  # read in the data
  data <- read.csv( filename, 
                    header = FALSE, 
                    col.names = c( "Name", "Gender", "Count" ) )
  
  # add a "Year" column by removing both "yob" and ".txt" from file name
  data$Year <- gsub( "yob|.txt", "", filename )
  
  return( data )
}

# execute that function across all files, outputting a data frame
doMC::registerDoMC( cores = 4 )
babynames <- plyr::ldply( .data = list.files(pattern="*.txt"),
                          .fun = readFun,
                          .parallel = TRUE )



#read in required packages
data_path = "C:/Users/Martin/Documents/Projekt_PTM_merge/results_PTM_melanoma/results"
require(data.table)
setDTthreads(threads = 0)

setwd(data_path)

#create a list of the files from your target directory
file_list <- list.files(path=data_path)

#initiate a blank data frame, each iteration of the loop will append the data from the given file to this variable
dataset <- data.frame()

#had to specify columns to get rid of the total column
for (i in 1:length(file_list)){
  temp_data <- fread(file_list[i], stringsAsFactors = F) #read in files using the fread function from the data.table package
  #temp_data2 = temp_data[,2:16]
  dataset <- rbindlist(list(dataset, temp_data), use.names = T, fill=TRUE) #for each iteration, bind the new data to the building dataset
}

files = lapply(files, basename)
for (i in files) {
  fl <- paste0(i)
  fl
}


ddat <- as.list(rep("", 20))
for(i in 1:20) {
  ddat[[i]] <- data.frame(ivec = 1:i)
  #other processing..
}


for(i in files) {
  assign(paste("d",i,sep="_"),data)
}


###Grouping and mean
airquality <- data.frame(City = c("CityA", "CityA","CityA",
                                  "CityB","CityB","CityB",
                                  "CityC", "CityC"),
                         year = c("1990", "1990", "2010", "1990", 
                                  "2000", "2010", "2000", "2010"),
                         month = c("June", "July", "August",
                                   "June", "July", "August",
                                   "June", "August"),
                         PM10 = c(runif(3), rnorm(5)),
                         PM25 = c(runif(3), rnorm(5)),
                         Ozone = c(runif(3), rnorm(5)),
                         CO2 = c(runif(3), rnorm(5)))
airquality

library(dplyr)
airquality %>%
  group_by(City, year) %>% 
  summarise_at(vars("PM25", "Ozone", "CO2"), mean)
(0.4513109+0.0416877)/2
