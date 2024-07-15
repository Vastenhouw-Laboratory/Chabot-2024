# Noémie Chabot, PhD, Vastenhouw lab, 2024/07/10
# R code for analysis and figures related to the Fig1d,e,f, Fig2a-e and FigS3a,b of the paper
# "Local DNA compaction creates TF-DNA clusters that enable transcription"
# As of submission to the journal NSMB

#Libraries ----
#Read all the librairies that might be required to run/plot the code below
library("ggplot2")
library("tidyverse")
library("dplyr")
library("viridisLite")
library("viridis")
library("data.table")

options(scipen=999)

#Change path to working directory
Input = "/Volumes/CIG/nvastenh/nanog_movement/D2c/Manuscript/240612_NCB_NSMB_submission/Material_make_figures/Fig1d,e,f_Fig2a,b,c,d,e_FigS3a,b/All_figures_Statistics/"

#=========================================================
# Read the files necessary for the analysis
#=========================================================


#Change path to working directory
setwd(Input)

############Pol II S5P
{
#Get to the good directory

#Get all the files in the working directory labeled "S5P"
List_directories_S5P <- grep("S5P_Statistics", list.dirs(path = ".", full.names = TRUE), value=TRUE)

#Make a list with all the statistics files to associate to the data frame
{
  List_files_S5P_number <- data.frame()
  S5P_statistics_total <- data.frame()
  dataset <- data.frame()
  Infos_file <- data.frame()
  Infinity <- data.frame()
  
  for(e in 1:length(List_directories_S5P)) {
    
    #Assemble all data in a unique dataframe
    print(length(list.files(List_directories_S5P[e])))
    print(List_directories_S5P[e])
    if(e == 1) {List_files_S5P_number <- data.frame(list.files(path = List_directories_S5P[e], pattern="S5P", recursive=TRUE)) } else {
      List_files_S5P_number <- List_files_S5P_number %>% data.frame(list.files(path = List_directories_S5P[e], pattern="S5P", recursive=TRUE))}
  }
}

#Read each statistic file and associate it in one unique data frame
{
  
  for(e in 1:length(List_files_S5P_number)) {
    setwd(Input)
    S5P_statistics_all <- data.frame()
    List_files_S5P_number <- list.files(path = List_directories_S5P[e], pattern="S5P", recursive=TRUE)
    setwd(List_directories_S5P[e])
    
    for (i in 1:length(List_files_S5P_number))
    {
      
      temp_data <- read.table(List_files_S5P_number[i],  sep=",", skip=3, header=TRUE) #each file will be read in, specify which columns you need read in to avoid any errors
      file_name <- as.character(List_files_S5P_number[i])
      
      #Select the info inside file name about date of acquisition and identity of the nucleus being analysing. These infos are stored in columns.
      #Adapt syntax based on the number of characters for different stages
      
      #Select the info inside file name about date of acquisition and identity of the nucleus being analysing. These infos are stored in columns.
      #Adapt syntax based on the number of characters.
      
      if (grepl("128|256|512", file_name)==TRUE){
        Infos_file[1,1] <- substr(file_name, 1,5)
        Infos_file[1,2] <- substr(file_name, 7,9)
        Infos_file[1,3] <- as.numeric(gsub(".*?([0-9]+).*", "\\1", substr(file_name, 10,12)))
        Infos_file[1,4] <- as.numeric(gsub(".*?([0-9]+).*", "\\1", substr(file_name, 13,14)))
      } else if (grepl("1k", file_name)==TRUE) {
        Infos_file[1,1] <- substr(file_name, 1,5)
        Infos_file[1,2] <- substr(file_name, 7,8)
        Infos_file[1,3] <- as.numeric(gsub(".*?([0-9]+).*", "\\1", substr(file_name, 9,11)))
        Infos_file[1,4] <- as.numeric(gsub(".*?([0-9]+).*", "\\1", substr(file_name, 12,13)))
      } else if (grepl("High", file_name)==TRUE) {
        Infos_file[1,1] <- substr(file_name, 1,5)
        Infos_file[1,2] <- substr(file_name, 7,10)
        Infos_file[1,3] <- as.numeric(gsub(".*?([0-9]+).*", "\\1", substr(file_name, 11,13)))
        Infos_file[1,4] <- as.numeric(gsub(".*?([0-9]+).*", "\\1", substr(file_name, 14,15)))
      }
      
      #Repeat n nows the infos about the nuclei, as much as they are segmented volumes in the file
      Infinity <- Infos_file[rep(seq_len(nrow(Infos_file)), each = nrow(temp_data)),]
      
      #Give a name to the column and make sure it looks like something
      colnames(Infinity) <- c("Date", "Stage", "Stack","Nucleus")
      print(List_files_S5P_number[i])
      
      temp_data <- cbind(Infinity, temp_data)
      
      #Assemble all data in a unique dataframe
      if(i == 1) {S5P_statistics_all <- temp_data } else {
        if(grepl("Time", file_name)==TRUE)
        {colnames(temp_data)[5] <- "Seconds"
        colnames(temp_data)[8] <- "Time"
        S5P_statistics_all <- S5P_statistics_all %>% full_join(temp_data, by=c("ID","Date","Stage","Nucleus","Stack","Time"))} 
        else {S5P_statistics_all <- S5P_statistics_all %>% full_join(temp_data, by=c("ID","Date","Stage","Nucleus","Stack","Time"))}
      }
      
      #Remove useless columns
      patterns <- c("Collection*", "Unit*","Image*","Channel*","X.x","X.y","Event*","Category*","TrackID*") 
      S5P_statistics_all <- S5P_statistics_all[,-grep(paste(patterns, collapse="|"), colnames(S5P_statistics_all))]
      
    }
    
    #Mix data frames with different number of columns by row
    if (e == 1){S5P_statistics_total <- S5P_statistics_all} else{S5P_statistics_total <- bind_rows(S5P_statistics_total, S5P_statistics_all)}
  }
  
  
  #Assign one unique nucleus ID for each nucleus
  {
    Unique_ID <- unique(S5P_statistics_total[,c(1:4)])
    S5P_statistics_total$Nucleus.number <- 0
    
    for(i in 1:nrow(S5P_statistics_total))
    {
      for(e in 1:nrow(Unique_ID))
      {
        if(S5P_statistics_total[i,1] == Unique_ID[e,1]&S5P_statistics_total[i,2] == Unique_ID[e,2]&S5P_statistics_total[i,3] == Unique_ID[e,3]&S5P_statistics_total[i,4] == Unique_ID[e,4])
          
          S5P_statistics_total[i,"Nucleus.number"] <- e
        
      }
    }
  }
}
}

#############Nanog clusters
{
#Get to the good directory
setwd(Input)

List_directories_Nanog <- grep("Nanog_Statistics", list.dirs(path = ".", full.names = TRUE, recursive = TRUE), value=TRUE)

#Associate all Nanog file together

List_files_Nanog_number <- data.frame()
Nanog_statistics_total <- data.frame()

dataset <- data.frame()
Infos_file <- data.frame()
Infinity <- data.frame()

for(e in 1:length(List_directories_Nanog)) {
  setwd(Input)
  
  #Assemble all data in a unique dataframe
  print(length(list.files(List_directories_Nanog[e])))
  print(List_directories_Nanog[e])
  if(e == 1) {List_files_Nanog_number <- data.frame(list.files(path = List_directories_Nanog[e], pattern="Nanog", recursive=TRUE)) } else {
    List_files_Nanog_number <- List_files_Nanog_number %>% data.frame(list.files(path = List_directories_Nanog[e], pattern="Nanog", recursive=TRUE))}
}

for(e in 1:length(List_files_Nanog_number)) {
  
  Nanog_statistics_all <- data.frame()
  setwd(Input)
  List_files_Nanog_number <- list.files(path = List_directories_Nanog[e], pattern="Nanog", recursive=TRUE)
  setwd(List_directories_Nanog[e])
  
  for (i in 1:length(List_files_Nanog_number))
  {
    
    temp_data <- read.table(List_files_Nanog_number[i],  sep=",", skip=3, header=TRUE) #each file will be read in, specify which columns you need read in to avoid any errors
    file_name <- as.character(List_files_Nanog_number[i])
    
    #Select the info inside file name about date of acquisition and identity of the nucleus being analysing. These infos are stored in columns.
    #Adapt syntax based on the number of characters
    
    #Select the info inside file name about date of acquisition and identity of the nucleus being analysing. These infos are stored in columns.
    #Adapt syntax based on the number of characters
    
    if (grepl("128|256|512", file_name)==TRUE){
      Infos_file[1,1] <- substr(file_name, 1,5)
      Infos_file[1,2] <- substr(file_name, 7,9)
      Infos_file[1,3] <- as.numeric(gsub(".*?([0-9]+).*", "\\1", substr(file_name, 10,12)))
      Infos_file[1,4] <- as.numeric(gsub(".*?([0-9]+).*", "\\1", substr(file_name, 13,14)))
    } else if (grepl("1k", file_name)==TRUE) {
      Infos_file[1,1] <- substr(file_name, 1,5)
      Infos_file[1,2] <- substr(file_name, 7,8)
      Infos_file[1,3] <- as.numeric(gsub(".*?([0-9]+).*", "\\1", substr(file_name, 9,11)))
      Infos_file[1,4] <- as.numeric(gsub(".*?([0-9]+).*", "\\1", substr(file_name, 12,13)))
    } else if (grepl("High", file_name)==TRUE) {
      Infos_file[1,1] <- substr(file_name, 1,5)
      Infos_file[1,2] <- substr(file_name, 7,10)
      Infos_file[1,3] <- as.numeric(gsub(".*?([0-9]+).*", "\\1", substr(file_name, 11,13)))
      Infos_file[1,4] <- as.numeric(gsub(".*?([0-9]+).*", "\\1", substr(file_name, 14,15)))
    }
    
    #Repeat n nows the infos about the nuclei, as much as they are segmented volumes in the file
    Infinity <- Infos_file[rep(seq_len(nrow(Infos_file)), each = nrow(temp_data)),]
    
    #Give a name to the column and make sure it looks like something
    colnames(Infinity) <- c("Date", "Stage", "Stack","Nucleus")
    print(List_files_Nanog_number[i])
    
    temp_data <- cbind(Infinity, temp_data)
    
    #Assemble all data in a unique dataframe
    if(i == 1) {Nanog_statistics_all <- temp_data } else {
      if(grepl("Time", file_name)==TRUE)
      {colnames(temp_data)[5] <- "Seconds"
      colnames(temp_data)[8] <- "Time"
      Nanog_statistics_all <- Nanog_statistics_all %>% full_join(temp_data, by=c("ID","Date","Stage","Nucleus","Stack","Time","TrackID"))} 
      else {Nanog_statistics_all <- Nanog_statistics_all %>% full_join(temp_data, by=c("ID","Date","Stage","Nucleus","Stack","Time","TrackID"))}
    }
    
    #Remove useless columns
    patterns <- c("Collection*", "Unit*","Image*","Channel*","X.x","X.y","Event*","Category*") 
    Nanog_statistics_all <- Nanog_statistics_all[,-grep(paste(patterns, collapse="|"), colnames(Nanog_statistics_all))]
    
  }
  
  #Mix data frames with different number of columns by row
  if (e == 1){
    print("A")
    Nanog_statistics_total <- Nanog_statistics_all
  } else{
    print("B")
    Nanog_statistics_total <- bind_rows(Nanog_statistics_total, Nanog_statistics_all)
  }
}

#Removing the Nanog that are not at the mir430 locus because they are not tracked and they don't have a TrackID
Nanog_statistics_total <-Nanog_statistics_total[complete.cases(Nanog_statistics_total[,c("TrackID")]),]

#Assign one unique nucleus ID for each nucleus

Unique_ID <- unique(unique(Nanog_statistics_total[,c(1:4)]))
Nanog_statistics_total$Nucleus.number <- 0

for(i in 1:nrow(Nanog_statistics_total))
{
  for(e in 1:nrow(Unique_ID))
  {
    if(Nanog_statistics_total[i,1] == Unique_ID[e,1]&Nanog_statistics_total[i,2] == Unique_ID[e,2]&Nanog_statistics_total[i,3] == Unique_ID[e,3]&Nanog_statistics_total[i,4] == Unique_ID[e,4])
      
      Nanog_statistics_total[i,"Nucleus.number"] <- e
    
  }
}

Unique_ID <- unique(Nanog_statistics_total[,c(1,3,4,7)])
Nanog_statistics_total$Unique_track_Nanog <- 0

for(i in 1:nrow(Nanog_statistics_total))
{
  for(e in 1:nrow(Unique_ID))
  {
    if(Nanog_statistics_total[i,1] == Unique_ID[e,1]&Nanog_statistics_total[i,3] == Unique_ID[e,2]&Nanog_statistics_total[i,4] == Unique_ID[e,3]&Nanog_statistics_total[i,7] == Unique_ID[e,4])
      
      Nanog_statistics_total[i,"Unique_track_Nanog"] <- e
    
  }
}



#Assign one unique TrackID for each nucleus

Nanog_statistics_total$Date[Nanog_statistics_total$Date == '06.10'] <- '06.11'

Nanog_statistics_total$Inhibition <- "No"
Nanog_statistics_total[Nanog_statistics_total$Date %in% c("10.08","17.01","27.07"),"Inhibition"] <- "Yes"
}

#############Associate Nanog to the corresponding Pol II S5P signal
{
#Calculate the distance between each pair of S5P-Nanog foci, and check if this distance if < 0.5 using the Shortest distance statistics from Imaris and < 1 for the distance calculated from the position of the Nanog spots. 
#If yes, join the two data frames at this specific line. If not, write just nothing.
#Add all columns from S5P file into Nanog file

colnames(S5P_statistics_total) <- paste0("S5P_",colnames(S5P_statistics_total))

Nanog_statistics_total[c(colnames(S5P_statistics_total))] <- NA
Nanog_statistics_total$Dist.Nanog.S5P <- NA

which(colnames(Nanog_statistics_total)=="S5P_Date")
max(col(Nanog_statistics_total))

for(r in 1:max(as.numeric((S5P_statistics_total$S5P_Nucleus.number)))) {
  
  if(nrow(S5P_statistics_total[(which(S5P_statistics_total$S5P_Nucleus.number==r)),])!=0) {
    
    for(w in min(S5P_statistics_total[(which(S5P_statistics_total$S5P_Nucleus.number==r)),]$S5P_Time):max(S5P_statistics_total[(which(S5P_statistics_total$S5P_Nucleus.number==r)),]$S5P_Time)) {
      
      Nanog_subset_time <- Nanog_statistics_total[(which(Nanog_statistics_total$Time==w&Nanog_statistics_total$Nucleus.number==r)),]
      S5P_subset_time <- S5P_statistics_total[(which(S5P_statistics_total$S5P_Time==w&S5P_statistics_total$S5P_Nucleus.number==r)),]
      #print(w)
      
      if(nrow(Nanog_subset_time)!=0) {
        
        for(a in 1:nrow(S5P_subset_time)) {
          
          S5P_x <- S5P_subset_time$S5P_Position.X[a]
          S5P_y <- S5P_subset_time$S5P_Position.Y[a]
          S5P_z <- S5P_subset_time$S5P_Position.Z[a]
          
          if(nrow(S5P_subset_time) != 0) {
            
            #print(S5P_statistics_total[a,])
            
            for(b in 1:nrow(Nanog_subset_time)) {
              
              Nanog_x <-  Nanog_subset_time$Position.X[b]
              Nanog_y <- Nanog_subset_time$Position.Y[b]
              Nanog_z <- Nanog_subset_time$Position.Z[b]
              
              Distance_Nanog_S5P <- sqrt((Nanog_x - S5P_x)^2 + (Nanog_y - S5P_y)^2 + (Nanog_z - S5P_z)^2)
              
              if(Distance_Nanog_S5P < 1 & Nanog_subset_time$Shortest.Distance.to.Surfaces[b] < 0.5) {
                
                Nanog_statistics_total[(which(Nanog_statistics_total$Time==w&Nanog_statistics_total$Nucleus.number==r)),34:60][b,] <- S5P_statistics_total[(which(S5P_statistics_total$S5P_Time==w&S5P_statistics_total$S5P_Nucleus.number==r)),][a,]
                Nanog_statistics_total[(which(Nanog_statistics_total$Time==w&Nanog_statistics_total$Nucleus.number==r)),61][b] <- Distance_Nanog_S5P
                
              }
            }
          }
        }
      }
    }
  }
}
}

#############Nanog subclusters (tracking for Nanog merging clusters for Fig.2e)
{
  #To measure the fluorescence of individual merging clusters, we tracked the merging clusters that are splitting for at least three consecutive time points. 
  #Get the statistics for these specific files, make a dataframe and fuse it to the main dataframe 
  #Get to the good directory
  setwd("/Volumes/CIG/nvastenh/nanog_movement/D2c/Manuscript/240612_NCB_NSMB_submission/Material_make_figures/Fig1d,e,f_Fig2a,b,c,d,e_FigS3a,b/Fig2e_Statistics_subclusters/b/")
  
  #Loading the name of the directories
  List_directories_Nanog <- grep("subclusters_Statistics", list.dirs(path = ".", full.names = TRUE, recursive = TRUE), value=TRUE)
  
  #Associate all Nanog file together
  
  List_files_Nanog_number <- data.frame()
  
  
  Nanog_statistics_total_sub <- data.frame()
  
  dataset <- data.frame()
  Infos_file <- data.frame()
  Infinity <- data.frame()
  
  for(e in 1:length(List_directories_Nanog)) {
    setwd("/Volumes/CIG/nvastenh/nanog_movement/D2c/Manuscript/240612_NCB_NSMB_submission/Material_make_figures/Fig1d,e,f_Fig2a,b,c,d,e_FigS3a,b/Fig2e_Statistics_subclusters/b/")
    
    #Assemble all data in a unique dataframe
    print(length(list.files(List_directories_Nanog[e])))
    print(List_directories_Nanog[e])
    if(e == 1) {List_files_Nanog_number <- data.frame(list.files(path = List_directories_Nanog[e], pattern="Nanog", recursive=TRUE)) } else {
      List_files_Nanog_number <- List_files_Nanog_number %>% data.frame(list.files(path = List_directories_Nanog[e], pattern="Nanog", recursive=TRUE))}
  }
  
  for(e in 1:length(List_files_Nanog_number)) {
    
    Nanog_statistics_all <- data.frame()
    setwd("/Volumes/CIG/nvastenh/nanog_movement/D2c/Manuscript/240612_NCB_NSMB_submission/Material_make_figures/Fig1d,e,f_Fig2a,b,c,d,e_FigS3a,b/Fig2e_Statistics_subclusters/b/")
    List_files_Nanog_number <- list.files(path = List_directories_Nanog[e], pattern="Nanog", recursive=TRUE)
    setwd(List_directories_Nanog[e])
    
    for (i in 1:length(List_files_Nanog_number))
    {
      
      temp_data <- read.table(List_files_Nanog_number[i],  sep=",", skip=3, header=TRUE) #each file will be read in, specify which columns you need read in to avoid any errors
      file_name <- as.character(List_files_Nanog_number[i])
      
      #Select the info inside file name about date of acquisition and identity of the nucleus being analysing. These infos are stored in columns.
      #Adapt syntax based on the number of characters
      
      #Select the info inside file name about date of acquisition and identity of the nucleus being analysing. These infos are stored in columns.
      #Adapt syntax based on the number of characters
      
      Infos_file[1,1] <- substr(file_name, 1,5)
      Infos_file[1,2] <- substr(file_name, 7,8)
      Infos_file[1,3] <- as.numeric(gsub(".*?([0-9]+).*", "\\1", substr(file_name, 9,11)))
      Infos_file[1,4] <- as.numeric(gsub(".*?([0-9]+).*", "\\1", substr(file_name, 12,13)))
      
      #Repeat n nows the infos about the nuclei, as much as they are segmented volumes in the file
      Infinity <- Infos_file[rep(seq_len(nrow(Infos_file)), each = nrow(temp_data)),]
      
      #Give a name to the column and make sure it looks like something
      colnames(Infinity) <- c("Date", "Stage", "Stack","Nucleus")
      print(List_files_Nanog_number[i])
      
      temp_data <- cbind(Infinity, temp_data)
      
      #Assemble all data in a unique dataframe
      if(i == 1) {Nanog_statistics_all <- temp_data } else {
        if(grepl("Time", file_name)==TRUE)
        {colnames(temp_data)[5] <- "Seconds"
        colnames(temp_data)[8] <- "Time"
        Nanog_statistics_all <- Nanog_statistics_all %>% full_join(temp_data, by=c("ID","Date","Stage","Nucleus","Stack","Time","TrackID"))} 
        else {Nanog_statistics_all <- Nanog_statistics_all %>% full_join(temp_data, by=c("ID","Date","Stage","Nucleus","Stack","Time","TrackID"))}
      }
      
      #Remove useless columns
      patterns <- c("Collection*", "Unit*","Image*","Channel*","X.x","X.y","Event*","Category*") 
      Nanog_statistics_all <- Nanog_statistics_all[,-grep(paste(patterns, collapse="|"), colnames(Nanog_statistics_all))]
      
    }
    
    #Mix data frames with different number of columns by row
    if (e == 1){
      print("A")
      Nanog_statistics_total_sub <- Nanog_statistics_all
    } else{
      print("B")
      Nanog_statistics_total_sub <- bind_rows(Nanog_statistics_total_sub, Nanog_statistics_all)
    }
  }
  
  Nanog_statistics_total_sub <- Nanog_statistics_total_sub[-which(is.na(Nanog_statistics_total_sub$TrackID)),]
  
  #Assemble the other Nanog dataframe with all the files to the one with only merging clusters
  Nanog_statistics_total <- Nanog_statistics_total %>% full_join(Nanog_statistics_total_sub, by=c("Date","Stage","Stack","Nucleus","Time","ID","Seconds","Volume"))
}

#############Nuclei
{
#All file names in one table
{
  setwd(Input)
  
  
  #Nucleus
  #Get a list of files containing the word "Nucleus" to analyze only the Nucleus foci
  
  setwd(Input)
  
  
  List_directories_Nuc <- grep(pattern="Nucleus", list.dirs(path = ".", full.names = TRUE), value=TRUE)
  List_directories_Nucleus <- grep(pattern="1k", List_directories_Nuc, value=TRUE)
  
  
  List_files_Nucleus_number <- data.frame()
  
  #Get all the Nucleus file in columns for each nucleus
  
  for(e in 1:length(List_directories_Nucleus))
  {
    setwd(Input)
    
    #print(length(list.files(List_directories_Nucleus[e])))
    if(e == 1) {List_files_Nucleus_number <- data.frame(list.files(path = List_directories_Nucleus[e], recursive=TRUE)) }
    else {
      #Assemble all data in a unique dataframe
      List_files_Nucleus_number <- List_files_Nucleus_number %>% data.frame(list.files(path = List_directories_Nucleus[e], pattern="Nucleus", recursive=TRUE))}
  }
}

#Assemble all data in one dataframe

dataset <- data.frame()
Infos_file <- data.frame()
Infinity <- data.frame()

for(e in 1:length(List_directories_Nucleus)) {
  
  setwd(Input)
  List_files_Nucleus_number <- list.files(path = List_directories_Nucleus[e], pattern="Nucleus", recursive=TRUE)
  setwd(List_directories_Nucleus[e])
  
  for (i in 1:length(List_files_Nucleus_number)){
    temp_data <- read.table(List_files_Nucleus_number[i],  sep=",", skip=3, header=TRUE) #each file will be read in, specify which columns you need read in to avoid any errors
    file_name <- as.character(List_files_Nucleus_number[i])
    
    #Select the info inside file name about date of acquisition and identity of the nucleus being analysing. These infos are stored in columns.
    #Adapt syntax based on the number of characters
    
    if (grepl("128|256|512", file_name)==TRUE){
      Infos_file[1,1] <- substr(file_name, 1,5)
      Infos_file[1,2] <- substr(file_name, 7,9)
      Infos_file[1,3] <- as.numeric(gsub(".*?([0-9]+).*", "\\1", substr(file_name, 10,12)))
      Infos_file[1,4] <- as.numeric(gsub(".*?([0-9]+).*", "\\1", substr(file_name, 13,14)))
    } else if (grepl("1k", file_name)==TRUE) {
      Infos_file[1,1] <- substr(file_name, 1,5)
      Infos_file[1,2] <- substr(file_name, 7,8)
      Infos_file[1,3] <- as.numeric(gsub(".*?([0-9]+).*", "\\1", substr(file_name, 9,11)))
      Infos_file[1,4] <- as.numeric(gsub(".*?([0-9]+).*", "\\1", substr(file_name, 12,13)))
    } else if (grepl("High", file_name)==TRUE) {
      Infos_file[1,1] <- substr(file_name, 1,5)
      Infos_file[1,2] <- substr(file_name, 7,10)
      Infos_file[1,3] <- as.numeric(gsub(".*?([0-9]+).*", "\\1", substr(file_name, 11,13)))
      Infos_file[1,4] <- as.numeric(gsub(".*?([0-9]+).*", "\\1", substr(file_name, 14,15)))
    }
    
    #Repeat n nows the infos about the nuclei, as much as they are segmented volumes in the file
    Infinity <- Infos_file[rep(seq_len(nrow(Infos_file)), each = nrow(temp_data)),]
    #Infinity$V5 <- sub("^\\D*(\\d+).*$", "\\1", List_directories_Nucleus[e])
    
    #Give a name to the column and make sure it looks like something
    colnames(Infinity) <- c("Date","Stage","Stack","Nucleus")
    #print(List_files_Nucleus_number[i])
    temp_data <- cbind(Infinity, temp_data)
    
    #Assemble all data in a unique dataframe
    if(i == 1) {Nucleus_statistics_all <- temp_data } else {
      if(grepl("Time", file_name)==TRUE)
      {colnames(temp_data)[5] <- "Seconds"
      colnames(temp_data)[8] <- "Time"
      Nucleus_statistics_all <- Nucleus_statistics_all %>% full_join(temp_data, by=c("ID","Date","Nucleus","Stack","Time","Stage"))} 
      else {Nucleus_statistics_all <- Nucleus_statistics_all %>% full_join(temp_data, by=c("ID","Date","Nucleus","Stack","Time","Stage"))}
    }
    
  }
  #Remove useless columns
  patterns <- c("Collection*", "Unit*","Image*","Channel*","X.x","X.y","Event*","Category*","Merging.foci.x*")
  Nucleus_statistics_all <- Nucleus_statistics_all[,-grep(paste(patterns, collapse="|"), colnames(Nucleus_statistics_all))]
  
  if(e == 1) {Nucleus_statistics_total <- Nucleus_statistics_all
  } else {
    #Mix data frames with different number of columns by row
    Nucleus_statistics_total <- bind_rows(Nucleus_statistics_total, Nucleus_statistics_all)
  }
  
}

#Assign a unique nucleus ID to each nucleus 
Unique_ID <- unique(Nucleus_statistics_total[,c(1,3:4)])
Nucleus_statistics_total$Nucleus.number <- 0

for(i in 1:nrow(Nucleus_statistics_total))
{
  for(e in 1:nrow(Unique_ID))
  {
    if(Nucleus_statistics_total[i,1] == Unique_ID[e,1]&Nucleus_statistics_total[i,3] == Unique_ID[e,2]&Nucleus_statistics_total[i,4] == Unique_ID[e,3])
      
      Nucleus_statistics_total[i,"Nucleus.number"] <- e
    
  }
}


#Normalization the level inside the nucleus for the mean intensity

for(i in 1:max(Nucleus_statistics_total$Nucleus.number))
{
  for(e in Nucleus_statistics_total[which(Nucleus_statistics_total$Nucleus.number==i),"Time"])
  {
    
    Nucleus_statistics_total[which(Nucleus_statistics_total$Nucleus.number==i&Nucleus_statistics_total$Time==e),"Mean_intensity_normalized"] <-  Nucleus_statistics_total[which(Nucleus_statistics_total$Nucleus.number==i&Nucleus_statistics_total$Time==e),"Intensity.Mean.y"]/max(Nucleus_statistics_total[which(Nucleus_statistics_total$Nucleus.number==i),"Intensity.Mean.y"]) 
    
  }
}

Nucleus_statistics_total$Stack <- as.numeric(Nucleus_statistics_total$Stack)
Nanog_statistics_total$Date <- as.numeric(Nanog_statistics_total$Date)
Nucleus_statistics_total$Date <- as.numeric(Nucleus_statistics_total$Date)

#Make the dataframe that will be used for plotting the figure
Nanog_statistics_total_nuc <- Nanog_statistics_total[which(Nanog_statistics_total$Stage=="1k"&Nanog_statistics_total$Inhibition=="No"),] %>% full_join(Nucleus_statistics_total, by=c("Date","Stage","Stack","Nucleus","Time"))
}

#=========================================================
# Centering the time on transcription initiation
#=========================================================

#Make the variable numeric
Nanog_statistics_total$Time <- as.numeric(Nanog_statistics_total$Time)
Nanog_statistics_total$Seconds <- as.numeric(Nanog_statistics_total$Seconds)
Nanog_statistics_total$Unique_track <- as.numeric(Nanog_statistics_total$Unique_track)

#Remove all the lines that are not tracked
Nanog_statistics_total <- Nanog_statistics_total[complete.cases(Nanog_statistics_total[,c("Unique_track")]),]


#Normalize the time of activation
i=0

#Find the time of activation for each track by selecting the first time point when the S5P signal is detected
for(i in sort(unique(Nanog_statistics_total$Unique_track)))
{
  Nanog_statistics_total[which(Nanog_statistics_total$Unique_track==i),"Time_activation"] <- min(na.omit(Nanog_statistics_total[!is.na(Nanog_statistics_total$S5P_Nucleus.number)&Nanog_statistics_total$Unique_track==i,"Time"]))
  print(i)
}

#Substract time of the imaging to the time of trx activation to get the time centered on trx activation
Nanog_statistics_total$Time_normalized <- Nanog_statistics_total$Time - Nanog_statistics_total$Time_activation

#Do the same as previously but with the time in seconds
Nanog_statistics_total$Unique_track <- as.numeric(Nanog_statistics_total$Unique_track)
for(i in sort(unique(Nanog_statistics_total$Unique_track)))
{
  Nanog_statistics_total[which(Nanog_statistics_total$Unique_track==i),"Seconds_activation"] <- min(Nanog_statistics_total[!is.na(Nanog_statistics_total$S5P_Nucleus.number)&Nanog_statistics_total$Unique_track==i&!is.na(Nanog_statistics_total$Seconds),"Seconds"])
}

Nanog_statistics_total$Seconds_normalized <- Nanog_statistics_total$Seconds - Nanog_statistics_total$Seconds_activation 

#=========================================================
# Determine if the track are merging or non-merging and calculate merging time
#=========================================================

#Rules:
#Non-merging = only one cluster detected at each time point prior to trx activation
#Merging = two or more clusters detected for at least one time point to trx activation
#Several_clusters_short = merging clusters with detection of one or more clusters for only one time point in a row
#Several_clusters_long = merging clusters with detection of one or more clusters for at least two time points in a row
#Several_clusters_never_merging = two or more clusters are always detected between five time points until three points around trx activation
#Time_merging = Time of merging 
#Length double = For how long maximum can two or more Nanog clusters be detected in a row?
#Phenotype = merging or non-merging?



#Define a function that count the number of clusters per time point in the track
seqle <- function(x,incr=1) { 
  if(!is.numeric(x)) x <- as.numeric(x) 
  n <- length(x)  
  y <- x[-1L] != x[-n] + incr 
  i <- c(which(y|is.na(y)),n) 
  list(lengths = diff(c(0L,i)),
       values = x[head(c(0L,i)+1L,-1L)]) 
} 

#Define the variables
Nanog_statistics_total$Phenotype <- "NA" 
Nanog_statistics_total$Fused_activation <- "NA" 

#1.Create a table with the number of spots per time point for each track
for(i in sort(unique(Nanog_statistics_total$Unique_track)))
{
  #Determine the time of activation
  Time_activation<-unique(Nanog_statistics_total[which(Nanog_statistics_total$Unique_track==i),"Time_activation"])
  
  if(!Time_activation=="Inf"){
    
    #Determine the time of the track before transcription start
    Nanog_Unique_track_activation <- Nanog_statistics_total[which(Nanog_statistics_total$Unique_track==i&Nanog_statistics_total$Time %in% seq(Time_activation-10,Time_activation)),]
    
    
    Table <- as.data.frame(matrix(nrow = 2,ncol = length(unique(Nanog_Unique_track_activation$Time))))
    
    k=1
    
    for(j in min(Nanog_Unique_track_activation$Time):max(Nanog_Unique_track_activation$Time))
    {
      Table[1,k] <- j
      Table[2,k] <- nrow(Nanog_Unique_track_activation[Nanog_Unique_track_activation$Time==j,])
      k=k+1
    }
    
    #2. Determine if there is one or two spots at the moment of activation
    if(nrow(Nanog_statistics_total[which(Nanog_statistics_total$Unique_track==i&Nanog_statistics_total$Time==Time_activation),])>=2) {
      Nanog_statistics_total[Nanog_statistics_total$Unique_track==i&!is.na(Nanog_statistics_total$Unique_track==i),"Fused_activation"] <- "No"
      print(nrow(Nanog_statistics_total[which(Nanog_statistics_total$Unique_track==i&Nanog_statistics_total$Time==Time_activation),]))
      
      Nanog_Unique_track_activation <- Nanog_statistics_total[which(Nanog_statistics_total$Unique_track==i&Nanog_statistics_total$Time %in% seq(Time_activation-5,Time_activation+5)),]
      Table <- as.data.frame(matrix(nrow = 2,ncol = length(unique(Nanog_Unique_track_activation$Time))))
      
      k=1
      for(j in min(Nanog_Unique_track_activation$Time):max(Nanog_Unique_track_activation$Time))
      {
        Table[1,k] <- j
        Table[2,k] <- nrow(Nanog_Unique_track_activation[Nanog_Unique_track_activation$Time==j,])
        k=k+1
      }
      
      if(all(Table[2,]>=2))
      {Nanog_statistics_total[Nanog_statistics_total$Unique_track==i&!is.na(Nanog_statistics_total$Unique_track==i),"Phenotype"] <- "Several_clusters_never_merging"
      Nanog_statistics_total[Nanog_statistics_total$Unique_track==i&!is.na(Nanog_statistics_total$Unique_track==i),"Time_merging"] <- "None"
      Nanog_statistics_total[Nanog_statistics_total$Unique_track==i&!is.na(Nanog_statistics_total$Unique_track==i),"Length_double"] <- NA
      }
      if(any(Table[2,]<2))
      {
        Nanog_statistics_total[Nanog_statistics_total$Unique_track==i&!is.na(Nanog_statistics_total$Unique_track==i),"Time_merging"] <- "To determine"
        Nanog_statistics_total[Nanog_statistics_total$Unique_track==i&!is.na(Nanog_statistics_total$Unique_track==i),"Length_double"] <- max(seqle(Table[1,which(Table[2,]>=2)])$lengths)
        
        if(any(seqle(Table[1,which(Table[2,]>=2)])$lengths>=2))
        {Nanog_statistics_total[Nanog_statistics_total$Unique_track==i&!is.na(Nanog_statistics_total$Unique_track==i),"Phenotype"] <- "Several_clusters_long"} else 
        {Nanog_statistics_total[Nanog_statistics_total$Unique_track==i&!is.na(Nanog_statistics_total$Unique_track==i),"Phenotype"] <- "Several_clusters_short"}
      }
    }
    
    #3.Determine how many time points in a row have 2 or more spots
    if(!nrow(Nanog_statistics_total[which(Nanog_statistics_total$Unique_track==i&Nanog_statistics_total$Time==Time_activation),])>=2)
    {Nanog_statistics_total[Nanog_statistics_total$Unique_track==i&!is.na(Nanog_statistics_total$Unique_track==i),"Fused_activation"] <- "Yes"
    if(!any(Table[2,]>2))
    {
      Nanog_statistics_total[Nanog_statistics_total$Unique_track==i&!is.na(Nanog_statistics_total$Unique_track==i),"Phenotype"] <- "One_cluster"
      Nanog_statistics_total[Nanog_statistics_total$Unique_track==i&!is.na(Nanog_statistics_total$Unique_track==i),"Length_double"] <- "0"
      Nanog_statistics_total[Nanog_statistics_total$Unique_track==i&!is.na(Nanog_statistics_total$Unique_track==i),"Time_merging"] <- "None"}
    if(any(Table[2,]>=2)) {
      if(any(seqle(Table[1,which(Table[2,]>=2)])$lengths>=2)) {
        Nanog_statistics_total[Nanog_statistics_total$Unique_track==i&!is.na(Nanog_statistics_total$Unique_track==i),"Phenotype"] <- "Several_clusters_long"
        Nanog_statistics_total[Nanog_statistics_total$Unique_track==i&!is.na(Nanog_statistics_total$Unique_track==i),"Length_double"] <- max(seqle(Table[1,which(Table[2,]>=2)])$lengths)
        Nanog_statistics_total[Nanog_statistics_total$Unique_track==i&!is.na(Nanog_statistics_total$Unique_track==i),"Time_merging"] <- 
          seqle(Table[1,which(Table[2,]>=2)])$values[which(seqle(Table[1,which(Table[2,]>=2)])$lengths==max(seqle(Table[1,which(Table[2,]>=2)])$lengths))]+max(seqle(Table[1,which(Table[2,]>=2)])$lengths)} else {Nanog_statistics_total[Nanog_statistics_total$Unique_track==i&!is.na(Nanog_statistics_total$Unique_track==i),"Phenotype"] <- "Several_clusters_short"
          Nanog_statistics_total[Nanog_statistics_total$Unique_track==i&!is.na(Nanog_statistics_total$Unique_track==i),"Length_double"] <- max(seqle(Table[1,which(Table[2,]>=2)])$lengths)
          Nanog_statistics_total[Nanog_statistics_total$Unique_track==i&!is.na(Nanog_statistics_total$Unique_track==i),"Time_merging"] <- (max(seqle(Table[1,which(Table[2,]>=2)])$values)) + (max(seqle(Table[1,which(Table[2,]>=2)])$lengths))}
      
    }}}}


#4. Normalize the timing for Time_merging

#Determine the time of merging for specific cases
#Make the list of all the tracks ID for which the time of merging needs to be determine manually 
unique(Nanog_statistics_total[Nanog_statistics_total$Time_merging=="To determine","Unique_track"])

#10 12 18 24 33 40 87 89 92

#Display the number of spots for each time point for chosen track
Nanog_statistics_total[Nanog_statistics_total$Unique_track==69,c("Time","Time_activation")]

#Rules:
#-If there are two Nanog clusters at the moment of trx activation, time of merging is the time point when we detect only one spot
#the closest to trx activation in a range from -3 time points to 3 time points after trx activation, with priority of of time points before trx activation.
#If in the previous case, there are no time point when only one spot is detected, the track is considered as never merging and there are no time of merging.



Nanog_statistics_total[Nanog_statistics_total$Unique_track==10&!is.na(Nanog_statistics_total$Unique_track==10),"Time_merging"] <- 11
Nanog_statistics_total[Nanog_statistics_total$Unique_track==12&!is.na(Nanog_statistics_total$Unique_track==10),"Time_merging"] <- 13
Nanog_statistics_total[Nanog_statistics_total$Unique_track==18&!is.na(Nanog_statistics_total$Unique_track==18),"Time_merging"] <- NA
Nanog_statistics_total[Nanog_statistics_total$Unique_track==18&!is.na(Nanog_statistics_total$Unique_track==18),"Phenotype"] <- "Several_clusters_never_merging"
Nanog_statistics_total[Nanog_statistics_total$Unique_track==24&!is.na(Nanog_statistics_total$Unique_track==24),"Time_merging"] <- 10
Nanog_statistics_total[Nanog_statistics_total$Unique_track==33&!is.na(Nanog_statistics_total$Unique_track==33),"Time_merging"] <- NA
Nanog_statistics_total[Nanog_statistics_total$Unique_track==33&!is.na(Nanog_statistics_total$Unique_track==33),"Phenotype"] <- "Several_clusters_never_merging"
Nanog_statistics_total[Nanog_statistics_total$Unique_track==40&!is.na(Nanog_statistics_total$Unique_track==40),"Time_merging"] <- 9
Nanog_statistics_total[Nanog_statistics_total$Unique_track==87&!is.na(Nanog_statistics_total$Unique_track==87),"Time_merging"] <- 15
Nanog_statistics_total[Nanog_statistics_total$Unique_track==89&!is.na(Nanog_statistics_total$Unique_track==89),"Time_merging"] <- 15
Nanog_statistics_total[Nanog_statistics_total$Unique_track==92&!is.na(Nanog_statistics_total$Unique_track==92),"Time_merging"] <- 8

#Calculating the time of merging by subtracting the times
Nanog_statistics_total$Time_merging <- as.numeric(Nanog_statistics_total$Time_merging)
Nanog_statistics_total$Time_merging_normalized <- Nanog_statistics_total$Time - Nanog_statistics_total$Time_merging
Nanog_statistics_total$Unique_track <- as.numeric(Nanog_statistics_total$Unique_track)
Nanog_statistics_total$Seconds_merging_normalized <- Nanog_statistics_total$Seconds - Nanog_statistics_total$Time_merging*15 +15
Nanog_statistics_total$Difference_activation_merging <- Nanog_statistics_total$Time_merging - Nanog_statistics_total$Time_activation 


#=========================================================
# Normalise the fluorescence
#=========================================================

#Normalise of fluorescence for Nanog clusters
{
  #Only select Nanog clusters that are maximum 10 time points before transcription, because mitotic Nanog cluster are very bright and we don't want to normalise with this maximum.
  Nanog_statistics_total_10 <- Nanog_statistics_total[Nanog_statistics_total$Time_normalized %in% c(-10:10),]
  
  #Normalising without summing the fluorescence
  #Normalization by setting max to 1 and divide all the rest by it.
  {
    
    Nanog_statistics_total_10$Unique_track<-as.numeric(Nanog_statistics_total_10$Unique_track)
    
    #Create the variable for the dataframe
    Nanog_statistics_total_10$Intensity.Mean.y.normalized <- 0
    Nanog_statistics_total_10$Volume.normalized <- 0
    Nanog_statistics_total_10$Intensity.Sum.y <- as.numeric(Nanog_statistics_total_10$Intensity.Sum.y)
    Nanog_statistics_total_10$Intensity.Sum.x.normalized <- 0
    
    for(i in sort(unique(Nanog_statistics_total_10$Unique_track)))
    {
      Nanog_statistics_total_10[which(Nanog_statistics_total_10$Unique_track==i),"Intensity.Mean.y.normalized"] <- Nanog_statistics_total_10[which(Nanog_statistics_total_10$Unique_track==i),"Intensity.Mean.y"]/max(Nanog_statistics_total_10[which(Nanog_statistics_total_10$Unique_track==i),"Intensity.Mean.y"])
      Nanog_statistics_total_10[which(Nanog_statistics_total_10$Unique_track==i),"Volume.normalized"] <- Nanog_statistics_total_10[which(Nanog_statistics_total_10$Unique_track==i),"Volume"]/max(Nanog_statistics_total_10[which(Nanog_statistics_total_10$Unique_track==i),"Volume"])
      Nanog_statistics_total_10[which(Nanog_statistics_total_10$Unique_track==i),"Intensity.Sum.y.normalized"] <- Nanog_statistics_total_10[which(Nanog_statistics_total_10$Unique_track==i),"Intensity.Sum.y"]/max(Nanog_statistics_total_10[which(Nanog_statistics_total_10$Unique_track==i),"Intensity.Sum.y"])
      Nanog_statistics_total_10[which(Nanog_statistics_total_10$Unique_track==i),"Intensity.Sum.x.normalized"] <- Nanog_statistics_total_10[which(Nanog_statistics_total_10$Unique_track==i),"Intensity.Sum.x"]/max(Nanog_statistics_total_10[which(Nanog_statistics_total_10$Unique_track==i),"Intensity.Sum.x"])
    }
    
  }
}

#S5P foci: one to one normalization
{
  Nanog_statistics_total_10$S5P_Intensity.Mean.x.normalized <- 0
  Nanog_statistics_total_10$S5P_Volume.normalized <- 0
  Nanog_statistics_total_10$S5P_Intensity.Sum.y.normalized <- 0
  
  
  for(i in sort(unique(Nanog_statistics_total_10$Unique_track)))
  {
    Nanog_statistics_total_10[which(Nanog_statistics_total_10$Unique_track==i),"S5P_Intensity.Mean.x.normalized"] <- Nanog_statistics_total_10[which(Nanog_statistics_total_10$Unique_track==i),"S5P_Intensity.Mean.y"]/min(Nanog_statistics_total_10[which(Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$S5P_Intensity.Mean.y!="NA"),"S5P_Intensity.Mean.y"])
    Nanog_statistics_total_10[which(Nanog_statistics_total_10$Unique_track==i),"S5P_Volume.normalized"] <- Nanog_statistics_total_10[which(Nanog_statistics_total_10$Unique_track==i),"S5P_Volume"]/min(Nanog_statistics_total_10[which(Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$S5P_Volume!="NA"),"S5P_Volume"])
    Nanog_statistics_total_10[which(Nanog_statistics_total_10$Unique_track==i),"S5P_Intensity.Sum.y.normalized"] <- Nanog_statistics_total_10[which(Nanog_statistics_total_10$Unique_track==i),"S5P_Intensity.Sum.y"]/min(Nanog_statistics_total_10[which(Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$S5P_Intensity.Sum.y!="NA"),"S5P_Intensity.Sum.y"])
  }
  
}  

#SUMMING ALL FOCI for volume and total intensity, make the mean for mean intensity
{
  #Create a new dataframe with one line per nucleus, track and time point with sum of volume or total intensity if there are more than two clusters per time point, and 
  {
    
    Sum_parameters_plotting <- data.frame()
    k=0  
    
    for(i in sort(unique(Nanog_statistics_total_10$Unique_track)))
    {
      
      if(nrow(Nanog_statistics_total_10[which(Nanog_statistics_total_10$Unique_track==i),])!=0)
      {
        for(e in 1:max(na.omit(Nanog_statistics_total_10[Nanog_statistics_total_10$Unique_track==i,'Time'])))
        {
          
          if(nrow(Nanog_statistics_total_10[Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==e,1:10])!=0)
          {
            
            k=k+1
            
            Sum_parameters_plotting[k,1:21] <- Nanog_statistics_total_10[Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==e,c("Date","Stage","Stack","Nucleus","Nucleus.number","Time","TrackID.x","Unique_track","Time_activation","Intensity.Mean.x","Intensity.Sum.x","Seconds","Time_normalized","Seconds_activation","Seconds_normalized","Inhibition","Phenotype","Seconds_merging_normalized","Length_double","TrackID.y")][1,]
            Sum_parameters_plotting[k,22] <- sum(Nanog_statistics_total_10[Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==e,"Intensity.Mean.y"])/length(Nanog_statistics_total_10[Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==e,"Intensity.Mean.y"])
            Sum_parameters_plotting[k,23] <- sum(Nanog_statistics_total_10[Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==e,"Volume"])
            Sum_parameters_plotting[k,24] <- sum(Nanog_statistics_total_10[Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==e,"Intensity.Sum.y"])
            
            
            if(length(Nanog_statistics_total_10[Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==e,"Position.X"])==1)
            {
              
              Sum_parameters_plotting[k,25] <- 0
              
            }else if(length(Nanog_statistics_total_10[Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==e,"Position.X"])==2) {
              Sum_parameters_plotting[k,25] <- sqrt((Nanog_statistics_total_10[Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==e,"Position.X"][1] - Nanog_statistics_total_10[Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==e,"Position.X"][2])^2 + (Nanog_statistics_total_10[Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==e,"Position.Y"][1] - Nanog_statistics_total_10[Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==e,"Position.Y"][2])^2 + (Nanog_statistics_total_10[Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==e,"Position.Z"][1] - Nanog_statistics_total_10[Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==e,"Position.Z"][2])^2)
              
            } else if(length(Nanog_statistics_total_10[Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==e,"Position.X"])==3)
            {
              
              Distance <- array()
              
              Distance[1] <- sqrt((Nanog_statistics_total_10[Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==e,"Position.X"][1] - Nanog_statistics_total_10[Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==e,"Position.X"][2])^2 + (Nanog_statistics_total_10[Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==e,"Position.Y"][1] - Nanog_statistics_total_10[Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==e,"Position.Y"][2])^2 + (Nanog_statistics_total_10[Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==e,"Position.Z"][1] - Nanog_statistics_total_10[Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==e,"Position.Z"][2])^2)
              
              Distance[2] <-   sqrt((Nanog_statistics_total_10[Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==e,"Position.X"][1] - Nanog_statistics_total_10[Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==e,"Position.X"][3])^2 + (Nanog_statistics_total_10[Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==e,"Position.Y"][1] - Nanog_statistics_total_10[Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==e,"Position.Y"][3])^2 + (Nanog_statistics_total_10[Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==e,"Position.Z"][1] - Nanog_statistics_total_10[Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==e,"Position.Z"][3])^2)
              
              Distance[3]  <- sqrt((Nanog_statistics_total_10[Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==e,"Position.X"][2] - Nanog_statistics_total_10[Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==e,"Position.X"][3])^2 + (Nanog_statistics_total_10[Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==e,"Position.Y"][2] - Nanog_statistics_total_10[Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==e,"Position.Y"][3])^2 + (Nanog_statistics_total_10[Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==e,"Position.Z"][2] - Nanog_statistics_total_10[Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==e,"Position.Z"][3])^2)
              
              Sum_parameters_plotting[k,25] <- max(Distance)
              
              
            } else if(length(Nanog_statistics_total_10[Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==e,"Position.X"])==4)
            {
              
              Distance[1] <- sqrt((Nanog_statistics_total_10[Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==e,"Position.X"][1] - Nanog_statistics_total_10[Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==e,"Position.X"][2])^2 + (Nanog_statistics_total_10[Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==e,"Position.Y"][1] - Nanog_statistics_total_10[Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==e,"Position.Y"][2])^2 + (Nanog_statistics_total_10[Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==e,"Position.Z"][1] - Nanog_statistics_total_10[Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==e,"Position.Z"][2])^2)
              
              Distance[2] <-   sqrt((Nanog_statistics_total_10[Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==e,"Position.X"][1] - Nanog_statistics_total_10[Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==e,"Position.X"][3])^2 + (Nanog_statistics_total_10[Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==e,"Position.Y"][1] - Nanog_statistics_total_10[Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==e,"Position.Y"][3])^2 + (Nanog_statistics_total_10[Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==e,"Position.Z"][1] - Nanog_statistics_total_10[Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==e,"Position.Z"][3])^2)
              
              Distance[3] <-   sqrt((Nanog_statistics_total_10[Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==e,"Position.X"][1] - Nanog_statistics_total_10[Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==e,"Position.X"][4])^2 + (Nanog_statistics_total_10[Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==e,"Position.Y"][1] - Nanog_statistics_total_10[Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==e,"Position.Y"][4])^2 + (Nanog_statistics_total_10[Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==e,"Position.Z"][1] - Nanog_statistics_total_10[Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==e,"Position.Z"][4])^2)
              
              Distance[4]  <- sqrt((Nanog_statistics_total_10[Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==e,"Position.X"][2] - Nanog_statistics_total_10[Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==e,"Position.X"][3])^2 + (Nanog_statistics_total_10[Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==e,"Position.Y"][2] - Nanog_statistics_total_10[Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==e,"Position.Y"][3])^2 + (Nanog_statistics_total_10[Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==e,"Position.Z"][2] - Nanog_statistics_total_10[Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==e,"Position.Z"][3])^2)
              
              Distance[5]  <- sqrt((Nanog_statistics_total_10[Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==e,"Position.X"][2] - Nanog_statistics_total_10[Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==e,"Position.X"][4])^2 + (Nanog_statistics_total_10[Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==e,"Position.Y"][2] - Nanog_statistics_total_10[Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==e,"Position.Y"][4])^2 + (Nanog_statistics_total_10[Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==e,"Position.Z"][2] - Nanog_statistics_total_10[Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==e,"Position.Z"][4])^2)
              
              Distance[6]  <- sqrt((Nanog_statistics_total_10[Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==e,"Position.X"][3] - Nanog_statistics_total_10[Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==e,"Position.X"][4])^2 + (Nanog_statistics_total_10[Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==e,"Position.Y"][3] - Nanog_statistics_total_10[Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==e,"Position.Y"][4])^2 + (Nanog_statistics_total_10[Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==e,"Position.Z"][3] - Nanog_statistics_total_10[Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==e,"Position.Z"][4])^2)
              
              
              Sum_parameters_plotting[k,25] <- max(Distance)
              
            }
          }
        }
      }
    }
  }
  
  #Naming the columns together
  colnames(Sum_parameters_plotting) <- c("Date","Stage","Stack","Nucleus","Nucleus.number","Time","TrackID","Unique_track","Time_activation","Intensity.Mean.x","Intensity.Sum.x","Seconds","Time_normalized","Seconds_activation","Seconds_normalized","Inhibition","Phenotype","Seconds_merging_normalized","Length_double","TrackID.y","Date","Mean_Intensity","Volume","Total_intensity","Distance")
  
  #Normalize the sum fluorescence to the max fluorescence 
  Sum_parameters_plotting$Mean_intensity.normalized <- 0
  Sum_parameters_plotting$Volume.normalized <- 0
  Sum_parameters_plotting$Total_intensity.normalized <- 0
  
  
  for(i in sort(unique(Sum_parameters_plotting$Unique_track)))
  {
    if(nrow(Sum_parameters_plotting[which(Sum_parameters_plotting$Unique_track==i),])!=0)
    {
      Sum_parameters_plotting[which(Sum_parameters_plotting$Unique_track==i&!is.na(Sum_parameters_plotting$Mean_Intensity)),"Mean_intensity.normalized"] <- Sum_parameters_plotting[which(Sum_parameters_plotting$Unique_track==i&!is.na(Sum_parameters_plotting$Mean_Intensity)),"Mean_Intensity"]/max(Sum_parameters_plotting[which(Sum_parameters_plotting$Unique_track==i&!is.na(Sum_parameters_plotting$Mean_Intensity)),"Mean_Intensity"])
      Sum_parameters_plotting[which(Sum_parameters_plotting$Unique_track==i&!is.na(Sum_parameters_plotting$Unique_track==i)),"Volume.normalized"] <- Sum_parameters_plotting[which(Sum_parameters_plotting$Unique_track==i&!is.na(Sum_parameters_plotting$Unique_track==i)),"Volume"]/max(Sum_parameters_plotting[which(Sum_parameters_plotting$Unique_track==i&!is.na(Sum_parameters_plotting$Unique_track==i)),"Volume"])
      Sum_parameters_plotting[which(Sum_parameters_plotting$Unique_track==i&!is.na(Sum_parameters_plotting$Unique_track==i)),"Total_intensity.normalized"] <- Sum_parameters_plotting[which(Sum_parameters_plotting$Unique_track==i&!is.na(Sum_parameters_plotting$Unique_track==i)),"Total_intensity"]/max(Sum_parameters_plotting[which(Sum_parameters_plotting$Unique_track==i&!is.na(Sum_parameters_plotting$Unique_track==i)),"Total_intensity"])
    }
  }
}

#Normalize to find the time when the maximum of fluorescence is reached for the different parameter
{
  Sum_parameters_plotting$Time_activation <- round(as.numeric(as.character(Sum_parameters_plotting$Time_activation)))
  
  for(i in  sort(unique(Sum_parameters_plotting$Unique_track)))
  {
    
    Sum_parameters_plotting[which(Sum_parameters_plotting$Unique_track==i),"Time_max_mean_intensity"] <- Sum_parameters_plotting[Sum_parameters_plotting$Unique_track==i&Sum_parameters_plotting$Mean_intensity.normalized==1,"Time"][1]
    
    Sum_parameters_plotting[which(Sum_parameters_plotting$Unique_track==i),"Time_max_volume"] <- Sum_parameters_plotting[Sum_parameters_plotting$Unique_track==i&Sum_parameters_plotting$Volume.normalized==1,"Time"][1]
    
    Sum_parameters_plotting[which(Sum_parameters_plotting$Unique_track==i),"Time_max_total_intensity"] <- Sum_parameters_plotting[Sum_parameters_plotting$Unique_track==i&Sum_parameters_plotting$Total_intensity.normalized==1,"Time"][1]
    
  }
  
  Sum_parameters_plotting$Time_difference_mean_intensity <- Sum_parameters_plotting$Time_max_mean_intensity - Sum_parameters_plotting$Time_activation
  Sum_parameters_plotting$Time_difference_volume <- Sum_parameters_plotting$Time_max_volume - Sum_parameters_plotting$Time_activation
  Sum_parameters_plotting$Time_difference_total_intensity <- Sum_parameters_plotting$Time_max_total_intensity - Sum_parameters_plotting$Time_activation
}


#Associate all the clusters from the same phenotype together
Nanog_statistics_total_10$Phenotype_simple <- "One_cluster"
Nanog_statistics_total_10[Nanog_statistics_total_10$Phenotype %in% c("Several_clusters_long","Several_clusters_short","Trx_before_merging","Several_clusters_never_merging"),"Phenotype_simple"] <- "More_cluster"

Sum_parameters_plotting$Phenotype_simple <- "One_cluster"
Sum_parameters_plotting[Sum_parameters_plotting$Phenotype %in% c("Several_clusters_long","Several_clusters_short","Trx_before_merging","Several_clusters_never_merging"),"Phenotype_simple"] <- "More_cluster"







#=========================================================
# Plots
#=========================================================

#For any plots, run these lines to round the time to the closest second  
Nanog_statistics_total_10$Seconds_normalized <- as.numeric(Nanog_statistics_total_10$Seconds_normalized)
Nanog_statistics_total_10$Seconds_normalized <- round(as.numeric(as.character(Nanog_statistics_total_10$Seconds_normalized)))

Sum_parameters_plotting$Seconds_normalized <- as.numeric(Sum_parameters_plotting$Seconds_normalized)
Sum_parameters_plotting$Seconds_normalized <- round(as.numeric(as.character(Sum_parameters_plotting$Seconds_normalized)))

Nanog_statistics_total_10$Seconds_merging_normalized <-round(as.numeric(as.character(Nanog_statistics_total_10$Seconds_merging_normalized)))
Sum_parameters_plotting$Seconds_merging_normalized <-round(as.numeric(as.character(Sum_parameters_plotting$Seconds_merging_normalized)))

colnames(Sum_parameters_plotting) <- make.unique(names(Sum_parameters_plotting))

############Fig 1d
{
  #Color_scheme <- unique(miR430DNA[,c("TrackID_unique","Time","Color")])
  #Color <- Color_scheme[which(Color_scheme$TrackID_unique %in% c(2)),][3]
  
  #Make the tracks align on one place by adding the number of Nanog clusters per time points and unique color, plus position number
  {
    
    Nanog_statistics_total_10 <- Nanog_statistics_total_10[with(Nanog_statistics_total_10,order(Unique_track,Time,TrackID.x)),]
    Nanog_statistics_total_10$Track_plotting <-0
    
    for(i in unique(as.numeric(Nanog_statistics_total_10$Unique_track)))
    {
      for(g in min(unique(Nanog_statistics_total_10[which(Nanog_statistics_total_10$Unique_track==i),"Time"])):max(unique(Nanog_statistics_total_10[which(Nanog_statistics_total_10$Unique_track==i),"Time"])))
      {
        print(nrow(Nanog_statistics_total_10[which(Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==g),]))
        
        if(nrow(Nanog_statistics_total_10[which(Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==g),])==0)
        {print("0")}
        
        else if(nrow(Nanog_statistics_total_10[which(Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==g),])==1)
        {Nanog_statistics_total_10[which(Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==g),"Track_plotting"][1] <- 0.5
        Nanog_statistics_total_10[which(Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==g),"Spot_number"][1] <- 1
        Nanog_statistics_total_10[which(Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==g),"Color"][1] <- "a"
        }
        
        else if(nrow(Nanog_statistics_total_10[which(Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==g),])==2)
        {Nanog_statistics_total_10[which(Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==g),"Track_plotting"][1] <- 0.25
        Nanog_statistics_total_10[which(Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==g),"Track_plotting"][2] <- 0.75
        Nanog_statistics_total_10[which(Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==g),"Spot_number"][1] <- 2
        Nanog_statistics_total_10[which(Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==g),"Spot_number"][2] <- 2
        Nanog_statistics_total_10[which(Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==g),"Color"][1] <- "b"
        Nanog_statistics_total_10[which(Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==g),"Color"][2] <- "b"
        
        }
        else if(nrow(Nanog_statistics_total_10[which(Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==g),])==3)
        {Nanog_statistics_total_10[which(Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==g),"Track_plotting"][1] <- 0
        Nanog_statistics_total_10[which(Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==g),"Track_plotting"][2] <- 0.5
        Nanog_statistics_total_10[which(Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==g),"Track_plotting"][3] <- 1
        Nanog_statistics_total_10[which(Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==g),"Spot_number"][1] <- 3
        Nanog_statistics_total_10[which(Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==g),"Spot_number"][2] <- 3
        Nanog_statistics_total_10[which(Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==g),"Spot_number"][3] <- 3
        Nanog_statistics_total_10[which(Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==g),"Color"][1] <- "c"
        Nanog_statistics_total_10[which(Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==g),"Color"][2] <- "c"
        Nanog_statistics_total_10[which(Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==g),"Color"][3] <- "c"
        
        
        
        
        }
        else if(nrow(Nanog_statistics_total_10[which(Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==g),])==4)
        {Nanog_statistics_total_10[which(Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==g),"Track_plotting"][1] <- 0
        Nanog_statistics_total_10[which(Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==g),"Track_plotting"][2] <- 0.5
        Nanog_statistics_total_10[which(Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==g),"Track_plotting"][3] <- 1
        Nanog_statistics_total_10[which(Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==g),"Track_plotting"][4] <- 1.5
        Nanog_statistics_total_10[which(Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==g),"Spot_number"][1] <- 4
        Nanog_statistics_total_10[which(Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==g),"Spot_number"][2] <- 4
        Nanog_statistics_total_10[which(Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==g),"Spot_number"][3] <- 4
        Nanog_statistics_total_10[which(Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==g),"Spot_number"][4] <- 4
        Nanog_statistics_total_10[which(Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==g),"Color"][1] <- "d"
        Nanog_statistics_total_10[which(Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==g),"Color"][2] <- "d"
        Nanog_statistics_total_10[which(Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==g),"Color"][3] <- "d"
        Nanog_statistics_total_10[which(Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==g),"Color"][4] <- "d"
        }
        else if(nrow(Nanog_statistics_total_10[which(Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==g),])==5)
        {Nanog_statistics_total_10[which(Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==g),"Track_plotting"][1] <- -0.5
        Nanog_statistics_total_10[which(Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==g),"Track_plotting"][2] <- 0
        Nanog_statistics_total_10[which(Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==g),"Track_plotting"][3] <- 0.5
        Nanog_statistics_total_10[which(Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==g),"Track_plotting"][4] <- 1
        Nanog_statistics_total_10[which(Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==g),"Track_plotting"][5] <- 1.5
        Nanog_statistics_total_10[which(Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==g),"Spot_number"][1] <- 5
        Nanog_statistics_total_10[which(Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==g),"Spot_number"][2] <- 5
        Nanog_statistics_total_10[which(Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==g),"Spot_number"][3] <- 5
        Nanog_statistics_total_10[which(Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==g),"Spot_number"][4] <- 5
        Nanog_statistics_total_10[which(Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==g),"Spot_number"][5] <- 5
        Nanog_statistics_total_10[which(Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==g),"Color"][1] <- "e"
        Nanog_statistics_total_10[which(Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==g),"Color"][2] <- "e"
        Nanog_statistics_total_10[which(Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==g),"Color"][3] <- "e"
        Nanog_statistics_total_10[which(Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==g),"Color"][4] <- "e"
        Nanog_statistics_total_10[which(Nanog_statistics_total_10$Unique_track==i&Nanog_statistics_total_10$Time==g),"Color"][5] <- "e"
        }
      }
    }
  }
  
  #Calculate the number of Nanog clusters at the moment of activation and the total number of clusters per track
  {
  Nanog_statistics_total_10$Unique_track <- as.numeric(Nanog_statistics_total_10$Unique_track)
  
  #Number of spots at the moment of activation
  for(i in 1:max(Nanog_statistics_total_10$Unique_track))
  {Nanog_statistics_total_10[which(Nanog_statistics_total_10$Unique_track==i),"Number_spot_activation"] <- Nanog_statistics_total_10[which(Nanog_statistics_total_10$Time_activation==Nanog_statistics_total_10$Time&Nanog_statistics_total_10$Unique_track==i),"Spot_number"][1]}
  
  for(i in 1:max(Nanog_statistics_total_10$Unique_track))
  {Nanog_statistics_total_10[which(Nanog_statistics_total_10$Unique_track == i),"Plot_spot_number"] <- sum(Nanog_statistics_total_10[Nanog_statistics_total_10$Unique_track==i,"Spot_number"])}
  
  Nanog_statistics_total_10$Unique_track <- factor(Nanog_statistics_total_10$Unique_track, levels =  unique(Nanog_statistics_total_10[order(Nanog_statistics_total_10$Number_spot_activation,Nanog_statistics_total_10$Plot_spot_number),"Unique_track"]))
  
  Nanog_statistics_total_10 <- Nanog_statistics_total_10[!is.na(Nanog_statistics_total_10$TrackID.x),]
  Nanog_statistics_total_10$group <- rleid(Nanog_statistics_total_10$Color)
  df_plot <- head(do.call(rbind, by(Nanog_statistics_total_10, Nanog_statistics_total_10$group, rbind, NA)), -1)
  df_plot[,c("group","Color")] <- lapply(df_plot[,c("group","Color")], na.locf)
  df_plot[] <- lapply(df_plot, na.locf, fromLast = TRUE)
  }
  
  #Order the tracks based on specific criterias
  {
    for(i in df_plot$Unique_track)
    {
      if(all(na.omit(df_plot[which(df_plot$Unique_track==i),"Phenotype_simple"]=="One_cluster")))
      {df_plot[which(df_plot$Unique_track==i),"Simple_phenotype"] <- 1
      print("True")
      } else {df_plot[which(df_plot$Unique_track==i),"Simple_phenotype"] <- 2
      print("False")
      }
    }
    
    #Get only the data for 1k and with no inhibition
    df_plot_1k_no_ihn <- df_plot[which(df_plot$Stage=="1k"&df_plot$Inhibition=="No"),]
    
    #Give a new unique track number to the unique track
    k=1
    for(i in na.omit(unique(df_plot_1k_no_ihn$Unique_track)))
    {
      df_plot_1k_no_ihn[which(df_plot_1k_no_ihn$Unique_track==i),"Unique_track_2"] <- k
      k=k+1
      print(k)
    }
    
    #Re-order the tracks in this order:
    #1. Clusters with more than two clusters at transcription start
    #2. Merging and then non-merging
    #3. Sum of the number of clusters detected in the tracks
    {
    df_plot_1k_no_ihn$Unique_track_2 <- as.numeric(as.character(df_plot_1k_no_ihn$Unique_track_2))     
    df_plot_1k_no_ihn$Simple_phenotype <- as.numeric(as.character(df_plot_1k_no_ihn$Simple_phenotype)) 
    df_plot_1k_no_ihn$Number_spot_activation <- as.numeric(as.character(df_plot_1k_no_ihn$Number_spot_activation))   
    df_plot_1k_no_ihn$Plot_spot_number <- as.numeric(as.character(df_plot_1k_no_ihn$Plot_spot_number))   
    
    df_plot_1k_no_ihn$Unique_track_2 <- factor(df_plot_1k_no_ihn$Unique_track_2, 
                                               levels =  unique(df_plot_1k_no_ihn[order(df_plot_1k_no_ihn$Simple_phenotype,
                                                                                        df_plot_1k_no_ihn$Number_spot_activation,
                                                                                        df_plot_1k_no_ihn$Plot_spot_number, 
                                                                                        decreasing=FALSE),"Unique_track_2"]))
    df_plot_1k_no_ihn$Unique_track_2 <- as.numeric(df_plot_1k_no_ihn$Unique_track_2)
    }
    
    #Add the number of spots equal to zero if no detection of nanog clusters at certain time points
    {
    df_plot_1k_no_ihn
    df_plot_1k_no_ihn$Unique_track_2 <- as.numeric(as.character(df_plot_1k_no_ihn$Unique_track_2))
    
    for(i in 1:max(unique(na.omit(df_plot_1k_no_ihn$Unique_track_2))))
    {
      Real <- df_plot_1k_no_ihn[df_plot_1k_no_ihn$Unique_track_2==i,"Time"]
      Ideal <- min(na.omit(df_plot_1k_no_ihn[df_plot_1k_no_ihn$Unique_track_2==i,"Time"])):max(na.omit(df_plot_1k_no_ihn[df_plot_1k_no_ihn$Unique_track_2==i,"Time"]))
      Ideal <- as.numeric(Ideal)
      
      Number <- Ideal[is.na(match(Ideal,Real))]
      
      for(y in Number)
      {
        Insert <- data.frame(matrix(ncol=ncol(df_plot_1k_no_ihn)))
        colnames(Insert) <- colnames(df_plot_1k_no_ihn)
        
        Insert[,c(1:4)] <- unique(df_plot_1k_no_ihn[df_plot_1k_no_ihn$Unique_track_2==i,c(1:4)])[1,]
        Insert$Unique_track_2 <- i
        Insert$Time <- y
        Insert$Spot_number <- 0
        Insert$Color <- "f"
        Insert$Seconds_normalized <- (y-df_plot_1k_no_ihn[df_plot_1k_no_ihn$Unique_track_2==i,"Time_activation"][1])*15
        df_plot_1k_no_ihn <- rbind(df_plot_1k_no_ihn,Insert)
      }
    }             
  }  
  }
  
  #Code for plotting Fig1d
  {
    ggplot() + geom_point(df_plot_1k_no_ihn,mapping=aes(y=Unique_track_2 , 
                                                        x=Seconds_normalized, 
                                                        fill=as.factor(Color)),
                          size=4,
                          shape=22, colour = "transparent")  +  
      #scale_fill_manual(values=c("grey80","grey50","grey30","grey30","grey30","grey95")) +
      scale_fill_manual(values=c("grey85","lightsteelblue3","lightsteelblue4","lightsteelblue4","grey95")) +
      labs(x= "Time Relative to Transcription (s)",y="Tracks") +
      theme_bw() +
      theme(legend.title=element_blank(), 
            panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank(),
            panel.grid.minor.x = element_blank(),
            #axis.text.y = element_blank(),
            #axis.ticks.y = element_blank(),
            aspect.ratio=10/3) + 
      #geom_rect(aes(xmin = -8, ymin = 0, xmax = 8, ymax = 91), color="red", alpha=0,size=0.3) + 
      geom_vline(aes(xintercept = 0,), size=0.5,color="red", linetype=2)+ 
      scale_x_continuous(breaks=seq(-270,270,30))+
      scale_y_continuous(breaks=seq(0,90,5))
  }
  
}

############Fig 1e
{
  
  #Take the following examples = merging
  {
  Plot_track_merging <- Nanog_statistics_total_10[which(Nanog_statistics_total_10$Unique_track %in% c("20","57","9")),]
  Plot_track_merging[which(Plot_track_merging$Unique_track=="57"),"Track_plotting"] <-  Plot_track_merging[which(Plot_track_merging$Unique_track=="57"),"Track_plotting"]+2
  Plot_track_merging[which(Plot_track_merging$Unique_track=="9"),"Track_plotting"] <-  Plot_track_merging[which(Plot_track_merging$Unique_track=="9"),"Track_plotting"]+4
  
  pdf(file = '/Users/nchabot/Desktop/Fig1c_mergings.pdf',height=4, width=8)
  ggplot() + 
    geom_vline(aes(xintercept = 0,), size=0.5,color="red", linetype=2)+
    geom_point(Plot_track_merging, 
               mapping = aes(y=Track_plotting, 
                             x=Seconds_normalized,
                             colour=as.numeric(Intensity.Sum.y.normalized)),
               size=5, shape=16) +
    scale_color_gradient(low="palegreen1", high="forestgreen") +
    theme_bw() + labs(x= "Time (s)")+ 
    ylim(0,5)+scale_x_continuous(breaks = seq(-90, 90, by = 15),limits = c(-90,90)) +
    theme(aspect.ratio=1, axis.text = element_text(size=10),axis.title = element_text(size=10),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          axis.ticks.y=element_blank(),axis.title.y=element_blank())
  dev.off()
  }
  #Take the following examples = non-merging
  {
  Plot_track_merging <- Nanog_statistics_total_10[which(Nanog_statistics_total_10$Unique_track %in% c("9","15","38")),]
  Plot_track_merging[which(Plot_track_merging$Unique_track=="15"),"Track_plotting"] <-  Plot_track_merging[which(Plot_track_merging$Unique_track=="15"),"Track_plotting"]+1
  Plot_track_merging[which(Plot_track_merging$Unique_track=="38"),"Track_plotting"] <-  Plot_track_merging[which(Plot_track_merging$Unique_track=="38"),"Track_plotting"]+2
  
  pdf(file = '/Users/nchabot/Desktop/Fig1c_non_merging.pdf',height=4, width=8)
  ggplot() + 
    geom_vline(aes(xintercept = 0,), size=0.5,color="red", linetype=2)+
    geom_point(Plot_track_merging, 
               mapping = aes(y=Track_plotting, 
                             x=Seconds_normalized,
                             colour=as.numeric(Intensity.Sum.y.normalized)),
               size=5, shape=16) +
    scale_color_gradient(low="palegreen1", high="forestgreen", limit=c(0.4,1)) +
    theme_bw() + labs(x= "Time (s)")+ 
    ylim(0,5)+scale_x_continuous(breaks = seq(-90, 90, by = 15),limits = c(-90,90)) +
    theme(aspect.ratio=1, axis.text = element_text(size=10),axis.title = element_text(size=10),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          axis.ticks.y=element_blank(),
          axis.title.y=element_blank())
  
  dev.off()
}
  
  
}

############Fig 1f = not the same plots!!
{
  
#No inhibition
  {
    #Calculate the values for the ribbons
    {
  Sum_parameters_plotting$Seconds_normalized <- as.numeric(Sum_parameters_plotting$Seconds_normalized)
  DT <- as.data.table(Sum_parameters_plotting[which(Sum_parameters_plotting$Stage=="1k"&Sum_parameters_plotting$Inhibition=="No"),])
  Summary <- DT[ , list(MEAN_DV = median(na.omit(Total_intensity.normalized)), SD_DV = sd(na.omit(Total_intensity.normalized)), Q10= quantile(Total_intensity.normalized, prob=c(0.10), na.rm=TRUE),Q25 = quantile(Total_intensity.normalized, prob=c(0.25),na.rm=TRUE),Q75 = quantile(Total_intensity.normalized, prob=c(0.75),na.rm=TRUE), Q90= quantile(Total_intensity.normalized, prob=c(0.90),na.rm=TRUE),N = .N),  by = c("Seconds_normalized")]
  Summary<- Summary[order(Summary$Seconds_normalized),]
  
  Median_circ <- Summary$MEAN_DV
  sd_circ <- Summary$SD_DV
  
  CI_infcirc <- Summary$Q25
  CI_supcirc <- Summary$Q75
  
  Circ_plot <- cbind(Median_circ,sd_circ,CI_infcirc,CI_supcirc)
  Circ_plot <-Circ_plot[complete.cases(Circ_plot[,c("CI_supcirc")]),]
  Circ_plot <- cbind(Circ_plot,unique(Summary[,c("Seconds_normalized")]))
  Circ_plot <- data.frame(Circ_plot)
  
  colnames(Circ_plot) <- c("Median","sd","CI_inf","CI_sup","Time")
  
  a <- Circ_plot %>%
    group_by(Time) %>%
    summarise_at(vars("Median","sd","CI_inf","CI_sup"), list(name = mean))
    }
    #Code for the plot
    {
  pdf(file = '/Users/nchabot/Desktop/Figure1f_wt',height=4, width=4)
  ggplot() + geom_jitter(Sum_parameters_plotting[Sum_parameters_plotting$Inhibition=="No"&Sum_parameters_plotting$Stage=="1k",], 
                         mapping = aes(y=as.numeric(as.character(Total_intensity.normalized)), 
                                       x=as.numeric(as.character(Seconds_normalized))), 
                         width=5,alpha=0.1, size = 2,show.legend = FALSE, color="black")+
    geom_line(a,mapping=aes(x=Time,y=Median_name), size=0.75,color="black")+
    geom_ribbon(a,mapping=aes(x=Time,y=Median_name,
                              ymin=CI_inf_name,ymax=CI_sup_name),
                alpha=0.15, fill="black")+
    scale_linetype_manual(values=c("dashed","solid"))+
    geom_vline(aes(xintercept = 0,), size=0.5,color="red", linetype=2)+ 
    theme_bw() + labs(x= "Time centered on transcription start (s)",y="Normalized Total Intensity")+ 
    ylim(-0.01,1.01)+scale_x_continuous(breaks = seq(-90, 150, by = 30),limits = c(-95,155)) +
    theme(aspect.ratio=1, axis.text = element_text(size=14),axis.title = element_text(size=14),
          panel.grid.minor.x = element_blank())
  dev.off()
    }
  }
#Inhibition 
  {
    #Calculate the values for the ribbons
    {
    DT <- as.data.table(Sum_parameters_plotting[which(Sum_parameters_plotting$Stage=="1k"&Sum_parameters_plotting$Inhibition=="Yes"),])
    Summary <- DT[ , list(MEAN_DV = median(na.omit(Total_intensity.normalized)), SD_DV = sd(na.omit(Total_intensity.normalized)), Q10= quantile(Total_intensity.normalized, prob=c(0.10), na.rm=TRUE),Q25 = quantile(Total_intensity.normalized, prob=c(0.25),na.rm=TRUE),Q75 = quantile(Total_intensity.normalized, prob=c(0.75),na.rm=TRUE), Q90= quantile(Total_intensity.normalized, prob=c(0.90),na.rm=TRUE),N = .N),  by = c("Seconds_normalized")]
    Summary<- Summary[order(Summary$Seconds_normalized),]
    
    Median_circ <- Summary$MEAN_DV
    sd_circ <- Summary$SD_DV
    
    CI_infcirc <- Summary$Q25
    CI_supcirc <- Summary$Q75
    
    Circ_plot <- cbind(Median_circ,sd_circ,CI_infcirc,CI_supcirc)
    Circ_plot <-Circ_plot[complete.cases(Circ_plot[,c("CI_supcirc")]),]
    Circ_plot <- cbind(Circ_plot,unique(Summary[,c("Seconds_normalized")]))
    Circ_plot <- data.frame(Circ_plot)
    
    colnames(Circ_plot) <- c("Median","sd","CI_inf","CI_sup","Time")
    
    a <- Circ_plot %>%
      group_by(Time) %>%
      summarise_at(vars("Median","sd","CI_inf","CI_sup"), list(name = mean))
    }
    #Code for the plot
    {
    pdf(file = '/Users/nchabot/Desktop/Figure_1f_inh.pdf',height=4, width=4)
    ggplot() + geom_jitter(Sum_parameters_plotting[which(Sum_parameters_plotting$Stage=="1k"&Sum_parameters_plotting$Inhibition=="Yes"),], 
                           mapping = aes(y=as.numeric(as.character(Total_intensity.normalized)), 
                                         x=as.numeric(as.character(Seconds_normalized))), 
                           width=5,alpha=0.1, size = 1,show.legend = FALSE, color="darkred")+
      geom_line(a,mapping=aes(x=Time,y=Median_name), size=0.75,color="firebrick3")+
      geom_ribbon(a,mapping=aes(x=Time,y=Median_name,
                                ymin=CI_inf_name,ymax=CI_sup_name),
                  alpha=0.15, fill="darkred")+
      scale_linetype_manual(values=c("dashed","solid"))+
      geom_vline(aes(xintercept = 0,), size=0.5,color="red", linetype=2)+ 
      theme_bw() + labs(x= "Time centered on transcription start (s)",y="Normalized Total Intensity")+ 
      ylim(-0.01,1.01)+scale_x_continuous(breaks = seq(-90, 150, by = 30),limits = c(-95,155)) +
      theme(aspect.ratio=1, axis.text = element_text(size=14),axis.title = element_text(size=14),
            panel.grid.minor.x = element_blank())
    dev.off()
    }
    }
    
  }
  
############Fig S3a = problem, does not work
{
  #Calculate the values for the ribbons
  {
    Nanog_statistics_total_nuc$Seconds_normalized <-round(as.numeric(as.character(Nanog_statistics_total_nuc$Seconds_normalized)))
    DT <- as.data.table(Nanog_statistics_total_nuc)
    Summary <- DT[ , list(MEAN_DV = median(Mean_intensity_normalized), SD_DV = sd(Mean_intensity_normalized), Q10= quantile(Mean_intensity_normalized, prob=c(0.10), na.rm=TRUE),Q25 = quantile(Mean_intensity_normalized, prob=c(0.25),na.rm=TRUE),Q75 = quantile(Mean_intensity_normalized, prob=c(0.75),na.rm=TRUE), Q90= quantile(Mean_intensity_normalized),N = .N),  by = c("Seconds_normalized")]
    Summary<- Summary[order(Summary$Seconds_normalized),]
    
    Median_circ <- Summary$MEAN_DV
    sd_circ <- Summary$SD_DV
    
    CI_infcirc <- Summary$Q25
    CI_supcirc <- Summary$Q75
    
    Circ_plot <- cbind(Median_circ,sd_circ,CI_infcirc,CI_supcirc)
    Circ_plot <-Circ_plot[complete.cases(Circ_plot[,c("CI_supcirc")]),]
    Circ_plot <- cbind(Circ_plot,unique(Summary[,c("Seconds_normalized")]))
    Circ_plot <- data.frame(Circ_plot)
    
    colnames(Circ_plot) <- c("Median","sd","CI_inf","CI_sup","Time")
    
    b <- Circ_plot %>%
      group_by(Time) %>%
      summarise_at(vars("Median","sd","CI_inf","CI_sup"), list(name = mean))
  }
  
  ggplot() + geom_jitter(DT, 
                         mapping = aes(y=as.numeric(as.character(Mean_intensity_normalized)), 
                                       x=as.numeric(as.character(Seconds_normalized))), 
                         width=15,shape=20, alpha=0.15, size=3) + 
    geom_ribbon(b,mapping=aes(x=Time,ymin=CI_inf_name,ymax=CI_sup_name),
                alpha=0.6, fill="grey70")+
    geom_line(b,mapping=aes(x=Time,y=Median_name),size=0.5, color="forestgreen") +
    geom_vline(data = Nanog_statistics_total_nuc,aes(xintercept = 0,), size=0.5,color="red", linetype=2)+ 
    theme_bw() + labs(x= "Time centered on transcription initiation (s)",y="Normalized mean intensity the nucleus") +   
    scale_y_continuous(breaks=seq(0.75,1,0.02),limits = c(0.75,1.01))+
    scale_x_continuous(breaks=seq(-150,150,30),limits = c(-150,150))+ 
    theme(aspect.ratio=1, axis.text = element_text(size=12),axis.title = element_text(size=12),
          panel.grid.minor.y = element_blank(),panel.grid.minor.x = element_blank())
}

############Fig S3b = not the same plots!!
{
  #No inhibition, non-merging clusters
  {
    DT <- as.data.table(Sum_parameters_plotting[which(Sum_parameters_plotting$Stage=="1k"&Sum_parameters_plotting$Inhibition=="No"&Sum_parameters_plotting$Phenotype_simple=="One_cluster"),])
    Summary <- DT[ , list(MEAN_DV = median(na.omit(Total_intensity.normalized)), SD_DV = sd(na.omit(Total_intensity.normalized)), Q10= quantile(Total_intensity.normalized, prob=c(0.10), na.rm=TRUE),Q25 = quantile(Total_intensity.normalized, prob=c(0.25),na.rm=TRUE),Q75 = quantile(Total_intensity.normalized, prob=c(0.75),na.rm=TRUE), Q90= quantile(Total_intensity.normalized, prob=c(0.90),na.rm=TRUE),N = .N),  by = c("Seconds_normalized")]
    Summary<- Summary[order(Summary$Seconds_normalized),]
    
    Median_circ <- Summary$MEAN_DV
    sd_circ <- Summary$SD_DV
    
    CI_infcirc <- Summary$Q25
    CI_supcirc <- Summary$Q75
    
    Circ_plot <- cbind(Median_circ,sd_circ,CI_infcirc,CI_supcirc)
    Circ_plot <-Circ_plot[complete.cases(Circ_plot[,c("CI_supcirc")]),]
    Circ_plot <- cbind(Circ_plot,unique(Summary[,c("Seconds_normalized")]))
    Circ_plot <- data.frame(Circ_plot)
    
    colnames(Circ_plot) <- c("Median","sd","CI_inf","CI_sup","Time")
    
    a <- Circ_plot %>%
      group_by(Time) %>%
      summarise_at(vars("Median","sd","CI_inf","CI_sup"), list(name = mean))
    
    pdf(file = '/Users/nchabot/Desktop/Figures3bf_wt_non_merging',height=4, width=4)
    ggplot() + geom_jitter(Sum_parameters_plotting[Sum_parameters_plotting$Inhibition=="No"&Sum_parameters_plotting$Stage=="1k"&Sum_parameters_plotting$Phenotype_simple=="One_cluster",], 
                           mapping = aes(y=as.numeric(as.character(Total_intensity.normalized)), 
                                         x=as.numeric(as.character(Seconds_normalized))), 
                           width=5,alpha=0.1, size = 1,show.legend = FALSE, color="black")+
      geom_line(a,mapping=aes(x=Time,y=Median_name), size=0.75,color="black")+
      geom_ribbon(a,mapping=aes(x=Time,y=Median_name,
                                ymin=CI_inf_name,ymax=CI_sup_name),
                  alpha=0.15, fill="black")+
      scale_linetype_manual(values=c("dashed","solid"))+
      geom_vline(aes(xintercept = 0,), size=0.5,color="red", linetype=2)+ 
      theme_bw() + labs(x= "Time centered on transcription start (s)",y="Normalized Total Intensity")+ 
      ylim(-0.01,1.01)+scale_x_continuous(breaks = seq(-90, 150, by = 30),limits = c(-95,155)) +
      theme(aspect.ratio=1, axis.text = element_text(size=14),axis.title = element_text(size=14),
            panel.grid.minor.x = element_blank())
    dev.off()
  }
  
  #No inhibition, merging clusters
  {
    DT <- as.data.table(Sum_parameters_plotting[which(Sum_parameters_plotting$Stage=="1k"&Sum_parameters_plotting$Inhibition=="No"&Sum_parameters_plotting$Phenotype_simple=="More_cluster"),])
    Summary <- DT[ , list(MEAN_DV = median(na.omit(Total_intensity.normalized)), SD_DV = sd(na.omit(Total_intensity.normalized)), Q10= quantile(Total_intensity.normalized, prob=c(0.10), na.rm=TRUE),Q25 = quantile(Total_intensity.normalized, prob=c(0.25),na.rm=TRUE),Q75 = quantile(Total_intensity.normalized, prob=c(0.75),na.rm=TRUE), Q90= quantile(Total_intensity.normalized, prob=c(0.90),na.rm=TRUE),N = .N),  by = c("Seconds_normalized")]
    Summary<- Summary[order(Summary$Seconds_normalized),]
    
    Median_circ <- Summary$MEAN_DV
    sd_circ <- Summary$SD_DV
    
    CI_infcirc <- Summary$Q25
    CI_supcirc <- Summary$Q75
    
    Circ_plot <- cbind(Median_circ,sd_circ,CI_infcirc,CI_supcirc)
    Circ_plot <-Circ_plot[complete.cases(Circ_plot[,c("CI_supcirc")]),]
    Circ_plot <- cbind(Circ_plot,unique(Summary[,c("Seconds_normalized")]))
    Circ_plot <- data.frame(Circ_plot)
    
    colnames(Circ_plot) <- c("Median","sd","CI_inf","CI_sup","Time")
    
    a <- Circ_plot %>%
      group_by(Time) %>%
      summarise_at(vars("Median","sd","CI_inf","CI_sup"), list(name = mean))
    
    pdf(file = '/Users/nchabot/Desktop/Figures3bf_wt_merging',height=4, width=4)
    ggplot() + geom_jitter(Sum_parameters_plotting[Sum_parameters_plotting$Inhibition=="No"&Sum_parameters_plotting$Stage=="1k"&Sum_parameters_plotting$Phenotype_simple=="More_cluster",], 
                           mapping = aes(y=as.numeric(as.character(Total_intensity.normalized)), 
                                         x=as.numeric(as.character(Seconds_normalized))), 
                           width=5,alpha=0.1, size = 1,show.legend = FALSE, color="black")+
      geom_line(a,mapping=aes(x=Time,y=Median_name), size=0.75,color="black")+
      geom_ribbon(a,mapping=aes(x=Time,y=Median_name,
                                ymin=CI_inf_name,ymax=CI_sup_name),
                  alpha=0.15, fill="black")+
      scale_linetype_manual(values=c("dashed","solid"))+
      geom_vline(aes(xintercept = 0,), size=0.5,color="red", linetype=2)+ 
      theme_bw() + labs(x= "Time centered on transcription start (s)",y="Normalized Total Intensity")+ 
      ylim(-0.01,1.01)+scale_x_continuous(breaks = seq(-90, 150, by = 30),limits = c(-95,155)) +
      theme(aspect.ratio=1, axis.text = element_text(size=14),axis.title = element_text(size=14),
            panel.grid.minor.x = element_blank())
    dev.off()
  }
  
  #Inhibition, non-merging clusters
  {
    DT <- as.data.table(Sum_parameters_plotting[which(Sum_parameters_plotting$Stage=="1k"&Sum_parameters_plotting$Inhibition=="Yes"&Sum_parameters_plotting$Phenotype_simple=="One_cluster"),])
    Summary <- DT[ , list(MEAN_DV = median(na.omit(Total_intensity.normalized)), SD_DV = sd(na.omit(Total_intensity.normalized)), Q10= quantile(Total_intensity.normalized, prob=c(0.10), na.rm=TRUE),Q25 = quantile(Total_intensity.normalized, prob=c(0.25),na.rm=TRUE),Q75 = quantile(Total_intensity.normalized, prob=c(0.75),na.rm=TRUE), Q90= quantile(Total_intensity.normalized, prob=c(0.90),na.rm=TRUE),N = .N),  by = c("Seconds_normalized")]
    Summary<- Summary[order(Summary$Seconds_normalized),]
    
    Median_circ <- Summary$MEAN_DV
    sd_circ <- Summary$SD_DV
    
    CI_infcirc <- Summary$Q25
    CI_supcirc <- Summary$Q75
    
    Circ_plot <- cbind(Median_circ,sd_circ,CI_infcirc,CI_supcirc)
    Circ_plot <-Circ_plot[complete.cases(Circ_plot[,c("CI_supcirc")]),]
    Circ_plot <- cbind(Circ_plot,unique(Summary[,c("Seconds_normalized")]))
    Circ_plot <- data.frame(Circ_plot)
    
    colnames(Circ_plot) <- c("Median","sd","CI_inf","CI_sup","Time")
    
    a <- Circ_plot %>%
      group_by(Time) %>%
      summarise_at(vars("Median","sd","CI_inf","CI_sup"), list(name = mean))
    
    pdf(file = '/Users/nchabot/Desktop/Figures3bf_inh_non_merging',height=4, width=4)
    ggplot() + geom_jitter(Sum_parameters_plotting[Sum_parameters_plotting$Inhibition=="Yes"&Sum_parameters_plotting$Stage=="1k"&Sum_parameters_plotting$Phenotype_simple=="One_cluster",], 
                           mapping = aes(y=as.numeric(as.character(Total_intensity.normalized)), 
                                         x=as.numeric(as.character(Seconds_normalized))), 
                           width=5,alpha=0.1, size = 1,show.legend = FALSE, color="darkred")+
      geom_line(a,mapping=aes(x=Time,y=Median_name), size=0.75,color="firebrick3")+
      geom_ribbon(a,mapping=aes(x=Time,y=Median_name,
                                ymin=CI_inf_name,ymax=CI_sup_name),
                  alpha=0.15, fill="darkred")+
      scale_linetype_manual(values=c("dashed","solid"))+
      geom_vline(aes(xintercept = 0,), size=0.5,color="red", linetype=2)+ 
      theme_bw() + labs(x= "Time centered on transcription start (s)",y="Normalized Total Intensity")+ 
      ylim(-0.01,1.01)+scale_x_continuous(breaks = seq(-90, 150, by = 30),limits = c(-95,155)) +
      theme(aspect.ratio=1, axis.text = element_text(size=14),axis.title = element_text(size=14),
            panel.grid.minor.x = element_blank())
    dev.off()
  }
  
  #Inhibition, merging clusters
  {
    DT <- as.data.table(Sum_parameters_plotting[which(Sum_parameters_plotting$Stage=="1k"&Sum_parameters_plotting$Inhibition=="Yes"&Sum_parameters_plotting$Phenotype_simple=="More_cluster"),])
    Summary <- DT[ , list(MEAN_DV = median(na.omit(Total_intensity.normalized)), SD_DV = sd(na.omit(Total_intensity.normalized)), Q10= quantile(Total_intensity.normalized, prob=c(0.10), na.rm=TRUE),Q25 = quantile(Total_intensity.normalized, prob=c(0.25),na.rm=TRUE),Q75 = quantile(Total_intensity.normalized, prob=c(0.75),na.rm=TRUE), Q90= quantile(Total_intensity.normalized, prob=c(0.90),na.rm=TRUE),N = .N),  by = c("Seconds_normalized")]
    Summary<- Summary[order(Summary$Seconds_normalized),]
    
    Median_circ <- Summary$MEAN_DV
    sd_circ <- Summary$SD_DV
    
    CI_infcirc <- Summary$Q25
    CI_supcirc <- Summary$Q75
    
    Circ_plot <- cbind(Median_circ,sd_circ,CI_infcirc,CI_supcirc)
    Circ_plot <-Circ_plot[complete.cases(Circ_plot[,c("CI_supcirc")]),]
    Circ_plot <- cbind(Circ_plot,unique(Summary[,c("Seconds_normalized")]))
    Circ_plot <- data.frame(Circ_plot)
    
    colnames(Circ_plot) <- c("Median","sd","CI_inf","CI_sup","Time")
    
    a <- Circ_plot %>%
      group_by(Time) %>%
      summarise_at(vars("Median","sd","CI_inf","CI_sup"), list(name = mean))
    
    pdf(file = '/Users/nchabot/Desktop/Figures3bf_inh_merging.pdf',height=4, width=4)
    ggplot() + geom_jitter(Sum_parameters_plotting[Sum_parameters_plotting$Inhibition=="Yes"&Sum_parameters_plotting$Stage=="1k"&Sum_parameters_plotting$Phenotype_simple=="More_cluster",], 
                           mapping = aes(y=as.numeric(as.character(Total_intensity.normalized)), 
                                         x=as.numeric(as.character(Seconds_normalized))), 
                           width=5,alpha=0.1, size = 1,show.legend = FALSE, color="darkred")+
      geom_line(a,mapping=aes(x=Time,y=Median_name), size=0.75,color="firebrick3")+
      geom_ribbon(a,mapping=aes(x=Time,y=Median_name,
                                ymin=CI_inf_name,ymax=CI_sup_name),
                  alpha=0.15, fill="darkred")+
      scale_linetype_manual(values=c("dashed","solid"))+
      geom_vline(aes(xintercept = 0,), size=0.5,color="red", linetype=2)+ 
      theme_bw() + labs(x= "Time centered on transcription start (s)",y="Normalized Total Intensity")+ 
      ylim(-0.01,1.01)+scale_x_continuous(breaks = seq(-90, 150, by = 30),limits = c(-95,155)) +
      theme(aspect.ratio=1, axis.text = element_text(size=14),axis.title = element_text(size=14),
            panel.grid.minor.x = element_blank())
    dev.off()
  }
}

############Fig 2a 
{
  #Merging clusters, centered on merging, mean intensity
  {
    DT <- as.data.table(Nanog_statistics_total_10[which(Nanog_statistics_total_10$Stage=="1k"&Nanog_statistics_total_10$Inhibition=="No"),])
    Summary <- DT[ , list(MEAN_DV = median(na.omit(Intensity.Mean.y.normalized)), SD_DV = sd(na.omit(Intensity.Mean.y.normalized)), Q10= quantile(Intensity.Mean.y.normalized, prob=c(0.10), na.rm=TRUE),Q25 = quantile(Intensity.Mean.y.normalized, prob=c(0.25),na.rm=TRUE),Q75 = quantile(Intensity.Mean.y.normalized, prob=c(0.75),na.rm=TRUE), Q90= quantile(Intensity.Mean.y.normalized, prob=c(0.90),na.rm=TRUE),N = .N),  by = c("Seconds_merging_normalized")]
    Summary<- Summary[order(Summary$Seconds_merging_normalized),]
    
    Median_circ <- Summary$MEAN_DV
    sd_circ <- Summary$SD_DV
    
    CI_infcirc <- Summary$Q25
    CI_supcirc <- Summary$Q75
    
    Circ_plot <- cbind(Median_circ,sd_circ,CI_infcirc,CI_supcirc)
    Circ_plot <-Circ_plot[complete.cases(Circ_plot[,c("CI_supcirc")]),]
    Circ_plot <- cbind(Circ_plot,unique(Summary[,c("Seconds_merging_normalized")]))
    Circ_plot <- data.frame(Circ_plot)
    
    colnames(Circ_plot) <- c("Median","sd","CI_inf","CI_sup","Time")
    
    a <- Circ_plot %>%
      group_by(Time) %>%
      summarise_at(vars("Median","sd","CI_inf","CI_sup"), list(name = mean))
    
    pdf(file = '/Users/nchabot/Desktop/Figure2a_mean_intensity',height=4, width=4)
    ggplot() + geom_jitter(Nanog_statistics_total_10[Nanog_statistics_total_10$Inhibition=="No"&Nanog_statistics_total_10$Stage=="1k",], 
                           mapping = aes(y=as.numeric(as.character(Intensity.Mean.y.normalized)), 
                                         x=as.numeric(as.character(Seconds_merging_normalized))), 
                           width=5,alpha=0.1, size = 1,show.legend = FALSE, color="black")+
      geom_line(a,mapping=aes(x=Time,y=Median_name), size=0.75,color="black")+
      geom_ribbon(a,mapping=aes(x=Time,y=Median_name,
                                ymin=CI_inf_name,ymax=CI_sup_name),
                  alpha=0.15, fill="black")+
      scale_linetype_manual(values=c("dashed","solid"))+
      geom_vline(aes(xintercept = 0,), size=0.5,color="red", linetype=2)+ 
      theme_bw() + labs(x= "Time centered on transcription start (s)",y="Normalized Mean Intensity")+ 
      ylim(-0.01,1.01)+scale_x_continuous(breaks = seq(-90, 150, by = 30),limits = c(-95,155)) +
      theme(aspect.ratio=1, axis.text = element_text(size=14),axis.title = element_text(size=14),
            panel.grid.minor.x = element_blank())
    dev.off()
  }
  
  #Merging clusters, centered on merging, total intensity
  {
    DT <- as.data.table(Nanog_statistics_total_10[which(Nanog_statistics_total_10$Stage=="1k"&Nanog_statistics_total_10$Inhibition=="No"),])
    Summary <- DT[ , list(MEAN_DV = median(na.omit(Intensity.Sum.y.normalized)), SD_DV = sd(na.omit(Intensity.Sum.y.normalized)), Q10= quantile(Intensity.Sum.y.normalized, prob=c(0.10), na.rm=TRUE),Q25 = quantile(Intensity.Sum.y.normalized, prob=c(0.25),na.rm=TRUE),Q75 = quantile(Intensity.Sum.y.normalized, prob=c(0.75),na.rm=TRUE), Q90= quantile(Intensity.Sum.y.normalized, prob=c(0.90),na.rm=TRUE),N = .N),  by = c("Seconds_merging_normalized")]
    Summary<- Summary[order(Summary$Seconds_merging_normalized),]
    
    Median_circ <- Summary$MEAN_DV
    sd_circ <- Summary$SD_DV
    
    CI_infcirc <- Summary$Q25
    CI_supcirc <- Summary$Q75
    
    Circ_plot <- cbind(Median_circ,sd_circ,CI_infcirc,CI_supcirc)
    Circ_plot <-Circ_plot[complete.cases(Circ_plot[,c("CI_supcirc")]),]
    Circ_plot <- cbind(Circ_plot,unique(Summary[,c("Seconds_merging_normalized")]))
    Circ_plot <- data.frame(Circ_plot)
    
    colnames(Circ_plot) <- c("Median","sd","CI_inf","CI_sup","Time")
    
    a <- Circ_plot %>%
      group_by(Time) %>%
      summarise_at(vars("Median","sd","CI_inf","CI_sup"), list(name = mean))
    
    pdf(file = '/Users/nchabot/Desktop/Figure2a_mean_intensity',height=4, width=4)
    ggplot() + geom_jitter(Nanog_statistics_total_10[Nanog_statistics_total_10$Inhibition=="No"&Nanog_statistics_total_10$Stage=="1k",], 
                           mapping = aes(y=as.numeric(as.character(Intensity.Sum.y.normalized)), 
                                         x=as.numeric(as.character(Seconds_merging_normalized))), 
                           width=5,alpha=0.1, size = 1,show.legend = FALSE, color="black")+
      geom_line(a,mapping=aes(x=Time,y=Median_name), size=0.75,color="black")+
      geom_ribbon(a,mapping=aes(x=Time,y=Median_name,
                                ymin=CI_inf_name,ymax=CI_sup_name),
                  alpha=0.15, fill="black")+
      scale_linetype_manual(values=c("dashed","solid"))+
      geom_vline(aes(xintercept = 0,), size=0.5,color="red", linetype=2)+ 
      theme_bw() + labs(x= "Time centered on transcription start (s)",y="Normalized Mean Intensity")+ 
      ylim(-0.01,1.01)+scale_x_continuous(breaks = seq(-90, 150, by = 30),limits = c(-95,155)) +
      theme(aspect.ratio=1, axis.text = element_text(size=14),axis.title = element_text(size=14),
            panel.grid.minor.x = element_blank())
    dev.off()
  }
  
  #Merging clusters, centered on merging, volume 
  {
    DT <- as.data.table(Nanog_statistics_total_10[which(Nanog_statistics_total_10$Stage=="1k"&Nanog_statistics_total_10$Inhibition=="No"),])
    Summary <- DT[ , list(MEAN_DV = median(na.omit(Volume.normalized)), SD_DV = sd(na.omit(Volume.normalized)), Q10= quantile(Volume.normalized, prob=c(0.10), na.rm=TRUE),Q25 = quantile(Volume.normalized, prob=c(0.25),na.rm=TRUE),Q75 = quantile(Volume.normalized, prob=c(0.75),na.rm=TRUE), Q90= quantile(Volume.normalized, prob=c(0.90),na.rm=TRUE),N = .N),  by = c("Seconds_merging_normalized")]
    Summary<- Summary[order(Summary$Seconds_merging_normalized),]
    
    Median_circ <- Summary$MEAN_DV
    sd_circ <- Summary$SD_DV
    
    CI_infcirc <- Summary$Q25
    CI_supcirc <- Summary$Q75
    
    Circ_plot <- cbind(Median_circ,sd_circ,CI_infcirc,CI_supcirc)
    Circ_plot <-Circ_plot[complete.cases(Circ_plot[,c("CI_supcirc")]),]
    Circ_plot <- cbind(Circ_plot,unique(Summary[,c("Seconds_merging_normalized")]))
    Circ_plot <- data.frame(Circ_plot)
    
    colnames(Circ_plot) <- c("Median","sd","CI_inf","CI_sup","Time")
    
    a <- Circ_plot %>%
      group_by(Time) %>%
      summarise_at(vars("Median","sd","CI_inf","CI_sup"), list(name = mean))
    
    pdf(file = '/Users/nchabot/Desktop/Figure2a_mean_intensity',height=4, width=4)
    ggplot() + geom_jitter(Nanog_statistics_total_10[Nanog_statistics_total_10$Inhibition=="No"&Nanog_statistics_total_10$Stage=="1k",], 
                           mapping = aes(y=as.numeric(as.character(Volume.normalized)), 
                                         x=as.numeric(as.character(Seconds_merging_normalized))), 
                           width=5,alpha=0.1, size = 1,show.legend = FALSE, color="black")+
      geom_line(a,mapping=aes(x=Time,y=Median_name), size=0.75,color="black")+
      geom_ribbon(a,mapping=aes(x=Time,y=Median_name,
                                ymin=CI_inf_name,ymax=CI_sup_name),
                  alpha=0.15, fill="black")+
      scale_linetype_manual(values=c("dashed","solid"))+
      geom_vline(aes(xintercept = 0,), size=0.5,color="red", linetype=2)+ 
      theme_bw() + labs(x= "Time centered on transcription start (s)",y="Normalized Mean Intensity")+ 
      ylim(-0.01,1.01)+scale_x_continuous(breaks = seq(-90, 150, by = 30),limits = c(-95,155)) +
      theme(aspect.ratio=1, axis.text = element_text(size=14),axis.title = element_text(size=14),
            panel.grid.minor.x = element_blank())
    dev.off()
  }
 
}

############Fig 2b
{
  Nanog_statistics_total_10$Difference_activation_merging <- Nanog_statistics_total_10$Time_merging - Nanog_statistics_total_10$Time_activation 
  Histo_fusion <- unique(Nanog_statistics_total_10[which(Nanog_statistics_total_10$Inhibition=="No"&Nanog_statistics_total_10$Stage=="1k"),c("Unique_track","Difference_activation_merging","Inhibition","Stage")])
  
  ggplot(Histo_fusion, aes(x=Difference_activation_merging)) + 
    geom_histogram(aes(y=(..count..)/sum(..count..)),color="black", fill="grey", binwidth=1)+
    theme(axis.text = element_text(size = 20),
          strip.text = element_text(size = 20), 
          legend.text = element_text(size = 20), 
          axis.title = element_text(size = 20),
          panel.background = element_blank(),
          axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) +
    scale_x_continuous(breaks = seq(-8, 8, by = 1))+
    scale_y_continuous(labels=percent)+
    labs(x= "Fusion time compared to activation",y="Number of occurences") 
}

############Fig 2c
{
  #Make the sum of volume of Nanog foci that
  {
  Nanog_statistics_total_10$Subtrack <- 0
  Nanog_statistics_total_10$Volume_difference <- 0
  Volume_activation <- data.frame()
  k=1
  
  for(i in 1:max(as.numeric(Nanog_statistics_total_10$Nucleus.number)))
  {
    if(nrow(Nanog_statistics_total_10[which(Nanog_statistics_total_10$Nucleus.number==i),])!=0)
      
    {
      if(length(unique(Nanog_statistics_total_10[which(Nanog_statistics_total_10$Nucleus.number==i),"Unique_track"]))==1)
        
      {Nanog_statistics_total_10[which(Nanog_statistics_total_10$Nucleus.number==i),"Subtrack"] <- 1}
      
      if(length(unique(Nanog_statistics_total_10[which(Nanog_statistics_total_10$Nucleus.number==i),"Unique_track"]))==2)
      {Nanog_statistics_total_10[which(Nanog_statistics_total_10$Nucleus.number==i&Nanog_statistics_total_10$Unique_track==unique(Nanog_statistics_total_10[which(Nanog_statistics_total_10$Nucleus.number==i),"Unique_track"])[1]),"Subtrack"] <- 1
      Nanog_statistics_total_10[which(Nanog_statistics_total_10$Nucleus.number==i&Nanog_statistics_total_10$Unique_track==unique(Nanog_statistics_total_10[which(Nanog_statistics_total_10$Nucleus.number==i),"Unique_track"])[2]),"Subtrack"] <- 2
      
      print(sum(Nanog_statistics_total_10[which(Nanog_statistics_total_10$Nucleus.number==i&Nanog_statistics_total_10$Time_normalized==0&Nanog_statistics_total_10$Subtrack==1),"Volume"]))
      print(sum(Nanog_statistics_total_10[which(Nanog_statistics_total_10$Nucleus.number==i&Nanog_statistics_total_10$Time_normalized==0&Nanog_statistics_total_10$Subtrack==2),"Volume"]))
      
      
      Nanog_statistics_total_10[which(Nanog_statistics_total_10$Nucleus.number==i),"Volume_difference"] <- sum(Nanog_statistics_total_10[which(Nanog_statistics_total_10$Nucleus.number==i&Nanog_statistics_total_10$Time_normalized==0&Nanog_statistics_total_10$Subtrack==1),"Volume"])-sum(Nanog_statistics_total_10[which(Nanog_statistics_total_10$Nucleus.number==i&Nanog_statistics_total_10$Time_normalized==0&Nanog_statistics_total_10$Subtrack==2),"Volume"])
      
      Substract_1 <-which(Nanog_statistics_total_10$Nucleus.number==i&Nanog_statistics_total_10$Subtrack==1)
      Substract_2 <- which(Nanog_statistics_total_10$Nucleus.number==i&Nanog_statistics_total_10$Subtrack==2)
      
      Volume_activation[k,1:4] <- Nanog_statistics_total_10[which(Nanog_statistics_total_10$Nucleus.number==i),1:4][1,]
      Volume_activation[k,"Unique_track_1"] <- unique(Nanog_statistics_total_10[Substract_1,"Unique_track"])
      Volume_activation[k,"Unique_track_2"] <- unique(Nanog_statistics_total_10[Substract_2,"Unique_track"])
      Volume_activation[k,"Mean_intensity_nucleus"] <- unique(Nanog_statistics_total_10[which(Nanog_statistics_total_10$Nucleus.number==i),"Intensity.Mean.y"])[2]
      Volume_activation[k,"Inhibition"] <- unique(Nanog_statistics_total_10[which(Nanog_statistics_total_10$Nucleus.number==i),"Inhibition"])[1]
      Volume_activation[k,"Phenotype_1"] <- unique(Nanog_statistics_total_10[Substract_1,"Phenotype_simple"])
      
      Volume_activation[k,"Phenotype_1"] <- unique(Nanog_statistics_total_10[Substract_1,"Phenotype_simple"])
      Volume_activation[k,"Phenotype_2"] <- unique(Nanog_statistics_total_10[Substract_2,"Phenotype_simple"])
      
      Volume_activation[k,"Sum_Volume_1"] <- sum(Nanog_statistics_total_10[which(Nanog_statistics_total_10$Nucleus.number==i&Nanog_statistics_total_10$Time_normalized==0&Nanog_statistics_total_10$Subtrack==1),"Volume"])
      Volume_activation[k,"Sum_Volume_2"] <- sum(Nanog_statistics_total_10[which(Nanog_statistics_total_10$Nucleus.number==i&Nanog_statistics_total_10$Time_normalized==0&Nanog_statistics_total_10$Subtrack==2),"Volume"])
      
      Volume_activation[k,"Volume_difference"] <- Volume_activation[k,"Sum_Volume_1"]- Volume_activation[k,"Sum_Volume_2"] 
      
      
      Volume_activation[k,"Sum_Total_intensity_1"] <- sum(Nanog_statistics_total_10[which(Nanog_statistics_total_10$Nucleus.number==i&Nanog_statistics_total_10$Time_normalized==0&Nanog_statistics_total_10$Subtrack==1),"Intensity.Sum.y"])
      Volume_activation[k,"Sum_Total_intensity_2"] <- sum(Nanog_statistics_total_10[which(Nanog_statistics_total_10$Nucleus.number==i&Nanog_statistics_total_10$Time_normalized==0&Nanog_statistics_total_10$Subtrack==2),"Intensity.Sum.y"])
      Volume_activation[k,"Sum_Total_difference"] <- Volume_activation[k,"Sum_Total_intensity_1"]- Volume_activation[k,"Sum_Total_intensity_2"] 
      
      Volume_activation[k,"Mean_intensity_1"] <- mean(Nanog_statistics_total_10[which(Nanog_statistics_total_10$Nucleus.number==i&Nanog_statistics_total_10$Time_normalized==0&Nanog_statistics_total_10$Subtrack==1),"Intensity.Mean.y"])
      Volume_activation[k,"Mean_intensity_2"] <- mean(Nanog_statistics_total_10[which(Nanog_statistics_total_10$Nucleus.number==i&Nanog_statistics_total_10$Time_normalized==0&Nanog_statistics_total_10$Subtrack==2),"Intensity.Mean.y"])
      
      Volume_activation[k,"Mean_intensity_difference"] <- Volume_activation[k,"Mean_intensity_1"] - Volume_activation[k,"Mean_intensity_2"]
      Volume_activation[k,"Time_activation_1"] <- unique(Nanog_statistics_total_10[which(Nanog_statistics_total_10$Nucleus.number==i&Nanog_statistics_total_10$Subtrack==1),"Time_activation"])
      Volume_activation[k,"Time_activation_2"] <- unique(Nanog_statistics_total_10[which(Nanog_statistics_total_10$Nucleus.number==i&Nanog_statistics_total_10$Subtrack==2),"Time_activation"])
      
      
      k=k+1
      
      }
      }
  }
  Volume_activation$Activation_difference <- Volume_activation$Time_activation_1-Volume_activation$Time_activation_2
  
  }
  #Create a column to indicate if the two tracks display the same phenotype
  {
    for(i in 1:nrow(Volume_activation)) {
      if(Volume_activation[i,"Phenotype_1"]==Volume_activation[i,"Phenotype_2"])
      {
        Volume_activation[i,"Same_phenotype"] <- "Same_phenotype"
      }
      if(Volume_activation[i,"Phenotype_1"]!=Volume_activation[i,"Phenotype_2"])
      {
        Volume_activation[i,"Same_phenotype"] <- "Different_phenotype"
      }
    }
  }
  #Code of the plot
  {
  pdf(file = '/Users/nchabot/Desktop/Fig2c.pdf',height=5, width=4)
  ggplot(Volume_activation[Volume_activation$Stage=="1k"&Volume_activation$Inhibition=="No",], aes(y=Activation_difference, x=Same_phenotype)) + geom_violin(Volume_activation[Volume_activation$Stage=="1k"&Volume_activation$Inhibition=="No",], mapping = aes(y=Activation_difference, x=Same_phenotype), size=0.5, width=0.5,scale="area",trim="FALSE") + 
    stat_summary(fun=median, fun.min = function(z) { quantile(z,0.25) },
                 fun.max = function(z) { quantile(z,0.75) },
                 geom="pointrange", color="black") + 
    labs(x= "Stage",y="Difference between time of activation with same or different phenotypes")    +
    theme(aspect.ratio=5/3, axis.text = element_text(size=14),axis.title = element_text(size=14))+  
    stat_compare_means(comparisons = list(c("Same_phenotype","Different_phenotype")), method="wilcox.test") +    
    theme_minimal()
  dev.off()
  }
}

############Fig 2d
{
  Sum_parameters_plotting$Length_double <- as.numeric(Sum_parameters_plotting$Length_double)
  Sum_parameters_plotting$Unique_track <- as.factor(Sum_parameters_plotting$Unique_track)
  Sum_parameters_plotting_1k_no_inh <- Sum_parameters_plotting[which(Sum_parameters_plotting$Stage=="1k"&Sum_parameters_plotting$Inhibition=="No"&Sum_parameters_plotting$Length_double>0),]
  
  ggplot(na.omit(unique(Sum_parameters_plotting_1k_no_inh[,c("Unique_track","Length_double"),])), aes(x=Length_double,width=0.7)) + 
    geom_histogram(aes(y = after_stat(count / sum(count))),color="black", fill="grey70", binwidth=1) +
    geom_text(aes(y = after_stat(count / sum(count)), label = scales::percent(after_stat(count / sum(count)))), stat = "count", vjust = -0.25) +
    scale_x_continuous(breaks = seq(0, 8, by = 1))+
    scale_y_continuous(labels=scales::percent, limits =c(0, 0.55))+
    labs(x= "Number of consecutive time points \n with two or more spots",y="Percentage") + theme_bw()+
    theme(aspect.ratio=2/4, axis.text = element_text(size=12),axis.title = element_text(size=14)) 
  
}

############Fig 2e
{
  Nanog_statistics_total_10_3 <- Nanog_statistics_total_10[which(!is.na(Nanog_statistics_total_10$TrackID.y)),]
  
  #Non-merging cluster
  Non_merging_1k_non_inh <- Nanog_statistics_total_10[Nanog_statistics_total_10$Inhibition=="No"&Nanog_statistics_total_10$Stage=="1k"&Nanog_statistics_total_10$Phenotype_simple=="One_cluster",]
  Non_merging_1k_non_inh_mean <- Nanog_statistics_total_10[Nanog_statistics_total_10$Inhibition=="No"&Nanog_statistics_total_10$Stage=="1k"&Nanog_statistics_total_10$Phenotype_simple=="One_cluster",c("Intensity.Sum.y.normalized","Seconds_normalized")] %>% 
    group_by(Seconds_normalized) %>%
    summarise_each(funs(mean,sd,se=sd(.)/sqrt(n())))  
  
  #Merging cluster >3 individuals
  Nanog_statistics_total_10_3
  Merging_1k_non_inh_mean_ind <- Nanog_statistics_total_10_3[,c("Intensity.Sum.y.normalized","Seconds_normalized")] %>% 
    group_by(Seconds_normalized) %>%
    summarise_each(funs(mean,sd,se=sd(.)/sqrt(n())))   
  
  #Merging cluster >3 sum
  Sum_parameters_plotting_3 <- Sum_parameters_plotting[which(!is.na(Sum_parameters_plotting$TrackID.y)),]
  Non_merging_1k_non_inh_mean_sum <- Sum_parameters_plotting_3[,c("Total_intensity.normalized","Seconds_normalized")] %>% 
    group_by(Seconds_normalized) %>%
    summarise_each(funs(mean,sd,se=sd(.)/sqrt(n())))  
  
  pdf(file = '/Users/nchabot/Desktop/210524/Pannels/240602_Fig2e.pdf',height=4, width=4)
  ggplot() + 
    geom_jitter(Non_merging_1k_non_inh, 
                mapping = aes(y=as.numeric(as.character(Intensity.Sum.y.normalized)), 
                              x=as.numeric(as.character(Seconds_normalized))), 
                width=2.5, alpha=0.2, size = 2,color="grey30",shape=16)+
    geom_line(Non_merging_1k_non_inh_mean,mapping=aes(x=Seconds_normalized,y=mean), size=0.5,color="grey30")+
    geom_ribbon(Non_merging_1k_non_inh_mean,mapping=aes(x=Seconds_normalized,y=mean,ymin=mean-se,ymax=mean+se),
                alpha=0.3, fill="grey30")+
    
    geom_jitter(Nanog_statistics_total_10_3,
                mapping = aes(y=Intensity.Sum.y.normalized, 
                              x=as.numeric(as.character(Seconds_normalized))), 
                width=2.5,alpha=0.2,size = 1.5,fill="skyblue",shape=24,stroke=0.2)+
    geom_line(Merging_1k_non_inh_mean_ind,mapping=aes(x=Seconds_normalized,y=mean), size=0.5,color="steelblue4",linetype="dashed")+
    geom_ribbon(Merging_1k_non_inh_mean_ind,mapping=aes(x=Seconds_normalized,y=mean,ymin=mean-se,ymax=mean+se),
                alpha=0.3, fill="skyblue")+
    
    geom_jitter(Sum_parameters_plotting_3, 
                mapping=aes(y=Total_intensity.normalized, 
                            x=as.numeric(as.character(Seconds_normalized))), 
                width=2.5,alpha=0.2, size =2,fill="steelblue4",shape=23,stroke=0.2)+
    geom_line(Non_merging_1k_non_inh_mean_sum,mapping=aes(x=Seconds_normalized,y=mean), size=0.5,color="steelblue4")+
    geom_ribbon(Non_merging_1k_non_inh_mean_sum,mapping=aes(x=Seconds_normalized,y=mean,ymin=mean-se,ymax=mean+se),
                alpha=0.3, fill="steelblue4")+
    theme_bw() + labs(x= "Time centered on transcription start (s)",y="Normalized Total Intensity")+ 
    ylim(-0.01,1.01)+scale_x_continuous(breaks = seq(-90, 0, by = 15),limits = c(-95,5)) +
    theme(aspect.ratio=1, axis.text = element_text(size=14),axis.title = element_text(size=14),
          panel.grid.minor.x = element_blank())
  dev.off()
  
  
}
