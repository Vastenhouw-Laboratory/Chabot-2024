# Noémie Chabot, PhD, Vastenhouw lab, 2024/07/10
# R code for analysis and figures related to the Fig.5a, S5a
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
Input = "/Volumes/CIG/nvastenh/nanog_movement/D2c/Manuscript/240612_NCB_NSMB_submission/Material_make_figures/Fig5a, S5a/3D_Imaris_files/"

#=========================================================
# Read the files necessary for the analysis
#=========================================================

#Change path to working directory
setwd(Input)

############MOVIE signal segmented with the shape algorithm
{
#All file names in one table
{
  #MOVIE
  #Get a list of files containing the word "MOVIE" to analyze only the MOVIE foci
  List_directories_MOVIE <- grep(pattern="MOVIE", list.dirs(path = ".", full.names = TRUE), value=TRUE)
  
  List_files_MOVIE_number <- data.frame()
  
  #Get all the MOVIE file in columns for each nucleus
  
  for(e in 1:length(List_directories_MOVIE))
  {
    setwd(Input)
    
    #print(length(list.files(List_directories_MOVIE[e])))
    if(e == 1) {List_files_MOVIE_number <- data.frame(list.files(path = List_directories_MOVIE[e], recursive=TRUE)) }
    else {
      #Assemble all data in a unique dataframe
      List_files_MOVIE_number <- List_files_MOVIE_number %>% data.frame(list.files(path = List_directories_MOVIE[e], pattern="MOVIE", recursive=TRUE))}
  }
}

#Assemble all data in one dataframe

dataset <- data.frame()
Infos_file <- data.frame()
Infinity <- data.frame()

for(e in 1:length(List_directories_MOVIE)) {
  
  setwd(Input)
  List_files_MOVIE_number <- list.files(path = List_directories_MOVIE[e], pattern="MOVIE", recursive=TRUE)
  setwd(List_directories_MOVIE[e])
  
  for (i in 1:length(List_files_MOVIE_number)){
    temp_data <- read.table(List_files_MOVIE_number[i],  sep=",", skip=3, header=TRUE) #each file will be read in, specify which columns you need read in to avoid any errors
    file_name <- as.character(List_files_MOVIE_number[i])
    
    #Select the info inside file name about date of acquisition and identity of the nucleus being analysing. These infos are stored in columns.
    #Adapt syntax based on the number of characters
    
    Infos_file[1,1] <- substr(file_name, 1,5)
    Infos_file[1,2] <- substr(file_name, 7,8)
    Infos_file[1,3] <- substr(file_name, 10,10)
    Infos_file[1,4] <- as.numeric(gsub(".*?([0-9]+).*", "\\1", substr(file_name, 12,13)))
    
    #Repeat n nows the infos about the nuclei, as much as they are segmented volumes in the file
    Infinity <- Infos_file[rep(seq_len(nrow(Infos_file)), each = nrow(temp_data)),]
    #Infinity$V5 <- sub("^\\D*(\\d+).*$", "\\1", List_directories_MOVIE[e])
    
    #Give a name to the column and make sure it looks like something
    colnames(Infinity) <- c("Date","Stage","Stack","Nucleus")
    #print(List_files_MOVIE_number[i])
    temp_data <- cbind(Infinity, temp_data)
    
    #Assemble all data in a unique dataframe
    if(i == 1) {MOVIE_statistics_all <- temp_data } else {
      if(grepl("Time", file_name)==TRUE)
      {colnames(temp_data)[5] <- "Seconds"
      colnames(temp_data)[8] <- "Time"
      MOVIE_statistics_all <- MOVIE_statistics_all %>% full_join(temp_data, by=c("ID","Date","Nucleus","Stack","Time","Stage","TrackID"))} 
      else {MOVIE_statistics_all <- MOVIE_statistics_all %>% full_join(temp_data, by=c("ID","Date","Nucleus","Stack","Time","Stage","TrackID"))}
    }
    
  }
  #Remove useless columns
  patterns <- c("Collection*", "Unit*","Image*","Channel*","X.x","X.y","Event*","Category*","Merging.foci.x*")
  MOVIE_statistics_all <- MOVIE_statistics_all[,-grep(paste(patterns, collapse="|"), colnames(MOVIE_statistics_all))]
  
  if(e == 1) {MOVIE_statistics_total <- MOVIE_statistics_all
  } else {
    #Mix data frames with different number of columns by row
    MOVIE_statistics_total <- bind_rows(MOVIE_statistics_total, MOVIE_statistics_all)
  }
  
}

Unique_ID <- unique(MOVIE_statistics_total[,c(1,3:4)])
MOVIE_statistics_total$Nucleus.number <- 0

for(i in 1:nrow(MOVIE_statistics_total))
{
  for(e in 1:nrow(Unique_ID))
  {
    if(MOVIE_statistics_total[i,1] == Unique_ID[e,1]&MOVIE_statistics_total[i,3] == Unique_ID[e,2]&MOVIE_statistics_total[i,4] == Unique_ID[e,3])
      
      MOVIE_statistics_total[i,"Nucleus.number"] <- e
    
  }
}

colnames(MOVIE_statistics_total)[which(names(MOVIE_statistics_total) == "TrackID")] <- "TrackID_MOVIE_3D"
}

###########Nanog clusters
{
  
  #All file names in one table
  {
    
    #Nanog
    #Get a list of files containing the word "Nanog" to analyze only the Nanog foci
    
    setwd(Input)
    
    List_directories_Nanog <- grep(pattern="Nanog", list.dirs(path = ".", full.names = TRUE), value=TRUE)
    List_files_Nanog_number <- data.frame()
    
    #Get all the Nanog file in columns for each nucleus
    
    for(e in 1:length(List_directories_Nanog))
    {
      setwd(Input)
      
      #print(length(list.files(List_directories_Nanog[e])))
      if(e == 1) {List_files_Nanog_number <- data.frame(list.files(path = List_directories_Nanog[e], recursive=TRUE)) }
      else {
        #Assemble all data in a unique dataframe
        List_files_Nanog_number <- List_files_Nanog_number %>% data.frame(list.files(path = List_directories_Nanog[e], pattern="Nanog", recursive=TRUE))}
    }
  }
  
  #Assemble all data in one dataframe
  
  dataset <- data.frame()
  Infos_file <- data.frame()
  Infinity <- data.frame()
  
  for(e in 1:length(List_directories_Nanog)) {
    
    setwd(Input)
    List_files_Nanog_number <- list.files(path = List_directories_Nanog[e], pattern="Nanog", recursive=TRUE)
    setwd(List_directories_Nanog[e])
    
    for (i in 1:length(List_files_Nanog_number)){
      temp_data <- read.table(List_files_Nanog_number[i],  sep=",", skip=3, header=TRUE) #each file will be read in, specify which columns you need read in to avoid any errors
      file_name <- as.character(List_files_Nanog_number[i])
      
      #Select the info inside file name about date of acquisition and identity of the nucleus being analysing. These infos are stored in columns.
      #Adapt syntax based on the number of characters
      
      Infos_file[1,1] <- substr(file_name, 1,5)
      Infos_file[1,2] <- substr(file_name, 7,8)
      Infos_file[1,3] <- substr(file_name, 10,10)
      Infos_file[1,4] <- as.numeric(gsub(".*?([0-9]+).*", "\\1", substr(file_name, 12,13)))
      
      #Repeat n nows the infos about the nuclei, as much as they are segmented volumes in the file
      Infinity <- Infos_file[rep(seq_len(nrow(Infos_file)), each = nrow(temp_data)),]
      #Infinity$V5 <- sub("^\\D*(\\d+).*$", "\\1", List_directories_Nanog[e])
      
      #Give a name to the column and make sure it looks like something
      colnames(Infinity) <- c("Date","Stage","Stack","Nucleus")
      
      #print(List_files_Nanog_number[i])
      temp_data <- cbind(Infinity, temp_data)
      
      #Assemble all data in a unique dataframe
      if(i == 1) {Nanog_statistics_all <- temp_data } else {
        if(grepl("Time", file_name)==TRUE)
        {colnames(temp_data)[5] <- "Seconds"
        colnames(temp_data)[8] <- "Time"
        Nanog_statistics_all <- Nanog_statistics_all %>% full_join(temp_data, by=c("ID","Date","Nucleus","Stack","Time","Stage","TrackID"))} 
        else {Nanog_statistics_all <- Nanog_statistics_all %>% full_join(temp_data, by=c("ID","Date","Nucleus","Stack","Time","Stage","TrackID"))}
      }
      
    }
    #Remove useless columns
    patterns <- c("Collection*", "Unit*","Image*","Channel*","X.x","X.y","Event*","Category*","Merging.foci.x*","Merging.foci.y*")
    Nanog_statistics_all <- Nanog_statistics_all[,-grep(paste(patterns, collapse="|"), colnames(Nanog_statistics_all))]
    
    if(e == 1) {Nanog_statistics_total <- Nanog_statistics_all
    } else {
      #Mix data frames with different number of columns by row
      Nanog_statistics_total <- bind_rows(Nanog_statistics_total, Nanog_statistics_all)
    }
    
  }
  
  
  #Assign one unique nucleus ID for each nucleus
  
  Unique_ID <- unique(Nanog_statistics_total[,c(1,3:4)])
  Nanog_statistics_total$Nucleus.number <- 0
  
  for(i in 1:nrow(Nanog_statistics_total))
  {
    for(e in 1:nrow(Unique_ID))
    {
      if(Nanog_statistics_total[i,1] == Unique_ID[e,1]&Nanog_statistics_total[i,3] == Unique_ID[e,2]&Nanog_statistics_total[i,4] == Unique_ID[e,3])
        
        Nanog_statistics_total[i,"Nucleus.number"] <- e
      
    }
  }
  
  #Assign one trackID for each track
  
  Unique_ID <- unique(Nanog_statistics_total[,c(1,3,4,7)])
  Nanog_statistics_total$Unique_track_Nanog_3D <- 0
  
  for(i in 1:nrow(Nanog_statistics_total))
  {
    for(e in 1:nrow(Unique_ID))
    {
      if(Nanog_statistics_total[i,1] == Unique_ID[e,1]&Nanog_statistics_total[i,3] == Unique_ID[e,2]&Nanog_statistics_total[i,4] == Unique_ID[e,3]&Nanog_statistics_total[i,7] == Unique_ID[e,4])
        
        Nanog_statistics_total[i,"Unique_track_Nanog_3D"] <- e
      
    }
  }
  
  colnames(Nanog_statistics_total)[which(names(Nanog_statistics_total) == "TrackID")] <- "TrackID_Nanog_3D"
  
}

###########Read MCP signal segmented in 2D by Damian 
{
  #MCP
  
  setwd("/Volumes/CIG/nvastenh/nanog_movement/D2c/Manuscript/240612_NCB_NSMB_submission/Material_make_figures/Fig5a, S5a/2D_Damian_Radial_distances/")
  
  #Get all the MCP file in columns for each nucleus
  
  List_files_MCP_number <- data.frame(list.files(path = ".", recursive=TRUE)) 
  
  #Assemble all data in one dataframe
  
  dataset <- data.frame()
  Infos_file <- data.frame()
  Infinity <- data.frame()
  MCP_statistics_total_2D_Damian <- data.frame()
  
  
  for (i in 1:nrow(List_files_MCP_number)) {
    
    setwd("/Volumes/CIG/nvastenh/nanog_movement/D2c/Manuscript/240612_NCB_NSMB_submission/Material_make_figures/Fig5a, S5a/2D_Damian_Radial_distances/")
    temp_data <- read.table(List_files_MCP_number[i,],  sep=",", header=TRUE) 
    #each file will be read in, specify which columns you need read in to avoid any errors
    file_name <- as.character(List_files_MCP_number[i,])
    
    if(grepl(pattern = "ims.csv", file_name)==FALSE)
    {
      if(ncol(temp_data)==18) {temp_data <- temp_data[,c(1,3,4,5,7,8,9,10,11,12,13,15,16,17,18)]}
      if(ncol(temp_data)==10) {temp_data <- temp_data[,c(1,3,4,5,7,8,9,10)]}
    }
    
    if(grepl(pattern = "ims.csv", file_name)==TRUE)
    {
      if(ncol(temp_data)==19)
      {temp_data <- temp_data[,c(1,2,5,6,7,8,9,10,11,14,15,16,17,18,19)]}
      if(ncol(temp_data)==10)
      {temp_data <- temp_data[,c(1,2,5,6,7,8,9,10)]}
    }
    
    if(ncol(temp_data)==8)
    {
      temp_data$Allele <- 1
      colnames(temp_data) <- c("Time","label_1_transcription","label_1_mean_edge_distace","label_1_std_edge_distace","label_1number_of_spots","label_1z_pos","label_1y_pos","label_1x_pos","Allele")
      
    } else if(ncol(temp_data)==15)
    {Allele_1 <- temp_data[,c(1,2:8)]
    Allele_1$Allele <- 1
    colnames(Allele_1) <- c("Time","label_1_transcription","label_1_mean_edge_distace","label_1_std_edge_distace","label_1number_of_spots","label_1z_pos","label_1y_pos","label_1x_pos","Allele")
    Allele_2 <- temp_data[,c(1,9:ncol(temp_data))]
    Allele_2$Allele <- 2
    colnames(Allele_2) <- colnames(Allele_1)
    temp_data <- rbind(Allele_1,Allele_2)
    }
    
    #Select the info inside file name about date of acquisition and identity of the nucleus being analysing. These infos are stored in columns.
    #Adapt syntax based on the number of characters
    
    Infos_file[1,1] <- substr(file_name, 1,5)
    Infos_file[1,2] <- substr(file_name, 7,8)
    Infos_file[1,3] <- as.numeric(gsub(".*?([0-9]+).*", "\\1", substr(file_name, 9,10)))
    Infos_file[1,4] <- as.numeric(gsub(".*?([0-9]+).*", "\\1", substr(file_name, 12,13)))
    
    #Repeat n nows the infos about the nuclei, as much as they are segmented volumes in the file
    Infinity <- Infos_file[rep(seq_len(nrow(Infos_file)), each = nrow(temp_data)),]
    #Infinity$V5 <- sub("^\\D*(\\d+).*$", "\\1", List_directories_MCP[e])
    
    #Give a name to the column and make sure it looks like something
    colnames(Infinity) <- c("Date","Stage","Stack","Nucleus")
    print(List_files_MCP_number[i,])
    temp_data <- cbind(Infinity, temp_data)
    
    if(i == 1) {MCP_statistics_total_2D_Damian <- temp_data} else {
      #Mix data frames with different number of columns by row
      MCP_statistics_total_2D_Damian <- rbind(MCP_statistics_total_2D_Damian, temp_data)
    }
  }
  
  #Assign one unique nucleus ID for each nucleus
  
  Unique_ID <- unique(MCP_statistics_total_2D_Damian[,c(1,3:4)])
  MCP_statistics_total_2D_Damian$Nucleus.number <- 0
  
  for(i in 1:nrow(MCP_statistics_total_2D_Damian))
  {
    for(e in 1:nrow(Unique_ID))
    {
      if(MCP_statistics_total_2D_Damian[i,1] == Unique_ID[e,1]&MCP_statistics_total_2D_Damian[i,3] == Unique_ID[e,2]&MCP_statistics_total_2D_Damian[i,4] == Unique_ID[e,3])
        
        MCP_statistics_total_2D_Damian[i,"Nucleus.number"] <- e
      
    }
  }
  
  #Assign one trackID for each track
  
  MCP_statistics_total_2D_Damian$Stack <- as.numeric(MCP_statistics_total_2D_Damian$Stack)
  MCP_statistics_total_2D_Damian$Nucleus <- as.numeric(MCP_statistics_total_2D_Damian$Nucleus)
  MCP_statistics_total_2D_Damian$Unique_track_2D_D <- 0
  
  Unique_ID <- unique(MCP_statistics_total_2D_Damian[,c("Date","Stack","Nucleus","Allele")])
  
  for(i in 1:nrow(MCP_statistics_total_2D_Damian))
  {
    for(e in 1:nrow(Unique_ID))
    {
      if(MCP_statistics_total_2D_Damian[i,"Date"] == Unique_ID[e,1]&MCP_statistics_total_2D_Damian[i,"Stack"] == Unique_ID[e,2]&MCP_statistics_total_2D_Damian[i,"Nucleus"] == Unique_ID[e,3]&MCP_statistics_total_2D_Damian[i,"Allele"] == Unique_ID[e,4])
        
        MCP_statistics_total_2D_Damian[i,"Unique_track_2D_D"] <- e
      
    }
  }
  
  colnames(MCP_statistics_total_2D_Damian)[which(names(MCP_statistics_total_2D_Damian) == "Allele")] <- "TrackID_MCP_2D_D"
  
  MCP_statistics_total_2D_Damian$Time <- MCP_statistics_total_2D_Damian$Time+1 
  
}

###########Nucleus
{
  #All file names in one table
  {
    #Nucleus
    #Get a list of files containing the word "Nucleus" to analyze only the Nucleus foci
    
    setwd(Input)
    
    List_directories_Nucleus <- grep(pattern="Nucleus", list.dirs(path = ".", full.names = TRUE), value=TRUE)
    
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
      
      Infos_file[1,1] <- substr(file_name, 1,5)
      Infos_file[1,2] <- substr(file_name, 7,8)
      Infos_file[1,3] <- substr(file_name, 10,10)
      Infos_file[1,4] <- as.numeric(gsub(".*?([0-9]+).*", "\\1", substr(file_name, 12,13)))
      
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
}

#=========================================================
# Assemble Nanog with MOVIE signal, center on transcription initiation and normalise fluorescence per time point
#=========================================================

############Associate Nanog with MOVIE
{
  #Calculate the distance between each pair of MOVIE-Nanog foci, and check if this distance if < 0.2. If yes, join the two data frames at this specific line. If not, write just nothing.
  
  #Add all columns from MOVIE file into Nanog file
  
  colnames(MOVIE_statistics_total) <- paste0("MOVIE_",colnames(MOVIE_statistics_total))
  
  Nanog_statistics_total[c(colnames(MOVIE_statistics_total))] <- NA
  Nanog_statistics_total$Dist.Nanog.MOVIE <- NA
  
  which(colnames(Nanog_statistics_total)=="MOVIE_Date")
  max(col(Nanog_statistics_total))
  
  for(r in 1:max(as.numeric((MOVIE_statistics_total$MOVIE_Nucleus.number)))) {
    
    if(nrow(MOVIE_statistics_total[(which(MOVIE_statistics_total$MOVIE_Nucleus.number==r)),])!=0) {
      
      for(w in min(MOVIE_statistics_total[(which(MOVIE_statistics_total$MOVIE_Nucleus.number==r)),]$MOVIE_Time):max(MOVIE_statistics_total[(which(MOVIE_statistics_total$MOVIE_Nucleus.number==r)),]$MOVIE_Time)) {
        
        Nanog_subset_time <- Nanog_statistics_total[(which(Nanog_statistics_total$Time==w&Nanog_statistics_total$Nucleus.number==r)),]
        MOVIE_subset_time <- MOVIE_statistics_total[(which(MOVIE_statistics_total$MOVIE_Time==w&MOVIE_statistics_total$MOVIE_Nucleus.number==r)),]
        #print(w)
        
        if(nrow(Nanog_subset_time)!=0) {
          
          for(a in 1:nrow(MOVIE_subset_time)) {
            
            MOVIE_x <- MOVIE_subset_time$MOVIE_Position.X[a]
            MOVIE_y <- MOVIE_subset_time$MOVIE_Position.Y[a]
            MOVIE_z <- MOVIE_subset_time$MOVIE_Position.Z[a]
            
            if(nrow(MOVIE_subset_time) != 0) {
              
              #print(MOVIE_statistics_total[a,])
              
              for(b in 1:nrow(Nanog_subset_time)) {
                
                Nanog_x <-  Nanog_subset_time$Position.X[b]
                Nanog_y <- Nanog_subset_time$Position.Y[b]
                Nanog_z <- Nanog_subset_time$Position.Z[b]
                
                Distance_Nanog_MOVIE <- sqrt((Nanog_x - MOVIE_x)^2 + (Nanog_y - MOVIE_y)^2 + (Nanog_z - MOVIE_z)^2)
                
                if(Distance_Nanog_MOVIE < 1.0 & Nanog_subset_time$Shortest.Distance.to.Surfaces.y[b] < 1.0) {
                  
                  Nanog_statistics_total[(which(Nanog_statistics_total$Time==w&Nanog_statistics_total$Nucleus.number==r)),55:123][b,] <- MOVIE_statistics_total[(which(MOVIE_statistics_total$MOVIE_Time==w&MOVIE_statistics_total$MOVIE_Nucleus.number==r)),][a,]
                  Nanog_statistics_total[(which(Nanog_statistics_total$Time==w&Nanog_statistics_total$Nucleus.number==r)),124][b] <- Distance_Nanog_MOVIE
                  
                }
              }
            }
          }
        }
      }
    }
  }
  
  #Check if it is correct
  
  MOVIE <- data.frame(matrix(nrow = 86,ncol=5))
  k=1
  for(i in 1:max(MOVIE_statistics_total$MOVIE_Nucleus.number)){
    for(y in unique(na.omit(MOVIE_statistics_total[MOVIE_statistics_total$MOVIE_Nucleus.number==i,"MOVIE_TrackID_MOVIE_3D"])))
    {MOVIE[k,] <- MOVIE_statistics_total[which(MOVIE_statistics_total$MOVIE_Nucleus.number==i&MOVIE_statistics_total$MOVIE_TrackID_MOVIE_3D==y),c("MOVIE_Date","MOVIE_Stack","MOVIE_Nucleus","MOVIE_TrackID_MOVIE_3D","MOVIE_Time")][1,]
    k=k+1}}
  
}

############Nanog find the minimum time for which transcription start
{
  Nanog_statistics_total$Time <- as.numeric(Nanog_statistics_total$Time)
  Nanog_statistics_total$Seconds <- as.numeric(Nanog_statistics_total$Seconds)
  Nanog_statistics_total$Unique_track <- as.numeric(Nanog_statistics_total$Unique_track)
  
  for(i in 1:max(na.omit(Nanog_statistics_total$Unique_track)))
  {
    Nanog_statistics_total[which(Nanog_statistics_total$Unique_track==i),"Time_activation"] <- min(Nanog_statistics_total[!is.na(Nanog_statistics_total$MOVIE_Nucleus.number)&Nanog_statistics_total$Unique_track==i,"Time"])
    #print(i)
  }
  
  Nanog_statistics_total$Time_normalized <- Nanog_statistics_total$Time - Nanog_statistics_total$Time_activation
  Nanog_statistics_total$Unique_track <- as.numeric(Nanog_statistics_total$Unique_track)
  
  for(i in 1:max(na.omit(Nanog_statistics_total$Unique_track)))
  {
    Nanog_statistics_total[which(Nanog_statistics_total$Unique_track==i),"Seconds_activation"] <- min(Nanog_statistics_total[!is.na(Nanog_statistics_total$MOVIE_Nucleus.number)&Nanog_statistics_total$Unique_track==i,"Seconds"])
  }
  
  Nanog_statistics_total$Seconds_normalized <- Nanog_statistics_total$Seconds - Nanog_statistics_total$Seconds_activation 
  
  
  
}

############Determine the phenotype for each track for Nanog
{
  
  Nanog_statistics_total$Time <- as.numeric(Nanog_statistics_total$Time)
  Nanog_statistics_total$Unique_track_Nanog_3D <- as.numeric(Nanog_statistics_total$Unique_track_Nanog_3D)
  
  #DETERMINE THE PHENOTYPE FOR THE DIFFERENT TRACKS AND THE MERGING TIME
  
  seqle <- function(x,incr=1) { 
    if(!is.numeric(x)) x <- as.numeric(x) 
    n <- length(x)  
    y <- x[-1L] != x[-n] + incr 
    i <- c(which(y|is.na(y)),n) 
    list(lengths = diff(c(0L,i)),
         values = x[head(c(0L,i)+1L,-1L)]) 
  } 
  
  
  Nanog_statistics_total$Phenotype <- "NA" 
  Nanog_statistics_total$Fused_activation <- "NA" 
  
  #1.Create a table with the number of spots per time point for each track
  for(i in sort(unique(Nanog_statistics_total$Unique_track_Nanog_3D)))
  {
    #Determine the time of activation
    Time_activation<-unique(Nanog_statistics_total[which(Nanog_statistics_total$Unique_track_Nanog_3D==i),"Time_activation"])
    
    if(!Time_activation=="Inf"){
      
      #Determine the time of the track before transcription start
      Nanog_Unique_track_Nanog_3D_activation <- Nanog_statistics_total[which(Nanog_statistics_total$Unique_track_Nanog_3D==i&Nanog_statistics_total$Time %in% seq(Time_activation-10,Time_activation)),]
      
      
      Table <- as.data.frame(matrix(nrow = 2,ncol = length(unique(Nanog_Unique_track_Nanog_3D_activation$Time))))
      
      k=1
      
      for(j in min(Nanog_Unique_track_Nanog_3D_activation$Time):max(Nanog_Unique_track_Nanog_3D_activation$Time))
      {
        Table[1,k] <- j
        Table[2,k] <- nrow(Nanog_Unique_track_Nanog_3D_activation[Nanog_Unique_track_Nanog_3D_activation$Time==j,])
        k=k+1
      }
      
      #2. Determine if there is one or two spots at the moment of activation
      if(nrow(Nanog_statistics_total[which(Nanog_statistics_total$Unique_track_Nanog_3D==i&Nanog_statistics_total$Time==Time_activation),])>=2) {
        Nanog_statistics_total[Nanog_statistics_total$Unique_track_Nanog_3D==i&!is.na(Nanog_statistics_total$Unique_track_Nanog_3D==i),"Fused_activation"] <- "No"
        print(nrow(Nanog_statistics_total[which(Nanog_statistics_total$Unique_track_Nanog_3D==i&Nanog_statistics_total$Time==Time_activation),]))
        
        Nanog_Unique_track_Nanog_3D_activation <- Nanog_statistics_total[which(Nanog_statistics_total$Unique_track_Nanog_3D==i&Nanog_statistics_total$Time %in% seq(Time_activation-3,Time_activation+3)),]
        Table1 <- as.data.frame(matrix(nrow = 2,ncol = length(unique(Nanog_Unique_track_Nanog_3D_activation$Time))))
        
        k=1
        for(j in min(Nanog_Unique_track_Nanog_3D_activation$Time):max(Nanog_Unique_track_Nanog_3D_activation$Time))
        {
          Table1[1,k] <- j
          Table1[2,k] <- nrow(Nanog_Unique_track_Nanog_3D_activation[Nanog_Unique_track_Nanog_3D_activation$Time==j,])
          k=k+1
        }
        
        if(all(Table1[2,]>=2))
        {Nanog_statistics_total[Nanog_statistics_total$Unique_track_Nanog_3D==i&!is.na(Nanog_statistics_total$Unique_track_Nanog_3D==i),"Phenotype"] <- "Several_clusters_never_merging"
        Nanog_statistics_total[Nanog_statistics_total$Unique_track_Nanog_3D==i&!is.na(Nanog_statistics_total$Unique_track_Nanog_3D==i),"Time_fusion"] <- "None"
        Nanog_statistics_total[Nanog_statistics_total$Unique_track_Nanog_3D==i&!is.na(Nanog_statistics_total$Unique_track_Nanog_3D==i),"Length_double"] <- max(seqle(Table[1,which(Table[2,]>=2)])$lengths)
        }
        if(any(Table1[2,]<2))
        {Nanog_statistics_total[Nanog_statistics_total$Unique_track_Nanog_3D==i&!is.na(Nanog_statistics_total$Unique_track_Nanog_3D==i),"Time_fusion"] <- "To determine"
        if(any(seqle(Table[1,which(Table[2,]>=2)])$lengths>=2))
        {Nanog_statistics_total[Nanog_statistics_total$Unique_track_Nanog_3D==i&!is.na(Nanog_statistics_total$Unique_track_Nanog_3D==i),"Phenotype"] <- "Several_clusters_long"
        Nanog_statistics_total[Nanog_statistics_total$Unique_track_Nanog_3D==i&!is.na(Nanog_statistics_total$Unique_track_Nanog_3D==i),"Length_double"] <- max(seqle(Table[1,which(Table[2,]>=2)])$lengths) 
        } else {Nanog_statistics_total[Nanog_statistics_total$Unique_track_Nanog_3D==i&!is.na(Nanog_statistics_total$Unique_track_Nanog_3D==i),"Phenotype"] <- "Several_clusters_short"
        Nanog_statistics_total[Nanog_statistics_total$Unique_track_Nanog_3D==i&!is.na(Nanog_statistics_total$Unique_track_Nanog_3D==i),"Length_double"] <- max(seqle(Table[1,which(Table[2,]>=2)])$lengths)}
        
        }
      }
    }
    
    #3.Determine how many time points in a row have 2 or more spots
    
    if(!nrow(Nanog_statistics_total[which(Nanog_statistics_total$Unique_track_Nanog_3D==i&Nanog_statistics_total$Time==Time_activation),])>=2)
    {Nanog_statistics_total[Nanog_statistics_total$Unique_track_Nanog_3D==i&!is.na(Nanog_statistics_total$Unique_track_Nanog_3D==i),"Fused_activation"] <- "Yes"
    if(!any(Table[2,]>2))
    {
      Nanog_statistics_total[Nanog_statistics_total$Unique_track_Nanog_3D==i&!is.na(Nanog_statistics_total$Unique_track_Nanog_3D==i),"Phenotype"] <- "One_cluster"
      Nanog_statistics_total[Nanog_statistics_total$Unique_track_Nanog_3D==i&!is.na(Nanog_statistics_total$Unique_track_Nanog_3D==i),"Length_double"] <- "0"
      Nanog_statistics_total[Nanog_statistics_total$Unique_track_Nanog_3D==i&!is.na(Nanog_statistics_total$Unique_track_Nanog_3D==i),"Time_fusion"] <- "None"}
    if(any(Table[2,]>=2)) {
      if(any(seqle(Table[1,which(Table[2,]>=2)])$lengths>=2)) {
        Nanog_statistics_total[Nanog_statistics_total$Unique_track_Nanog_3D==i&!is.na(Nanog_statistics_total$Unique_track_Nanog_3D==i),"Phenotype"] <- "Several_clusters_long"
        Nanog_statistics_total[Nanog_statistics_total$Unique_track_Nanog_3D==i&!is.na(Nanog_statistics_total$Unique_track_Nanog_3D==i),"Length_double"] <- max(seqle(Table[1,which(Table[2,]>=2)])$lengths)
        Nanog_statistics_total[Nanog_statistics_total$Unique_track_Nanog_3D==i&!is.na(Nanog_statistics_total$Unique_track_Nanog_3D==i),"Time_fusion"] <- 
          max(Table[1,which(Table[2,]>=2)])+1
      } else {Nanog_statistics_total[Nanog_statistics_total$Unique_track_Nanog_3D==i&!is.na(Nanog_statistics_total$Unique_track_Nanog_3D==i),"Phenotype"] <- "Several_clusters_short"
      Nanog_statistics_total[Nanog_statistics_total$Unique_track_Nanog_3D==i&!is.na(Nanog_statistics_total$Unique_track_Nanog_3D==i),"Length_double"] <- max(seqle(Table[1,which(Table[2,]>=2)])$lengths)
      Nanog_statistics_total[Nanog_statistics_total$Unique_track_Nanog_3D==i&!is.na(Nanog_statistics_total$Unique_track_Nanog_3D==i),"Time_fusion"] <- max(Table[1,which(Table[2,]>=2)])+1}
      
    }}}
  
  
  #4. Normalize the timing for Time_fusion
  #Determine the time of fusion for specific cases
  unique(Nanog_statistics_total[Nanog_statistics_total$Time_fusion=="To determine","Unique_track_Nanog_3D"])
  Nanog_statistics_total[Nanog_statistics_total$Unique_track_Nanog_3D==69,c("Time","Time_activation")]
  #2 14 24 35 45 52 59 61 67 69 72
  
  Nanog_statistics_total[Nanog_statistics_total$Unique_track_Nanog_3D==2&!is.na(Nanog_statistics_total$Unique_track_Nanog_3D==2),"Time_fusion"] <- 11
  Nanog_statistics_total[Nanog_statistics_total$Unique_track_Nanog_3D==14&!is.na(Nanog_statistics_total$Unique_track_Nanog_3D==14),"Time_fusion"] <- 13
  Nanog_statistics_total[Nanog_statistics_total$Unique_track_Nanog_3D==24&!is.na(Nanog_statistics_total$Unique_track_Nanog_3D==24),"Time_fusion"] <- 13
  Nanog_statistics_total[Nanog_statistics_total$Unique_track_Nanog_3D==35&!is.na(Nanog_statistics_total$Unique_track_Nanog_3D==35),"Time_fusion"] <- 18
  Nanog_statistics_total[Nanog_statistics_total$Unique_track_Nanog_3D==45&!is.na(Nanog_statistics_total$Unique_track_Nanog_3D==45),"Time_fusion"] <- 8
  Nanog_statistics_total[Nanog_statistics_total$Unique_track_Nanog_3D==52&!is.na(Nanog_statistics_total$Unique_track_Nanog_3D==52),"Time_fusion"] <- 14
  Nanog_statistics_total[Nanog_statistics_total$Unique_track_Nanog_3D==59&!is.na(Nanog_statistics_total$Unique_track_Nanog_3D==59),"Time_fusion"] <- 3
  Nanog_statistics_total[Nanog_statistics_total$Unique_track_Nanog_3D==61&!is.na(Nanog_statistics_total$Unique_track_Nanog_3D==61),"Time_fusion"] <- NA
  Nanog_statistics_total[Nanog_statistics_total$Unique_track_Nanog_3D==67&!is.na(Nanog_statistics_total$Unique_track_Nanog_3D==67),"Time_fusion"] <- 11
  Nanog_statistics_total[Nanog_statistics_total$Unique_track_Nanog_3D==69&!is.na(Nanog_statistics_total$Unique_track_Nanog_3D==69),"Time_fusion"] <- NA
  Nanog_statistics_total[Nanog_statistics_total$Unique_track_Nanog_3D==72&!is.na(Nanog_statistics_total$Unique_track_Nanog_3D==72),"Time_fusion"] <- 13
  
  Nanog_statistics_total$Time_fusion <- as.numeric(Nanog_statistics_total$Time_fusion)
  Nanog_statistics_total$Time_fusion_normalized <- Nanog_statistics_total$Time - Nanog_statistics_total$Time_fusion
  Nanog_statistics_total$Unique_track_Nanog_3D <- as.numeric(Nanog_statistics_total$Unique_track_Nanog_3D)
  Nanog_statistics_total$Seconds_fusion_normalized <- Nanog_statistics_total$Seconds - Nanog_statistics_total$Time_fusion*15 
  
  #Associate all the clusters from the same phenotype together
  Nanog_statistics_total$Phenotype_simple <- "One_cluster"
  Nanog_statistics_total[Nanog_statistics_total$Phenotype %in% c("Several_clusters_long","Several_clusters_short","Several_clusters_never_merging"),"Phenotype_simple"] <- "More_cluster"
  
  
}

############Count Nanog cluster per time point + generate position for aligned position in graph
{
  
  #Make the tracks align on one place
  {
    
    
    Nanog_statistics_total <- Nanog_statistics_total[with(Nanog_statistics_total,order(Unique_track_Nanog_3D,Time,TrackID_Nanog_3D)),]
    Nanog_statistics_total$Track_plotting <-0
    
    for(i in 1:max(as.numeric(Nanog_statistics_total$Unique_track_Nanog_3D)))
    {
      for(g in 1:max(Nanog_statistics_total[which(Nanog_statistics_total$Unique_track_Nanog_3D==i),"Time"]))
      {
        print(nrow(Nanog_statistics_total[which(Nanog_statistics_total$Unique_track_Nanog_3D==i&Nanog_statistics_total$Time==g),]))
        
        if(nrow(Nanog_statistics_total[which(Nanog_statistics_total$Unique_track_Nanog_3D==i&Nanog_statistics_total$Time==g),])==0)
        {print(0)}
        else if(nrow(Nanog_statistics_total[which(Nanog_statistics_total$Unique_track_Nanog_3D==i&Nanog_statistics_total$Time==g),])==1)
        {Nanog_statistics_total[which(Nanog_statistics_total$Unique_track_Nanog_3D==i&Nanog_statistics_total$Time==g),"Track_plotting"][1] <- 0.5
        Nanog_statistics_total[which(Nanog_statistics_total$Unique_track_Nanog_3D==i&Nanog_statistics_total$Time==g),"Spot_number"][1] <- 1
        Nanog_statistics_total[which(Nanog_statistics_total$Unique_track_Nanog_3D==i&Nanog_statistics_total$Time==g),"Color"][1] <- "a"
        }
        
        else if(nrow(Nanog_statistics_total[which(Nanog_statistics_total$Unique_track_Nanog_3D==i&Nanog_statistics_total$Time==g),])==2)
        {Nanog_statistics_total[which(Nanog_statistics_total$Unique_track_Nanog_3D==i&Nanog_statistics_total$Time==g),"Track_plotting"][1] <- 0.25
        Nanog_statistics_total[which(Nanog_statistics_total$Unique_track_Nanog_3D==i&Nanog_statistics_total$Time==g),"Track_plotting"][2] <- 0.75
        Nanog_statistics_total[which(Nanog_statistics_total$Unique_track_Nanog_3D==i&Nanog_statistics_total$Time==g),"Spot_number"][1] <- 2
        Nanog_statistics_total[which(Nanog_statistics_total$Unique_track_Nanog_3D==i&Nanog_statistics_total$Time==g),"Spot_number"][2] <- 2
        Nanog_statistics_total[which(Nanog_statistics_total$Unique_track_Nanog_3D==i&Nanog_statistics_total$Time==g),"Color"][1] <- "b"
        Nanog_statistics_total[which(Nanog_statistics_total$Unique_track_Nanog_3D==i&Nanog_statistics_total$Time==g),"Color"][2] <- "b"
        
        
        
        }
        else if(nrow(Nanog_statistics_total[which(Nanog_statistics_total$Unique_track_Nanog_3D==i&Nanog_statistics_total$Time==g),])==3)
        {Nanog_statistics_total[which(Nanog_statistics_total$Unique_track_Nanog_3D==i&Nanog_statistics_total$Time==g),"Track_plotting"][1] <- 0
        Nanog_statistics_total[which(Nanog_statistics_total$Unique_track_Nanog_3D==i&Nanog_statistics_total$Time==g),"Track_plotting"][2] <- 0.5
        Nanog_statistics_total[which(Nanog_statistics_total$Unique_track_Nanog_3D==i&Nanog_statistics_total$Time==g),"Track_plotting"][3] <- 1
        Nanog_statistics_total[which(Nanog_statistics_total$Unique_track_Nanog_3D==i&Nanog_statistics_total$Time==g),"Spot_number"][1] <- 3
        Nanog_statistics_total[which(Nanog_statistics_total$Unique_track_Nanog_3D==i&Nanog_statistics_total$Time==g),"Spot_number"][2] <- 3
        Nanog_statistics_total[which(Nanog_statistics_total$Unique_track_Nanog_3D==i&Nanog_statistics_total$Time==g),"Spot_number"][3] <- 3
        Nanog_statistics_total[which(Nanog_statistics_total$Unique_track_Nanog_3D==i&Nanog_statistics_total$Time==g),"Color"][1] <- "c"
        Nanog_statistics_total[which(Nanog_statistics_total$Unique_track_Nanog_3D==i&Nanog_statistics_total$Time==g),"Color"][2] <- "c"
        Nanog_statistics_total[which(Nanog_statistics_total$Unique_track_Nanog_3D==i&Nanog_statistics_total$Time==g),"Color"][3] <- "c"
        
        
        
        
        }
        else if(nrow(Nanog_statistics_total[which(Nanog_statistics_total$Unique_track_Nanog_3D==i&Nanog_statistics_total$Time==g),])==4)
        {Nanog_statistics_total[which(Nanog_statistics_total$Unique_track_Nanog_3D==i&Nanog_statistics_total$Time==g),"Track_plotting"][1] <- 0
        Nanog_statistics_total[which(Nanog_statistics_total$Unique_track_Nanog_3D==i&Nanog_statistics_total$Time==g),"Track_plotting"][2] <- 0.5
        Nanog_statistics_total[which(Nanog_statistics_total$Unique_track_Nanog_3D==i&Nanog_statistics_total$Time==g),"Track_plotting"][3] <- 1
        Nanog_statistics_total[which(Nanog_statistics_total$Unique_track_Nanog_3D==i&Nanog_statistics_total$Time==g),"Track_plotting"][4] <- 1.5
        Nanog_statistics_total[which(Nanog_statistics_total$Unique_track_Nanog_3D==i&Nanog_statistics_total$Time==g),"Spot_number"][1] <- 4
        Nanog_statistics_total[which(Nanog_statistics_total$Unique_track_Nanog_3D==i&Nanog_statistics_total$Time==g),"Spot_number"][2] <- 4
        Nanog_statistics_total[which(Nanog_statistics_total$Unique_track_Nanog_3D==i&Nanog_statistics_total$Time==g),"Spot_number"][3] <- 4
        Nanog_statistics_total[which(Nanog_statistics_total$Unique_track_Nanog_3D==i&Nanog_statistics_total$Time==g),"Spot_number"][4] <- 4
        Nanog_statistics_total[which(Nanog_statistics_total$Unique_track_Nanog_3D==i&Nanog_statistics_total$Time==g),"Color"][1] <- "d"
        Nanog_statistics_total[which(Nanog_statistics_total$Unique_track_Nanog_3D==i&Nanog_statistics_total$Time==g),"Color"][2] <- "d"
        Nanog_statistics_total[which(Nanog_statistics_total$Unique_track_Nanog_3D==i&Nanog_statistics_total$Time==g),"Color"][3] <- "d"
        Nanog_statistics_total[which(Nanog_statistics_total$Unique_track_Nanog_3D==i&Nanog_statistics_total$Time==g),"Color"][4] <- "d"
        }
        else if(nrow(Nanog_statistics_total[which(Nanog_statistics_total$Unique_track_Nanog_3D==i&Nanog_statistics_total$Time==g),])==5)
        {Nanog_statistics_total[which(Nanog_statistics_total$Unique_track_Nanog_3D==i&Nanog_statistics_total$Time==g),"Track_plotting"][1] <- -0.5
        Nanog_statistics_total[which(Nanog_statistics_total$Unique_track_Nanog_3D==i&Nanog_statistics_total$Time==g),"Track_plotting"][2] <- 0
        Nanog_statistics_total[which(Nanog_statistics_total$Unique_track_Nanog_3D==i&Nanog_statistics_total$Time==g),"Track_plotting"][3] <- 0.5
        Nanog_statistics_total[which(Nanog_statistics_total$Unique_track_Nanog_3D==i&Nanog_statistics_total$Time==g),"Track_plotting"][4] <- 1
        Nanog_statistics_total[which(Nanog_statistics_total$Unique_track_Nanog_3D==i&Nanog_statistics_total$Time==g),"Track_plotting"][5] <- 1.5
        Nanog_statistics_total[which(Nanog_statistics_total$Unique_track_Nanog_3D==i&Nanog_statistics_total$Time==g),"Spot_number"][1] <- 5
        Nanog_statistics_total[which(Nanog_statistics_total$Unique_track_Nanog_3D==i&Nanog_statistics_total$Time==g),"Spot_number"][2] <- 5
        Nanog_statistics_total[which(Nanog_statistics_total$Unique_track_Nanog_3D==i&Nanog_statistics_total$Time==g),"Spot_number"][3] <- 5
        Nanog_statistics_total[which(Nanog_statistics_total$Unique_track_Nanog_3D==i&Nanog_statistics_total$Time==g),"Spot_number"][4] <- 5
        Nanog_statistics_total[which(Nanog_statistics_total$Unique_track_Nanog_3D==i&Nanog_statistics_total$Time==g),"Spot_number"][5] <- 5
        Nanog_statistics_total[which(Nanog_statistics_total$Unique_track_Nanog_3D==i&Nanog_statistics_total$Time==g),"Color"][1] <- "e"
        Nanog_statistics_total[which(Nanog_statistics_total$Unique_track_Nanog_3D==i&Nanog_statistics_total$Time==g),"Color"][2] <- "e"
        Nanog_statistics_total[which(Nanog_statistics_total$Unique_track_Nanog_3D==i&Nanog_statistics_total$Time==g),"Color"][3] <- "e"
        Nanog_statistics_total[which(Nanog_statistics_total$Unique_track_Nanog_3D==i&Nanog_statistics_total$Time==g),"Color"][4] <- "e"
        Nanog_statistics_total[which(Nanog_statistics_total$Unique_track_Nanog_3D==i&Nanog_statistics_total$Time==g),"Color"][5] <- "e"
        }
        else if(nrow(Nanog_statistics_total[which(Nanog_statistics_total$Unique_track_Nanog_3D==i&Nanog_statistics_total$Time==g),])==6)
        {Nanog_statistics_total[which(Nanog_statistics_total$Unique_track_Nanog_3D==i&Nanog_statistics_total$Time==g),"Track_plotting"][1] <- -1
        Nanog_statistics_total[which(Nanog_statistics_total$Unique_track_Nanog_3D==i&Nanog_statistics_total$Time==g),"Track_plotting"][2] <- -0.5
        Nanog_statistics_total[which(Nanog_statistics_total$Unique_track_Nanog_3D==i&Nanog_statistics_total$Time==g),"Track_plotting"][3] <- 0
        Nanog_statistics_total[which(Nanog_statistics_total$Unique_track_Nanog_3D==i&Nanog_statistics_total$Time==g),"Track_plotting"][4] <- 0.5
        Nanog_statistics_total[which(Nanog_statistics_total$Unique_track_Nanog_3D==i&Nanog_statistics_total$Time==g),"Track_plotting"][5] <- 1
        Nanog_statistics_total[which(Nanog_statistics_total$Unique_track_Nanog_3D==i&Nanog_statistics_total$Time==g),"Track_plotting"][6] <- 1.5
        Nanog_statistics_total[which(Nanog_statistics_total$Unique_track_Nanog_3D==i&Nanog_statistics_total$Time==g),"Spot_number"][1] <- 6
        Nanog_statistics_total[which(Nanog_statistics_total$Unique_track_Nanog_3D==i&Nanog_statistics_total$Time==g),"Spot_number"][2] <- 6
        Nanog_statistics_total[which(Nanog_statistics_total$Unique_track_Nanog_3D==i&Nanog_statistics_total$Time==g),"Spot_number"][3] <- 6
        Nanog_statistics_total[which(Nanog_statistics_total$Unique_track_Nanog_3D==i&Nanog_statistics_total$Time==g),"Spot_number"][4] <- 6
        Nanog_statistics_total[which(Nanog_statistics_total$Unique_track_Nanog_3D==i&Nanog_statistics_total$Time==g),"Spot_number"][5] <- 6
        Nanog_statistics_total[which(Nanog_statistics_total$Unique_track_Nanog_3D==i&Nanog_statistics_total$Time==g),"Spot_number"][6] <- 6
        Nanog_statistics_total[which(Nanog_statistics_total$Unique_track_Nanog_3D==i&Nanog_statistics_total$Time==g),"Color"][1] <- "f"
        Nanog_statistics_total[which(Nanog_statistics_total$Unique_track_Nanog_3D==i&Nanog_statistics_total$Time==g),"Color"][2] <- "f"
        Nanog_statistics_total[which(Nanog_statistics_total$Unique_track_Nanog_3D==i&Nanog_statistics_total$Time==g),"Color"][3] <- "f"
        Nanog_statistics_total[which(Nanog_statistics_total$Unique_track_Nanog_3D==i&Nanog_statistics_total$Time==g),"Color"][4] <- "f"
        Nanog_statistics_total[which(Nanog_statistics_total$Unique_track_Nanog_3D==i&Nanog_statistics_total$Time==g),"Color"][5] <- "f"
        Nanog_statistics_total[which(Nanog_statistics_total$Unique_track_Nanog_3D==i&Nanog_statistics_total$Time==g),"Color"][6] <- "f"
        }
        else if(nrow(Nanog_statistics_total[which(Nanog_statistics_total$Unique_track_Nanog_3D==i&Nanog_statistics_total$Time==g),])==7)
        {Nanog_statistics_total[which(Nanog_statistics_total$Unique_track_Nanog_3D==i&Nanog_statistics_total$Time==g),"Track_plotting"][1] <- -1.5
        Nanog_statistics_total[which(Nanog_statistics_total$Unique_track_Nanog_3D==i&Nanog_statistics_total$Time==g),"Track_plotting"][2] <- -1
        Nanog_statistics_total[which(Nanog_statistics_total$Unique_track_Nanog_3D==i&Nanog_statistics_total$Time==g),"Track_plotting"][3] <- -0.5
        Nanog_statistics_total[which(Nanog_statistics_total$Unique_track_Nanog_3D==i&Nanog_statistics_total$Time==g),"Track_plotting"][4] <- 0
        Nanog_statistics_total[which(Nanog_statistics_total$Unique_track_Nanog_3D==i&Nanog_statistics_total$Time==g),"Track_plotting"][5] <- 1
        Nanog_statistics_total[which(Nanog_statistics_total$Unique_track_Nanog_3D==i&Nanog_statistics_total$Time==g),"Track_plotting"][6] <- 1.5
        Nanog_statistics_total[which(Nanog_statistics_total$Unique_track_Nanog_3D==i&Nanog_statistics_total$Time==g),"Track_plotting"][7] <- 2
        Nanog_statistics_total[which(Nanog_statistics_total$Unique_track_Nanog_3D==i&Nanog_statistics_total$Time==g),"Spot_number"][1] <- 7
        Nanog_statistics_total[which(Nanog_statistics_total$Unique_track_Nanog_3D==i&Nanog_statistics_total$Time==g),"Spot_number"][2] <- 7
        Nanog_statistics_total[which(Nanog_statistics_total$Unique_track_Nanog_3D==i&Nanog_statistics_total$Time==g),"Spot_number"][3] <- 7
        Nanog_statistics_total[which(Nanog_statistics_total$Unique_track_Nanog_3D==i&Nanog_statistics_total$Time==g),"Spot_number"][4] <- 7
        Nanog_statistics_total[which(Nanog_statistics_total$Unique_track_Nanog_3D==i&Nanog_statistics_total$Time==g),"Spot_number"][5] <- 7
        Nanog_statistics_total[which(Nanog_statistics_total$Unique_track_Nanog_3D==i&Nanog_statistics_total$Time==g),"Spot_number"][6] <- 7
        Nanog_statistics_total[which(Nanog_statistics_total$Unique_track_Nanog_3D==i&Nanog_statistics_total$Time==g),"Spot_number"][7] <- 7
        Nanog_statistics_total[which(Nanog_statistics_total$Unique_track_Nanog_3D==i&Nanog_statistics_total$Time==g),"Color"][1] <- "f"
        Nanog_statistics_total[which(Nanog_statistics_total$Unique_track_Nanog_3D==i&Nanog_statistics_total$Time==g),"Color"][2] <- "g"
        Nanog_statistics_total[which(Nanog_statistics_total$Unique_track_Nanog_3D==i&Nanog_statistics_total$Time==g),"Color"][3] <- "g"
        Nanog_statistics_total[which(Nanog_statistics_total$Unique_track_Nanog_3D==i&Nanog_statistics_total$Time==g),"Color"][4] <- "g"
        Nanog_statistics_total[which(Nanog_statistics_total$Unique_track_Nanog_3D==i&Nanog_statistics_total$Time==g),"Color"][5] <- "g"
        Nanog_statistics_total[which(Nanog_statistics_total$Unique_track_Nanog_3D==i&Nanog_statistics_total$Time==g),"Color"][6] <- "g"
        Nanog_statistics_total[which(Nanog_statistics_total$Unique_track_Nanog_3D==i&Nanog_statistics_total$Time==g),"Color"][7] <- "g"
        }
      }
    }
    
  }
  
  
}

############Calculate the distance between Nanog clusters + normalisation of fluorescence
{
  
  Distance_Nanog_dataframe <- data.frame(matrix(ncol=ncol(Nanog_statistics_total)+1))
  colnames(Distance_Nanog_dataframe) <- colnames(Nanog_statistics_total)
  
  for(i in 1:max(Nanog_statistics_total$Unique_track_Nanog_3D))
  {
    for(j in unique(Nanog_statistics_total[Nanog_statistics_total$Unique_track_Nanog_3D==i,"Time"]))
    {
      
      Nanog_statistics_total[Nanog_statistics_total$Unique_track_Nanog_3D==i&Nanog_statistics_total$Time==j,"Spot_number_Nanog"] <- nrow(Nanog_statistics_total[Nanog_statistics_total$Unique_track_Nanog_3D==i&Nanog_statistics_total$Time==j,])
      
    }
  }
  
  k=1
  
  for(r in 1:max(as.numeric((Nanog_statistics_total$Unique_track_Nanog_3D)))) {
    
    if(nrow(Nanog_statistics_total[(which(Nanog_statistics_total$Unique_track_Nanog_3D==r)),])!=0) {
      
      for(w in min(Nanog_statistics_total[(which(Nanog_statistics_total$Unique_track_Nanog_3D==r)),]$Time):max(Nanog_statistics_total[(which(Nanog_statistics_total$Unique_track==r)),]$Time)) {
        
        Nanog_subset_time <- Nanog_statistics_total[(which(Nanog_statistics_total$Time==w&Nanog_statistics_total$Unique_track_Nanog_3D==r)),]
        
        if(nrow(Nanog_subset_time)!=0) {
          
          if(nrow(Nanog_subset_time)==1)
          {Distance_Nanog_dataframe[k,] <- Nanog_subset_time[1,]
          Distance_Nanog_dataframe$Distance_clusters[k] <- 0
          Distance_Nanog_dataframe[k,"Intensity.Sum_sum"] <- Nanog_subset_time[1,"Intensity.Sum"]
          Distance_Nanog_dataframe[k,"Volume_sum"] <- Nanog_subset_time[1,"Volume"]
          
          k=k+1}
          else if(nrow(Nanog_subset_time)==2)
          {Distance_Nanog_dataframe[k,] <- Nanog_subset_time[1,]
          
          Distance_Nanog_dataframe[k,"Intensity.Sum_sum"] <- ((Nanog_subset_time[1,"Intensity.Sum"] + Nanog_subset_time[2,"Intensity.Sum"]))
          Distance_Nanog_dataframe[k,"Volume_sum"] <- ((Nanog_subset_time[1,"Volume"] + Nanog_subset_time[2,"Volume"]))
          
          Distance_Nanog_dataframe$Distance_clusters[k]  <- sqrt((Nanog_subset_time$Position.X[1] - Nanog_subset_time$Position.X[2])^2 + (Nanog_subset_time$Position.Y[1] - Nanog_subset_time$Position.Y[2])^2 + (Nanog_subset_time$Position.Z[1] - Nanog_subset_time$Position.Z[2])^2)
          k=k+1}
          else if(nrow(Nanog_subset_time)==3)
          {Distance_Nanog_dataframe[k,] <- Nanog_subset_time[1,]
          
          Distance_Nanog_dataframe[k,"Intensity.Sum_sum"] <- ((Nanog_subset_time[1,"Intensity.Sum"] + Nanog_subset_time[2,"Intensity.Sum"]+Nanog_subset_time[3,"Intensity.Sum"]))
          Distance_Nanog_dataframe[k,"Volume_sum"] <- ((Nanog_subset_time[1,"Volume"] + Nanog_subset_time[2,"Volume"]+Nanog_subset_time[3,"Volume"]))
          
          Distance_Nanog_MOVIE_1 <- sqrt((Nanog_subset_time$Position.X[1] - Nanog_subset_time$Position.X[2])^2 + (Nanog_subset_time$Position.Y[1] - Nanog_subset_time$Position.Y[2])^2 + (Nanog_subset_time$Position.Z[1] - Nanog_subset_time$Position.Z[2])^2) 
          Distance_Nanog_MOVIE_2 <- sqrt((Nanog_subset_time$Position.X[1] - Nanog_subset_time$Position.X[3])^2 + (Nanog_subset_time$Position.Y[1] - Nanog_subset_time$Position.Y[3])^2 + (Nanog_subset_time$Position.Z[1] - Nanog_subset_time$Position.Z[3])^2)
          Distance_Nanog_MOVIE_3 <- sqrt((Nanog_subset_time$Position.X[2] - Nanog_subset_time$Position.X[3])^2 + (Nanog_subset_time$Position.Y[2] - Nanog_subset_time$Position.Y[3])^2 + (Nanog_subset_time$Position.Z[2] - Nanog_subset_time$Position.Z[3])^2) 
          
          Distance_Nanog_dataframe$Distance_clusters[k] <- max(c(Distance_Nanog_MOVIE_1,Distance_Nanog_MOVIE_2,Distance_Nanog_MOVIE_3))
          k=k+1}
          
          else if(nrow(Nanog_subset_time)==4)
          {Distance_Nanog_dataframe[k,] <- Nanog_subset_time[1,]
          
          Distance_Nanog_dataframe[k,"Intensity.Sum_sum"] <- ((Nanog_subset_time[1,"Intensity.Sum"] + Nanog_subset_time[2,"Intensity.Sum"]+Nanog_subset_time[3,"Intensity.Sum"]+Nanog_subset_time[4,"Intensity.Sum"]))
          Distance_Nanog_dataframe[k,"Volume_sum"] <- ((Nanog_subset_time[1,"Volume"] + Nanog_subset_time[2,"Volume"]+Nanog_subset_time[3,"Volume"]+Nanog_subset_time[4,"Volume"]))
          
          Distance_Nanog_MOVIE_1 <- sqrt((Nanog_subset_time$Position.X[1] - Nanog_subset_time$Position.X[2])^2 + (Nanog_subset_time$Position.Y[1] - Nanog_subset_time$Position.Y[2])^2 + (Nanog_subset_time$Position.Z[1] - Nanog_subset_time$Position.Z[2])^2) 
          Distance_Nanog_MOVIE_2 <- sqrt((Nanog_subset_time$Position.X[1] - Nanog_subset_time$Position.X[3])^2 + (Nanog_subset_time$Position.Y[1] - Nanog_subset_time$Position.Y[3])^2 + (Nanog_subset_time$Position.Z[1] - Nanog_subset_time$Position.Z[3])^2) 
          Distance_Nanog_MOVIE_3 <- sqrt((Nanog_subset_time$Position.X[1] - Nanog_subset_time$Position.X[3])^2 + (Nanog_subset_time$Position.Y[1] - Nanog_subset_time$Position.Y[3])^2 + (Nanog_subset_time$Position.Z[1] - Nanog_subset_time$Position.Z[3])^2) 
          Distance_Nanog_MOVIE_4 <- sqrt((Nanog_subset_time$Position.X[1] - Nanog_subset_time$Position.X[4])^2 + (Nanog_subset_time$Position.Y[1] - Nanog_subset_time$Position.Y[4])^2 + (Nanog_subset_time$Position.Z[1] - Nanog_subset_time$Position.Z[4])^2) 
          Distance_Nanog_MOVIE_5 <- sqrt((Nanog_subset_time$Position.X[2] - Nanog_subset_time$Position.X[3])^2 + (Nanog_subset_time$Position.Y[2] - Nanog_subset_time$Position.Y[3])^2 + (Nanog_subset_time$Position.Z[2] - Nanog_subset_time$Position.Z[3])^2) 
          Distance_Nanog_MOVIE_6 <- sqrt((Nanog_subset_time$Position.X[2] - Nanog_subset_time$Position.X[4])^2 + (Nanog_subset_time$Position.Y[2] - Nanog_subset_time$Position.Y[4])^2 + (Nanog_subset_time$Position.Z[2] - Nanog_subset_time$Position.Z[4])^2) 
          Distance_Nanog_MOVIE_7 <- sqrt((Nanog_subset_time$Position.X[3] - Nanog_subset_time$Position.X[4])^2 + (Nanog_subset_time$Position.Y[3] - Nanog_subset_time$Position.Y[4])^2 + (Nanog_subset_time$Position.Z[3] - Nanog_subset_time$Position.Z[4])^2) 
          
          Distance_Nanog_dataframe$Distance_clusters[k] <- max(c(Distance_Nanog_MOVIE_1,Distance_Nanog_MOVIE_2,Distance_Nanog_MOVIE_3,Distance_Nanog_MOVIE_4,Distance_Nanog_MOVIE_5,Distance_Nanog_MOVIE_6,Distance_Nanog_MOVIE_7))
          
          k=k+1
          }
          
          else if(nrow(Nanog_subset_time)==5)
          {Distance_Nanog_dataframe[k,] <- Nanog_subset_time[1,]
          
          
          Distance_Nanog_dataframe[k,"Intensity.Sum_sum"] <- ((Nanog_subset_time[1,"Intensity.Sum"] + Nanog_subset_time[2,"Intensity.Sum"]+Nanog_subset_time[3,"Intensity.Sum"]+Nanog_subset_time[4,"Intensity.Sum"]+Nanog_subset_time[5,"Intensity.Sum"]))
          Distance_Nanog_dataframe[k,"Volume_sum"] <- ((Nanog_subset_time[1,"Volume"] + Nanog_subset_time[2,"Volume"]+Nanog_subset_time[3,"Volume"]+Nanog_subset_time[4,"Volume"]+Nanog_subset_time[5,"Volume"]))
          
          
          
          Distance_Nanog_MOVIE_1 <- sqrt((Nanog_subset_time$Position.X[1] - Nanog_subset_time$Position.X[2])^2 + (Nanog_subset_time$Position.Y[1] - Nanog_subset_time$Position.Y[2])^2 + (Nanog_subset_time$Position.Z[1] - Nanog_subset_time$Position.Z[2])^2) 
          
          Distance_Nanog_MOVIE_2 <- sqrt((Nanog_subset_time$Position.X[1] - Nanog_subset_time$Position.X[3])^2 + (Nanog_subset_time$Position.Y[1] - Nanog_subset_time$Position.Y[3])^2 + (Nanog_subset_time$Position.Z[1] - Nanog_subset_time$Position.Z[3])^2) 
          
          Distance_Nanog_MOVIE_3 <- sqrt((Nanog_subset_time$Position.X[1] - Nanog_subset_time$Position.X[3])^2 + (Nanog_subset_time$Position.Y[1] - Nanog_subset_time$Position.Y[3])^2 + (Nanog_subset_time$Position.Z[1] - Nanog_subset_time$Position.Z[3])^2) 
          
          Distance_Nanog_MOVIE_4 <- sqrt((Nanog_subset_time$Position.X[1] - Nanog_subset_time$Position.X[4])^2 + (Nanog_subset_time$Position.Y[1] - Nanog_subset_time$Position.Y[4])^2 + (Nanog_subset_time$Position.Z[1] - Nanog_subset_time$Position.Z[4])^2) 
          
          Distance_Nanog_MOVIE_5 <- sqrt((Nanog_subset_time$Position.X[1] - Nanog_subset_time$Position.X[5])^2 + (Nanog_subset_time$Position.Y[1] - Nanog_subset_time$Position.Y[5])^2 + (Nanog_subset_time$Position.Z[1] - Nanog_subset_time$Position.Z[5])^2) 
          
          Distance_Nanog_MOVIE_6 <- sqrt((Nanog_subset_time$Position.X[2] - Nanog_subset_time$Position.X[3])^2 + (Nanog_subset_time$Position.Y[2] - Nanog_subset_time$Position.Y[3])^2 + (Nanog_subset_time$Position.Z[2] - Nanog_subset_time$Position.Z[3])^2) 
          
          Distance_Nanog_MOVIE_7 <- sqrt((Nanog_subset_time$Position.X[2] - Nanog_subset_time$Position.X[4])^2 + (Nanog_subset_time$Position.Y[2] - Nanog_subset_time$Position.Y[4])^2 + (Nanog_subset_time$Position.Z[2] - Nanog_subset_time$Position.Z[4])^2) 
          
          Distance_Nanog_MOVIE_8 <- sqrt((Nanog_subset_time$Position.X[2] - Nanog_subset_time$Position.X[5])^2 + (Nanog_subset_time$Position.Y[2] - Nanog_subset_time$Position.Y[5])^2 + (Nanog_subset_time$Position.Z[2] - Nanog_subset_time$Position.Z[5])^2) 
          
          Distance_Nanog_MOVIE_9 <- sqrt((Nanog_subset_time$Position.X[3] - Nanog_subset_time$Position.X[4])^2 + (Nanog_subset_time$Position.Y[3] - Nanog_subset_time$Position.Y[4])^2 + (Nanog_subset_time$Position.Z[3] - Nanog_subset_time$Position.Z[4])^2) 
          
          Distance_Nanog_MOVIE_10 <- sqrt((Nanog_subset_time$Position.X[3] - Nanog_subset_time$Position.X[5])^2 + (Nanog_subset_time$Position.Y[3] - Nanog_subset_time$Position.Y[5])^2 + (Nanog_subset_time$Position.Z[3] - Nanog_subset_time$Position.Z[5])^2) 
          
          Distance_Nanog_MOVIE_11 <- sqrt((Nanog_subset_time$Position.X[4] - Nanog_subset_time$Position.X[5])^2 + (Nanog_subset_time$Position.Y[4] - Nanog_subset_time$Position.Y[5])^2 + (Nanog_subset_time$Position.Z[4] - Nanog_subset_time$Position.Z[5])^2) 
          
          Distance_Nanog_dataframe$Distance_clusters[k] <- max(c(Distance_Nanog_MOVIE_1,Distance_Nanog_MOVIE_2,Distance_Nanog_MOVIE_3,Distance_Nanog_MOVIE_4,Distance_Nanog_MOVIE_5,Distance_Nanog_MOVIE_6,Distance_Nanog_MOVIE_7,Distance_Nanog_MOVIE_8,Distance_Nanog_MOVIE_9,Distance_Nanog_MOVIE_10,Distance_Nanog_MOVIE_11))
          
          k=k+1
          }
          
          else if(nrow(Nanog_subset_time)==6)
          {Distance_Nanog_dataframe[k,] <- Nanog_subset_time[1,]
          Distance_Nanog_MOVIE_1 <- sqrt((Nanog_subset_time$Position.X[1] - Nanog_subset_time$Position.X[2])^2 + (Nanog_subset_time$Position.Y[1] - Nanog_subset_time$Position.Y[2])^2 + (Nanog_subset_time$Position.Z[1] - Nanog_subset_time$Position.Z[2])^2) 
          
          Distance_Nanog_MOVIE_2 <- sqrt((Nanog_subset_time$Position.X[1] - Nanog_subset_time$Position.X[3])^2 + (Nanog_subset_time$Position.Y[1] - Nanog_subset_time$Position.Y[3])^2 + (Nanog_subset_time$Position.Z[1] - Nanog_subset_time$Position.Z[3])^2) 
          
          Distance_Nanog_MOVIE_3 <- sqrt((Nanog_subset_time$Position.X[1] - Nanog_subset_time$Position.X[3])^2 + (Nanog_subset_time$Position.Y[1] - Nanog_subset_time$Position.Y[3])^2 + (Nanog_subset_time$Position.Z[1] - Nanog_subset_time$Position.Z[3])^2) 
          
          Distance_Nanog_MOVIE_4 <- sqrt((Nanog_subset_time$Position.X[1] - Nanog_subset_time$Position.X[4])^2 + (Nanog_subset_time$Position.Y[1] - Nanog_subset_time$Position.Y[4])^2 + (Nanog_subset_time$Position.Z[1] - Nanog_subset_time$Position.Z[4])^2) 
          
          Distance_Nanog_MOVIE_5 <- sqrt((Nanog_subset_time$Position.X[1] - Nanog_subset_time$Position.X[5])^2 + (Nanog_subset_time$Position.Y[1] - Nanog_subset_time$Position.Y[5])^2 + (Nanog_subset_time$Position.Z[1] - Nanog_subset_time$Position.Z[5])^2) 
          
          Distance_Nanog_MOVIE_6 <- sqrt((Nanog_subset_time$Position.X[1] - Nanog_subset_time$Position.X[6])^2 + (Nanog_subset_time$Position.Y[1] - Nanog_subset_time$Position.Y[6])^2 + (Nanog_subset_time$Position.Z[1] - Nanog_subset_time$Position.Z[6])^2) 
          
          Distance_Nanog_MOVIE_7 <- sqrt((Nanog_subset_time$Position.X[2] - Nanog_subset_time$Position.X[3])^2 + (Nanog_subset_time$Position.Y[2] - Nanog_subset_time$Position.Y[3])^2 + (Nanog_subset_time$Position.Z[2] - Nanog_subset_time$Position.Z[3])^2) 
          
          Distance_Nanog_MOVIE_8 <- sqrt((Nanog_subset_time$Position.X[2] - Nanog_subset_time$Position.X[4])^2 + (Nanog_subset_time$Position.Y[2] - Nanog_subset_time$Position.Y[4])^2 + (Nanog_subset_time$Position.Z[2] - Nanog_subset_time$Position.Z[4])^2) 
          
          Distance_Nanog_MOVIE_9 <- sqrt((Nanog_subset_time$Position.X[2] - Nanog_subset_time$Position.X[5])^2 + (Nanog_subset_time$Position.Y[2] - Nanog_subset_time$Position.Y[5])^2 + (Nanog_subset_time$Position.Z[2] - Nanog_subset_time$Position.Z[5])^2) 
          
          Distance_Nanog_MOVIE_10 <- sqrt((Nanog_subset_time$Position.X[2] - Nanog_subset_time$Position.X[6])^2 + (Nanog_subset_time$Position.Y[2] - Nanog_subset_time$Position.Y[6])^2 + (Nanog_subset_time$Position.Z[2] - Nanog_subset_time$Position.Z[6])^2) 
          
          Distance_Nanog_MOVIE_11 <- sqrt((Nanog_subset_time$Position.X[3] - Nanog_subset_time$Position.X[4])^2 + (Nanog_subset_time$Position.Y[3] - Nanog_subset_time$Position.Y[4])^2 + (Nanog_subset_time$Position.Z[3] - Nanog_subset_time$Position.Z[4])^2) 
          
          Distance_Nanog_MOVIE_12 <- sqrt((Nanog_subset_time$Position.X[3] - Nanog_subset_time$Position.X[5])^2 + (Nanog_subset_time$Position.Y[3] - Nanog_subset_time$Position.Y[5])^2 + (Nanog_subset_time$Position.Z[3] - Nanog_subset_time$Position.Z[5])^2) 
          
          Distance_Nanog_MOVIE_13 <- sqrt((Nanog_subset_time$Position.X[3] - Nanog_subset_time$Position.X[6])^2 + (Nanog_subset_time$Position.Y[3] - Nanog_subset_time$Position.Y[6])^2 + (Nanog_subset_time$Position.Z[3] - Nanog_subset_time$Position.Z[6])^2) 
          
          Distance_Nanog_MOVIE_14 <- sqrt((Nanog_subset_time$Position.X[4] - Nanog_subset_time$Position.X[5])^2 + (Nanog_subset_time$Position.Y[4] - Nanog_subset_time$Position.Y[5])^2 + (Nanog_subset_time$Position.Z[4] - Nanog_subset_time$Position.Z[5])^2) 
          
          Distance_Nanog_MOVIE_15 <- sqrt((Nanog_subset_time$Position.X[4] - Nanog_subset_time$Position.X[6])^2 + (Nanog_subset_time$Position.Y[4] - Nanog_subset_time$Position.Y[6])^2 + (Nanog_subset_time$Position.Z[4] - Nanog_subset_time$Position.Z[6])^2) 
          
          Distance_Nanog_MOVIE_16 <- sqrt((Nanog_subset_time$Position.X[5] - Nanog_subset_time$Position.X[6])^2 + (Nanog_subset_time$Position.Y[5] - Nanog_subset_time$Position.Y[6])^2 + (Nanog_subset_time$Position.Z[5] - Nanog_subset_time$Position.Z[6])^2) 
          
          Distance_Nanog_dataframe$Distance_clusters[k] <- max(c(Distance_Nanog_MOVIE_1,Distance_Nanog_MOVIE_2,Distance_Nanog_MOVIE_3,Distance_Nanog_MOVIE_4,Distance_Nanog_MOVIE_5,Distance_Nanog_MOVIE_6,Distance_Nanog_MOVIE_7,Distance_Nanog_MOVIE_8,Distance_Nanog_MOVIE_9,Distance_Nanog_MOVIE_10,Distance_Nanog_MOVIE_11,Distance_Nanog_MOVIE_12,Distance_Nanog_MOVIE_13,Distance_Nanog_MOVIE_14,Distance_Nanog_MOVIE_15,Distance_Nanog_MOVIE_16))
          
          k=k+1
          }
          
        }
      }
    }
  }
  
  
}

#=========================================================
# Associate MCP 2D segmentation (Damian's data) with Nanog 
#=========================================================

############Associate MCP/Nanog/MOVIE shape in 2D with Damian data (X, Y position)
{
  #Assign one trackID for each track
  
  MCP_statistics_total_2D_Damian$Nucleus_ID <- paste(MCP_statistics_total_2D_Damian$Date,MCP_statistics_total_2D_Damian$Stage,MCP_statistics_total_2D_Damian$Stack,MCP_statistics_total_2D_Damian$Nucleus, sep="_")
  Distance_Nanog_dataframe
  
  #TO CHANGE
  
  #Normalization of the position in X and Y
  
  for(i in 1:max(na.omit(Distance_Nanog_dataframe$Unique_track_Nanog_3D)))
  {Distance_Nanog_dataframe[which(Distance_Nanog_dataframe$Unique_track_Nanog_3D==i),"Position.X"] <- Distance_Nanog_dataframe[which(Distance_Nanog_dataframe$Unique_track_Nanog_3D==i),"Position.X"] -min(na.omit(Distance_Nanog_dataframe[which(Distance_Nanog_dataframe$Unique_track_Nanog_3D==i),"Position.X"]))
  Distance_Nanog_dataframe[which(Distance_Nanog_dataframe$Unique_track_Nanog_3D==i),"Position.Y"] <- Distance_Nanog_dataframe[which(Distance_Nanog_dataframe$Unique_track_Nanog_3D==i),"Position.Y"] -min(na.omit(Distance_Nanog_dataframe[which(Distance_Nanog_dataframe$Unique_track_Nanog_3D==i),"Position.Y"]))}
  
  MCP_statistics_total_2D_Damian$label_1x_pos_um <-  MCP_statistics_total_2D_Damian$label_1x_pos*0.1087
  MCP_statistics_total_2D_Damian$label_1y_pos_um <- MCP_statistics_total_2D_Damian$label_1y_pos*0.1113
  
  for(i in 1:max(MCP_statistics_total_2D_Damian$Unique_track_2D_D))
  {MCP_statistics_total_2D_Damian[MCP_statistics_total_2D_Damian$Unique_track_2D_D==i,"label_1x_pos_um"] <- MCP_statistics_total_2D_Damian[MCP_statistics_total_2D_Damian$Unique_track_2D_D==i,"label_1x_pos_um"] -min(na.omit(MCP_statistics_total_2D_Damian[MCP_statistics_total_2D_Damian$Unique_track_2D_D==i,"label_1x_pos_um"]))
  MCP_statistics_total_2D_Damian[MCP_statistics_total_2D_Damian$Unique_track_2D_D==i,"label_1x_pos_um"] <- MCP_statistics_total_2D_Damian[MCP_statistics_total_2D_Damian$Unique_track_2D_D==i,"label_1x_pos_um"] -min(na.omit(MCP_statistics_total_2D_Damian[MCP_statistics_total_2D_Damian$Unique_track_2D_D==i,"label_1x_pos_um"]))}
  
  MCP_statistics_total_2D_Damian <- MCP_statistics_total_2D_Damian[!is.na(MCP_statistics_total_2D_Damian$label_1_mean_edge_distace),]
  
  
  Distance_Nanog_dataframe_subset <- unique(Distance_Nanog_dataframe[,c("Date","Stage","Stack","Nucleus","TrackID_Nanog_3D","Time","Unique_track_Nanog_3D","Intensity.Sum_sum","Time_activation","Seconds_activation","Seconds_normalized","Phenotype","Fused_activation","Length_double","Time_fusion","Time_fusion_normalized","Seconds_fusion_normalized","Position.X","Position.Y","Position.Z","Phenotype_simple","Spot_number","Volume_sum"),])
  
  MCP_statistics_total_2D_Damian$Nucleus_ID <- paste(MCP_statistics_total_2D_Damian$Date,MCP_statistics_total_2D_Damian$Stage,MCP_statistics_total_2D_Damian$Stack,MCP_statistics_total_2D_Damian$Nucleus, sep="_")
  Distance_Nanog_dataframe_subset$Nucleus_ID <- paste(Distance_Nanog_dataframe_subset$Date,Distance_Nanog_dataframe_subset$Stage,Distance_Nanog_dataframe_subset$Stack,Distance_Nanog_dataframe_subset$Nucleus, sep="_")
  
  for(i in 1:length(unique(MCP_statistics_total_2D_Damian$Nucleus_ID)))
  {
    
    Nucleus_ID <- unique(MCP_statistics_total_2D_Damian$Nucleus_ID)[i]
    
    #X and Y position each allele in MCP 2D 
    
    a <- unique(Distance_Nanog_dataframe_subset[which(Distance_Nanog_dataframe_subset$Nucleus_ID == Nucleus_ID),"TrackID_Nanog_3D"])[1]
    b <- unique(Distance_Nanog_dataframe_subset[which(Distance_Nanog_dataframe_subset$Nucleus_ID == Nucleus_ID),"TrackID_Nanog_3D"])[2]
    
    if(length(unique(Distance_Nanog_dataframe_subset[Distance_Nanog_dataframe_subset$Nucleus_ID == Nucleus_ID,"TrackID_Nanog_3D"])) >1){
      
      MCP_X_1 <- MCP_statistics_total_2D_Damian[which(MCP_statistics_total_2D_Damian$Nucleus_ID ==Nucleus_ID&MCP_statistics_total_2D_Damian$TrackID_MCP_2D_D=="1"),c("label_1x_pos_um","Time")]
      MCP_X_2 <- MCP_statistics_total_2D_Damian[which(MCP_statistics_total_2D_Damian$Nucleus_ID ==Nucleus_ID&MCP_statistics_total_2D_Damian$TrackID_MCP_2D_D=="2"),c("label_1x_pos_um","Time")]
      
      Nanog_X_1 <- Distance_Nanog_dataframe_subset[which(Distance_Nanog_dataframe_subset$Nucleus_ID ==Nucleus_ID&Distance_Nanog_dataframe_subset$TrackID_Nanog_3D==a),c("Position.X","Time")]
      Nanog_X_2 <- Distance_Nanog_dataframe_subset[which(Distance_Nanog_dataframe_subset$Nucleus_ID ==Nucleus_ID&Distance_Nanog_dataframe_subset$TrackID_Nanog_3D==b),c("Position.X","Time")]
      
      X_1_1 <- merge(MCP_X_1,Nanog_X_1,by="Time")
      X_1_1_mean <- mean(abs(X_1_1$label_1x_pos_um-X_1_1$Position.X))
      X_1_2 <- merge(MCP_X_1,Nanog_X_2,by="Time")
      X_1_2_mean <-mean(abs(X_1_2$label_1x_pos_um-X_1_2$Position.X))
      
      if(X_1_1_mean<X_1_2_mean)
      {
        MCP_statistics_total_2D_Damian[which(MCP_statistics_total_2D_Damian$Nucleus_ID ==Nucleus_ID&MCP_statistics_total_2D_Damian$TrackID_MCP_2D_D=="1"),"TrackID_Nanog_3D"] <- a
        MCP_statistics_total_2D_Damian[which(MCP_statistics_total_2D_Damian$Nucleus_ID ==Nucleus_ID&MCP_statistics_total_2D_Damian$TrackID_MCP_2D_D=="2"),"TrackID_Nanog_3D"] <- b
        Distance_Nanog_dataframe_subset[Distance_Nanog_dataframe_subset$Nucleus_ID ==Nucleus_ID&Distance_Nanog_dataframe_subset$TrackID_Nanog==a,"TrackID_MCP_2D_D"] <- 1
        Distance_Nanog_dataframe_subset[Distance_Nanog_dataframe_subset$Nucleus_ID ==Nucleus_ID&Distance_Nanog_dataframe_subset$TrackID_Nanog==b,"TrackID_MCP_2D_D"] <- 2
        
      } else if(X_1_1_mean>X_1_2_mean){MCP_statistics_total_2D_Damian[which(MCP_statistics_total_2D_Damian$Nucleus_ID ==Nucleus_ID&MCP_statistics_total_2D_Damian$TrackID_MCP_2D_D=="1"),"TrackID_Nanog_3D"] <- b
      MCP_statistics_total_2D_Damian[which(MCP_statistics_total_2D_Damian$Nucleus_ID ==Nucleus_ID&MCP_statistics_total_2D_Damian$TrackID_MCP_2D_D=="2"),"TrackID_Nanog_3D"] <- a
      Distance_Nanog_dataframe_subset[Distance_Nanog_dataframe_subset$Nucleus_ID ==Nucleus_ID&Distance_Nanog_dataframe_subset$TrackID_Nanog==a,"TrackID_MCP_2D_D"] <- 2
      Distance_Nanog_dataframe_subset[Distance_Nanog_dataframe_subset$Nucleus_ID ==Nucleus_ID&Distance_Nanog_dataframe_subset$TrackID_Nanog==b,"TrackID_MCP_2D_D"] <- 1
      }
    }
  }
  
  
  MCP_statistics_total_2D_Damian$TrackID_Nanog_3D <- as.numeric(as.character(MCP_statistics_total_2D_Damian$TrackID_Nanog_3D))
  Distance_Nanog_dataframe_subset$TrackID_Nanog_3D <- as.numeric(as.character(Distance_Nanog_dataframe_subset$TrackID_Nanog_3D))
  MCP_statistics_total_2D_Damian$Time <- as.numeric(as.character(MCP_statistics_total_2D_Damian$Time))
  Distance_Nanog_dataframe_subset$Time <- as.numeric(as.character(Distance_Nanog_dataframe_subset$Time))
  MCP_statistics_total_2D_Damian$Stack <- as.numeric(as.character(MCP_statistics_total_2D_Damian$Stack))
  Distance_Nanog_dataframe_subset$Stack <- as.numeric(as.character(Distance_Nanog_dataframe_subset$Stack))
  
  MCP_statistics_total_2D_Damian[MCP_statistics_total_2D_Damian$Nucleus_ID=="09.01_1k_3_9","TrackID_Nanog_3D"] <- 1000000000
  MCP_statistics_total_2D_Damian[MCP_statistics_total_2D_Damian$Nucleus_ID=="12.01_1k_9_4","TrackID_Nanog_3D"] <- 1000000000
  Distance_Nanog_dataframe_subset[Distance_Nanog_dataframe_subset$Nucleus_ID=="09.01_1k_3_9","TrackID_MCP_2D_D"] <- 1
  Distance_Nanog_dataframe_subset[Distance_Nanog_dataframe_subset$Nucleus_ID=="12.01_1k_9_4","TrackID_MCP_2D_D"] <- 1
  
  MCP_statistics_shape_total_2D_Nanog_D <- Distance_Nanog_dataframe_subset %>% full_join(MCP_statistics_total_2D_Damian, by=c("Date","Stage","Stack","Nucleus","Nucleus_ID","TrackID_Nanog_3D","Time","TrackID_MCP_2D_D"))
  
  
  #Allow all the lines with MCP shape detected but no Nanog to have a "Time_activation" parameter as well
  for(i in 1:max(na.omit(MCP_statistics_shape_total_2D_Nanog_D$Unique_track_Nanog_3D)))
  {H <- unique(na.omit(MCP_statistics_shape_total_2D_Nanog_D[which(MCP_statistics_shape_total_2D_Nanog_D$Unique_track_Nanog_3D==i),"TrackID_MCP_2D_D"]))
  N <- unique(na.omit(MCP_statistics_shape_total_2D_Nanog_D[which(MCP_statistics_shape_total_2D_Nanog_D$Unique_track_Nanog_3D==i),"Nucleus_ID"]))
  
  MCP_statistics_shape_total_2D_Nanog_D[MCP_statistics_shape_total_2D_Nanog_D$Nucleus_ID==N&MCP_statistics_shape_total_2D_Nanog_D$TrackID_MCP_2D_D==H,"Time_activation"]<- unique(na.omit(MCP_statistics_shape_total_2D_Nanog_D[which(MCP_statistics_shape_total_2D_Nanog_D$Unique_track_Nanog_3D==i),"Time_activation"]))
  MCP_statistics_shape_total_2D_Nanog_D[MCP_statistics_shape_total_2D_Nanog_D$Nucleus_ID==N&MCP_statistics_shape_total_2D_Nanog_D$TrackID_MCP_2D_D==H,"Phenotype"]<- unique(na.omit(MCP_statistics_shape_total_2D_Nanog_D[which(MCP_statistics_shape_total_2D_Nanog_D$Unique_track_Nanog_3D==i),"Phenotype"]))
  MCP_statistics_shape_total_2D_Nanog_D[MCP_statistics_shape_total_2D_Nanog_D$Nucleus_ID==N&MCP_statistics_shape_total_2D_Nanog_D$TrackID_MCP_2D_D==H,"Phenotype_simple"]<- unique(na.omit(MCP_statistics_shape_total_2D_Nanog_D[which(MCP_statistics_shape_total_2D_Nanog_D$Unique_track_Nanog_3D==i),"Phenotype_simple"]))}
  
  MCP_statistics_shape_total_2D_Nanog_D$Time_normalized <- MCP_statistics_shape_total_2D_Nanog_D$Time - MCP_statistics_shape_total_2D_Nanog_D$Time_activation
  MCP_statistics_shape_total_2D_Nanog_D$Unique_track <- as.numeric(MCP_statistics_shape_total_2D_Nanog_D$Unique_track_2D_D)
  
  MCP_statistics_shape_total_2D_Nanog_D$Seconds_normalized <-   MCP_statistics_shape_total_2D_Nanog_D$Time_normalized*15
  
  MCP_statistics_shape_total_2D_Nanog_D <- MCP_statistics_shape_total_2D_Nanog_D[-which(MCP_statistics_shape_total_2D_Nanog_D$Unique_track_Nanog_3D==81),]
  
  
  
  
}

############Calculate other parameters in 2D, my data + Damian's data
{
  #NORMALIZATION NANOG TOTAL INTENSITY
  MCP_statistics_shape_total_2D_Nanog_D$Unique_track_2D_D <- as.numeric(MCP_statistics_shape_total_2D_Nanog_D$Unique_track_2D_D)
  MCP_statistics_shape_total_2D_Nanog_D$Intensity.Sum_sum <- as.numeric(MCP_statistics_shape_total_2D_Nanog_D$Intensity.Sum_sum)
  
  
  for(i in 1:max(na.omit(MCP_statistics_shape_total_2D_Nanog_D$Unique_track_2D_D)))
  { MCP_statistics_shape_total_2D_Nanog_D[which(MCP_statistics_shape_total_2D_Nanog_D$Unique_track_2D_D==i),"Intensity.Sum_sum_norm"] <-
    MCP_statistics_shape_total_2D_Nanog_D[which(MCP_statistics_shape_total_2D_Nanog_D$Unique_track_2D_D==i),"Intensity.Sum_sum"]/max(na.omit(MCP_statistics_shape_total_2D_Nanog_D[which(MCP_statistics_shape_total_2D_Nanog_D$Unique_track_2D_D==i),"Intensity.Sum_sum"]))}
  
  MCP_statistics_shape_total_2D_Nanog_D$Seconds_normalized<-plyr::round_any(MCP_statistics_shape_total_2D_Nanog_D$Seconds_normalized, 5, f = round)
  MCP_statistics_shape_total_2D_Nanog_D_subset <- MCP_statistics_shape_total_2D_Nanog_D[MCP_statistics_shape_total_2D_Nanog_D$Seconds_normalized %in% c(-105:105),]
  MCP_statistics_shape_total_2D_Nanog_D_subset$Intensity.Sum_sum <- as.numeric(MCP_statistics_shape_total_2D_Nanog_D_subset$Intensity.Sum_sum)
  
  
  #Calculate the average value for AR during mitosis
  mean(na.omit(MCP_statistics_shape_total_2D_Nanog_D[which(MCP_statistics_shape_total_2D_Nanog_D$Nucleus_Sphericity<0.6),c("Radial_distance_cv")]))
  0.1183266
  
  #Calculate and normalize Radial distance
  
  MCP_statistics_shape_total_2D_Nanog_D_subset$Radial_distance_cv <- MCP_statistics_shape_total_2D_Nanog_D_subset$label_1_std_edge_distace/MCP_statistics_shape_total_2D_Nanog_D_subset$label_1_mean_edge_distace
  
  MCP_statistics_shape_total_2D_Nanog_D$Radial_distance_cv <- MCP_statistics_shape_total_2D_Nanog_D$label_1_std_edge_distace/MCP_statistics_shape_total_2D_Nanog_D$label_1_mean_edge_distace
}

#=========================================================
# Plots
#=========================================================

#Fig. 5a
{
#Determine values for ribbon
{
  DT <- as.data.table(MCP_statistics_shape_total_2D_Nanog_D_subset)
  Summary <- DT[ , list(MEAN_DV = mean(na.omit(Radial_distance_cv)), SD_DV = sd(na.omit(Radial_distance_cv)), Q10= quantile(Radial_distance_cv, prob=c(0.10), na.rm=TRUE),Q25 = quantile(Radial_distance_cv, prob=c(0.25),na.rm=TRUE),Q75 = quantile(Radial_distance_cv, prob=c(0.75),na.rm=TRUE), Q90= quantile(Radial_distance_cv, prob=c(0.90),na.rm=TRUE),N = .N),  by = c("Seconds_normalized")]
  Summary<- Summary[order(Summary$Seconds_normalized),]
  
  Mean_circ <- Summary$MEAN_DV
  sd_circ <- Summary$SD_DV
  
  CI_infcirc <- Summary$Q25
  CI_supcirc <- Summary$Q75
  
  Circ_plot <- cbind(Mean_circ,sd_circ,CI_infcirc,CI_supcirc)
  Circ_plot <-Circ_plot[complete.cases(Circ_plot[,c("CI_supcirc")]),]
  Circ_plot <- cbind(Circ_plot,unique(Summary$Seconds_normalized))
  Circ_plot <- data.frame(Circ_plot)
  
  colnames(Circ_plot) <- c("Mean","sd","CI_inf","CI_sup","Time")
  
  a <- Circ_plot %>%
    group_by(Time) %>%
    summarise_at(vars("Mean","sd","CI_inf","CI_sup"), list(name = mean)) }

#Plot
{
ggplot() + geom_jitter(MCP_statistics_shape_total_2D_Nanog_D_subset,
                       mapping = aes(y=Radial_distance_cv, x=as.numeric(as.character(Seconds_normalized))),
                       width=2.5,alpha=0.15, size = 2)+
  stat_summary_bin(MCP_statistics_shape_total_2D_Nanog_D_subset, 
                   mapping=aes(y=Radial_distance_cv, x=as.numeric(as.character(Seconds_normalized))),
                   fun="median",  geom="line", size=0.5, color="black") +
  geom_ribbon(a,mapping=aes(x=Time,ymin=CI_inf_name,ymax=CI_sup_name),
              alpha=0.5, fill="grey70")+
  geom_hline(data = MCP_statistics_shape_total_2D_Nanog_D_subset,aes(yintercept = 0.12,), size=0.5,color="black", linetype=2)+
  #geom_rect(mapping=aes(xmin=-3,xmax=3,ymin=0, ymax=0.5), alpha=0,color="red",size=0.3)+ 
  geom_vline(data = MCP_statistics_shape_total_2D_Nanog_D_subset,aes(xintercept = 0,), size=0.5,color="red", linetype=2)+ 
  theme_bw() + labs(x= "Time centered on transcription start (s)",y="Normalized Radial distance")+ 
  #ggtitle('Merging Tracks') + 
  ylim(0.00,0.5)+scale_x_continuous(breaks=seq(-90,90,30))+ 
  theme(aspect.ratio=3/2, axis.text = element_text(size=14),axis.title = element_text(size=14),
        panel.grid.minor.y = element_blank(),panel.grid.minor.x = element_blank())
dev.off()
}
}
  
#Fig. S4a, non-merging tracks
{
  #Determine values for ribbon
  {
    DT <- as.data.table(MCP_statistics_shape_total_2D_Nanog_D_subset[MCP_statistics_shape_total_2D_Nanog_D_subset$Phenotype_simple=="One_cluster",])
    Summary <- DT[ , list(MEAN_DV = mean(na.omit(Radial_distance_cv)), SD_DV = sd(na.omit(Radial_distance_cv)), Q10= quantile(Radial_distance_cv, prob=c(0.10), na.rm=TRUE),Q25 = quantile(Radial_distance_cv, prob=c(0.25),na.rm=TRUE),Q75 = quantile(Radial_distance_cv, prob=c(0.75),na.rm=TRUE), Q90= quantile(Radial_distance_cv, prob=c(0.90),na.rm=TRUE),N = .N),  by = c("Seconds_normalized")]
    Summary<- Summary[order(Summary$Seconds_normalized),]
    
    Mean_circ <- Summary$MEAN_DV
    sd_circ <- Summary$SD_DV
    
    CI_infcirc <- Summary$Q25
    CI_supcirc <- Summary$Q75
    
    Circ_plot <- cbind(Mean_circ,sd_circ,CI_infcirc,CI_supcirc)
    Circ_plot <-Circ_plot[complete.cases(Circ_plot[,c("CI_supcirc")]),]
    Circ_plot <- cbind(Circ_plot,unique(Summary$Seconds_normalized))
    Circ_plot <- data.frame(Circ_plot)
    
    colnames(Circ_plot) <- c("Mean","sd","CI_inf","CI_sup","Time")
    
    a <- Circ_plot %>%
      group_by(Time) %>%
      summarise_at(vars("Mean","sd","CI_inf","CI_sup"), list(name = mean)) }
  
  #Plot
  {
    ggplot() + geom_jitter(MCP_statistics_shape_total_2D_Nanog_D_subset[MCP_statistics_shape_total_2D_Nanog_D_subset$Phenotype_simple=="One_cluster",],
                           mapping = aes(y=Radial_distance_cv, x=as.numeric(as.character(Seconds_normalized))),
                           width=2.5,alpha=0.15, size = 2)+
      stat_summary_bin(MCP_statistics_shape_total_2D_Nanog_D_subset[MCP_statistics_shape_total_2D_Nanog_D_subset$Phenotype_simple=="One_cluster",], 
                       mapping=aes(y=Radial_distance_cv, x=as.numeric(as.character(Seconds_normalized))),
                       fun="median",  geom="line", size=0.5, color="black") +
      geom_ribbon(a,mapping=aes(x=Time,ymin=CI_inf_name,ymax=CI_sup_name),
                  alpha=0.5, fill="grey70")+
      geom_hline(data = MCP_statistics_shape_total_2D_Nanog_D_subset[MCP_statistics_shape_total_2D_Nanog_D_subset$Phenotype_simple=="One_cluster",],aes(yintercept = 0.12,), size=0.5,color="black", linetype=2)+
      #geom_rect(mapping=aes(xmin=-3,xmax=3,ymin=0, ymax=0.5), alpha=0,color="red",size=0.3)+ 
      geom_vline(data = MCP_statistics_shape_total_2D_Nanog_D_subset[MCP_statistics_shape_total_2D_Nanog_D_subset$Phenotype_simple=="One_cluster",],aes(xintercept = 0,), size=0.5,color="red", linetype=2)+ 
      theme_bw() + labs(x= "Time centered on transcription start (s)",y="Normalized Radial distance")+ 
      ggtitle('Non-merging clusters') + 
      ylim(0.00,0.5)+scale_x_continuous(breaks=seq(-90,90,30))+ 
      theme(aspect.ratio=3/2, axis.text = element_text(size=14),axis.title = element_text(size=14),
            panel.grid.minor.y = element_blank(),panel.grid.minor.x = element_blank())
  }
}

#Fig. S4a, merging tracks
{
  #Determine values for ribbon
  {
    DT <- as.data.table(MCP_statistics_shape_total_2D_Nanog_D_subset[MCP_statistics_shape_total_2D_Nanog_D_subset$Phenotype_simple=="More_cluster",])
    Summary <- DT[ , list(MEAN_DV = mean(na.omit(Radial_distance_cv)), SD_DV = sd(na.omit(Radial_distance_cv)), Q10= quantile(Radial_distance_cv, prob=c(0.10), na.rm=TRUE),Q25 = quantile(Radial_distance_cv, prob=c(0.25),na.rm=TRUE),Q75 = quantile(Radial_distance_cv, prob=c(0.75),na.rm=TRUE), Q90= quantile(Radial_distance_cv, prob=c(0.90),na.rm=TRUE),N = .N),  by = c("Seconds_normalized")]
    Summary<- Summary[order(Summary$Seconds_normalized),]
    
    Mean_circ <- Summary$MEAN_DV
    sd_circ <- Summary$SD_DV
    
    CI_infcirc <- Summary$Q25
    CI_supcirc <- Summary$Q75
    
    Circ_plot <- cbind(Mean_circ,sd_circ,CI_infcirc,CI_supcirc)
    Circ_plot <-Circ_plot[complete.cases(Circ_plot[,c("CI_supcirc")]),]
    Circ_plot <- cbind(Circ_plot,unique(Summary$Seconds_normalized))
    Circ_plot <- data.frame(Circ_plot)
    
    colnames(Circ_plot) <- c("Mean","sd","CI_inf","CI_sup","Time")
    
    a <- Circ_plot %>%
      group_by(Time) %>%
      summarise_at(vars("Mean","sd","CI_inf","CI_sup"), list(name = mean)) }
  
  #Plot
  {
    ggplot() + geom_jitter(MCP_statistics_shape_total_2D_Nanog_D_subset[MCP_statistics_shape_total_2D_Nanog_D_subset$Phenotype_simple=="More_cluster",],
                           mapping = aes(y=Radial_distance_cv, x=as.numeric(as.character(Seconds_normalized))),
                           width=2.5,alpha=0.15, size = 2)+
      stat_summary_bin(MCP_statistics_shape_total_2D_Nanog_D_subset[MCP_statistics_shape_total_2D_Nanog_D_subset$Phenotype_simple=="More_cluster",], 
                       mapping=aes(y=Radial_distance_cv, x=as.numeric(as.character(Seconds_normalized))),
                       fun="median",  geom="line", size=0.5, color="black") +
      geom_ribbon(a,mapping=aes(x=Time,ymin=CI_inf_name,ymax=CI_sup_name),
                  alpha=0.5, fill="grey70")+
      geom_hline(data = MCP_statistics_shape_total_2D_Nanog_D_subset[MCP_statistics_shape_total_2D_Nanog_D_subset$Phenotype_simple=="More_cluster",],aes(yintercept = 0.12,), size=0.5,color="black", linetype=2)+
      #geom_rect(mapping=aes(xmin=-3,xmax=3,ymin=0, ymax=0.5), alpha=0,color="red",size=0.3)+ 
      geom_vline(data = MCP_statistics_shape_total_2D_Nanog_D_subset[MCP_statistics_shape_total_2D_Nanog_D_subset$Phenotype_simple=="More_cluster",],aes(xintercept = 0,), size=0.5,color="red", linetype=2)+ 
      theme_bw() + labs(x= "Time centered on transcription start (s)",y="Normalized Radial distance")+ 
      ggtitle('Non-merging clusters') + 
      ylim(0.00,0.5)+scale_x_continuous(breaks=seq(-90,90,30))+ 
      theme(aspect.ratio=3/2, axis.text = element_text(size=14),axis.title = element_text(size=14),
            panel.grid.minor.y = element_blank(),panel.grid.minor.x = element_blank())
  }
}

#Fig. 5e, S6
{

  #Analysis
  {
    
    MCP_statistics_shape_total_2D_Nanog_D_subset$Unique_track_2D_D <- as.numeric(as.character(MCP_statistics_shape_total_2D_Nanog_D_subset$Unique_track_2D_D))
    MCP_statistics_shape_total_2D_Nanog_D_subset <- MCP_statistics_shape_total_2D_Nanog_D_subset[complete.cases(MCP_statistics_shape_total_2D_Nanog_D_subset[,c("Unique_track_2D_D")]),]
    levels(MCP_statistics_shape_total_2D_Nanog_D_subset$Phenotype_simple) <- c("One_cluster","More_cluster")
    
    
    #Connect the lines only when there is no breaking
    MCP_statistics_shape_total_2D_Nanog_D_subset <- MCP_statistics_shape_total_2D_Nanog_D_subset[order(MCP_statistics_shape_total_2D_Nanog_D_subset$Date,MCP_statistics_shape_total_2D_Nanog_D_subset$Stack,MCP_statistics_shape_total_2D_Nanog_D_subset$Nucleus,MCP_statistics_shape_total_2D_Nanog_D_subset$TrackID_MCP_2D,MCP_statistics_shape_total_2D_Nanog_D_subset$Time),]
    MCP_statistics_shape_total_2D_Nanog_D_subset$Group <- 0
    k=1
    for(i in 1:max(na.omit(as.numeric(MCP_statistics_shape_total_2D_Nanog_D_subset$Unique_track_2D_D))))
    {
      if(!max(MCP_statistics_shape_total_2D_Nanog_D_subset[which(MCP_statistics_shape_total_2D_Nanog_D_subset$Unique_track_2D_D==i),"Time"])=="-Inf")
      {
        for(j in 1:max(MCP_statistics_shape_total_2D_Nanog_D_subset[which(MCP_statistics_shape_total_2D_Nanog_D_subset$Unique_track_2D_D==i),"Time"]))
        {
          if(!is.na(MCP_statistics_shape_total_2D_Nanog_D_subset[which(MCP_statistics_shape_total_2D_Nanog_D_subset$Unique_track_2D_D==i),"Time"][j])) {
            if(MCP_statistics_shape_total_2D_Nanog_D_subset[which(MCP_statistics_shape_total_2D_Nanog_D_subset$Unique_track_2D_D==i),"Time"][j]==max(MCP_statistics_shape_total_2D_Nanog_D_subset[which(MCP_statistics_shape_total_2D_Nanog_D_subset$Unique_track_2D_D==i),"Time"]))
            {MCP_statistics_shape_total_2D_Nanog_D_subset[which(MCP_statistics_shape_total_2D_Nanog_D_subset$Unique_track_2D_D==i),"Group"][j] <- k
            k=k+1}
            else if(MCP_statistics_shape_total_2D_Nanog_D_subset[which(MCP_statistics_shape_total_2D_Nanog_D_subset$Unique_track_2D_D==i),"Time"][j]==(MCP_statistics_shape_total_2D_Nanog_D_subset[which(MCP_statistics_shape_total_2D_Nanog_D_subset$Unique_track_2D_D==i),"Time"][j+1]-1))
            {MCP_statistics_shape_total_2D_Nanog_D_subset[which(MCP_statistics_shape_total_2D_Nanog_D_subset$Unique_track_2D_D==i),"Group"][j] <- k}
            else if(MCP_statistics_shape_total_2D_Nanog_D_subset[which(MCP_statistics_shape_total_2D_Nanog_D_subset$Unique_track_2D_D==i),"Time"][j]!=(MCP_statistics_shape_total_2D_Nanog_D_subset[which(MCP_statistics_shape_total_2D_Nanog_D_subset$Unique_track_2D_D==i),"Time"][j+1]-1))
            {MCP_statistics_shape_total_2D_Nanog_D_subset[which(MCP_statistics_shape_total_2D_Nanog_D_subset$Unique_track_2D_D==i),"Group"][j] <- k
            k=k+1} else{k=k+1}
          }
        }
      }
    }
    
    
    for(i in 1:max(MCP_statistics_shape_total_2D_Nanog_D_subset$Unique_track_2D_D))
    {MCP_statistics_shape_total_2D_Nanog_D_subset[which(MCP_statistics_shape_total_2D_Nanog_D_subset$Unique_track_2D_D == i),"Plot_spot_number"] <- sum(na.omit(MCP_statistics_shape_total_2D_Nanog_D_subset[MCP_statistics_shape_total_2D_Nanog_D_subset$Unique_track_2D_D==i,"Spot_number"]))}
    
    
    #Maximum value for the coefficient of variation
    for(i in 1:max(MCP_statistics_shape_total_2D_Nanog_D_subset$Unique_track_2D_D))
    {MCP_statistics_shape_total_2D_Nanog_D_subset[MCP_statistics_shape_total_2D_Nanog_D_subset$Unique_track_2D_D==i,"Radial_distance_cv_max"] <- max(na.omit(MCP_statistics_shape_total_2D_Nanog_D_subset[MCP_statistics_shape_total_2D_Nanog_D_subset$Unique_track_2D_D==i,"Radial_distance_cv"]))
    }
    
    MCP_statistics_shape_total_2D_Nanog_D_subset$Unique_track_2D_D <- factor(MCP_statistics_shape_total_2D_Nanog_D_subset$Unique_track_2D_D, levels =  unique(MCP_statistics_shape_total_2D_Nanog_D_subset[order(MCP_statistics_shape_total_2D_Nanog_D_subset$Phenotype_simple,MCP_statistics_shape_total_2D_Nanog_D_subset$Radial_distance_cv_max),"Unique_track_2D_D"]))

    #Radial_distance_cv 
    
    #Make the Nanog fluorescence normalized to Radial_distance_cv levels
    MCP_statistics_shape_total_2D_Nanog_D_subset$Unique_track_2D_D <- as.numeric(as.character(MCP_statistics_shape_total_2D_Nanog_D_subset$Unique_track_2D_D))
  
    for(i in 1:max(MCP_statistics_shape_total_2D_Nanog_D_subset$Unique_track_2D_D))
    {
      MCP_statistics_shape_total_2D_Nanog_D_subset[MCP_statistics_shape_total_2D_Nanog_D_subset$Unique_track_2D_D==i,"Intensity.Sum_sum_norm_RD"] <- (0.5-(0))*
        ((MCP_statistics_shape_total_2D_Nanog_D_subset[MCP_statistics_shape_total_2D_Nanog_D_subset$Unique_track_2D_D==i,"Intensity.Sum_sum_norm"] -   
            min(na.omit(MCP_statistics_shape_total_2D_Nanog_D_subset[MCP_statistics_shape_total_2D_Nanog_D_subset$Unique_track_2D_D==i,"Intensity.Sum_sum_norm"])))/ 
           (max(na.omit(MCP_statistics_shape_total_2D_Nanog_D_subset[MCP_statistics_shape_total_2D_Nanog_D_subset$Unique_track_2D_D==i,"Intensity.Sum_sum_norm"])) -  
              min(na.omit(MCP_statistics_shape_total_2D_Nanog_D_subset[MCP_statistics_shape_total_2D_Nanog_D_subset$Unique_track_2D_D==i,"Intensity.Sum_sum_norm"])))) + 0
    }
  
  #Plot tracks
  {
   
    ggplot() + 
      geom_line(MCP_statistics_shape_total_2D_Nanog_D_subset,mapping=aes(y=as.numeric(Radial_distance_cv), 
                                  x=Seconds_normalized,group=Group), 
                color="black",size=0.25) + 
      facet_wrap(MCP_statistics_shape_total_2D_Nanog_D_subset$Unique_track_2D_D, ncol=7)+
      geom_vline(MCP_statistics_shape_total_2D_Nanog_D_subset,mapping=aes(xintercept = 0,), size=0.25,color="red", linetype=2)+
      geom_hline(MCP_statistics_shape_total_2D_Nanog_D_subset,mapping=aes(yintercept = 0.24,), size=0.25,color="forestgreen", linetype=2)+
      geom_hline(MCP_statistics_shape_total_2D_Nanog_D_subset,mapping=aes(yintercept = 0.16,), size=0.25,color="black", linetype=2)+
      
      geom_line(MCP_statistics_shape_total_2D_Nanog_D_subset,mapping=aes(y=as.numeric(Intensity.Sum_sum_norm_RD), 
                                  x=Seconds_normalized), 
                color="chartreuse3",size=0.25) +
      theme_bw() + xlab("Time (s)")+ #ggtitle('Merging Tracks')+
      scale_x_continuous(breaks = seq(-90, 90, by = 90),limits = c(-105,105)) +
      scale_y_continuous(breaks = seq(0, 0.5, by = 0.25),limits = c(0,0.5),labels=c(0,0.25,0.5),
                         name= "Coefficient of variation of the radial distance",
                         sec.axis=sec_axis(~.,name="Total Normalized Nanog Intensity",
                                           breaks = seq(0, 1, by = 0.5))) +
      theme(aspect.ratio=2/3, panel.grid.minor.x = element_blank(),
            panel.grid.minor.y = element_blank(),strip.text.x=element_blank(),
            axis.title.x = element_text(color = "black", size=10),
            axis.text.x = element_text(color = "black", size=7),
            axis.title.y = element_text(color = "black", size=10),
            axis.text.y = element_text(color = "black", size=7),
            axis.text.y.right = element_text(color = "chartreuse3"),
            axis.title.y.right = element_text(color = "chartreuse3", size=10))
  }

  }
  
}