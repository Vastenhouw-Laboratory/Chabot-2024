# Noémie Chabot, PhD, Vastenhouw lab, 2024/07/10
# R code for analysis and figures related to the Fig3d,f
# "Local DNA compaction creates TF-DNA clusters that enable transcription"
# As of submission to the journal NSMB

#Libraries ----
#Read all the librairies that might be required to run/plot the code below
library("ggplot2")
library("viridisLite")
library("viridis")
library("tidyverse")
library("ggpubr")
library("cowplot")
library("ggpubr")
library("grid")
library("gridExtra")
library("dplyr")


options(scipen=999)

#=========================================================
# Read the files necessary for the analysis
#=========================================================
# Loading the MCP signal files for WT samples
 {
  
   Input = "/Volumes/CIG/nvastenh/nanog_movement/D2c/Manuscript/240612_NCB_NSMB_submission/Material_make_figures/Fig3d,f_Volume_mir430DNA_RNA_WT_inh/Stat/Statistics_WT/"  
   setwd(Input)
   
  #All file names in one table
  {
    #Change path to working directory
    
    #MCP
    #Get a list of files containing the word "MCP" to analyze only the MCP foci
  
    
    List_directories_MCP <- grep(pattern="MCP", list.dirs(path = ".", full.names = TRUE), value=TRUE)
    
    List_files_MCP_number <- data.frame()
    
    #Get all the MCP file in columns for each nucleus
    for(e in 1:length(List_directories_MCP))
    {
      print(length(list.files(List_directories_MCP[e])))
      if(e == 1) {List_files_MCP_number <- data.frame(list.files(path = List_directories_MCP[e], recursive=TRUE)) }
      else {
        #Assemble all data in a unique dataframe
        List_files_MCP_number <- List_files_MCP_number %>% data.frame(list.files(path = List_directories_MCP[e], pattern="MCP", recursive=TRUE))}
    }
  }
  
  #Assemble all data in one dataframe
  
  dataset <- data.frame()
  Infos_file <- data.frame()
  Infinity <- data.frame()
  
  for(e in 1:length(List_directories_MCP)) {
    
    setwd(Input)
    List_files_MCP_number <- list.files(path = List_directories_MCP[e], pattern="MCP", recursive=TRUE)
    setwd(List_directories_MCP[e])
    
    for (i in 1:length(List_files_MCP_number)){
      temp_data <- read.table(List_files_MCP_number[i],  sep=",", skip=3, header=TRUE) #each file will be read in, specify which columns you need read in to avoid any errors
      file_name <- as.character(List_files_MCP_number[i])
      
      #Select the info inside file name about date of acquisition and identity of the nucleus being analysing. These infos are stored in columns.
      #Adapt syntax based on the number of characters
      
      Infos_file[1,1] <- substr(file_name, 1,5)
      Infos_file[1,2] <- substr(file_name, 7,7)
      Infos_file[1,3] <- substr(file_name, 9,9)
      
      #Repeat n nows the infos about the nuclei, as much as they are segmented volumes in the file
      Infinity <- Infos_file[rep(seq_len(nrow(Infos_file)), each = nrow(temp_data)),]
      #Infinity$V5 <- sub("^\\D*(\\d+).*$", "\\1", List_directories_MCP[e])
      
      #Give a name to the column and make sure it looks like something
      colnames(Infinity) <- c("Date","Stack","Nucleus")
      print(List_files_MCP_number[i])
      temp_data <- cbind(Infinity, temp_data)
      
      #Assemble all data in a unique dataframe
      if(i == 1) {Spots_statistics_all <- temp_data } else {
        if(grepl("Time", file_name)==TRUE)
        {colnames(temp_data)[4] <- "Seconds"
        colnames(temp_data)[7] <- "Time"
        Spots_statistics_all <- Spots_statistics_all %>% full_join(temp_data, by=c("ID","Date","Nucleus","Stack","Time","Set.1"))} 
        else {Spots_statistics_all <- Spots_statistics_all %>% full_join(temp_data, by=c("ID","Date","Nucleus","Stack","Time","Set.1"))}
      }
      
    }
    #Remove useless columns
    patterns <- c("Collection*", "Unit*","Image*","Channel*","X.x","X.y","Event*","Category*","Merging.foci.x*","Merging.foci.y*")
    Spots_statistics_all <- Spots_statistics_all[,-grep(paste(patterns, collapse="|"), colnames(Spots_statistics_all))]
    
    if(e == 1) {Spots_statistics_total <- Spots_statistics_all
    } else {
      #Mix data frames with different number of columns by row
      Spots_statistics_total <- bind_rows(Spots_statistics_total, Spots_statistics_all)
    }
    
  }
  
  #Adressing one unique nucleus ID for each nucleus
  Unique_ID <- unique(Spots_statistics_total[,1:3])
  Spots_statistics_total$Nucleus.number <- 0
  
  for(i in 1:nrow(Spots_statistics_total))
  {
    for(e in 1:nrow(Unique_ID))
    {
      if(Spots_statistics_total[i,1] == Unique_ID[e,1]&Spots_statistics_total[i,2] == Unique_ID[e,2]&Spots_statistics_total[i,3] == Unique_ID[e,3])
        
        Spots_statistics_total[i,"Nucleus.number"] <- e
      
    }
  }
  
  #Clean dataframe with data not associated with a class
  Spots_statistics_total <- Spots_statistics_total[Spots_statistics_total$Set.1 %in% c("Class A","Class B"),]
  Spots_statistics_total[Spots_statistics_total$Set.1 == "Class A","Set.1"] <- 0
  Spots_statistics_total[Spots_statistics_total$Set.1 == "Class B","Set.1"] <- 1
  
  Spots_statistics_total[which(Spots_statistics_total$Set.1=="0"),"Unique_track"] <- Spots_statistics_total[which(Spots_statistics_total$Set.1=="0"),"Nucleus.number"]
  Spots_statistics_total[which(Spots_statistics_total$Set.1=="1"),"Unique_track"] <- Spots_statistics_total[which(Spots_statistics_total$Set.1=="1"),"Nucleus.number"]*30

  #Sum-up the MOVIE surfaces that have been detected in two different surfaces
  for(i in 1:max(Spots_statistics_total$Unique_track))
  {
    if(nrow(Spots_statistics_total[Spots_statistics_total$Unique_track==i&Spots_statistics_total$Time==e,])!=0)
    {
      for(e in 1:max(Spots_statistics_total[Spots_statistics_total$Unique_track==i,"Time"]))
      {
        for(k in 1:max(as.numeric(Spots_statistics_total[Spots_statistics_total$Unique_track==i&Spots_statistics_total$Time==e,"Set.1"])))
        {
          if(nrow(Spots_statistics_total[which(Spots_statistics_total$Unique_track==i&Spots_statistics_total$Time==e&Spots_statistics_total$Set.1==k),])==2)
          {
            Spots_statistics_total[Spots_statistics_total$Unique_track==i&Spots_statistics_total$Time==e&Spots_statistics_total$Set.1==k,"Intensity.Mean.x"][1] <- mean(Spots_statistics_total[Spots_statistics_total$Unique_track==i&Spots_statistics_total$Time==e&Spots_statistics_total$Set.1==k,"Intensity.Mean.x"])
            Spots_statistics_total[Spots_statistics_total$Unique_track==i&Spots_statistics_total$Time==e&Spots_statistics_total$Set.1==k,"Intensity.Mean.y"][1] <-   mean(Spots_statistics_total[Spots_statistics_total$Unique_track==i&Spots_statistics_total$Time==e&Spots_statistics_total$Set.1==k,"Intensity.Mean.y"])
            Spots_statistics_total[Spots_statistics_total$Unique_track==i&Spots_statistics_total$Time==e&Spots_statistics_total$Set.1==k,"Intensity.Sum.x"][1] <-       sum(Spots_statistics_total[Spots_statistics_total$Unique_track==i&Spots_statistics_total$Time==e&Spots_statistics_total$Set.1==k,"Intensity.Sum.x"])
            Spots_statistics_total[Spots_statistics_total$Unique_track==i&Spots_statistics_total$Time==e&Spots_statistics_total$Set.1==k,"Intensity.Sum.y"][1] <- sum(Spots_statistics_total[Spots_statistics_total$Unique_track==i&Spots_statistics_total$Time==e&Spots_statistics_total$Set.1==k,"Intensity.Sum.y"])
            Spots_statistics_total[Spots_statistics_total$Unique_track==i&Spots_statistics_total$Time==e&Spots_statistics_total$Set.1==k,"Intensity.Sum.of.Square.x"][1] <- sum(Spots_statistics_total[Spots_statistics_total$Unique_track==i&Spots_statistics_total$Time==e&Spots_statistics_total$Set.1==k,"Intensity.Sum.of.Square.x"])
            Spots_statistics_total[Spots_statistics_total$Unique_track==i&Spots_statistics_total$Time==e&Spots_statistics_total$Set.1==k,"Intensity.Sum.of.Square.y"][1] <-   sum(Spots_statistics_total[Spots_statistics_total$Unique_track==i&Spots_statistics_total$Time==e&Spots_statistics_total$Set.1==k,"Intensity.Sum.of.Square.y"])
            Spots_statistics_total[Spots_statistics_total$Unique_track==i&Spots_statistics_total$Time==e&Spots_statistics_total$Set.1==k,"Volume"][1]<- sum(Spots_statistics_total[Spots_statistics_total$Unique_track==i&Spots_statistics_total$Time==e&Spots_statistics_total$Set.1==k,"Volume"])
          }
        }
      }
    }
  }
  
  Spots_statistics_total <- Spots_statistics_total %>% distinct(Date,Stack,Nucleus,Time,Set.1, .keep_all = TRUE)
}

# Loading the MCP signal files for inhibited samples
 {
  #All file names in one table
  {
    Input = "/Volumes/CIG/nvastenh/nanog_movement/D2c/Manuscript/240612_NCB_NSMB_submission/Material_make_figures/Fig3d,f_Volume_mir430DNA_RNA_WT_inh/Stat/Statistics_Inh/"  
    setwd(Input)
    
    
    #MCP
    #Get a list of files containing the word "MCP" to analyze only the MCP foci
    List_directories_MCP_Inh <- grep(pattern="MCP", list.dirs(path = ".", full.names = TRUE), value=TRUE)
    List_files_MCP_number_Inh <- data.frame()
    
    #Get all the MCP file in columns for each nucleus
    
    for(e in 1:length(List_directories_MCP_Inh))
    {
      setwd(Input)      
      print(length(list.files(List_directories_MCP_Inh[e])))
      if(e == 1) {List_files_MCP_number_Inh <- data.frame(list.files(path = List_directories_MCP_Inh[e], recursive=TRUE)) }
      else {
        #Assemble all data in a unique dataframe
        List_files_MCP_number_Inh <- List_files_MCP_number_Inh %>% data.frame(list.files(path = List_directories_MCP_Inh[e], pattern="MCP", recursive=TRUE))}
    }
    
  }
  
  #Assemble all data in one dataframe
  
  dataset <- data.frame()
  Infos_file <- data.frame()
  Infinity <- data.frame()
  
  for(e in 1:length(List_directories_MCP_Inh)) {
    
    setwd(Input) 
    List_files_MCP_number_Inh <- list.files(path = List_directories_MCP_Inh[e], pattern="MCP", recursive=TRUE)
    setwd(List_directories_MCP_Inh[e])
    
    for (i in 1:length(List_files_MCP_number_Inh)){
      temp_data <- read.table(List_files_MCP_number_Inh[i],  sep=",", skip=3, header=TRUE) #each file will be read in, specify which columns you need read in to avoid any errors
      file_name <- as.character(List_files_MCP_number_Inh[i])
      
      #Select the info inside file name about date of acquisition and identity of the nucleus being analysing. These infos are stored in columns.
      #Adapt syntax based on the number of characters
      
      Infos_file[1,1] <- substr(file_name, 1,5)
      Infos_file[1,2] <- substr(file_name, 7,7)
      Infos_file[1,3] <- substr(file_name, 9,9)
      
      #Repeat n nows the infos about the nuclei, as much as they are segmented volumes in the file
      Infinity <- Infos_file[rep(seq_len(nrow(Infos_file)), each = nrow(temp_data)),]
      #Infinity$V5 <- sub("^\\D*(\\d+).*$", "\\1", List_directories_MCP[e])
      
      #Give a name to the column and make sure it looks like something
      colnames(Infinity) <- c("Date","Stack","Nucleus")
      print(List_files_MCP_number_Inh[i])
      temp_data <- cbind(Infinity, temp_data)
      
      #Assemble all data in a unique dataframe
      if(i == 1) {MCP_inh_statistics_all <- temp_data } else {
        if(grepl("Time", file_name)==TRUE)
        {colnames(temp_data)[4] <- "Seconds"
        colnames(temp_data)[7] <- "Time"
        MCP_inh_statistics_all <- MCP_inh_statistics_all %>% full_join(temp_data, by=c("ID","Date","Nucleus","Stack","Time","Set.1"))} 
        else {MCP_inh_statistics_all <- MCP_inh_statistics_all %>% full_join(temp_data, by=c("ID","Date","Nucleus","Stack","Time","Set.1"))}
      }
      
    }
    #Remove useless columns
    patterns <- c("Collection*", "Unit*","Image*","Channel*","X.x","X.y","Event*","Category*","Merging.foci.x*","Merging.foci.y*","TrackID*")
    MCP_inh_statistics_all <- MCP_inh_statistics_all[,-grep(paste(patterns, collapse="|"), colnames(MCP_inh_statistics_all))]
    
    if(e == 1) {MCP_inh_statistics_total <- MCP_inh_statistics_all
    } else {
      #Mix data frames with different number of columns by row
      MCP_inh_statistics_total <- bind_rows(MCP_inh_statistics_total, MCP_inh_statistics_all)
    }
    
  }
  
  #Adressing one unique nucleus ID for each nucleus
  
  Unique_ID <- unique(MCP_inh_statistics_total[,1:3])
  MCP_inh_statistics_total$Nucleus.number <- 0
  
  for(i in 1:nrow(MCP_inh_statistics_total))
  {
    for(e in 1:nrow(Unique_ID))
    {
      if(MCP_inh_statistics_total[i,1] == Unique_ID[e,1]&MCP_inh_statistics_total[i,2] == Unique_ID[e,2]&MCP_inh_statistics_total[i,3] == Unique_ID[e,3])
        
        MCP_inh_statistics_total[i,"Nucleus.number"] <- e
      
    }
  }
  
  #Clean dataframe with data not associated with a class
  
  MCP_inh_statistics_total <- MCP_inh_statistics_total[MCP_inh_statistics_total$Set.1 %in% c("Class A","Class B"),]
  MCP_inh_statistics_total[MCP_inh_statistics_total$Set.1 == "Class A","Set.1"] <- 1
  MCP_inh_statistics_total[MCP_inh_statistics_total$Set.1 == "Class B","Set.1"] <- 2
  
  MCP_inh_statistics_total[which(MCP_inh_statistics_total$Set.1=="1"),"Unique_track"] <- MCP_inh_statistics_total[which(MCP_inh_statistics_total$Set.1=="1"),"Nucleus.number"]
  MCP_inh_statistics_total[which(MCP_inh_statistics_total$Set.1=="2"),"Unique_track"] <- MCP_inh_statistics_total[which(MCP_inh_statistics_total$Set.1=="2"),"Nucleus.number"]*30
  

  #Sum-up the MOVIE surfaces that have been detected in two different surfaces
  
  for(i in 1:max(MCP_inh_statistics_total$Unique_track))
  {
    
    if(nrow(MCP_inh_statistics_total[MCP_inh_statistics_total$Unique_track==i&MCP_inh_statistics_total$Time==e,])!=0)
    {
      
      for(e in 1:max(MCP_inh_statistics_total[MCP_inh_statistics_total$Unique_track==i,"Time"]))
      {
        for(k in 1:max(as.numeric(MCP_inh_statistics_total[MCP_inh_statistics_total$Unique_track==i&MCP_inh_statistics_total$Time==e,"Set.1"])))
        {
          
          if(nrow(MCP_inh_statistics_total[which(MCP_inh_statistics_total$Unique_track==i&MCP_inh_statistics_total$Time==e&MCP_inh_statistics_total$Set.1==k),])==2)
          {
            MCP_inh_statistics_total[MCP_inh_statistics_total$Unique_track==i&MCP_inh_statistics_total$Time==e&MCP_inh_statistics_total$Set.1==k,"Intensity.Mean.x"][1] <- mean(MCP_inh_statistics_total[MCP_inh_statistics_total$Unique_track==i&MCP_inh_statistics_total$Time==e&MCP_inh_statistics_total$Set.1==k,"Intensity.Mean.x"])
            MCP_inh_statistics_total[MCP_inh_statistics_total$Unique_track==i&MCP_inh_statistics_total$Time==e&MCP_inh_statistics_total$Set.1==k,"Intensity.Mean.y"][1] <-   mean(MCP_inh_statistics_total[MCP_inh_statistics_total$Unique_track==i&MCP_inh_statistics_total$Time==e&MCP_inh_statistics_total$Set.1==k,"Intensity.Mean.y"])
            MCP_inh_statistics_total[MCP_inh_statistics_total$Unique_track==i&MCP_inh_statistics_total$Time==e&MCP_inh_statistics_total$Set.1==k,"Intensity.Sum.x"][1] <-       sum(MCP_inh_statistics_total[MCP_inh_statistics_total$Unique_track==i&MCP_inh_statistics_total$Time==e&MCP_inh_statistics_total$Set.1==k,"Intensity.Sum.x"])
            MCP_inh_statistics_total[MCP_inh_statistics_total$Unique_track==i&MCP_inh_statistics_total$Time==e&MCP_inh_statistics_total$Set.1==k,"Intensity.Sum.y"][1] <- sum(MCP_inh_statistics_total[MCP_inh_statistics_total$Unique_track==i&MCP_inh_statistics_total$Time==e&MCP_inh_statistics_total$Set.1==k,"Intensity.Sum.y"])
            MCP_inh_statistics_total[MCP_inh_statistics_total$Unique_track==i&MCP_inh_statistics_total$Time==e&MCP_inh_statistics_total$Set.1==k,"Intensity.Sum.of.Square.x"][1] <- sum(MCP_inh_statistics_total[MCP_inh_statistics_total$Unique_track==i&MCP_inh_statistics_total$Time==e&MCP_inh_statistics_total$Set.1==k,"Intensity.Sum.of.Square.x"])
            MCP_inh_statistics_total[MCP_inh_statistics_total$Unique_track==i&MCP_inh_statistics_total$Time==e&MCP_inh_statistics_total$Set.1==k,"Intensity.Sum.of.Square.y"][1] <-   sum(MCP_inh_statistics_total[MCP_inh_statistics_total$Unique_track==i&MCP_inh_statistics_total$Time==e&MCP_inh_statistics_total$Set.1==k,"Intensity.Sum.of.Square.y"])
            MCP_inh_statistics_total[MCP_inh_statistics_total$Unique_track==i&MCP_inh_statistics_total$Time==e&MCP_inh_statistics_total$Set.1==k,"Volume"][1]<- sum(MCP_inh_statistics_total[MCP_inh_statistics_total$Unique_track==i&MCP_inh_statistics_total$Time==e&MCP_inh_statistics_total$Set.1==k,"Volume"])
          }
        }
      }
    }
    
    
  }
  
  MCP_inh_statistics_total <- MCP_inh_statistics_total %>% distinct(Date,Stack,Nucleus,Time,Set.1, .keep_all = TRUE)
  
}

# Loading the MOVIE signal in WT samples
 {
  
  #All file names in one table
  {
    Input = "/Volumes/CIG/nvastenh/nanog_movement/D2c/Manuscript/240612_NCB_NSMB_submission/Material_make_figures/Fig3d,f_Volume_mir430DNA_RNA_WT_inh/Stat/Statistics_WT/"  
    setwd(Input)    
    #MOVIE
    #Get a list of files containing the word "MOVIE" to analyze only the MOVIE foci
  
    List_directories_MOVIE <- grep(pattern="MOVIE", list.dirs(path = ".", full.names = TRUE), value=TRUE)
    
    List_files_MOVIE_number <- data.frame()
    
    #Get all the MOVIE file in columns for each nucleus
    
    for(e in 1:length(List_directories_MOVIE))
    {
      setwd(Input) 
      
      print(length(list.files(List_directories_MOVIE[e])))
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
      Infos_file[1,2] <- substr(file_name, 7,7)
      Infos_file[1,3] <- substr(file_name, 9,9)
      
      #Repeat n nows the infos about the nuclei, as much as they are segmented volumes in the file
      Infinity <- Infos_file[rep(seq_len(nrow(Infos_file)), each = nrow(temp_data)),]
      #Infinity$V5 <- sub("^\\D*(\\d+).*$", "\\1", List_directories_MOVIE[e])
      
      #Give a name to the column and make sure it looks like something
      colnames(Infinity) <- c("Date","Stack","Nucleus")
      print(List_files_MOVIE_number[i])
      temp_data <- cbind(Infinity, temp_data)
      
      #Assemble all data in a unique dataframe
      if(i == 1) {MOVIE_statistics_all <- temp_data } else {
        if(grepl("Time", file_name)==TRUE)
        {colnames(temp_data)[4] <- "Seconds"
        colnames(temp_data)[7] <- "Time"
        MOVIE_statistics_all <- MOVIE_statistics_all %>% full_join(temp_data, by=c("ID","Date","Nucleus","Stack","Time","TrackID"))} 
        else {MOVIE_statistics_all <- MOVIE_statistics_all %>% full_join(temp_data, by=c("ID","Date","Nucleus","Stack","Time","TrackID"))}
      }
      
    }
    #Remove useless columns
    patterns <- c("Collection*", "Unit*","Image*","Channel*","X.x","X.y","Event*","Category*","Merging.foci.x*","Merging.foci.y*")
    MOVIE_statistics_all <- MOVIE_statistics_all[,-grep(paste(patterns, collapse="|"), colnames(MOVIE_statistics_all))]
    
    if(e == 1) {MOVIE_statistics_total <- MOVIE_statistics_all
    } else {
      #Mix data frames with different number of columns by row
      MOVIE_statistics_total <- bind_rows(MOVIE_statistics_total, MOVIE_statistics_all)
    }
    
  }
  
  #Adressing one unique nucleus ID for each nucleus
  
  Unique_ID <- unique(MOVIE_statistics_total[,1:3])
  MOVIE_statistics_total$Nucleus.number <- 0
  
  for(i in 1:nrow(MOVIE_statistics_total))
  {
    for(e in 1:nrow(Unique_ID))
    {
      if(MOVIE_statistics_total[i,1] == Unique_ID[e,1]&MOVIE_statistics_total[i,2] == Unique_ID[e,2]&MOVIE_statistics_total[i,3] == Unique_ID[e,3])
        
        MOVIE_statistics_total[i,"Nucleus.number"] <- e
      
    }
  }
  
  #Clean dataframe with data not associated with a class
  
  MOVIE_statistics_total <- MOVIE_statistics_total[complete.cases(MOVIE_statistics_total), ]
  
  #Adressing one unique TrackID for each nucleus
  
  Unique_ID <- unique(MOVIE_statistics_total[,c(1:3,6)])
  MOVIE_statistics_total$Unique_track <- 0
  
  for(i in 1:nrow(MOVIE_statistics_total))
  {
    for(e in 1:nrow(Unique_ID))
    {
      if(MOVIE_statistics_total[i,1] == Unique_ID[e,1]&MOVIE_statistics_total[i,2] == Unique_ID[e,2]&MOVIE_statistics_total[i,3] == Unique_ID[e,3]&MOVIE_statistics_total[i,6] == Unique_ID[e,4])
        
        MOVIE_statistics_total[i,"Unique_track"] <- e
      
    }
  }
}

#=========================================================
# Plots
#=========================================================

# Fig. 3d - Plot MOVIE + MCP volume in average in the same time
{
  ggplot() + 
    geom_line(data=Spots_statistics_total, 
              aes(x=Time*2, y=Volume, group=as.factor(Unique_track)), 
              color="forestgreen", 
              size=0.7, 
              alpha=0.2) + 
    geom_line(data=MOVIE_statistics_total, aes(x=Time*2, y=Volume, group=as.factor(Unique_track)), 
              color="magenta3",
              size=0.7, 
              alpha=0.2) +
    stat_summary(Spots_statistics_total, mapping=aes(y=Volume, x=Time*2), fun=mean, geom="line", size=0.7, color="darkgreen") +
    stat_summary(MOVIE_statistics_total, mapping=aes(y=Volume, x=Time*2), fun=mean, geom="line", size=0.7, color="darkmagenta") +
    labs(x= "Time (min)",y="Volume (um3)") +  ylim(c(0,17.5))+xlim(c(0,41))+theme_bw()+
    theme(aspect.ratio=1, axis.text = element_text(size=14),axis.title = element_text(size=14))
  
}

# Fig. 3f - Plot  MCP volume inhibition in average in the same time
{
  ggplot() + 
    geom_line(data=MCP_inh_statistics_total, aes(x=Time*2, y=Volume, group=as.factor(Unique_track)), 
              color="forestgreen", 
              size=0.7, 
              alpha=0.2) + 
    stat_summary(MCP_inh_statistics_total, mapping=aes(y=Volume, x=Time*2), 
                 fun=mean, geom="line", 
                 size=0.7, 
                 color="darkgreen") +
    labs(x= "Time (min)",y="Volume (um3)")  + ylim(c(0,17.5))+ xlim(c(0,41))+theme_bw()+
    theme(aspect.ratio=1, axis.text = element_text(size=14),axis.title = element_text(size=14))
}