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
Input_data = "Path_directory"

#############Associate Nanog to the corresponding Pol II S5P signal
{
#Calculate the distance between each pair of S5P-Nanog foci, and check if this distance if < 0.5 using the Shortest distance statistics from Imaris and < 1 for the distance calculated from the position of the Nanog spots. 
#If yes, join the two data frames at this specific line. If not, write just nothing.
#Add all columns from S5P file into Nanog file

setwd(Input_data)
Nanog_statistics_total <- read.csv('Nanog_statistics.csv')
S5P_statistics_total <- read.csv('S5P_statistics.csv')

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
  
  #Assemble the other Nanog dataframe with all the files to the one with only merging clusters
  
setwd(Input_data)
Nanog_statistics_total_sub <- read.csv('Nanog_statistics_subclusters.csv')
Nanog_statistics_total <- Nanog_statistics_total %>% full_join(Nanog_statistics_total_sub, by=c("Date","Stage","Stack","Nucleus","Time","ID","Seconds","Volume"))


#Normalization the level inside the nucleus for the mean intensity

Nucleus_statistics_total <- read.csv('Nucleus_statistics.csv')

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

#Export data frames used to plot figures
Nanog_statistics_total_10
Sum_parameters_plotting

  
}
