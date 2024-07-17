# Overview

The Figure 5b-d displays information about the percentage of overlap between the Nanog and DNA masks. Briefly, Nanog and mir430 DNA signal have been segmented using Imaris, and the resulting masks resulting were isolated in two channels in the same image file. The following script reads already available information about the data, as well as the images, and calculate:
- The percentage of overlap between Nanog and mir430 DNA masks in comparison to the size of the Nanog mask
- The percentage of overlap between Nanog and mir430 DNA masks in comparison to the size of the mir430 DNA mask
- The Pearson's correlation score between Nanog and mir430 DNA raw pixels inside the mir430 DNA mask.

# Folder content

Input files:
- `Filelist_ShortestDist_MCP.txt`
- `ImageInfo_List.txt`
- `MasterTable_AllNuclei_AlleleNum_TrackID_TxnTime_30Oct2023.csv`
- Two demo images : `11.01_1k_6_5.ims.tif` and `12.01_1k_9_1.ims.tif`

Scripts for data management and analysis:
- `1_Parse_FusionTime_NumClusterOvlp.pl`
- `2_Calculate_Pearson_OverlapVolume.m`
- `3_Align_OverlapVol_Pearson_TxnStart_Merging.pl`

Output files:
- `AllNuclei_Allele_All_Values_Pooled.txt`
- `AllNuclei_TrackID_FusionTimes_NumNanogClust-Ovlp-Cutoff_RadiusZ.txt`
- `OverlapVolume_Matlab.txt`
- `PearsonCoeff_Matlab_actual_scramble.txt`

# System requirements
## Hardware requirements

No specific hardware requirements. 

## Software dependencies and operating systems

Operating system: The system used for the analysis is MacOS Sonoma 14.4. However, any other system that can support MatLab and Perl can be used.
Softwares and dependencies: The following softwares were used for these analysis:

- 1.Matlab 2023a with following toolboxes:
  - Image Processing Toolbox
  - Statistics and Machine Learning Toolbox
- 2.Perl v5.34.0

#  Installation guide

1. Place the input files in the same folder: Input files: `Filelist_ShortestDist_MCP.txt`, `ImageInfo_List.txt`, `MasterTable_AllNuclei_AlleleNum_TrackID_TxnTime_30Oct2023.csv`,`11.01_1k_6_5.ims.tif` and `12.01_1k_9_1.ims.tif`.
2. Open a command line window
3. Run the commands as described below.

# Demo
## Intructions to run on data (typical runtime on test data, 5 seconds, not counting image loading time)

There are three pieces of code that are numbered from 1 to 3 and have to run in this specific order.
The Perl code 1_Parse_FusionTime_NumClusterOvlp.pl calculates (a) the number of nanog clusters at each time point (b) number of nanog clusters overlapping with miR430 DNA at each time point, and (c) the time slices when nanog clusters merge for each miR430 allele.
 
To launch the different scripts for the analysis:
```
	cd /path/
	perl 1_Parse_FusionTime_NumClusterOvlp.pl
```

This file is doing the following function:

###1.Reading Input File List:
The script starts by reading a list of filenames from `Filelist_ShortestDist_MCP.txt`. These filenames correspond to CSV files that contain shortest distance data to miR430 DNA surfaces as exported from Imaris software.

###2.Extracting Nuclei IDs and Preparing File Paths:
For each filename, the script extracts a nuclei ID and constructs another filename for diameter data by modifying the path and filename.

###3.Processing Diameter Data: he script reads the diameter data from the constructed filename. It extracts the smallest diameter values and stores them in a hash %cutoff keyed by a certain identifier.

###4.Processing Shortest Distance Data: It opens the original shortest distance file again and reads the data, extracting track IDs and time points. To judge whether the clusters are overlapping, it compares the shortest distance values to the corresponding cutoff values (as dictated by the shortest diameter of the cluster).
Overlap information and counts of spots (entities being tracked) are stored in nested hashes %ovlp_MCP and %spotnum.

###5.Generating Output File: The script writes the processed data to an output file `AllNuclei_TrackID_FusionTimes_NumNanogClust-Ovlp-Cutoff_RadiusZ.txt`.

The MATLAB code `2_Calculate_Pearson_OverlapVolume.m calculates` (a) Volumetric overlap between Nanog and miR430-DNA signals (b) Pearson correlation between Nanog and mir430-DNA signal intensities within miR430-DNA mask. It also scrambles Nanog signal values within miR430-DNA mask and re-calculates the Pearson correlation between scrambled Nanog and miR430-DNA signal intenstities as negative control 

In MATLAB interface, 
```
	code 2_Calculate_Pearson_OverlapVolume.m
```

This MATLAB script performs the following tasks:

###1.Read Input Data:
It reads data from `ImageInfo_List.txt` using tdfread, which contains information such as image file names, number of frames, slices, channels, and alleles for each nuclei. 
Each image file contains three original channels corresponding to Nanog, miR430-DNA and MOVIE. Then C4 and C5 correspond to masks for allele1 and allele2 in Nanog channel, C6 and C7 correspond to masks for allele1 and allele2 in miR430_DNA channel whereas C8 and C9 correspond to masks for allele1 and allele2 in MoViE channel. For single allele cases the channels C4,C5 and C6 correspond to masks in Nanog, miR430-DNA and MOVIE channel respectively.

###2.Open Output Files:
It opens two files for writing:
`PearsonCoeff_Matlab_actual_scramble.txt` to store Pearson correlation coefficients.
OverlapVolume_Matlab.txt to store overlap volume data.

###3.Write Headers to Output Files:
It writes headers to both output files for better organization and readability.

###4.Iterate Over Each File:
For each image file in the list, it initializes a counter and an empty matrix to store image data.

###5.Load Image Data:
For each frame, slice, and channel, it loads the image data into the matrix.
Demo images can be found at the following adress:
`12.01_1k_9_1.ims.tif`: `https://tinyurl.com/3htdmh4j` and `11.01_1k_6_5.ims.tif`: `https://tinyurl.com/y3kce7b8`.

###6.Process Each Frame:
For each frame, it initializes variables for storing Nanog and DNA data and overlap volumes for each allele.

###7.Process Each Slice:
For each slice, it processes the image data to:
Extract masks for Nanog, DNA, and RNA.
Calculate volumes for Nanog, DNA, RNA, and overlap between Nanog and DNA.
Accumulate this data for all slices in the current frame.

###8.Calculate Pearson Correlation Coefficients:
It calculates the Pearson correlation coefficient between Nanog and DNA signal intensities within DNA mask for the first allele. For negative control, it scrambles the Nanog intensity values within DNA mask and recalculates the correlation coefficient. It also performs bootstrap ten times and reoports the average scrambled correlation. In the figure, the non-bootstrapped values are represented.
If there is a second allele, it repeats these calculations for the second allele.

###9.Write Data to Output Files:
It writes the calculated overlap volumes and Pearson correlation coefficients (both actual and scrambled) to the respective output files.

###10. Log Processed File:
It prints the filename to the console for tracking progress.

The Perl script `3_Align_OverlapVol_Pearson_TxnStart_Merging.pl` parses the output files as generated in previous two codes along with additional information about transcription initiation, consolidating it into a master table. 

In bash,
```
	perl 3_Align_OverlapVol_Pearson_TxnStart_Merging.pl
```

This Perl script processes data from four input files and generates a comprehensive output file with aggregated and computed data. Let's break down the purpose of the script and the roles of each input and output file:

Function of the Script

###1. The script reads data from four input files:
`OverlapVolume_Matlab.txt`: Contains volume data for Nanog, DNA, RNA, and their overlaps.
`PearsonCoeff_Matlab_actual_scramble.txt`: Contains Pearson correlation coefficients for actual and scrambled data.
`AllNuclei_TrackID_FusionTimes_NumNanogClust-Ovlp-Cutoff_RadiusZ.txt`: Contains information about nuclei ID, track ID, fusion times, and the number of Nanog clusters over time.
`MasterTable_AllNuclei_AlleleNum_TrackID_TxnTime_30Oct2023.csv`: Contains master table information with nuclei ID, allele number, track ID, and time of activation (txn time).

###2. The script processes these files to:
Extract relevant data.
Compute relationships between different metrics (e.g., merging times, relative times).
Aggregate all the information into a single output file.

###3. Output File:
`AllNuclei_Allele_All_Values_Pooled.txt`:

## Running the software on other data than this paper

This pieces of code was written to operate on the data generated for this publication. If one wants to run these analysis on new dataset, one needs:
1. A .csv file containing all the information about the nuclei to be analyzed named date, ID, time points, TrackID
2. Images for which signal of interest has been already segmented and isolated as another new channel, in addition to the raw pixel values in the same file. These images needs to be 3D time-lapse and contain at least two channels with objects segmented in both channels and some overlap between the objects of the two channels.


