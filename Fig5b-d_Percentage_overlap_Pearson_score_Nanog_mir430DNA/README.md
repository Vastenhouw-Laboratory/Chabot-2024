# Overview

The Figure 5a displays the coefficient of the radial distances of the mir430 DNA shapes for the whole dataset at 1k-cell stage. 
Briefly, the mir430 DNA signal was segmented in 3D using Imaris and the pixels contained in the segmented DNA mask were isolated.
Using this mask, distances between the center of the gravity and the center of edge pixels was retrieved. 

# Folder content

Input files:
- `Filelist_ShortestDist_MCP.txt`
- `ImageInfo_List.txt`
- `MasterTable_AllNuclei_AlleleNum_TrackID_TxnTime_30Oct2023.csv`
- Two demo images : `11.01_1k_6_5.ims.tif` and `12.01_1k_9_1.ims.tif`

Scripts for data management and analysis:
- 1_Parse_FusionTime_NumClusterOvlp.pl
- 2_Calculate_Pearson_OverlapVolume.m
- 3_Align_OverlapVol_Pearson_TxnStart_Merging.pl

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



# Demo
## Intructions to run on data
## Expected output and run time


# Citation

