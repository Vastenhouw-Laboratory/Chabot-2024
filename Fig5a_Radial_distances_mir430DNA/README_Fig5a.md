# Overview

- The Figure 5a displays the coefficient of the radial distances of the mir430 DNA shapes for the whole dataset at 1k-cell stage. 
Briefly, the mir430 DNA signal was segmented in 3D using Imaris and the pixels contained in the segmented DNA mask were isolated.
Using this mask, distances between the center of the gravity and the center of edge pixels was retrieved and plotted centered on transcription initiation.
- Script developed by Damian Dalle Nogare at the Human Technopole BioImage Analysis Infrastructure Unit, National Facility for Data Handling and Analysis. Licensed under the BSD-3 license.

# Folder content

- `calculate_circle_distances.ipynb`
- `environment.yml`
- `LICENSE.txt`
- `radius_5.csv`
- `single_image_analysis.ipynb`
- `test_ouput.csv`
- `https://tinyurl.com/chabot2024testdata`

# System requirements (Hardware requirements, software dependencies and operating system)

Jupyter notebooks tested on Red Hat Enterprise Linux 8.6 on a AMD Epyc 7763 processor under Python 3.8.19 using jupyter-core 5.7.2

#  Installation guide

1. Place the scripts `single_image_analysis.ipynb`, `calculate_circle_distances.ipynb` and the `environment.yml` file into a folder.
2. Activate a terminal with conda installed and install the environment and dependencies using 
`conda env create -f environment.yml`
3. Activate the environment by typing `conda activate chabot_2024`
4. Launch a jupyter instance by typing `jupyter notebook`

# Demo
## Intructions to run on data (typical runtime on test data, 2.5 seconds, not counting image loading time)``

1. Download the example data from [`https://tinyurl.com/chabot2024testdata`](https://tinyurl.com/chabot2024testdata) and place the two files (`test_raw.ims` and `test_tracked.tif` in a folder. Note that this data will be moved after publication to a more stable repository (ie the bioimage archive).
2. In a jupyter notebook, open the `single_image_analysis.ipynb` notebook. 
3. Execute the first two cells to load the dependencies (1) and initialize functions (2)
4. Update the `raw_image` and `tracked_image` variables to contain the path to the raw image and tracked image respectively (downloaded above).
5. Run the analysis by executing the line `df = process_single_image(tracked_image, raw_image)`. The results will be stored in the resulting dataframe
6. The output of this script on the supplied test data should reproduce the results in `test_output.csv`

## Running the software on other data than this paper

This software was written to operate on the data generated for the current publication. In order to run this analysis on your own data, you must first ensure that the data is formatted in a manner consistent with the test data. In this case, you need two files

1.  A tracked label file, in the shape `TZYX` containing 3D segmented labels for each locus
2. A raw data file, in the shape `TCZYX`, containing the raw data. In this case, it is important that the 6th channel contains the segmented (binary) RNA transcription signal which will determine when transcription is initiated
3. Run the script, updating where appropriate the file paths

# Citation



