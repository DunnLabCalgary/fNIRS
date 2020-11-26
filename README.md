# fNIRS Pipeline
Matlab-based user friendly Functional Near Infrared Spectroscopy Brain Imaging Data Processing Pipeline

## Software requirements 
* MATLAB R2019a or newer
* Curve Fitting Toolbox
* Signal Processing Toolbox
## New requirements: 
* Image processing Toolbox


## Sections: 

##### Quick Test
##### Overview
##### Required Folder Structure
 ##### Initialization
  ##### Using the Stim Builder
##### Preprocessing
##### Postprocessing
  ##### Channel Selection
##### Analysis
##### Appendix and Sample Figures

## Quick Test with the Data Provided:

1. Open the fNIRSPipeline.mlapp file the start the pipeline

2. Initialization:

* The selected output folder should be the **fNIRS2** folder
* The Raw Data folder should be the **Data** folder

Some subjects in this data set are already initialized, they will appear in the Initialized subjects list.

You can then run Preprocessing, Postprocessing or Analysis individually or all at once via the Run All Pane. Some configuration of each process is required within their respective panes before running.

**To run with other data:** Simply replace the relevant files: A Data folder with a folder containing your raw data and a demographics file. *Important:* Keep the folder struture the same. eg: your Data folder should still be contained one level deep within the study output folder (fNIRS2).

## Overview: 

The pipeline consists of three processing phases, pre-processing, post-processing and analysis with an initialization phase preceding them. Each phase can be configured and run separately in their respective panes or they can be configured all at once and then run consecutively in the ‘Run All’ pane. 

## Required Folder Structure: 

To ensure the pipeline can locate the necessary data the following folder structure should be followed:

<div align="center">
<img src="https://github.com/carterrandall/fNIRS/blob/master/images/folderstructure.png" height = 400>
</div>

### ~~Demographics File:   ~~
#### Age is now entered in the Initialization table of the GUI.
An optional excel file with subject demographics information in the following format: format coming soon - for now just view provided demographics file
Age is used in the calculation of the lights differential path length through the brain. See preprocessing - converting the raw signal to concentration for more info. *If no demographics file is provided an age of 30 is assumed.*



### Probe File: 
A .sd file describing the source detector setup (probe) that was used in recording the raw data. This is typically automatically generated when recording data. If you do not have a .sd file (see Appendix **Creating a probe (.sd) file ** for instructions on how to generate one.)

### Raw Data Folder: 
This folder should only contain the raw data .nirs or .nirx files of each subject.

### Processed: 
All data output from the pipeline will be found here in subject-wise output folders. This folder will be created for you if it doesn’t exist already and should not be renamed. 

## Initialization:
The initialization phase contains the majority of setup work required for running a study. Here you will choose the desired folder outputs and build stim designs for each subject.

<div align="center">
<img src="https://github.com/carterrandall/fNIRS/blob/master/images/Initialization_new.png" height = 400>
</div>


#### A: Select Data Type and Folder
Select your data type and folders in accordance with the Required Folder Structure section. This step must be completed before anything else can be done.

#### B: Initialization Table
After selecting your raw data folder, subjects will appear in the initialization table. Here you can indicate the age of each subject under the 'Age' column. If none is specified a default age of 30 is used. Subjects must also first be initialized before any data processing. Clicking the Initialize button under the 'Action' column will open the stim builder. Once subjects have been initialized using the stim builder, they will be marked as such under the 'Stims Built' column.

## Using the Stim Builder
The stim view will populate with the stim-time relation according to the subjects raw data file. Vertical blue lines represent either the start of stop of stims. Since these rarely perfectly match data collection the stim builder allows you to trim the data and or shift the start/end of stim times.


<div align="center">
<img src="https://github.com/carterrandall/fNIRS/blob/master/images/stims.png" height = 400>
</div>
  
 ### Designing a stim from scratch:

1. Mark the desired start and end of the data. These can be adjusted by either dragging left and right or more precisely via the edit fields next to their respective button.
2. Mark a stim by clicking the *Mark Stim* button. This will appear as a green box in the stim view. Adjust the start and end times by dragging or more precisely in the edit fields. The names of each stim can be edited and will be used in figures and file names.
3. (Optional) Once your stim design is complete, save it as a template by pressing “Save as Template” for future use.
4. Click Finish Subject. If there are more subjects to be initialized, the stim builder will reopen automatically ready for the next subject. If you saved your design as a template, read below to see how to initialize subjects faster.

### Designing a stim from an existing template: 

If you have already built a template, click Load Template and open it. You can then adjust each component to the subject as needed and then click finish subject to move to the next subject or finish Initialization.

## Preprocessing

This is the first of the data processing phases and will prepare your data for the next two. Preprocessing can be performed on a subject before it is initialized. It’s main objectives are:

### 1 Pruning poor-quality data:

Remove data from channels that does not meet criteria specified by the user. Channel pruning is performed using the homer2 “hmrDataQualityCheck”.

### 2 Converting raw data to optical density:

The subjects age provided in the demographics file is used (otherwise an age of 30 is assumed) in calculating the length of the light path or the differential pathlength factor (DPF). The DPF is a key factor in the Modified Beer Lambert equation which is used to calculate hemoglobin concentration.

### 3 Motion and artifact detection/correction: 

Automated motion correction is done via the hmrMotionCorrectSplineSG function from the homer2 package. Motion is detected with hmrtInc_baselineshift_Ch which detects motion based on std variations, gradient outliers, baseline shift and spikes, in consideration of heart rate amplitude, and then corrected using Spline Interpolation.

Manual motion correction is performed with MARA (Currently not fully implemented, use automated for now). How to detect and reduce movement artifacts in near-infrared imaging using moving standard deviation and spline interpolation."Physiological Measurement, 31, 649-662. Scholkmann et al. 2010. and requires the signal processing toolbox ‘MATLAB_ToolsForUCL_Felix_2014'

### 4 Convert Optical Density to Hemoglobin Concentration:
Performed using the hmrOD2Conc function.

### 5 De-trend, z-score normalize, and Downsample (optional)

<div align="center">
<img src="https://github.com/carterrandall/fNIRS/blob/master/images/preprocessing.png" height = 400>
</div>

#### A: Check Data Quality:

Determine the criteria used to prune channels.

* **Amplitude:** A dB range that specifies the min and max allowable dB value for a raw signal.
* **S-D Distance:** The min and max allowable separation between a source and detector. 	
* **Noise:** The maximum value permitted for snr (signal to noise ratio).
* **hb Match:** Prune the oxy and deoxy channels if the corresponding oxy and deoxy channel is bad.

#### B: Data Processing:

* **Detrend:** Remove any trends in the raw data.
* **Normalize:** Normalize the raw data signal.
* **Downsample:** If selected, specify the target sample rate in Hz. The raw data is only	downsampled if the target frequency is lower than the raw data sampling frequency.

#### C: Motion Artifact Correction:

Decide whether to perform manual or automated motion correction. If automated, additionally specify whether to use a smoothing window and specify the length in seconds in the SD window field.

#### D: Data Visualization:

Select which figures you want for output.

#### E: Run Preprocessing:

To run only preprocessing, fill out the desired preprocessing configuration, select the subjects in the selection box and click Run Preprocessing.

## Postprocessing

Post-processing performs a channel pair-wise continuous wavelet analysis on the preprocessed data. The resulting wavelet is then operated on to obtain the following vectors.

#### Time Varying Coherence:

The wavelet coherence matrix used in this calculated is from taking the transform over the entire time series. To obtain a time varying vector we collapse the fz dimension across the fz band specified by the user.

#### Time Varying Phase:

For interregional comparisons, the phase vs time of the wavelet cross-spectrum matrix is generated by taking the phase angle in the interval [-PI, PI] for each element. This matrix is then collapsed across the fz band specified by the user.

#### Time Varying Power:

For each distinct region, a time varying power vector is calculated by taking a continuous 1-D wavelet transform on each the regions channels. The output is collapsed across the fz band specified by the user. Each channels vector within a region is then combined and a mean is taken to determine the regions time varying power.

#### Frequency Varying Power:

The signal for each channel is separated based on task. Wavelet coherence is the generated for the time series of each task individually. Collapsing across the time series we obtain a frequency varying coherence vector for each task. 

<div align="center">
<img src="https://github.com/carterrandall/fNIRS/blob/master/images/postprocessing.png" height = 400>
</div>

#### Channel Selection: 
After selecting regions for processing, the channels contained within that region will appear in each region box. You can edit the names given to each region for use in figures and file names. *See below for instructions on how to select channels*

#### FZ Band of Interest:

Min and max fz of the signal to be considered for calculations.

#### Hemoglobin Parameters:

The coherence analysis will be performed for each hemoglobin parameter selected.

### Run Postprocessing:

To run only postprocessing, fill out the desired postprocessing configuration, select the subjects in the selection box and click Run Postprocessing. *Only subjects who have been initialized and preprocessed will be available for postprocessing.*

## Channel Selection

The channel selection interface allows you to select to regions for comparison and how to compare them. After launching the Channels Selection Interface from the postprocessing pane select your channel configuration by:

**1.** Click the SD File button and select your probe (.sd) file.

**2.** Draw regions by clicking either of the blue of purple buttons. Regions must be closed polygons and can be edited once drawn. Channels that appear green are considered to be within the region and will do so if both endpoints (a source and detector) are located within the polygon.

**3.** The region overview displays a list of the channel numbers selected in each regions.

**4.** The region comparison section allows you to choose one or more comparison options when determining coherence: Intra-hemispheric for each region and/or Inter-hemispheric.

<div align="center">
<img src="https://github.com/carterrandall/fNIRS/blob/master/images/channelselection.png" height = 400>
</div>

## Analysis

This phase is to visualize the results from postprocessing and output data in tables. 

<div align="center">
<img src="https://github.com/carterrandall/fNIRS/blob/master/images/Analysis.png" height = 400>
</div>


#### Figure Selection:

Determine which figures you want to be output from the pipeline. View the appendix for examples of each. *Generating figures for large studies can greatly slow run time.*

#### Run Analysis: 

To run only analysis, fill out the desired analysis configuration, select the subjects in the selection box and click Run Analysis. *Only subjects who have been post-processed will be available for analysis.*

<interface>

## Appendix:

### Creating a probe (.sd) file:

1. In Matlab under the “Home” pane select the “Import Data" button.
2. In the file browser enable all file types to be availible for selection. 
3. Open one of your raw data files. (.nirs or .nirx)
4. This will bring up a window title “Import Wizard”. Click “Finish”
5. In the command line type : save(‘myfilename.sd’, ‘SD’). Important: Keep the quotes, only replace the myfilename with your desired filename. Press Enter. The probe file will now be saved in your current directory.

### Preprocessing Sample Figures:

#### Corrected Signal: 
<div align="center">
<img src="https://github.com/carterrandall/fNIRS/blob/master/images/preprocesssingCleanUp.jpg" height = 400>
</div>

#### Signal Quality: 
<div align="center">
<img src="https://github.com/carterrandall/fNIRS/blob/master/images/subject3_day1_SignalQuality.jpg" height = 400>
</div>

#### Cardiac Check:
<div align="center">
<img src="https://github.com/carterrandall/fNIRS/blob/master/images/subject3_day1CardiacPlot.jpg" height = 400>
</div>

#### Quality Check:
<div align="center">
<img src="https://github.com/carterrandall/fNIRS/blob/master/images/subject3_day1QualityCheck.jpg" height = 400>
</div>

### Postprocessing Sample Figures:

#### Channel Wise Static Coherene Plot:

<div align="center">
<img src="https://github.com/carterrandall/fNIRS/blob/master/images/WCohMatrix_0p001_0p15_Hz_Region 1 Region 2 Inter.jpg" height = 400>
</div>

### Analysis Sample Figures:

#### Time Varying Wavelet Transform (Power)
<div align="center">
<img src="https://github.com/carterrandall/fNIRS/blob/master/images/R1R2_wt.png" height = 400>
</div>

#### Time Varying Coherence
<div align="center">
<img src="https://github.com/carterrandall/fNIRS/blob/master/images/R1R2_tvCoh.png" height = 400>
</div>

#### Time Varying Phase
<div align="center">
<img src="https://github.com/carterrandall/fNIRS/blob/master/images/R1R2_tvPhase.png" height = 400>
</div>


#### FZ Varying Coherence (By Task)
<div align="center">
<img src="https://github.com/carterrandall/fNIRS/blob/master/images/R1R2_fzByTask.png" height = 400>
</div>
