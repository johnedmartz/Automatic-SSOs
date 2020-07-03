# Automatic-SSOs
Automatic data calibration of Solar System Objects using the ssos pipeline, Filabres and GLS.


## Objectives:
The main purpose of this pipeline is to process a set of raw images to find and characterize candidates to Solar System Objects (both known and unknown). This includes the image reduction (bias subtraction and flat-field correction), astrometric calibration using Gaia DR2 as a reference catalogue, and photometric calibration, also using Gaia DR2 photometry. The pipeline is implemented to automatically detect bad astrometry and photometry (as defined below). Finally, the rotation period of the asteroid can be obtained from its cleaned light curve.

The pipeline uses for the image preparation (reduction and astrometric calibration) the Filabres pipeline and for the SSOs identification, the ssos pipeline

## Pipeline:
These are the main steps that the pipeline will follow in a normal execution:
-	With Filabres, images from each night will be reduced with their corresponding masterflat and masterbias, considering if there was a meridian flip at some point (characteristic of LSSS survey). This will also write a WCS header on those images using Astrometry in a first place and SCAMP afterwards. For more details on the Filabres performance, we refer the user to the specific documentation.
-	If Astrometry can’t astrometrically solve solve an image, it won’t be used in the following steps.
-	For images with astrometric solution, Filabres saves a detection catalogue of extrated sources with their corresponding counterparts in the reference catalogue (Gaia DR2). For each catalogue, a boxplot with median separations between the coordinates of extracted and reference sources will be used to detect outliers among the calibrated images. These outliers would be images with median separations between extracted and reference sources above Q3 (third quartile) + 1.5 IQR (interquartile range). 
-	The user can decide whether automatically recalibrate these outliers in the configuration file (config.autossos).
-	The “good” images (astrometrically calibrated and within the accepted threshold) will be copied to an SSOS folder where the main software (ssos) to detect objects will be executed.
-	A check for the “FILTER” keyword in the header will be done. If this keyword doesn’t exist it will be set to “NONE”, as the ssos pipeline requires a “FILTER” keyword.
-	A file with statistics (mean, median, standard deviation and MAD) of the background of these images will be created. This will help to identify images with bad photometry using another boxplot. Outliers will be removed following the same criteria as before (Q3 + 1.5 IQR).
-	The ssos pipeline will be executed on the SSOS folder, containing all astrometrically calibrated images. The key in the identification of SSOs with this pipeline lies on the identification of a linear motion of an object over the sky. It performs a crossmatch with the SkyBoT (Sky Body Tracker) VO service to identify known SSOs and differentiate the unknown (new candidates) objects.
-	A visual inspection of candidates is needed. The ones marked as ASTEROID but without SkyBot name will be named AST_X and the ones marked as UFOs will be named as UNK_X.
-	After ssos execution, each full catalogue (this is, detection catalogues) will be matched with its corresponding reference sources from Gaia DR2 to carry out the photometric calibration.
-	A linear regression between instrumental and reference (Gmag) magnitudes with sigma clipping will be performed on these matched catalogues in order to perform the photometric calibration.
-	The parameters of these fits will be used to convert instrumental magnitudes to physical magnitudes the SSOs candidates.
-	Calibrated magnitudes given by a regression with R2 lower than certain threshold will be removed.
-	Finally, a file named validcalib_ssos_xxxxxx.csv will be created with the positions and calibrated magnitudes of the SSOs candidates, together with the SkyBoTs recovered information, when available.
-	If desired, the user can call a function (explained at the end of this document) to extract rotation periods using GLS software.

## Requirements: (Python 3.8)
###	Main software:
+ [Filabres](https://filabres.readthedocs.io/en/latest/index.html)
+ [SSOs](https://ssos.readthedocs.io/en/latest/)
+ [GLS](https://github.com/mzechmeister/GLS)

###	Based on:
+ SExtractor 2.25.0
+ SCAMP 2.7.8
+ SWarp 2.38.1

###	Aditional Python packages:
+ sbpy
+ astroML
+ photutils

## Installation:
First part is like the Filabres installation.

#### 1. Conda installation:
Visit the Miniconda webpage and download the installer corresponding to your operative system.
If you have updated the $PATH system variable during the miniconda or conda installation, you can call conda commands directly in the shell, like this:
$ conda info	
If not, you will need the add the path to the command, like:
$ /path/to/conda/bin/conda info
In this guide we will write the commands without the full path, for simplicity.

#### 2. Create a conda environment:
```bash
$ conda create --name autossos python=3 \
astropy \
ipython \
matplotlib \
numpy \
pandas \
python-dateutil \
PyYaml \
scipy \
setuptools
```
and answer "y".

#### 3. Activate environment:
```bash
$ conda activate autossos
```
which yields a different system prompt to the user:
```bash
(autossos) $
```
#### 4. Installing filabres:
```bash
(autossos) $ git clone https://github.com/nicocardiel/filabres.git
```

A folder named filabres will be created with the setup files.
```bash
(autossos) $ cd filabres
(autossos) $ python setup.py build
(autossos) $ python setup.py install
```

If you have filabres already installed in your system, but want to update the code with the latest version, you need to move to the same directory where you previously cloned the repository, pull the latest changes of the code, and reinstall it:
```bash
(autossos) $ cd filabres
(autossos) $ git pull
(autossos) $ python setup.py build
(autossos) $ python setup.py install
(autossos) $ cd ..
```
#### 5. Installing additional packages:
```bash
(autossos) $ conda install -c conda-forge pyvo
(autossos) $ conda install -c conda-forge astrometry
(autossos) $ conda install -c conda-forge astromatic-source-extractor
(autossos) $ conda install -c conda-forge astromatic-scamp
(autossos) $ conda install -c conda-forge astromatic-swarp 

(autossos) $ pip install sbpy
(autossos) $ pip install photutils
(autossos) $ pip install astroml
(autossos) $ pip install tqdm
(autossos) $ pip install statsmodels
```
Note: These packages can also be installed from the conda repository.

#### 6. Download ssos and Automatic-SSOs:
* ssos:
```bash
(autossos) $ git clone https://github.com/maxmahlke/ssos.git
```
A folder named ssos will be created with the setup files.

* Automatic-SSOs:
```bash
(autossos) $ git clone https://github.com/johnedmartz/Automatic-SSOs.git
```
A folder named Automatic-SSOs will be created with the setup files.

#### 7.	Installing SSOs:

* Replace utils.py inside ssos with with the one included in this pipeline.
```bash
(autossos) $ cp Automatic-SSOs/ssos/utils.py ssos/ssos/
```
* Proceed with SSOS installation.
```bash
(autossos) $ cd ssos
(autossos) $ python setup.py build
(autossos) $ python setup.py install
(autossos) $ cd ..
```
#### 8.	 Installing autossos:
```bash
(autossos) $ cd Automatic-SSOs 
(autossos) $ python setup.py build
(autossos) $ python setup.py install
(autossos) $ cd ..
```

## Step-by-step workflow:

#### 1.	Create a main directory and go into it:
```bash
(autossos) $ mkdir test_20190411
(autossos) $ cd test_20190411
```

#### 2.	Generate the default configuration files and folders in the main directory:
```bash
(autossos) $ autossos -d
```
Four configuration files and two folders will be created:
```bash
SSOS/
multiple_nights/
config.scamp
config.sex
config.autossos
default.param
```

-	Multiple_nights/ contains an example of the required structure:
```bash
Multiple_nights/night1/images
Multiple_nights/night2/images
...
```
masterbias and masterflat fits must be included in each nightx/ folder.
This folder and all included can be named differently. There should be only one folder with these characteristics to avoid conflicts between folders.
The pipeline can process all the images at the same time, but it’s recommended to run each night separately or, if they are consecutive nights, adjust the parameter CROSSMATCH_SKYBOT as explained in the SSOS webpage.


-	SSOS/ contains the default.ssos file and semp/ folder required by SSOs.
You can edit its configuration files. By default, the scamp execution used by SSOs will have the parameters MATCH and SOLVE_ASTROM disabled. Usually, Filabres returns good results and enabling these parameters is time consuming and doesn’t improve the astrometric solution.
 
-	default.param, config.scamp and config.sex are the configuration files used by Filabres.

-	config.autossos is the main configuration file used by this pipeline.

+ Filabres parameters
INSTRUMENT:		Name of the instrument used for the observations.
FLIPSTAT:			True/False if there is a meridian flip.
RECALIB:	True/False to attempt automatically an individual recalibration of outliers.
+ Sigma clipping parameters
MAX_MAG	Highest magnitude of the linear region.
MIN_MAG	Lowest magnitude of the linear region.
MAX_MAGERR	Magnitudes with higher error will be ignored.
SIGM	Sigma value used for sigma clipping.

Catalogues with R2 lower than this threshold will not be used.


#### 3. Execution:
In the terminal write
```bash
(autossos) $ autossos path/to/directory
```
Where the directory is the path to the main directory we have created before. Or, if in the terminal we are in this directory
```bash
(autossos) $ autossos .
```
#### 4.	User input through the pipeline execution:

-	The pipeline has a conflict with files with duplicated names, so if the pipeline finds any images with the same name inside the nights directories it will warn you and give you the option in terminal to exit the pipeline and change them or automatically append the corresponding night directory to their names.

-	If RECALIB is set to True in the configuration file (config.autossos), it will attempt an individual recalibration of the outliers. If there are still outliers it will ask you to choose between use the outliers in the next steps, attempt another recalibration (you can edit Filabres configuration files while the pipeline is paused), continue the pipeline ignoring the outliers or exit.

-	After these first steps, the pipeline will ask you if you want to enable/disable the SCAMP astrometric solution (by default is disabled) or continue with your current ssos parameters.

-	After ssos execution you will need to inspect the candidates.


## Output:
Filabres will return the calibrated images inside the science-imaging folder(for the structure of the output of Filabres, we refer the used to the proper documentation).
The list of images which cannot be solved or are outliers will be saved in files inside each night folder (list_ouliers.csv and list_notsolved.csv).
The main output will be inside SSOS folder:

./SSOS/Stats.dat	These are the values of the statistics of each image with a column identifying if they are outliers

./SSOS/checkplots/	This folder contains SSOs plots, a plot of the statistics, R2 evolution over the night and the light curve of each asteroid.

./SSOS/cats/match_[…].cat	These are the full catalogues matched with the reference sources.

./SSOS/cats/ssos_[…].csv	SSOs output.

./SSOS/cats/calib_[…].csv	SSOs output including a column with the calibrated magnitudes.

./SSOS/cats/validcalib_[…].csv	SSOs output with calibrated magnitudes after removing bad values (given by bad photometry or R2 less than threshold).

## Optional executions:
-	If the Filabres has already been executed you can run the pipeline normally, Filabres will skip automatically images already calibrated. But if you want to skip this step anyway you can use:
```bash
(autossos) $ autossos path/to/directory -skip_filabres
```
-	If you want to skip the recalibration:
```bash
(autossos) $ autossos path/to/directory -skip_recalibration
```
Using both arguments at the same times also works.

-	If you want to check again the calibrations, you can use:
```bash
(autossos) $ autossos -check
(autossos) $ autossos -–check_calibrations
```
-	In case you want to run the pipeline from the ssos execution until the end use this from the main directory:
```bash
(autossos) $ autossos -SSOS
```
-	If you want to generate a plot of a single linear regression with sigma clipping:
```bash
(autossos) $ autossos -plot_sc [number_full] [catalogue] [sigma] [min_mag] [max_mag] [max_magerr]
(autossos) $ autossos –-plot_sigmaclipping [number_full] [catalogue] [sigma] [min_mag] [max_mag] [max_magerr]
```
For example, using the catalogue 50 from the full_2.cat, with 2.5 sigma, lower magnitude 12, higher magnitude 17 and maximum error 0.01:
```bash
(autossos) $ autossos -plot_sc 2 50 2.5 12 17 0.01
```
-	If you want to get the output again with different values for sigma clipping:
```bash
(autossos) $ autossos -re_sc [sigma] [min_mag] [max_mag] [max_magerr] [R2]
(autossos) $ autossos --redo_sigmaclipping [sigma] [min_mag] [max_mag] [max_magerr] [R2]
```

-	To obtain periods you can use:
```bash
(autossos) $ autossos -period [file] [name] [iterations]
```
where file is the output file named “validcalib_[…].csv” and name is the name provided by SkyBoT or the one assigned in the column SKYBOT_NAME.

For example, to find periods of the asteroid Polyxo using the output “validcalib_[…].csv” and with 2 iterations:
```bash
(autossos) $ autossos -period ./SSOS/cats/validcalib_[…].csv Polyxo 2
```
This will create a file, Polyxo_period.txt in this case, in the directory where this command was executed. This file contains the best period of each iteration and its false alarm probability.
