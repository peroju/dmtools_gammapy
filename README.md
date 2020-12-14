# DMTools - Gammapy pipeline
This is the github repo fo the CTA - DM tools package in Gammapy (python based):
https://docs.gammapy.org/0.18.2/index.html

## Installation
1 - First you need to install the last stable version of Gammapy in your computer. For this, you need to follow the steps in https://docs.gammapy.org/0.18.2/install/index.html
We recommend you to install Gammapy via conda. If you are not familiar with Gammapy, maybe it will be useful to you to download some of the tutorial notebooks and go through them. The developer version will in principle not be needed for working and testing this package.

2 - Second, you need to do a second installation of the last stable version of Gammapy in your computer, just follow again the same steps in the fisrt paragraph and when naming the conda environment, choose a different name so later you will be able to distinguish between the plain Gammapy enviroment and the one with the dm-module. 

3 - Third, our dm-module is written as a part of of the `astro` package in the `Dark Matter` sub-package. You need to find in your second Gammapy installation this package, for this, follow this path:
`/Users/your_home_user/opt/anaconda3/envs/name_second_installation/lib/python3.7/site-packages/gammapy/astro/darkmatter`

4 - Replace the `spectra.py` and `utils.py` in your folder for the ones that are in this repo in the folder Installation. The `spectra.py` contains the `DMAnnihilationModel` and `DMDatasetOnOff`, and `utils.py` contains `SigmaVEstimator`. 

5 - To perform the first tests and to check if your installation works, execute the jupyter notebook `DarkMatterUseCaseSigmaVEstimator.ipynb`.

