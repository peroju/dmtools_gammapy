# DMTools - Gammapy pipeline
This is the github repo fo the CTA - DM tools package in Gammapy (python based):
https://docs.gammapy.org/0.18.2/index.html

Task list: https://docs.google.com/document/d/1d9CnTxpa384FKHirjk3OHzYNMG0no7EWq4QU_5nix1Y/edit

Main Task Force page: https://portal.cta-observatory.org/WG/PHYS/SitePages/DM%20Tools%20Task%20Force.aspx

If you want to be a collaborator of this repo and contribute, contact Judit Pérez-Romero (judit.perez@uam.es).

## Option 1: dm-module and SigmaVEstimator
Available notebooks:
- DarkMatterUseCaseSigmaVEstimator.ipynb (point-like case, DM annihilation, nuissance parameters)

### Installation
1 - First you need to install the last stable version of Gammapy in your computer. For this, you need to follow the steps in https://docs.gammapy.org/0.18.2/install/index.html
We recommend you to install Gammapy via conda. If you are not familiar with Gammapy, maybe it will be useful to you to download some of the tutorial notebooks and go through them. The developer version will in principle not be needed for working and testing this package.

2 - Second, you need to do a second installation of the last stable version of Gammapy in your computer, just follow again the same steps in the fisrt paragraph and when naming the conda environment, choose a different name so later you will be able to distinguish between the plain Gammapy enviroment and the one with the dm-module. 

3 - Third, our dm-module is written as a part of of the `astro` package in the `Dark Matter` sub-package. You need to find in your second Gammapy installation this package, for this, follow this path:
`<your_path_to_anacodna>/anaconda3/envs/name_second_installation/lib/python3.7/site-packages/gammapy/astro/darkmatter`

4 - Replace the `spectra.py` and `utils.py` in your folder for the ones that are in this repo in the folder `Installation`. The `spectra.py` contains the `DMAnnihilationModel` and `DMDatasetOnOff`, and `utils.py` contains `SigmaVEstimator`. 

5 - To perform the first tests and to check if your installation works, execute the jupyter notebook `DarkMatterUseCaseSigmaVEstimator.ipynb`.
**Note**: It should work for the latest version of Gammapy, if not, please contact us. 

## Option 2: Gammapy with some addendums
Available notebooks/scripts:
- fluxpointestimator_dmannihilation.py (point-like case, DM annihilation)
- fluxpointestimator_dmdecay.py (point-like case, DM decay)
- reding_likelihoodprofile_files.ipynb, tau_vs_mass_comps.ipynb, sigma_vs_mass_comps.ipynb (reading varius output files)
- Gammapy_DM_3DSim_Limits-ExtendedSource_notParallel.py (extended source, 3D analysis, DM annihilation)

We strongly encourage you to use the python scripts if you intend to run more than 10 simulations, and use the notebooks as a safe test place.
To run the python scripts:
- use the flag: -W ignore
- BEFORE RUNNING IT: the parameters of the observation, the original spectrum, the numer of realizations, the masses and the channels are instantiated in the script, so please take a look to it to customize everything you may need.

### Installation
Just install gammapy as described in step 1) of Option 1.

## Which option should I use?
The Option 1 (dm-module and the SigmaVEstimator) was the first DMTool based in Gammapy. It allows to look carefully to the likelihoods and also has implemented the J-factor uncertainty as a nuissance parameter. It should also work for the latest Gammapy version. However, it needs lot of testing still, since the estimator does not seem pretty stable until reaching O(1000), which usually is very computationally expensive. The notebook is also not up-to-date, so it doesn't take profit of the latest changes in Gammapy.

The Option 2 (Gammapy with some addendums) is the more recent path the DMTools has taken. It allows to perform a complete standard DM analysis all within Gammapy, just with some definitions (included in the code) that have to be instantiated every time you run it. It's completely stable, since it only make use of Gammapy internal functionalities, and way more fast. We will work to add the j-factor uncertainty as a nuissance parameter. 

So in conclusion, use the one you feel more confortable with, but the one that would keep in manteinance would be Option 2.


## Gammapy tutorials/notebooks/documentation to take into account
- https://docs.gammapy.org/0.18.2/tutorials/simulate_3d.html
- https://docs.gammapy.org/0.18.2/tutorials/spectrum_analysis.html
- https://docs.gammapy.org/0.18.2/tutorials/sed_fitting.html
- https://docs.gammapy.org/0.18.2/tutorials/spectrum_simulation.html
- https://docs.gammapy.org/0.18.2/tutorials/extended_source_spectral_analysis.html
- https://docs.gammapy.org/0.18.2/tutorials/exclusion_mask.html
- https://docs.gammapy.org/0.18.2/tutorials/maps.html
- https://docs.gammapy.org/0.18.2/tutorials/modeling.html
- https://docs.gammapy.org/0.18.2/makers/reflected.html
- https://docs.gammapy.org/0.18.2/estimators/index.html

