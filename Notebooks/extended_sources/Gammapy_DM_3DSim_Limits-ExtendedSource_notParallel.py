#!/usr/bin/env python

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=Warning)

from scipy.optimize import brentq
from scipy.interpolate import interp1d
import scipy.stats as st
import time,sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import multiprocessing as mp

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.units import Quantity
from astropy.io import fits

from gammapy.irf import load_cta_irfs
from gammapy.maps import WcsGeom, MapAxis, WcsNDMap
from gammapy.modeling.models import (
    PowerLawSpectralModel,
    TemplateSpatialModel,
    PointSpatialModel,
    SkyModel,
    Models,
    FoVBackgroundModel
    
)
from gammapy.makers import MapDatasetMaker, SafeMaskMaker
from gammapy.modeling import Fit
from gammapy.data import Observation
from gammapy.datasets import MapDataset
from gammapy.astro.darkmatter import DarkMatterAnnihilationSpectralModel

# ## Load IRFs

#filename = "$CALDB/data/cta/prod3b-v1/bcf/North_z20_average_50h/irf_file.fits"
filename = "$CALDB/data/cta/prod3b-v1/bcf/South_z20_average_50h/irf_file.fits"
irfs = load_cta_irfs(filename)

jfactor_filename = 'jfactor_maps/Einasto/annihil_RetII2D_FOVdiameter10.0deg_nside1024-JFACTOR-Jtot-image.fits'
hdul = fits.open(jfactor_filename)

# 
# # Define map geometry

GLON = hdul[0].header['PSI_0']   * u.Unit("deg")
GLAT = hdul[0].header['THETA_0'] * u.Unit("deg")
src_pos = SkyCoord(GLON, GLAT, frame="galactic")
emin = 0.03 
emax = 100
unit = "TeV"
lg_emin = np.log10(emin)
lg_emax = np.log10(emax)
ENERGY_BINS = 31

axis = MapAxis.from_edges(
    np.logspace(lg_emin, lg_emax, ENERGY_BINS),
    unit=unit,
    name="energy",
    interp="log",
)
geom = WcsGeom.create(
    skydir=src_pos, 
    binsz=0.02, 
    width=(2,2), 
    frame="galactic", 
    axes=[axis]
)


# ## Build model and create dataset

# **Declare constants and parameters for our DM model**

#JFAC = 2.00e19 * u.Unit("GeV2 cm-5") # <--- Reticulum II Point Source
JFAC = 3.27e19 * u.Unit("GeV2 cm-5") # <--- Reticulum II Extended
#JFAC = 1.26e20 * u.Unit("GeV2 cm-5") # <--- Draco I Extended

mDM = 40000*u.Unit("GeV")
channel = "b"
xsection = 1e-26 * u.Unit("cm3 s-1")
redshift = 0
RATIO=2.71


# **Define 3D Sky Model**
DarkMatterAnnihilationSpectralModel.THERMAL_RELIC_CROSS_SECTION = xsection
flux_model = DarkMatterAnnihilationSpectralModel(
    mass=mDM, 
    channel=channel, 
    jfactor=JFAC
)
    
spatial_model = TemplateSpatialModel.read(jfactor_filename)

sky_model = SkyModel(
    spatial_model=spatial_model, spectral_model=flux_model,
    name="model-simu"
)

bkg_model = FoVBackgroundModel(dataset_name="dataset-simu")

models = Models([sky_model, bkg_model])

# ## Declare observation values

pointing = src_pos
livetime = 100 * u.hour
offset = 2.0 * u.deg
#offset = 0.5 * u.deg

# Create an in-memory observation
obs = Observation.create(pointing=pointing, livetime=livetime, irfs=irfs)

# ## Start the simulations and get the limits

#masses = [70, 200, 500, 800, 1000, 5000, 8000, 10000, 30000, 50000, 60000, 100000]*u.GeV
masses = [200, 1000, 50000]*u.GeV
#channels = ["b", "mu", "tau", "Z"] 
channels = ["b"] 


columns = ['ch','mass', 'sigma_v']
table = pd.DataFrame(columns=columns)
    
# Make the MapDataset
empty = MapDataset.create(geom, name="dataset-simu")
maker = MapDatasetMaker(selection=["exposure", "background", "psf", "edisp"])
maker_safe_mask = SafeMaskMaker(methods=["offset-max"], offset_max=offset)

    


NumberOfRuns = 2
for run in range(NumberOfRuns):
    
    # Create simulated dataset
    dataset = maker.run(empty, obs)
    dataset = maker_safe_mask.run(dataset, obs)
    dataset.models = models
    dataset.fake(int(time.time()))

    for ch in channels:
        for mass in masses:
    
            flux_model_fit = DarkMatterAnnihilationSpectralModel(
                mass=mass, 
                channel=ch, 
                jfactor=JFAC
            )

            #spatial_model_fit = TemplateSpatialModel.read(jfactor_filename)

            model_fit = SkyModel(
                spatial_model=spatial_model, spectral_model=flux_model_fit, 
                name="model-fit"
            )

    
            bkg_fit = FoVBackgroundModel(dataset_name="dataset")
            models_fit = Models([model_fit, bkg_fit])
    
    
            # We do not want to fit the background in this case, so we will freeze the parameters
            models_fit["dataset-bkg"].spectral_model.norm.frozen = False
            models_fit["dataset-bkg"].spectral_model.tilt.frozen = False
            models_fit["dataset-bkg"].spectral_model.norm.value = 0.95

            models_fit.parameters["scale"].value=1
            models_fit.parameters["scale"].min=1e-6
            models_fit.parameters["scale"].max=1e10
            models_fit.parameters["norm"].frozen =True
    
            dataset.models = models_fit

            fit = Fit([dataset])
            result = fit.run()
        
            total_stat  = result.total_stat
            scale_min   = result.parameters["scale"].value
            scale_error = result.parameters["scale"].error  
        
            if scale_min<1e-4:
                scale_min=1e-4
        
            scale_max = 10*scale_min
            #print('0: ',run, mp_mass, scale_max)
            max_iter = 100
            counter = 0
            while (np.max(fit.stat_profile(parameter='scale',nvalues=1,
                                   bounds=(scale_max,scale_max))['stat_scan']-total_stat)<3 
                   and counter<max_iter):
        
                scale_max = scale_max * 5
                counter = counter + 1
    
            profile = fit.stat_profile(parameter='scale',nvalues=100, bounds=(scale_min,scale_max))
                    
            xvals = profile["scale_scan"]
            yvals = profile["stat_scan"] - total_stat - RATIO
    
            try:
                scale_found = brentq(
                        interp1d(xvals, yvals, kind="quadratic"),
                        scale_min,
                        scale_max,
                        maxiter=100,
                        rtol=1e-5,
                        )

                sigma_v = scale_found * xsection
                sigma_v = sigma_v.value
                
            except Exception as ex:
                sigma_v = None
                scale_min = None
                profile = None
                likemin = None
                #print(ex)

           
            del models_fit
            del model_fit
            del bkg_fit
            del profile

            print('1: ',run, ch, mass, sigma_v)
            table=table.append({'ch':ch, 'mass':mass.value, 'sigma_v':sigma_v},
                         ignore_index=True)
    
    del dataset
    

table[['mass']]    = table[['mass']].astype(int)
table[['sigma_v']] = table[['sigma_v']].astype(float)

#r = table.groupby(['ch','mass']).mean()
#r.reset_index(inplace=True)

fileout=sys.argv[1]
table.to_csv(fileout, index=False)
#table.to_csv('outputRaw_gammapy_RetII_ExtendedSource_T100_N1000.csv', index=False)
#table.to_csv('test.csv', index=False)

#r.to_csv('output_gammapy_RetII_ExtendedSource_T100_N1000.csv')






