#!/usr/bin/env python
# coding: utf-8

# # Importations and set-up checking


# Check package versions
import gammapy
import numpy as np
import astropy
import regions
import math

print("gammapy:", gammapy.__version__)
print("numpy:", np.__version__)
print("astropy", astropy.__version__)
print("regions", regions.__version__)


import matplotlib.pyplot as plt
from matplotlib.pyplot import gca

import astropy.units as u
from astropy.coordinates import SkyCoord, Angle
from regions import CircleSkyRegion

from gammapy.datasets import SpectrumDatasetOnOff, SpectrumDataset, Datasets, FluxPointsDataset
from gammapy.makers import SpectrumDatasetMaker
from gammapy.modeling import Fit
from gammapy.modeling.models import (
    PowerLawSpectralModel,
    SkyModel,
    PointSpatialModel,
    EBLAbsorptionNormSpectralModel
)
from gammapy.astro.darkmatter import DarkMatterAnnihilationSpectralModel
from gammapy.irf import load_cta_irfs
from gammapy.data import Observation
from gammapy.maps import MapAxis
from gammapy.estimators import FluxPointsEstimator


# # Simulate the observation
print("Simulating the observation: spectral+spatial")

# DEFINE OBSERVATION PARAMETERS AND OBSERVATION STRATEGY
livetime = 300 * u.h

pointing = SkyCoord(150.57, -13.26, unit="deg", frame="galactic")
offset = 1.5 * u.deg
on_region_radius = Angle("0.5 deg")

on_region = CircleSkyRegion(center=pointing, radius=on_region_radius)

# Energy axis in TeV
emin = 50/1000
emax = 50

energy_axis = MapAxis.from_energy_bounds(
    emin, emax, 10, unit="TeV", name="energy"
)


# Define spectral model, in our case DM
JFAC = 3.03e18 * u.Unit("GeV2 cm-5") # Perseus c-m moliné SRD VL-II
mDM = 10000*u.Unit("GeV")
channel = "b"
redshift = 0.017284
spectral_model = DarkMatterAnnihilationSpectralModel(
    mass=mDM, 
    channel=channel, 
    jfactor=JFAC, 
    z=redshift
)
absorption = EBLAbsorptionNormSpectralModel.read_builtin("dominguez", redshift=redshift)

model_simu = spectral_model * absorption

# Plot the simulated model
fig_1 = plt.figure()
model_simu.plot([(emin*1000)/mDM.value, (emax*1000)/mDM.value], energy_power=1)
plt.title("Simulated model")
form = plt.FormatStrFormatter('$%g$')
gca().xaxis.set_major_formatter(form)
plt.close(fig_1)
fig_1.savefig('X/original_spectra_t.png', quality=95, dpi=1000)



# Set the sky model used for the observation
model = SkyModel(spectral_model=model_simu, name="perseus")
print("This is the simulated model")
print(model)


# Load the IRFs
irfs = load_cta_irfs(
    "$GAMMAPY_DATA/prod3b-v2/bcf/North_z20_50h/irf_file.fits"
)

# Create the observation
obs = Observation.create(pointing=pointing, livetime=livetime, irfs=irfs)
print("Characteristics of the simulated observation")
print(obs)


# Make the SpectrumDataset
# NOTE: Even we don't set different energy ranges for recovered and true, if edisp is not considered then the 
# FluxPointEstimator breaks
dataset_empty = SpectrumDataset.create(
    e_reco=energy_axis, region=on_region, name="obs-0"
)
maker = SpectrumDatasetMaker(selection=["exposure", "edisp", "background"])

dataset = maker.run(dataset_empty, obs)


# Set the model on the dataset, and fake the first one to create the rest from here
dataset.models = model
dataset.fake(random_state=42)


# # Create the On/Off simulations

# Set off regions
dataset_on_off = SpectrumDatasetOnOff.from_spectrum_dataset(
    dataset=dataset, acceptance=1, acceptance_off=3
)

n_obs = 50

print("Creating the", n_obs, "On/Off simulations")

datasets = Datasets()

for idx in range(n_obs):
    dataset_on_off.fake(
        random_state=idx, npred_background=dataset.npred_background()
    )
    dataset_fake = dataset_on_off.copy(name=f"obs-{idx}")
    dataset_fake.meta_table["OBS_ID"] = [idx]
    datasets.append(dataset_fake)


table = datasets.info_table()
table.write('X/observations.txt', format='ascii')



# Check counts in one realization
fig_2 = plt.figure(1)
datasets[0].npred().plot_hist(label='Predicted S+B')
datasets[0].npred_signal().plot_hist(label='Predicted S')
datasets[0].npred_background().plot_hist(label='Predicted B')
plt.legend()
form = plt.FormatStrFormatter('$%g$')
gca().xaxis.set_major_formatter(form)
plt.close(fig_2)
fig_2.savefig('X/obs_counts.png', quality=95, dpi=1000)



## Check consistency in the sample of observaitions


mean_counts = table["counts"].mean()
mean_error = table["counts"].std()
mean_counts_off = table["counts_off"].mean()
mean_off_error = table["counts_off"].std()


fig_3, axes = plt.subplots(1, 2, figsize=(12, 4))
axes[0].hist(table["counts"], label=f"{mean_counts} +- {mean_error}")
axes[0].set_xlabel("Counts")
axes[0].axvline(x=table["counts"].mean(), color="red")
axes[0].axvspan(table["counts"].mean()-table["counts"].std(), table["counts"].mean()+table["counts"].std(), facecolor='r', alpha=0.2)
axes[0].legend()
axes[1].hist(table["counts_off"], label=f"{mean_counts_off} +- {mean_off_error}")
axes[1].set_xlabel("Counts Off")
axes[1].axvline(x=table["counts_off"].mean(), color="red")
axes[1].axvspan(table["counts_off"].mean()-table["counts_off"].std(), table["counts_off"].mean()+table["counts_off"].std(), facecolor='r', alpha=0.2)
axes[1].legend()

form = plt.FormatStrFormatter('$%g$')
gca().xaxis.set_major_formatter(form)
plt.close(fig_3)
fig_3.savefig('X/distr_counts.png', quality=95, dpi=1000)


# # Flux point estimator

# DEFINE LIST OF CHANNELS AND LIST OF MASSES
channels = ["b", "tau", "W"]
masses = [100, 250, 500, 750, 1000, 2500, 5000, 7500, 10000] #mas puntos, tipo log-scale 5 ptos por orden de magnitud 100-1000, 1000-10000

results = dict(mean={}, runs={})

print("Entering in the loops, this will take some minutes, so be patient!")

for ch in channels:
    results["runs"][ch] = {}
    results["mean"][ch] = {}

    for m in masses:

        # Fix to which DM model we want to fit the data
        flux_model_fit = DarkMatterAnnihilationSpectralModel(
            mass=m * u.Unit("GeV"),
            channel=ch,
            jfactor=JFAC,
            z=redshift
        )
        flux_fit = flux_model_fit * absorption
        model_fit = SkyModel(spectral_model=flux_fit, name="model-fit")

        # Set the energy bins to use in the flux point estimator
        energy_edges = np.array([emin, emax]) * u.TeV

        # Instantiate the estimator
        fpe = FluxPointsEstimator(energy_edges=energy_edges)

        upper_container = []
        results["runs"][ch][m] = {}
        results["mean"][ch][m] = {}

        for i in range(n_obs):
            # Run the estimator
            datasets[i].models = model_fit
            flux_points = fpe.run(datasets=datasets[i])

            # Clean the possible NaNs in the ULs
            for j in range(len(flux_points.table['norm_ul'])):
                x = math.isnan(flux_points.table['norm_ul'][j])
                if x == True:
                    flux_points.table['norm_ul'][j] = 0

            # Save the data
            upper_container.append(np.sum(flux_points.table['norm_ul']))
            results["runs"][ch][m][i] = flux_points.table_formatted

        results["mean"][ch][m][0] = np.mean(upper_container) * 3e-26
        results["mean"][ch][m][1] = np.std(upper_container) * 3e-26

sigmav = dict(ul={}, one_sigma={})
for ch in channels:
    sigmav["ul"][ch] = []
    sigmav["one_sigma"][ch] = []
    for m in masses:
        sigmav["ul"][ch].append(results["mean"][ch][m][0])
        sigmav["one_sigma"][ch].append(results["mean"][ch][m][1])



for ch in channels:
    sigmav["ul"][ch] = np.asarray(sigmav["ul"][ch])
    sigmav["one_sigma"][ch] = np.asarray(sigmav["one_sigma"][ch])
masses = np.asarray(masses)



fig_4 = plt.figure(figsize=(9,7))

for ch in channels:
   plt.plot(masses*1e-3, sigmav["ul"][ch], label='{}'.format(ch))
   plt.fill_between(masses*1e-3, sigmav["ul"][ch] - sigmav["one_sigma"][ch], sigmav["ul"][ch] + sigmav["one_sigma"][ch], alpha=0.2)

plt.yscale('log')
plt.xscale('log')
plt.xlabel('mass [TeV]')
plt.ylabel(r'$<\sigma v$> [cm$^3$s$^{-1}$]')
plt.legend()
plt.title(r'$m_{DM}$= 10 TeV, $b\bar{b}$, N$_{runs}$ = 50')

plt.hlines(3e-26, 0.1, 10, ls="--")

plt.show()

fig_4.savefig('X/sigmav_vs_mass_one_bin_10.png', quality=95, dpi=1000)





