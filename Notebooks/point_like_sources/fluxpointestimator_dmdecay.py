#!/usr/bin/env python
# coding: utf-8

# # Importations and set-up checking
# Check package versions
import gammapy
import numpy as np
import astropy
import regions
import math

print("Versions of the packages in usage")
print("gammapy:", gammapy.__version__)
print("numpy:", np.__version__)
print("astropy:", astropy.__version__)
print("regions:", regions.__version__)

import matplotlib
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
    EBLAbsorptionNormSpectralModel,
    SpectralModel,
    Parameter,
    PointSpatialModel
)
from gammapy.astro.darkmatter import DarkMatterAnnihilationSpectralModel
from gammapy.irf import load_cta_irfs
from gammapy.data import Observation
from gammapy.maps import MapAxis
from gammapy.estimators import FluxPointsEstimator
from gammapy.astro.darkmatter import PrimaryFlux


# # DEFINITION OF DECAY DM MODEL
class DarkMatterDecaySpectralModel(SpectralModel):

    r"""Dark matter annihilation spectral model.

    The gamma-ray flux is computed as follows:

    .. math::
        \frac{\mathrm d \phi}{\mathrm d E} =
        \frac{1}{4\pi \tau m_{\mathrm{DM}}}
        \frac{\mathrm d N}{\mathrm dE} \times D(\Delta\Omega)

    Parameters
    ----------
    mass : `~astropy.units.Quantity`
        Dark matter mass
    channel : str
        Annihilation channel for `~gammapy.astro.darkmatter.PrimaryFlux`
    scale : float
        Scale parameter for model fitting
    dfactor : `~astropy.units.Quantity`
        Integrated D-Factor needed when `~gammapy.modeling.models.PointSpatialModel` spatial model is used
    z: float
        Redshift value

    Examples
    --------
    This is how to instantiate a `DarkMatterDecaySpectralModel` model::

        channel = "b"
        massDM = 5000*u.Unit("GeV")
        dfactor = 3.41e19 * u.Unit("GeV cm-2")
        modelDM = DarkMatterDecaySpectralModel(mass=massDM, channel=channel, dfactor=dfactor)

    References
    ----------
    * `2011JCAP...03..051 <https://ui.adsabs.harvard.edu/abs/2011JCAP...03..051>`_
    """
    
    tag = "DarkMatterDecaySpectralModel"
    
    MEAN_LIFETIME = 1e+27 * u.Unit("s")
    """Mean lifetime of the particle"""

    scale = Parameter("scale", 1, min=0)

    def __init__(self, mass, channel, scale=scale.quantity, dfactor=1e19, z=0):
        self.z = z
        self.mass = mass
        self.channel = channel
        self.dfactor = dfactor
        self.primary_flux = PrimaryFlux(mass, channel=self.channel).table_model
        super().__init__(scale=scale)

    def evaluate(self, energy, scale):
        """Evaluate dark matter decay model."""
        flux = (
            scale
            * self.dfactor
            * self.primary_flux(energy=energy * (1 + self.z))
            / self.MEAN_LIFETIME
            / self.mass
            / (4 * np.pi)
        )
        return flux


# # START OF THE ANALYSIS


# # Simulate the observation

print("Simulating the observation: spectral+spatial")

# Define simulation parameters parameters
livetime = 300 * u.h
l = 150.57
b = -13.26

pointing = SkyCoord(150.57, -13.26, unit="deg", frame="galactic")
offset = 1.0 * u.deg
on_region_radius = Angle("1.0 deg")
on_region = CircleSkyRegion(center=pointing, radius=on_region_radius)

# Energy axis in TeV
emin = 30/1000
emax = 100
energy_axis = MapAxis.from_energy_bounds(
    emin, emax, 10, unit="TeV", name="energy"
)

# Define spectral model 
DFAC = 1.59e19 * u.Unit("GeV cm-2") # Perseus
mDM = 10000*u.Unit("GeV")
channel = "b"
redshift = 0.017284
spectral_model = DarkMatterDecaySpectralModel(
    mass=mDM, 
    channel=channel, 
    dfactor=DFAC, 
    z=redshift
)

# If want to consider EBL, uncomment following to lines and in the following change spectral_model->model_simu
#absorption = EBLAbsorptionNormSpectralModel.read_builtin("dominguez", redshift=redshift)
#model_simu = spectral_model * absorption

# Plot spectral model
fig_1 = plt.figure()
spectral_model.plot(energy_range=[emin, emax]* u.TeV)
form = plt.FormatStrFormatter('$%g$')
gca().xaxis.set_major_formatter(form)
plt.close(fig_1)
fig_1.savefig('original_spectra_10tev_b.png', quality=95, dpi=1000)

#  Define spatial model
spatial_model = PointSpatialModel(lon_0=l*u.Unit("deg"), lat_0=b*u.Unit("deg"), frame="galactic")

# Set the sky model used in the dataset
model = SkyModel(spectral_model=model_simu, spatial_model=spatial_model, name="perseus")
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
dataset_on_off = SpectrumDatasetOnOff.from_spectrum_dataset(dataset=dataset, acceptance=1, acceptance_off=3)

# Set the number of observations we want to create
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

# Save the data from the simulations
table = datasets.info_table()
table.write('observations.txt', format='ascii')

# Check counts in one selected realization
fig_2 = plt.figure(1)
datasets[0].npred().plot_hist(label='Predicted S+B')
datasets[0].npred_signal().plot_hist(label='Predicted S')
datasets[0].npred_background().plot_hist(label='Predicted B')
plt.legend()
form = plt.FormatStrFormatter('$%g$')
gca().xaxis.set_major_formatter(form)
plt.close(fig_2)
fig_2.savefig('obs_counts_one_simu.png', quality=95, dpi=1000)


# # Check consistency in the sample of observations

mean_counts = table["counts"].mean()
mean_error = table["counts"].std()
mean_counts_off = table["counts_off"].mean()
mean_off_error = table["counts_off"].std()

# Consistency in On counts and Off counts
fig_3, axes = plt.subplots(1, 2, figsize=(12, 4))
axes[0].hist(table["counts"], label=f"{mean_counts} +- {mean_error}")
axes[0].set_xlabel("Counts")
axes[0].axvline(x=table["counts"].mean(), color="red")
axes[0].axvspan(table["counts"].mean()-table["counts"].std(),
                table["counts"].mean()+table["counts"].std(), facecolor='r', alpha=0.2)
axes[0].legend()
axes[1].hist(table["counts_off"], label=f"{mean_counts_off} +- {mean_off_error}")
axes[1].set_xlabel("Counts Off")
axes[1].axvline(x=table["counts_off"].mean(), color="red")
axes[1].axvspan(table["counts_off"].mean()-table["counts_off"].std(),
                table["counts_off"].mean()+table["counts_off"].std(), facecolor='r', alpha=0.2)
axes[1].legend()

form = plt.FormatStrFormatter('$%g$')
gca().xaxis.set_major_formatter(form)
plt.close(fig_3)
fig_3.savefig('distr_counts.png', quality=95, dpi=1000)


# # Flux point estimator 

# Set list of channels and masses we want to fit
channels = ["b", "tau", "W"]
masses = [100, 250, 500, 750, 1000, 2500, 5000, 7500, 10000, 25000, 50000, 75000, 100000]

print("Entering in the loops, this will take some minutes, so be patient!")

results = dict(mean={}, runs={})

for ch in channels:
    results["runs"][ch] = {}
    results["mean"][ch] = {}

    for m in masses:

        # Fix to which DM model we want to fit the data
        flux_model_fit = DarkMatterDecaySpectralModel(
            mass=m * u.Unit("GeV"),
            channel=ch,
            dfactor=DFAC,
            z=redshift
        )
        # Same here for the EBL
        #flux_fit = flux_model_fit * absorption
        model_fit = SkyModel(spectral_model=flux_model_fit, spatial_model=spatial_model, name="model-fit")

        # Set the energy bins to use in the flux point estimator
        # IMPORTANT: Don't try to fit further than the m, since the spectra drops abruptly
        energy_edges = np.array([emin, m/1000]) * u.TeV

        # Instantiate the estimator, here you can specify the range of values to search for norm, we restrict norm > 0
        fpe = FluxPointsEstimator(energy_edges=energy_edges, norm_min=1e-6, norm_max=1e3, norm_n_values=50)

        upper_container = []
        results["runs"][ch][m] = {}
        results["mean"][ch][m] = {}

        print("Computing", ch, "channel and", m, "GeV mass")

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

        results["mean"][ch][m][0] = 1e27/np.mean(upper_container)
        results["mean"][ch][m][1] = np.sqrt(
            ((1e27) ** 2) * (np.std(upper_container) ** 2) / (np.mean(upper_container) ** 4))

# Gammapy has some troubles writing directly the table of the runs
# if you want to save the data as txt file uncomment the next few lines of code
# these results tables contain information about the likelihood and the fit, but not the complete likelihood profile
#res = np.zeros([n_obs, 23])
#for ch in channels:
#    for m in masses:
#        for i in range(n_obs):
#            for j in range(23):
#                if j==16:
#                    res[i, j] = -90
#                elif j==17:
#                    res[i, j] = results["runs"][ch][m][i][0][j][0]
#                else:
#                    res[i, j] = results["runs"][ch][m][i][0][j]
#        np.savetxt('results_{}_{}.txt'.format(ch, m), res, header='Columns corresponding to flux_points.table_formated in gammapy')

# Arrange the data from the fits so we can plot it easily later
livetime = dict(ul={}, one_sigma={})
for ch in channels:
    livetime["ul"][ch] = []
    livetime["one_sigma"][ch] = []
    for m in masses:
        livetime["ul"][ch].append(results["mean"][ch][m][0])
        livetime["one_sigma"][ch].append(results["mean"][ch][m][1])

for ch in channels:
    livetime["ul"][ch] = np.asarray(livetime["ul"][ch])
    livetime["one_sigma"][ch] = np.asarray(livetime["one_sigma"][ch])
masses = np.asarray(masses)


# # Plot the constraints

# This is to have the latex font
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
fig_4 = plt.figure(figsize=(9,7))

for ch in channels:
    plt.plot(masses*1e-3, sigmav["ul"][ch], label='{}'.format(ch))
    plt.fill_between(masses*1e-3, sigmav["ul"][ch] - sigmav["one_sigma"][ch], sigmav["ul"][ch] + sigmav["one_sigma"][ch], alpha=0.2)

plt.yscale('log')
plt.xscale('log')
plt.xlim(1e-1, 100)
plt.xlabel('Mass [TeV]', fontsize=15)
plt.ylabel(r'$\tau_{DM}$ [s]', fontsize=15)
plt.legend()

form = plt.FormatStrFormatter('$%g$')
gca().xaxis.set_major_formatter(form)
for tick in gca().xaxis.get_major_ticks():
    tick.label.set_fontsize(15)
for tick in gca().yaxis.get_major_ticks():
    tick.label.set_fontsize(15)

plt.show()
fig_4.savefig('tau_vs_mass.png', quality=95, dpi=1000)





