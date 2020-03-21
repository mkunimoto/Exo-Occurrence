# Exo-Occurrence

Grid-based exoplanet occurrence rates described in Kunimoto et al. (2020). Currently, inputs are specific to the planet search and vetting pipeline described in [Kunimoto et al. (2020)](https://iopscience.iop.org/article/10.3847/1538-3881/ab6cf8), but a more general form that works for any pipeline/catalogue is in progress.

The implementation of approximate Bayesian computation (ABC) is done with [CosmoABC](https://github.com/COINtoolbox/CosmoABC). See that link for more information on the form of the input and function files required for CosmoABC to run.

If you use the following, please cite the following:

Kunimoto, M. & Matthews, J. 2020, Searching the Entirety of Kepler Data. II. Occurrence Rate Estimates for FGK Stars, submitted to AAS journals.

## Provided Files

`make_input_file.py`: creates an input file in the format required by CosmoABC. Takes arbitrary period and radius bin limits as inputs. For example, if you want to fit 5 bins simultaneously using period limits of P = [50, 100] days and radius limits of R = [1.0, 1.414, 2.0, 2.828, 4.0, 5.657] Re, you would run `python make_input_file.py --pbin 50 100 --rbin 1 1.414 2 2.828 4 5.657`.

`func_ExoOccurrence.py`: function file used as an input to CosmoABC. Within the function file, the exoplanet population simulator and distance function are defined. The provided example fits 5 bins simultaneously using the bin limits above.

`input_ExoOccurrence.txt`: input file used as an input to CosmoABC. Within the input file, the priors on the occurrence rates and ABC settings (e.g. size of each population) are defined. The provided example fits 5 bins simultaneously using the bin limits above.

`FGK_properties.dat`: stellar properties for all 96,280 FGK stars. Columns are: KIC ID, radius (Rsun), lower uncertainty on radius, upper uncertainty on r_adius, average of lower/upper uncertainties, mass (Msun), effective temperature (K), logg (cm/s^2), metallicity (dex), quadratic limb darkening coefficient 1, quadratic limb darkening coefficient 2, standard devation of flux measurements, total length of observations (days), duty cycle

`FGK_planets.dat`: planet properties for the planets corresponding to the stellar sample. Columns are: orbital period (days), radius (Re)

The algorithm can be run using `python run_ABC.py -i input_ExoOccurrence.txt -f func_ExoOccurrence.py`, with all of the above files in the same directory. `run_ABC.py` is provided by CosmoABC.

