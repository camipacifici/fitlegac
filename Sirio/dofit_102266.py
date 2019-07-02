from specfit import *
from spectrum import *
from sfhdust import *
from mcmc import *
import util, sys, os

# Use a piecewise SFH and a Calzetti law for the attenuation
csp = CSP(sfh=PiecewiseSFH(), dust=CalzettiDust())
age = 6.0e9
lookback_nodes = [3.0e9, 1.5e9, 0.70e9, 0.35e9, 0.0e9]
csp.param['timenodes'] = [ age - item for item in lookback_nodes ]
csp.param['age'] = age

# Use a gaussian broadening for the velocity dispersion
model = SpectrumModel(lsf=GaussianVelocityLSF(300.), csp=[csp])

# Priors on the model
priors = [{'z':     GaussianPrior(0.7788, 0.005),
           'sigma': UniformPrior(50, 500)},
          {'amplitude1':    LogGaussianPrior(1.0, 3.0),
           'amplitude2':    LogGaussianPrior(1.0, 3.0),
           'amplitude3':    LogGaussianPrior(1.0, 3.0),
           'amplitude4':    LogGaussianPrior(1.0, 3.0),
           'amplitude5':    LogGaussianPrior(1.0, 3.0),
#           'amplitude6':    LogGaussianPrior(1.0, 3.0),
           'A_V':           UniformPrior(0, 4),
           'metal':         GaussianPrior(0.02, 0.005)}]

# Various options
options = util.options()
options.evidence_tolerance = 0.50
options.nlive = 200
options.poly_order = 10              # let the spectrum shape be modulated
                                    # by a multiplicative polynomial
options.filterfile = '/data/jcosmos/data/catalogs/uvista/uvista_photometry/FILTER.RES.v7.R300'
options.grid_path = '/data/lib/stellarlib/bc03/models/Padova1994/chabrier/bc2003_hr_stelib_m[metal]_chab_ssp.ised'
options.grid_param = ['metal']
options.grid_paramlabels = [['52', '62', '72']]
options.grid_paramvalues = [[0.008, 0.02, 0.05]]
options.grid_paramlogint = [True]   # interpolate in log space over metallicity
options.air2vac = True       # Convert to vacuum wavelengths, since the data are in vacuum
options.output_prefix = 'fit/pysf_102266'
options.rescale_spectrum_errors = 1.3712 # multiply specerr.fits by this
options.add_to_lnprob = 5000.
options.age_lt_universe = False

# Set up the data
photdata = 'phot_102266.dat'
specdata = ['spec_102266.fits']
specerror = ['specerr_102266.fits']
specmaskfiles = ['specmask_102266.fits']
inst_lsf = None              # ignoring the instrumental resolution for now

specfit(specdata, specerror, photdata, priors, model, options,
        inst_lsf, specmaskfiles=specmaskfiles)
