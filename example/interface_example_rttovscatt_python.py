
# Example program calling RTTOV from Python directly through wrapper
# Requires the file lib/rttov_wrapper_f2py.so to be in your $PYTHONPATH or the current directory

# Example profile data is contained in example_data.py


from rttov_wrapper_f2py import *
import sys
from example.example_data_rttovscatt import *


# =================================================================
# Specify the profile data

# This example demonstrates how to run a simulation for HIRS for two
# profiles with variable CO2 and a very simple cloud profile

# example_data.py contains p, T, q, co2 for a single profile
# It also contains surface variables and other data for two profiles


# Define number of profiles and number of levels
nprofiles = 2
nlevels = len(p_ex)
ngases = 6 

# See wrapper user guide for gas IDs
gas_id_q                = 1
gas_id_scatt_hydro_frac = 60
gas_id_scatt_clw        = 61
gas_id_scatt_ciw        = 62
gas_id_scatt_rain       = 63
gas_id_scatt_snow       = 64

# Specify gas list (includes RTTOV-SCATT cloud inputs)
gas_id = np.array([gas_id_q, gas_id_scatt_hydro_frac, gas_id_scatt_clw, gas_id_scatt_rain, gas_id_scatt_ciw, gas_id_scatt_snow], dtype=np.int32)


# Define arrays for pressure, temperature and gases/clouds/aerosols;
# specify Fortran ('F') order for array storage to be more efficient
p = np.empty((nlevels, nprofiles), order='F', dtype=np.float64)
t = np.empty((nlevels, nprofiles), order='F', dtype=np.float64)
gases = np.empty((nlevels, nprofiles, len(gas_id)), order='F', dtype=np.float64)
ph = np.empty((nlevels+1, nprofiles), order='F', dtype=np.float64)
usercfrac = np.empty((nprofiles,), order='F', dtype=np.float64)


# The non-vertical-profile data is specified in example_data_rttovscatt.py
# Note that most arrays in example_data_rttovscatt.py need to be transposed
# for use with the direct wrapper interface rather than pyrttov:
datetimes   = datetimes.transpose()
angles      = angles.transpose()
surftype    = surftype.transpose()
surfgeom    = surfgeom.transpose()
s2m         = s2m.transpose()
skin        = skin.transpose()
zeeman      = zeeman.transpose()


# Populate the pressure, temperature, q and cloud arrays: these are the same for both profiles
for i in range(nprofiles):
    p[:, i] = p_ex[:]
    t[:, i] = t_ex[:]
    gases[:, i, 0] = q_ex[:]              # index 0 in gas_id array above is water vapour
    gases[:, i, 1] = cc_ex[:]             # similarly for cloud inputs...
    gases[:, i, 2] = clw_ex[:]
    gases[:, i, 3] = rain_ex[:]
    gases[:, i, 4] = ciw_ex[:]
    gases[:, i, 5] = snow_ex[:]
    ph[:nlevels, i] = ph_ex[:]
    ph[nlevels, i] = s2m[0, i]            # Bottom pressure half-level set to 2m pressure
    usercfrac[i] = usercfrac_ex

# =================================================================



# =================================================================
# Load the instrument

# Specify RTTOV and wrapper options. In this case:
# - provide access to the full radiance structure after calling RTTOV
# - turn on the verbose wrapper option
# - specify interpolation mode
# NB the spaces in the string between option names and values are important!
opts_str = 'store_rad 1 ' \
           'verbose_wrapper 1 ' \
           'opts%interpolation%interp_mode 1 ' \
           'nprofs_per_call 2 '

# Specify instrument and channel list and add coefficient files to the options string
rtcoef_dir = '../rtcoef_rttov13/'

rtcoef_file = rtcoef_dir + 'rttov13pred54L/rtcoef_noaa_15_amsua.dat'
hydrotable_file = rtcoef_dir + 'hydrotable/hydrotable_noaa_amsua.dat'

nchannels = 15
channel_list = np.arange(1, nchannels+1, 1, dtype=np.int32)

opts_str += ' file_coef ' + rtcoef_file + \
            ' file_hydrotable ' + hydrotable_file

calc_zef = False          # Turn off radar simulations
multi_hydro_frac = False  # Single hydro_frac profile

# Call the wrapper subroutine to load the instrument and check we obtained a valid instrument ID
inst_id = rttov_load_inst(opts_str, np.array((0,), dtype=np.int32))
if inst_id < 1:
    print('Error loading instrument')
    sys.exit(1)
# =================================================================


# =================================================================
# Initialise emissivity atlas

emis_atlas_path = '../../emis_data/'
month = datetimes[1, 0]            # Month is taken from the profile date

# Call the wrapper subroutine to set up the IR emissivity atlas
# NB we specify inst_id here so the atlas is initialised for this specific instrument for faster access;
#    to initialise the atlas for use with multiple instruments pass 0 as the inst_id
#    (see wrapper user guide for more information)
#err = rttov_mw_emis_atlas_setup(emis_atlas_path, month, -1, inst_id, 0)
#if err != 0: print('Error initialising IR emissivity atlas: atlas will not be used')
# =================================================================


# =================================================================
# Declare arrays for other inputs and outputs

# Define array for input/output surface emissivity and BRDF
surfemis = np.empty((nchannels, nprofiles), order='F', dtype=np.float64)

# Define direct model outputs
bt = np.empty((nchannels, nprofiles), order='F', dtype=np.float64)
# =================================================================


# =================================================================
# Call RTTOV

# Initialise the surface emissivity and reflectance before every call to RTTOV:
# in this case we specify a negative number to use the IR atlas over land
# (because we initialised it above) and to use RTTOV's emissivity models over sea surfaces
surfemis[:,:] = -1.

# Call the wrapper subroutine to run RTTOV direct
err = rttov_scatt_call_direct(inst_id, channel_list, datetimes, angles, surfgeom, surftype, skin, s2m, \
                              zeeman, p, t, gas_units, gas_id, gases, ph, usercfrac, \
                              multi_hydro_frac, calc_zef, surfemis, bt)
if err != 0:
    print('Error running RTTOV-SCATT direct')
    sys.exit(1)
# =================================================================


# =================================================================
# Examine outputs

# Outputs available are:
# - surfemis array contains surface emissivities used by RTTOV
# - rad array contains RTTOV radiance%total array
# - bt array contains RTTOV radiance%bt array
# - it is also possible to access the whole radiance structure because we set the store_rad option above

print('Surface emissivity used by RTTOV')
print(surfemis[:,:].transpose())

print('Total cloudy BT')    # This example has no visible/near-IR channels so this array contains BTs only
print(bt.transpose())

# To obtain data from RTTOV output structures, declare an array and call the relevant wrapper subroutine.
# For example for the clear-sky BTs:
btclear = np.empty((nchannels, nprofiles), order='F', dtype=np.float64)
err = rttov_get_bt_clear(inst_id, btclear)
print('Clear-sky BT')
print(btclear.transpose())

# =================================================================


# =================================================================
# Deallocate memory for all instruments and atlases

err = rttov_drop_all()
if err != 0: print('Error deallocating wrapper')
# =================================================================


