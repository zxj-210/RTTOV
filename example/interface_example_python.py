
# Example program calling RTTOV from Python directly through wrapper
# Requires the file lib/rttov_wrapper_f2py.so to be in your $PYTHONPATH or the current directory

# Example profile data is contained in example_data.py


from rttov_wrapper_f2py import *
import sys
from example.example_data import *


# =================================================================
# Specify the profile data

# This example demonstrates how to run a simulation for HIRS for two
# profiles with variable CO2 and a very simple cloud profile

# example_data.py contains p, T, q, co2 for a single profile
# It also contains surface variables and other data for two profiles


# Define number of profiles and number of levels
nprofiles = 2
nlevels = len(p_ex)


# See wrapper user guide for gas IDs
gas_id_q = 1
gas_id_o3 = 2
gas_id_co2 = 3

gas_id_cfrac = 20
gas_id_lwc1  = 21
gas_id_iwc   = 30


# The gas ID array tells the wrapper which gases, aerosol and cloud profiles are being supplied:
# it must specify water vapour in all cases plus any other optional items;
# the order of the elements is not important, but it must be consistent with the data in the gases array
gas_id = np.array([gas_id_q, gas_id_co2, gas_id_cfrac, gas_id_lwc1, gas_id_iwc], dtype=np.int32)

# Define arrays for pressure, temperature and gases/clouds/aerosols;
# specify Fortran ('F') order for array storage to be more efficient
p = np.empty((nlevels, nprofiles), order='F', dtype=np.float64)
t = np.empty((nlevels, nprofiles), order='F', dtype=np.float64)
gases = np.empty((nlevels, nprofiles, len(gas_id)), order='F', dtype=np.float64)

# Populate the pressure, temperature, q and co2 arrays: these are the same for both profiles
for i in range(nprofiles):
    p[:, i] = p_ex[:]
    t[:, i] = t_ex[:]
    gases[:, i, 0] = q_ex[:]      # index 0 in gas_id array above is water vapour
    gases[:, i, 1] = co2_ex[:]    # index 1 in gas_id array above is co2

# Initialise cloud inputs to zero
gases[:, :, 2:] = 0.

# Specify some very simple cloud inputs (in layers 50 and 60) in profile 1; profile 2 contains no cloud
gases[49, 0, 2] = 0.5  # cfrac in layer 50, profile 1
gases[49, 0, 4] = 1.0  # IWC in layer 50, profile 1
gases[59, 0, 2] = 0.8  # cfrac in layer 60, profile 1
gases[59, 0, 3] = 1.0  # LWC in layer 60, profile 1

mmr_cldaer = 0

# The remaining profile data is specified in example_data.py
# Note that most arrays in example_data.py need to be transposed
# for use with the direct wrapper interface rather than pyrttov:
datetimes   = datetimes.transpose()
angles      = angles.transpose()
surftype    = surftype.transpose()
surfgeom    = surfgeom.transpose()
s2m         = s2m.transpose()
skin        = skin.transpose()
simplecloud = simplecloud.transpose()
clwscheme   = clwscheme.transpose()
icecloud    = icecloud.transpose()
zeeman      = zeeman.transpose()

# =================================================================



# =================================================================
# Load the instrument

# Specify RTTOV and wrapper options. In this case:
# - turn interpolation on
# - supply CO2 as a variable gas
# - turn cloudy IR simulations on
# - provide access to the full radiance structure after calling RTTOV
# - turn on the verbose wrapper option
# NB the spaces in the string between option names and values are important!
opts_str = 'opts%interpolation%addinterp 1 ' \
           'opts%rt_all%co2_data 1 '         \
           'opts%rt_ir%addclouds 1 '         \
           'store_rad 1 '                    \
           'verbose_wrapper 1 '

# Specify instrument and channel list and add coefficient files to the options string
rtcoef_dir = '../rtcoef_rttov13/'

rtcoef_file = rtcoef_dir + 'rttov13pred54L/rtcoef_metop_1_hirs_o3co2.dat'
sccldcoef_file = rtcoef_dir + 'cldaer_visir/sccldcoef_metop_1_hirs.dat'

nchannels = 19
channel_list = np.arange(1, nchannels+1, 1, dtype=np.int32)

opts_str += ' file_coef ' + rtcoef_file + \
            ' file_sccld ' + sccldcoef_file


# Call the wrapper subroutine to load the instrument and check we obtained a valid instrument ID
inst_id = rttov_load_inst(opts_str, channel_list)
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
atlas_wrap_id = rttov_load_ir_emis_atlas(emis_atlas_path, month, -1, inst_id, 0)
if atlas_wrap_id < 1: print('Error loading IR emissivity atlas: atlas will not be used')
# =================================================================


# =================================================================
# Declare arrays for other inputs and outputs

# Define array for input/output surface emissivity and BRDF
surfemisrefl = np.empty((nchannels, nprofiles, 4), order='F', dtype=np.float64)

# Define direct model outputs
btrefl  = np.empty((nchannels, nprofiles), order='F', dtype=np.float64)
rad     = np.empty((nchannels, nprofiles), order='F', dtype=np.float64)
# =================================================================


# =================================================================
# Call RTTOV

# Initialise the surface emissivity and reflectance before every call to RTTOV:
# in this case we specify a negative number to use the IR atlas over land
# (because we initialised it above) and to use RTTOV's emissivity models over sea surfaces
surfemisrefl[:,:,:] = -1.

# Use atlas
if atlas_wrap_id > 0:
    err = rttov_get_emisbrdf(atlas_wrap_id, surfgeom[0,:], surfgeom[1,:], surftype[0,:], surftype[1,:], \
                             angles[0,:], angles[1,:], angles[2,:], angles[3,:], skin[2,:], \
                             inst_id, channel_list, surfemisrefl[:,:,0])
    if err != 0:
        print('Error returning atlas emissivities: not using atlas')
        surfemisrefl[:,:,:] = -1.

# Call the wrapper subroutine to run RTTOV direct
err = rttov_call_direct(inst_id, channel_list, datetimes, angles, surfgeom, surftype, skin, s2m, \
                        simplecloud, clwscheme, icecloud, zeeman, p, t, gas_units, mmr_cldaer, \
                        gas_id, gases, surfemisrefl, btrefl, rad)
if err != 0:
    print('Error running RTTOV direct')
    sys.exit(1)
# =================================================================


# =================================================================
# Examine outputs

# Outputs available are:
# - surfemisrefl array contains surface emissivities (and reflectances) used by RTTOV
# - rad array contains RTTOV radiance%total array
# - btrefl array contains RTTOV radiance%bt and radiance%refl arrays (depending on channel wavelength)
# - it is also possible to access the whole radiance structure because we set the store_rad option above

print('Surface emissivity used by RTTOV')
print(surfemisrefl[:,:,0].transpose())

print('Total cloudy BT')    # This example has no visible/near-IR channels so this array contains BTs only
print(btrefl.transpose())

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


