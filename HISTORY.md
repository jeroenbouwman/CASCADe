## History of versions

**version 1.0.6**
- First online release

**version 1.0.7**
- Fixed bug in limb darkening correction
- Added disance to object configuration file
- Added derivation of stellar spectrum

**version 1.0.8**
- Improved parallel calcuation of regression modeling
- Added regularization to correction of averige TD subtraction
- Fixed a bug in the reading in of spectral fits files
- Added additional verbose plot of the calibrated stellar spectrum
- Updated processing_exceptions.ini

**version 1.0.9**
- Added sensitivity curves for Spitzer IRS data
- Added a fix for the low resolution Atlas 2000 models at longer wavelengths
- Updated Verbose plots for stellar spectrum

**version 1.1.0**
- Added asymmetric frame drop for spatial scans
- Updated processing_exceptions.ini
- Added passbands and wavelength definitions for the long wavelength modules
  of the Spitzer speectrograph

**version 1.1.1**
- Fixed the LL1 passband file for Spitzer IRS
- Added aditional output file containing systematics, residual and lc model
- Added aditional functionality to spectrally rebin timeseries data and to save data files.

**version 1.1.1**
- Added input stellar model as aditional output file to save_results step.
- Fixed issue with SL2 bandpass file.
- Added a higher spectral resolution wavelength bins file for Spizer IRS.
- Added the posibility of using higher orders of the  sources position as
  an aditional regressor .

**version 1.1.2**
- Changed the directory structure for the examples coming with the distribution.
- Added the option to get data files directly from git repository.

**version 1.1.3**
- Cleaned up the initilization module.
- Fixed bug in unpacking of zip file from git.
- Resolved problem with versioning of numpy, numba and batman in conda environment
  Currently only tested for python 3.8.

**version 1.1.4**
- Made a minor bug fix to the initialize module to make sure the warnings flag was
  not ignored.
- Small fixed to the example scrips and .ini files.
- Small bug fix to the Generic instrument module for the time unit.
- Bug fix in the extract_1d pipeline step when calling the
  correct_initial_wavelength_shift function.

**version 1.1.5**
 - Updated the setup.py script.
 - Updated the README file to reflect all changes
 - Fixed some spelling mistakes in the plot captions and variable names.

**version 1.1.6**
 - Fix bug in HST module for non constant integration times / samples op the ramp.
 - Updated the HST database file.

**version 1.1.7**
 - Added JWST functionality
 - Updated requirements

**version 1.1.8**
 - Solved issue #100, added normalization=None to phase cross-correlation.
 - Updated requirements for scikit-image to >= 0.19

**version 1.1.9**
 - Updated the requirements for numpy and nummba.
 - Updated the JWST functionality for MIRI/LRS

**version 1.1.10**
 - Made the number of data servers configurable

**version 1.1.11**
 - Update to JWST MIRI LRS module
 - Addedconfigurable regularization parameters for TD correction

**version 1.1.12**
 - Added the posibility to read in position file for JWST observations.
 - Added a parameter to regulate the number of data chuncks to be loaded at
   the same time

**version 1.1.15**
 - Updated database file for HST observations

**version 1.1.16**
 - Changed orbital phase calculation for MIRI LRS
 - Added option to remove integrations at the end of the time series
 - Added FWHM as addition regressor

**version 1.1.17**
 - Added lower resolution rebin option for combining observations.
 - Updated build_archive module to allow for skipping of archive incase of no connection.
 - Fixed depreciation warnings for methods used in HST and exoplanet_tools modules.
 - Removed regression matrix from data needed by post_proceccing step in the cpm_model module.

**version 1.2.0**
 - Major update to the cpm_model module, with the processing step now completely parallel.
 - Updated the responce of JWST/LRS for the exotethys package.

**version 1.2.1**
 - Fixed a memory issue with the parallel proseeeing of the fit results.

**version 1.2.2**
 - Fixed further memory issues with the parallel fitting and parameter storage.