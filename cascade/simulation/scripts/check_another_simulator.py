# IPython log file
#
from cascade.utilities.readwrite_spectra import read_spectra_dir
from cascade.utilities.readwrite_spectra import plot_spectrum_2d
from another_simulator import another_simulator
#
file_model='init_files/cascade_WASP43b_calibrate_Atlasstar.ini'
file_object='init_files/cascade_WASP43b_object.ini'
#
another_simulator(file_object, file_model, verbose=False, plot=True)
#
in_dir='data/WASP43b/SPECTRA/'
#
mx1, ferror, waves, times = read_spectra_dir(in_dir, pattern='*fake_atlas*.fits', verbose=False)
flux1 = mx1.data
#
mx2, ferror, waves, times = read_spectra_dir(in_dir, pattern='giuseppe*.fits', verbose=False)
flux2 = mx2.data
#
diff = flux1-flux2
print(diff.min(), diff.max())

plot_spectrum_2d(diff.value, waves, times, 'diff noisy flux electons', normalise=False)
