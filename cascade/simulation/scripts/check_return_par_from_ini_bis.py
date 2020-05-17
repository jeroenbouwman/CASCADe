# check return_par_from_ini_bis

from compute_planet_flux_evolution import *

file_model='init_files/cascade_WASP43b_calibrate_Atlasstar.ini'
file_object='init_files/cascade_WASP43b_object.ini'
planet_radius, planet_semi_axis, albedo = return_par_from_ini_bis(file_object, file_model, verbose=False)
print(planet_radius, planet_semi_axis, albedo )
