# -*- coding: utf-8 -*-
import cascade


# test downloading catalog
ct = cascade.exoplanet_tools.parse_database('EXOPLANETS.ORG', update=True)
# test extracting data record for single system
dr = cascade.exoplanet_tools.extract_exoplanet_data(ct, 'HD 189733 b')


# test generation of ligthcurve model
# first generate standard ini file and initialize cascade
cascade.initialize.generate_default_initialization()
path = cascade.initialize.default_initialization_path
cascade_param = cascade.initialize.configurator(path+"cascade_default.ini")

# define ligthcurve model
lc_model = cascade.exoplanet_tools.lightcuve()
