# -*- coding: utf-8 -*-
import cascade
import matplotlib.pyplot as plt

# test downloading catalog
ct = cascade.exoplanet_tools.parse_database('EXOPLANETS.ORG', update=True)
# test extracting data record for single system
dr = cascade.exoplanet_tools.extract_exoplanet_data(ct, 'HD 189733 b')

print(dr)

# test generation of ligthcurve model
# first generate standard ini file and initialize cascade
cascade.initialize.generate_default_initialization()
path = cascade.initialize.default_initialization_path
cascade_param = cascade.initialize.configurator(path+"cascade_default.ini")

# define ligthcurve model
lc_model = cascade.exoplanet_tools.lightcuve()
print(lc_model.valid_models)
print(lc_model.par)

fig, axs = plt.subplots(1, 1, figsize=(12, 10))
axs.plot(lc_model.lc[0], lc_model.lc[1])
axs.set_ylabel(r'Normalized Signal')
axs.set_xlabel(r'Phase')
axes = plt.gca()
axes.set_xlim([0, 1])
axes.set_ylim([-1.1, 0.1])
plt.show()
