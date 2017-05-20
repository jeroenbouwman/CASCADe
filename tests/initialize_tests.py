# -*- coding: utf-8 -*-
import cascade

print(cascade.initialize.cascade_configuration)
print(cascade.initialize.cascade_configuration.__dict__)

# generate standard ini file
cascade.initialize.generate_default_initialization()
path = cascade.initialize.default_initialization_path
cascade_param = cascade.initialize.configurator(path+"cascade_default.ini")

# one can split the init files
cascade_param = cascade.initialize.configurator(path+"cascade_cpm.ini",
                                                path+"cascade_object.ini",
                                                path+"cascade_data.ini")

print(cascade.initialize.cascade_configuration)
print(cascade.initialize.cascade_configuration.__dict__)

cascade_param.reset()
print(cascade.initialize.cascade_configuration)
print(cascade.initialize.cascade_configuration.__dict__)

