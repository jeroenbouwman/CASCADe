# -*- coding: utf-8 -*-

import cascade
import numpy as np
from scipy.stats import norm as norm_stats

import pandas as pd

import matplotlib.pyplot as plt

from matplotlib.pylab import rcParams

# create linear system with RV
n = 1000

x1 = norm_stats.rvs(0, 1, size=n)
x2 = -x1 + norm_stats.rvs(0, 10**-3, size=n)
x3 = norm_stats.rvs(0, 1, size=n)

A = np.column_stack([x1, x2, x3])
b = 10 * x1 + 10 * x2 + 0.1 * x3

# Solve system

regularization_pararameters = {"lam0": 1.e-20, "lam1": 1.e4, "nlam": 150}

parameters, error_parameters, optimal_regularization_parameter = \
     cascade.cpm_model.solve_linear_equation(A, b, cv_method='gcv',
                                        reg_par=regularization_pararameters)

print(parameters, error_parameters, optimal_regularization_parameter)

parameters, error_parameters, optimal_regularization_parameter = \
     cascade.cpm_model.solve_linear_equation(A, b, cv_method='b100',
                                        reg_par=regularization_pararameters)

print(parameters, error_parameters, optimal_regularization_parameter)

parameters, error_parameters, optimal_regularization_parameter = \
     cascade.cpm_model.solve_linear_equation(A, b, cv_method='b95',
                                        reg_par=regularization_pararameters)

print(parameters, error_parameters, optimal_regularization_parameter)

A = np.array([[1, 0, -1], [0, 1, 0], [1, 0, 1], [1, 1, 0], [-1, 1, 0]])
coef = np.array([4, 2, 7])
b = np.dot(A, coef)
b = b + np.random.normal(0.0, 0.01, size=b.size)


parameters, error_parameters, optimal_regularization_parameter = \
     cascade.cpm_model.solve_linear_equation(A, b, cv_method='gcv',
                                        reg_par=regularization_pararameters)

print(parameters, error_parameters, optimal_regularization_parameter)

parameters, error_parameters, optimal_regularization_parameter = \
     cascade.cpm_model.solve_linear_equation(A, b, cv_method='b100',
                                        reg_par=regularization_pararameters)

print(parameters, error_parameters, optimal_regularization_parameter)

########################################
# FIT SIN(X) with polynomial regression
########################################
rcParams['figure.figsize'] = 12, 10

# Define input array with angles from 60deg to 300deg converted to radians
x = np.array([i*np.pi/180 for i in range(60, 300, 4)])
np.random.seed(10)   # Setting seed for reproducability
y = np.sin(x) + np.random.normal(0, 0.15, len(x))
data = pd.DataFrame(np.column_stack([x, y]), columns=['x', 'y'])
plt.plot(data['x'], data['y'], '.')

for i in range(2, 16):  # power of 1 is already there
    colname = 'x_%d' % i      # new var will be x_power
    data[colname] = data['x']**i
print(data.head())


def linear_regression(data, power, models_to_plot, cv_method='gcv',
                      reg_par={"lam0": 1.e-6, "lam1": 1.e2, "nlam": 60}):
    # initialize predictors:
    predictors = ['x']
    if power >= 2:
        predictors.extend(['x_%d' % i for i in range(2, power+1)])

    # Fit the model
    A = np.vstack((np.ones_like(data['x']), data[predictors].T)).T
    par, error_par, optimal_reg_par = \
        cascade.cpm_model.solve_linear_equation(A, data['y'],
                                                cv_method=cv_method,
                                                reg_par=reg_par)
    y_pred = np.dot(A, par)

    # Check if a plot is to be made for the entered power
    if power in models_to_plot:
        plt.subplot(models_to_plot[power])
        plt.tight_layout()
        plt.plot(data['x'], y_pred)
        plt.plot(data['x'], data['y'], '.')
        plt.title('Plot for power: %d' % power)

    # Return the result in pre-defined format
    rss = sum((y_pred-data['y'])**2)
    ret = [rss]
    ret.extend([par[0]])
    ret.extend(par[1:])
    return ret

# Initialize a dataframe to store the results:
col = ['rss', 'intercept'] + ['coef_x_%d' % i for i in range(1,16)]
ind = ['model_pow_%d' % i for i in range(1, 16)]
coef_matrix_simple = pd.DataFrame(index=ind, columns=col)

# Define the powers for which a plot is required:
models_to_plot = {1:231,3:232,6:233,9:234,12:235,15:236}

regularization_pararameters = {"lam0": 1.e-20, "lam1": 1.e5, "nlam": 200}

# Iterate through all powers and assimilate results
for i in range(1, 16):
    coef_matrix_simple.iloc[i-1,0:i+2] = linear_regression(data, power=i, models_to_plot=models_to_plot, reg_par=regularization_pararameters)

# Iterate through all powers and assimilate results
for i in range(1, 16):
    coef_matrix_simple.iloc[i-1,0:i+2] = linear_regression(data, power=i, models_to_plot=models_to_plot, cv_method='b100', reg_par=regularization_pararameters)

# Iterate through all powers and assimilate results
for i in range(1, 16):
    coef_matrix_simple.iloc[i-1,0:i+2] = linear_regression(data, power=i, models_to_plot=models_to_plot, cv_method='b95', reg_par=regularization_pararameters)



