#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 16 10:17:51 2019

@author: bouwman
"""

import numpy as np
import pandas as pd
from causality.inference.search import IC
from causality.inference.independence_tests import RobustRegressionTest
from statsmodels.regression.linear_model import OLS
import networkx as nx
import matplotlib.pyplot as plt
# from networkx.drawing.nx_agraph import graphviz_layout
import seaborn as sns

sns.set_context("talk", font_scale=1.0, rc={"lines.linewidth": 2.5})
sns.set_style("white", {"xtick.bottom": True, "ytick.left": True})

size = 3000
source = np.random.normal(size=size)
noise = np.random.normal(size=size)
detector_y_ws = source + np.random.normal(size=size)
detector_y = detector_y_ws + noise
# detector_x1 = noise
detector_x1 = noise + np.random.normal(size=size)
detector_x2 = noise + np.random.normal(size=size)
detector_x3 = noise + np.random.normal(size=size)
detector_x4 = noise + np.random.normal(size=size)
detector_x5 = noise + np.random.normal(size=size)
detector_x6 = noise + np.random.normal(size=size)
detector_x7 = noise + np.random.normal(size=size)
detector_x8 = noise + np.random.normal(size=size)
detector_x9 = noise + np.random.normal(size=size)
detector_x10 = noise + np.random.normal(size=size)
detector_x11 = noise + np.random.normal(size=size)
detector_x12 = noise + np.random.normal(size=size)
detector_x13 = noise + np.random.normal(size=size)
detector_x14 = noise + np.random.normal(size=size)
detector_x15 = noise + np.random.normal(size=size)
detector_x16 = noise + np.random.normal(size=size)
detector_x17 = noise + np.random.normal(size=size)
detector_x18 = noise + np.random.normal(size=size)
detector_x19 = noise + np.random.normal(size=size)
detector_x20 = noise + np.random.normal(size=size)
detector_x21 = noise + np.random.normal(size=size)
detector_x22 = noise + np.random.normal(size=size)
detector_x23 = noise + np.random.normal(size=size)
detector_x24 = noise + np.random.normal(size=size)
detector_x25 = noise + np.random.normal(size=size)
detector_x26 = noise + np.random.normal(size=size)
detector_x27 = noise + np.random.normal(size=size)
detector_x28 = noise + np.random.normal(size=size)
detector_x29 = noise + np.random.normal(size=size)
detector_x30 = noise + np.random.normal(size=size)
detector_x31 = noise + np.random.normal(size=size)
detector_x32 = noise + np.random.normal(size=size)
detector_x33 = noise + np.random.normal(size=size)
detector_x34 = noise + np.random.normal(size=size)
detector_x35 = noise + np.random.normal(size=size)
detector_x36 = noise + np.random.normal(size=size)
detector_x37 = noise + np.random.normal(size=size)
detector_x38 = noise + np.random.normal(size=size)
detector_x39 = noise + np.random.normal(size=size)
detector_x40 = noise + np.random.normal(size=size)

X = pd.DataFrame({'source': source, 'noise': noise, 'det_y': detector_y,
                  'det_x1': detector_x1,
                  'det_x2': detector_x2,
                  'det_x3': detector_x3,
                  'det_x4': detector_x4,
                  'det_x5': detector_x5,
                  'det_x6': detector_x6,
                  'det_x7': detector_x7,
                  'det_x8': detector_x8,
                  'det_x9': detector_x9,
                  'det_x10': detector_x10,
                  'det_x11': detector_x11,
                  'det_x12': detector_x12,
                  'det_x13': detector_x13,
                  'det_x14': detector_x14,
                  'det_x15': detector_x15,
                  'det_x16': detector_x16,
                  'det_x17': detector_x17,
                  'det_x18': detector_x18,
                  'det_x19': detector_x19,
                  'det_x20': detector_x20,
                  'det_x21': detector_x21,
                  'det_x22': detector_x22,
                  'det_x23': detector_x23,
                  'det_x24': detector_x24,
                  'det_x25': detector_x25,
                  'det_x26': detector_x26,
                  'det_x27': detector_x27,
                  'det_x28': detector_x28,
                  'det_x29': detector_x29,
                  'det_x30': detector_x30,
                  'det_x31': detector_x31,
                  'det_x32': detector_x32,
                  'det_x33': detector_x33,
                  'det_x34': detector_x34,
                  'det_x35': detector_x35,
                  'det_x36': detector_x36,
                  'det_x37': detector_x37,
                  'det_x38': detector_x38,
                  'det_x39': detector_x39,
                  'det_x40': detector_x40,
                  'det_y_ws': detector_y_ws})

print(X.head())

plt.rcParams["patch.force_edgecolor"] = True
pd.plotting.scatter_matrix(X[['source', 'noise', 'det_y',
                              'det_x1']],
                           alpha=0.9,
                           diagonal='hist', grid=True)
plt.show()

# list_use = ['det_x1', 'det_x2', 'det_x3', 'det_x4', 'det_x5']
# list_use = ['det_x1', 'det_x2', 'det_x3', 'det_x4', 'det_x5',
#            'det_x6', 'det_x7', 'det_x8', 'det_x9', 'det_x10']
# list_use = ['det_x1', 'det_x2', 'det_x3', 'det_x4', 'det_x5',
#            'det_x6', 'det_x7', 'det_x8', 'det_x9', 'det_x10',
#            'det_x11', 'det_x12', 'det_x13', 'det_x14', 'det_x15']
list_use = ['det_x1', 'det_x2', 'det_x3', 'det_x4', 'det_x5',
            'det_x6', 'det_x7', 'det_x8', 'det_x9', 'det_x10',
            'det_x11', 'det_x12', 'det_x13', 'det_x14', 'det_x15',
            'det_x16', 'det_x17', 'det_x18', 'det_x19', 'det_x20',
            'det_x21', 'det_x22', 'det_x23', 'det_x24', 'det_x25',
            'det_x26', 'det_x27', 'det_x28', 'det_x29', 'det_x30',
            'det_x31', 'det_x32', 'det_x33', 'det_x34', 'det_x35',
            'det_x36', 'det_x37', 'det_x38', 'det_x39', 'det_x40']
# list_use = ['det_x1', 'det_x2']
# list_use = ['det_x1']

model = OLS(X['det_y'], X[list_use])
result = model.fit()
print(result.summary())

p = result.params.values
plt.plot(X['det_y'].values, color='blue')
# plt.plot(p*X['det_x1', 'det_x2', 'det_x3'].values)
plt.plot(np.dot(p, X[list_use].values.T), color='red')
plt.ylim([-6, 6])
plt.show()

# plt.plot(X['det_y'].values - p*X['det_x1', 'det_x2', 'det_x3'].values)
plt.plot(X['det_y'].values - np.dot(p, X[list_use].values.T), color='blue')
# plt.plot(X['det_y_ws'].values)
plt.ylim([-6, 6])
plt.show()

#plt.plot((X['det_y'].values - p*X['det_x1', 'det_x2', 'det_x3'].values) -
#         X['det_y_ws'].values)
plt.plot((X['det_y'].values - np.dot(p, X[list_use].values.T)) -
         X['det_y_ws'].values, color='blue')
plt.ylim([-6, 6])
plt.show()

print(np.std((X['det_y'].values - np.dot(p, X[list_use].values.T)) -
             X['det_y_ws'].values))

variable_types = {'noise': 'c', 'det_y': 'c', 'source': 'c',
                  'det_x1': 'c'}

# run the search
ic_algorithm = IC(RobustRegressionTest)
graph = ic_algorithm.search(X, variable_types)

for i in graph.edges(data=True):
    print(i)

pos = nx.spring_layout(graph)
nx.draw_networkx_nodes(graph, pos, cmap=plt.get_cmap('jet'), node_size=500)
nx.draw_networkx_labels(graph, pos)
nx.draw_networkx_edges(graph, pos, arrowstyle='->', cmap=plt.get_cmap('jet'),
                       arrowsize=10, node_size=500)
plt.title("Graph Title")
plt.axis('off')
plt.show()