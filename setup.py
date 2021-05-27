import setuptools
from setuptools import setup

config = {
    'name': 'CASCADe',
    'description': 'CASCADe : Calibration of trAnsit Spectroscopy using CAusal Data',
    'author': 'Jeroen Bouwman',
    'url': 'https://jbouwman.gitlab.io/CASCADe/',
    'download_url': 'https://gitlab.com/jbouwman/CASCADe',
    'author_email': 'bouwman@mpia.de',
    'version': '1.0.6',
    'python_requires': '>=3.7',
    'install_requires': ['batman-package', 'astropy', 'jplephem', 'scipy',
                         'numpy', 'configparser', 'photutils', 'pandas',
                         'scikit-learn', 'matplotlib', 'tqdm', 'seaborn',
                         'pytest', 'scikit-image', 'sphinx', 'alabaster',
                         'networkx', 'cython', 'astroquery', 'numba', 'ray',
                         'pyfiglet', 'termcolor', 'exotethys', 'statsmodels'],
    'packages': setuptools.find_packages(),
    'classifiers': ["Programming Language :: Python :: 3",
                    "License :: OSI Approved :: GNU GPLv3 License",
                    "Operating System :: OS Independent",
                    ],
    'scripts': []
}

setup(**config)
