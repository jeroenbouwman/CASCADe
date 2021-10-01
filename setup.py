import setuptools
from setuptools import setup

config = {
    'name': 'CASCADe',
    'description': 'CASCADe : Calibration of trAnsit Spectroscopy using CAusal Data',
    'author': 'Jeroen Bouwman',
    'url': 'https://jbouwman.gitlab.io/CASCADe/',
    'download_url': 'https://gitlab.com/jbouwman/CASCADe',
    'author_email': 'bouwman@mpia.de',
    'version': '1.0.9',
    'python_requires': '>=3.7',
    'license': 'GNU General Public License v3 (GPLv3)',
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
                    'Intended Audience :: Science/Research',
                    'Topic :: Scientific/Engineering :: Astronomy',
                    ],
    'package_data': {"data":["exoplanet_data/wavelength_bins/*.txt",
                             "exotethys/passbands/*.pass",
                             "exotethys/wavelength_bins/*.txt",
                             "calibration/SPITZER/IRS/S18.18.0/*.fits",
                             "calibration/SPITZER/IRS/S18.18.0/*.tbl",
                             "calibration/HST/WFC3/*.fits",
                             "calibration/HST/WFC3/*.conf",
                             "archive_databases/HST/WFC3/*.pickle",
                             "archive_databases/HST/WFC3/*.ini"]},
    'scripts': [],
    'data_files': [('', ['README.md'])]
}

setup(**config)
