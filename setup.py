import setuptools
import os
import pathlib

def package_files(directory):
    exclude = ['data', 'examples']
    paths = []
    for (path, directories, filenames) in os.walk(directory, topdown=True):
        directories[:] = [d for d in directories if d not in exclude]
        for filename in filenames:
            paths.append(os.path.join('..', path, filename))
    return paths

# The directory containing this file, MUST be RELATIVE
HERE = pathlib.Path(__file__)
HERE = HERE.relative_to(HERE.parent).parent

# The text of the README file
README = (HERE / "README.md").read_text()

#extra_files = package_files((HERE/'data').as_posix())

scripts =  [(HERE / 'scripts/build_local_hst_archive.py').as_posix(),
            (HERE / 'scripts/setup_cascade.py').as_posix(),
            (HERE / 'scripts/run_cascade.py').as_posix(),
            (HERE / 'scripts/run_cascade.sh').as_posix(),
            (HERE / 'scripts/run_stats.sh').as_posix()]

config = {
    'name': 'CASCADe-spectroscopy',
    'description': 'CASCADe : Calibration of trAnsit Spectroscopy using CAusal Data',
    'long_description': README,
    'long_description_content_type': "text/markdown",
    'author': 'Jeroen Bouwman',
    'url': 'https://jbouwman.gitlab.io/CASCADe/',
    'download_url': 'https://gitlab.com/jbouwman/CASCADe',
    'author_email': 'bouwman@mpia.de',
    'version': '1.1.5',
    'python_requires': '>=3.7, <=3.9',
    'license': 'GNU General Public License v3 (GPLv3)',
    'classifiers': ["Programming Language :: Python :: 3",
                    "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
                    "Operating System :: OS Independent",
                    'Intended Audience :: Science/Research',
                    'Topic :: Scientific/Engineering :: Astronomy',
                    ],
    'packages': setuptools.find_packages(exclude=("tests_private", "docs", "examples",)),
    'include_package_data': True,
    'package_data': {"": scripts},
    'install_requires': ['batman-package', 'astropy', 'jplephem', 'scipy',
                         'numpy==1.22.1', 'configparser', 'photutils', 'pandas',
                         'scikit-learn', 'matplotlib', 'tqdm', 'seaborn',
                         'pytest', 'scikit-image', 'sphinx', 'alabaster',
                         'networkx', 'cython', 'astroquery', 'numba==0.53.1',
                         'ray[default]', 'pyfiglet', 'termcolor',
                         'statsmodels',
                         'exotethys==2.0.2'],
    'scripts': scripts,
    'data_files': [('', [(HERE / 'README.md').as_posix(),
                         (HERE / 'LICENSE.txt').as_posix()])]
}

setuptools.setup(**config)
