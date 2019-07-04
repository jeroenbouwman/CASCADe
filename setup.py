try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

config = {
    'description': 'CASCADe : Calibration of trAnsit Spectroscopy using CAusal Data',
    'author': 'Jeroen Bouwman',
    'url': 'https://gitlab.com/jbouwman/CASCADe/wikis/home',
    'download_url': 'https://gitlab.com/jbouwman/CASCADe',
    'author_email': 'bouwman@mpia.de',
    'version': '1.0',
    'install_requires': ['batman-package', 'astropy', 'jplephem', 'scipy',
                         'numpy', 'configparser', 'photutils', 'pandas',
			             'scikit-learn', 'matplotlib', 'tqdm', 'seaborn',
                         'pytest', 'scikit-image', 'sphinx', 'alabaster',
                         'networkx', 'cython', 'astroquery'],
    'packages': ['cascade'],
    'scripts': [],
    'name': 'CASCADe'
}

setup(**config)

