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
    'install_requires': ['nose','batman','astropy','scipy','numpy','configparser',
			 'os', 'weakref', 'sklearn', 'matplotlib', 'tqdm', 'seaborn',
			 'skimage', 'warnings', 'abc', 'sphinx'],
    'packages': ['cascade'],
    'scripts': [],
    'name': 'CASCADe'
}

setup(**config)
