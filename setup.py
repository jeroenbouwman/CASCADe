try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

config = {
    'description': 'CASCADe : Calibration of trAnsit Spectroscopy using CAusal Data',
    'author': 'Jeroen Bouwman',
    'url': 'http://www.mpia.de',
    'download_url': 'http://www.mpia.de',
    'author_email': 'bouwman.mpia.de',
    'version': '0.1',
    'install_requires': ['nose','batman','astropy','scipy','numpy','configparser'],
    'packages': ['cascade'],
    'scripts': [],
    'name': 'CASCADe'
}

setup(**config)
