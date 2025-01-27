import distutils.core
from distutils.core import setup

version_number = '0.2'
setup(name = 'qvalue',
      version = version_number,
      description = 'Converts p-values in q-values in order to account for multiple hypotheses testing, see (Storey and Tibshirani, 2003)',
      long_description = open('README.md').read(),
      author = 'Nicolo Fusi, Neo Christopher Chung',
      author_email = 'nicolo.fusi@sheffield.ac.uk, nchchung@gmail.com',
      packages = ['qvalue'],
      requires = ['numpy (>=1.24.1)',
                  'scipy (>=1.13.1)',
                  'matplotlib (>=3.9.2)'],
      license = '3-clause BSD',
      )
