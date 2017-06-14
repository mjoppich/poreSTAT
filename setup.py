from setuptools import setup

setup(name='porestat',
      version='1.0',
      description='Collection of python tools for examining sequencing experiments from FAST5 files',
      url='http://github.com/mjoppich/porestat',
      author='Markus Joppich',
      author_email='joppich@bio.ifi.lmu.de',
      license='GNU GPL v3',
      packages=['.', 'porestat.hdf5tool', 'porestat.tools', 'porestat.plots', 'porestat.utils', 'porestat.analysis'],
      package_data = {'': ['data/*.*']},
      zip_safe=False,
      install_requires=['numpy', 'scipy', 'matplotlib', 'pathos', 'openpyxl', 'h5py']
      )