from setuptools import setup

setup(name='porestat',
      version='1.0',
      description='Collection of python tools for examining sequencing experiments from FAST5 files',
      url='http://github.com/mjoppich/porestat',
      author='Markus Joppich',
      author_email='joppich@bio.ifi.lmu.de',
      license='GNU GPL v3',
      packages=['.', 'porestat.hdf5tool', 'porestat.tools', 'porestat.plots', 'porestat.utils', 'porestat.analysis'],
      package_data = {'': ['porestat/data/*.*']},
      zip_safe=False,
      install_requires=['UpSetPlot', 'h5py', 'htseq', 'matplotlib', 'matplotlib-venn', 'mpld3', 'natsort', 'networkx', 'numpy', 'openpyxl', 'pandas', 'pathos', 'scikit-learn', 'scipy', 'seaborn', 'statsmodels', 'upsetplot', 'xlsxwriter']
      )
