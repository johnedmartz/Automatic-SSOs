import os
from setuptools import setup, find_packages

def read(fname):
      return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
      name='ssos_calib',
      version='1.0',
      author='John Martinez',
      author_email='johnedmf@gmail.com',
      description='Automatic data calibration of SSOS pipeline using Fialbres and GLS',
      packages=find_packages(),
      package_data={'':['config.ssos_calib', 'config.scamp', 'config.sex', 'default.param']},
      license='BSD',
      keywords='filabres ssos gls',
      long_description=read('README.md'),
      entry_points={'console_scripts': ['ssos_calib = ssos_calib.__main__:main']},
      include_package_data=True
      )
