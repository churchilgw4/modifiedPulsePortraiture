# from distutils.core import setup
from setuptools import setup, find_packages

setup(name='modifiedPulsePortraiture',
      version='1.0',
      description='Data analysis package for wideband and multiband wideband pulsar timing with MLA and MLAN technique for InPTA',
      author='Avinash Kumar Paladi',
      author_email='avinashkumarpaladi@gmail.com',
      url='https://github.com/AvinashKumarPaladi/modifiedPulsePortraiture',
      py_modules=[ 'mpplib', 'mpptoas', 'mpptoaslib', 'mpp_telescope_codes', 'mpplib_b35', 'mpptoas_b35', 'mpptoaslib_b35', 'mpptoaslib_MLAN', 'mpptoaslib_b35_MLAN', 'mpp_gaussian_test'],
      scripts=['mpptoas.py', 'mpptoas_b35.py', 'metafile_maker_cc_b35.py']
     )
