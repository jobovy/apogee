from setuptools import setup #, Extension

with open('README.md') as dfile:
    long_description = dfile.read()


setup(name='apogee',
      version='1.',
      description='APOGEE data tools',
      author='Jo Bovy',
      author_email='bovy@ias.edu',
      license='New BSD',
      long_description=long_description,
      url='https://github.com/jobovy/apogee',
      package_dir = {'apogee/': ''},
      packages=['apogee','apogee/tools','apogee/select','apogee/test',
                'apogee/util','apogee/samples'],
      package_data={'apogee/samples':['data/rcmodel_mode_jkz_ks_parsec_newlogg.sav']},
      dependency_links = ['https://github.com/jobovy/galpy/tarball/master#egg=galpy',
                          'https://github.com/jobovy/isodist/tarball/master#egg=isodist'],
      install_requires=['numpy','scipy','matplotlib',
                        'fitsio','esutil','galpy',
                        'isodist']
      )
