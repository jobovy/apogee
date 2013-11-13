from setuptools import setup #, Extension

longDescription= ""


setup(name='apogee',
      version='1.',
      description='APOGEE data tools',
      author='Jo Bovy',
      author_email='bovy@ias.edu',
      license='New BSD',
      long_description=longDescription,
      url='https://github.com/jobovy/apogee',
      package_dir = {'apogee/': ''},
      packages=['apogee'],
      dependency_links = ['https://github.com/jobovy/galpy/tarball/master#egg=galpy'],
      install_requires=['numpy','scipy','fitsio','esutil','galpy']
      )
