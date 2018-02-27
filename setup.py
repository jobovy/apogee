import os
from setuptools import setup #, Extension
import sys
import shutil
import subprocess
import tempfile

long_description = "Tools for APOGEE data analysis; see `here <https://github.com/jobovy/apogee>`__ for further documentation"

# Install FERRE when specifying --install-ferre; needs a FORTRAN compiler, e.g., http://hpc.sourceforge.net/
try:
    ferre_pos= sys.argv.index('--install-ferre')
except ValueError:
    _INSTALL_FERRE= False
else:
    del sys.argv[ferre_pos]
    _INSTALL_FERRE= True

try:
    ferre_openmp= sys.argv.index('--ferre-noopenmp')
except ValueError:
    _FERRE_NO_OPENMP= False
else:
    del sys.argv[ferre_openmp]
    _FERRE_NO_OPENMP= True

try:
    ferre_flen= sys.argv.index('--ferre-flen')
except ValueError:
    _FERRE_FLEN= 180
else:
    _FERRE_FLEN= int(sys.argv[ferre_flen+1])
    del sys.argv[ferre_flen]
    del sys.argv[ferre_flen]

if _INSTALL_FERRE:
    # Code to determine the binary install directory, from http://jasonstitt.com/setuptools-bin-directory
    from setuptools import Distribution
    from setuptools.command.install import install
    class OnlyGetScriptPath(install):
        def run(self):
            self.distribution.install_scripts = self.install_scripts
            
    def get_setuptools_script_dir():
        " Get the directory setuptools installs scripts to for current python "
        dist = Distribution({'cmdclass': {'install': OnlyGetScriptPath}})
        dist.dry_run = True  # not sure if necessary
        dist.parse_config_files()
        command = dist.get_command_obj('install')
        command.ensure_finalized()
        command.run()
        return dist.install_scripts

if _INSTALL_FERRE:
    # Download the code
    #_FERRE_FILE= 'ferre_4.5.6.tar.gz'
    _FERRE_FILE= 'ferre_4.6.6.tar.gz'
    #_FERRE_URL= 'http://leda.as.utexas.edu/ferre/%s' % _FERRE_FILE
    _FERRE_URL= 'http://www.as.utexas.edu/~hebe/ferre/%s' % _FERRE_FILE
    print('\033[1m'+"Downloading and installing FERRE from %s ..." % _FERRE_URL +'\033[0m')
    # Create temporary directory
    tmpdir= tempfile.mkdtemp(dir='./')
    os.mkdir(os.path.join(tmpdir,'ferre'))
    try:
        subprocess.check_call(['wget',_FERRE_URL,'-O',
                               os.path.join(tmpdir,'ferre',_FERRE_FILE)])
    except subprocess.CalledProcessError:
        print('\033[1m'+"Downloading FERRE from %s failed ..." % _FERRE_URL +'\033[0m')
    # Unpack and install
    os.chdir(os.path.join(tmpdir,'ferre'))
    try:
        subprocess.check_call(['tar','xvzf',_FERRE_FILE])
    except subprocess.CalledProcessError:
        print('\033[1m'+"Untarring/gunzipping FERRE failed ..." % _FERRE_URL +'\033[0m')
    os.chdir('src')
    # Change flen in share.f90
    with open("tmp.f90", "w") as fout:
        with open("share.f90", "r") as fin:
            for line in fin:
                fout.write(line.replace('flen=120','flen=%i' % _FERRE_FLEN))
    os.rename('tmp.f90','share.f90')
    # Change output format in ferre.f90
    with open("tmp.f90", "w") as fout:
        with open("ferre.f90", "r") as fin:
            for line in fin:
                fout.write(line.replace("write(3,'(1x,a30,100(1x,F9.3))')",
                                        "write(3,'(1x,a40,100(1x,F9.4))')"))
    os.rename('tmp.f90','ferre.f90')
    try:
        if _FERRE_NO_OPENMP:
            subprocess.check_call(['make','OPT=-O2'])
        else:
            subprocess.check_call(['make'])
    except subprocess.CalledProcessError:
        print('\033[1m'+"Compiling FERRE failed ..." % _FERRE_URL +'\033[0m')
    os.rename('a.out','../../../ferre')
    os.rename('ascii2bin','../../../ascii2bin')
    # Remove everything
    os.chdir('../../../')
    try:
        subprocess.check_call(['rm','-rf',tmpdir])
    except subprocess.CalledProcessError:
        print('\033[1m'+"Removing FERRE temporary files failed ..." % _FERRE_URL +'\033[0m')
    shutil.copy('ferre',get_setuptools_script_dir())
    shutil.copy('ascii2bin',get_setuptools_script_dir())

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
                'apogee/util','apogee/samples','apogee/spec','apogee/modelatm',
                'apogee/modelspec'],
      package_data={'apogee/samples':['data/rcmodel_mode_jkz_ks_parsec_newlogg.sav',
                                      'data/rcmodel_mode_jkz_h_parsec_newlogg.sav',
                                      'data/rcmodel_mass_agez.sav',
                                      'data/rcmodel_mass_agez_coarseage.sav',
                                      'data/rcmodel_omega_agez.sav'],
                    'apogee/spec':['filter/dr12/*.filt',
                                   'cannon/training/*.txt',
                                   'cannon/trained/*.txt'],
                    'apogee/modelspec':['scripts/makemoogmodel.awk'],},
      dependency_links = ['https://github.com/jobovy/galpy/tarball/master#egg=galpy',
                          'https://github.com/jobovy/isodist/tarball/master#egg=isodist'],
      install_requires=['numpy','scipy','matplotlib',
                        'astropy','galpy',
                        'isodist','periodictable','tqdm']
      )
