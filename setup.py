import os
import sys
import subprocess

from setuptools import Command, find_packages, setup
from setuptools.command.install import install
from setuptools.command.develop import develop
from setuptools.command.egg_info import egg_info
from LineageTracker import __VERSION__


MIN_PY_VER = (3, 6)
if sys.version_info[:2] < MIN_PY_VER:
    sys.stderr.write(
        ("Biopython requires Python %i.%i or later. " % MIN_PY_VER)
        + ("Python %d.%d detected.\n" % sys.version_info[:2])
    )
    sys.exit(1)


required = ['pandas',
            'numpy',
            'scikit-learn',
            'scipy',
            'matplotlib',
            'ete3',
            'networkx',
            'BioPython',
            'pysam',
            'pyqt5',
            'adjustText']
request = ['pygraphviz']


def get_virtualenv_path():
    '''Used to work out path to install compiled binaries to.'''

    if hasattr(sys, 'real_prefix'):
        return sys.prefix

    if hasattr(sys, 'base_prefix') and sys.base_prefix != sys.prefix:
        return sys.prefix

    if 'conda' in sys.prefix:
        return sys.prefix

    return None


def compile_and_install_software():

    '''Used the subprocess module to compile/install the C software.'''
    src_path = './LineageTracker/src/'

    # install the software (into the virtualenv bin dir if present)
    try:
        subprocess.check_call('make', cwd=src_path, shell=True)
    except:
        print('make command is required for installation')

    try:
        import platform
        system = platform.system()
        if platform == 'Windows':
            subprocess.check_call('pip install pygraphviz', shell=True)
        else:
            subprocess.check_call('pip install --install-option="--include-path=/usr/local/include/" --install-option="--library-path=/usr/local/lib/" pygraphviz', shell=True)
    except:
        print('Failed to install pygraphviz')

class CustomInstall(install):
    '''Custom handler for the 'install' command.'''

    def run(self):
        install.run(self)
        compile_and_install_software()
        super().run()

class CustomDevelopCommand(develop):

    def run(self):
        develop.run(self)
        compile_and_install_software()


class CustomEggInfoCommand(egg_info):

    def run(self):
        egg_info.run(self)
        compile_and_install_software()


class TestProgram(Command):

    description = 'run to test Y-LineageTracker'
    user_options = [('cmd', None, 'test all commands'),]

    def initialize_options(self):
        '''Set default values for options'''
        self.cmd = 'all'

    def finalize_options(self):
        '''Post-process options'''
        print('[Y-LineageTracker] [Test] Testing started')

    def run(self):
      '''Run command'''
      import importlib
      for i in required:
          found = importlib.util.find_spec(i)
      if not found:
          print('[Y-LineageTracker] [Test] Package [%s] is required for testing')
          sys.exit()

      current_dir = os.path.split(os.path.realpath(__file__))[0]
      sys.path.insert(0, current_dir+ '/Test')
      from Test import TestProgram
      if self.cmd == 'all':
          TestProgram.main()


setup(name='Y-LineageTracker',
      version=__VERSION__,
      author='Hao Chen',
      author_email='chenhao@picb.ac.cn',
      description='Y-LineageTracker: A framework to fully analyze human Y-chromosome sequencing data',
      license='GPL-3.0',
      keywords=['Bioinformatics', 'Computational biology', 'Population genetics',
                'Genomics', 'Y chromosome', 'haplogroup', 'Y-STR'],
      url='http://www.picb.ac.cn/PGG/resource.php',
      packages=['LineageTracker',
                'LineageTracker.ProcessData',
                'LineageTracker.GetLineage',
                'LineageTracker.Test'],
      package_data={'LineageTracker': ['Data/*', 'sans-serif/*', 'src/*', 'Test/TestData/*', 'Config.ini']},
      entry_points={'console_scripts':
                    ['LineageTracker = LineageTracker.RunLineagerTracker:main']},
      install_requires=required,
      classifiers=['Development Status :: 3 - Alpha',
                   'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
                   'Intended Audience :: Developers',
                   'Intended Audience :: Science/Research',
                   'Programming Language :: Python :: 3',
                   'Programming Language :: Python :: 3.6',
                   'Programming Language :: C',
                   'Topic :: Scientific/Engineering :: Bio-Informatics',
                   'Topic :: Scientific/Engineering :: Visualization'],
      zip_safe=False,
      cmdclass={'test': TestProgram,
                'install': CustomInstall,
                'develop': CustomDevelopCommand,
                'egg_info': CustomEggInfoCommand},
      python_requires=">=%i.%i" % MIN_PY_VER)
