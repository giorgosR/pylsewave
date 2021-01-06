# setup.py
# build from python setup.py build_ext --build-lib <build directory>
import os
from Cython.Build import cythonize
import numpy
try:
    from setuptools import setup
    from setuptools import find_packages
    from setuptools import Extension
except ImportError:
    from distutils.core import setup
    from distutils.extension import Extension
    from distutils.extension import find_packages
from distutils.command.build_ext import build_ext

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

ext_cynum = Extension("pylsewave.cynum",
                      sources=["pylsewave/cynum.pyx",
                               "pylsewave/cypwfdm.cpp",
                               "pylsewave/cypwmesh.cpp",
                               "pylsewave/cypwfuns.cpp"],
                      language='c++',
                      include_dirs=["pylsewave/include/", numpy.get_include()])

class build_ext_compiler_check(build_ext):
    def build_extensions(self):
        compiler = self.compiler.compiler_type
        print('\n\ncompiler', compiler)
        if 'msvc' in compiler:
            print("compiling on Windows")
            for extension in self.extensions:
                if extension.name == "pylsewave.cynum":
                    extension.extra_compile_args.append('/O2')
                    extension.extra_compile_args.append('/openmp')
                    extension.extra_compile_args.append('/EHsc')
                    extension.extra_compile_args.append('/GA')
                    extension.extra_compile_args.append('/favor:INTEL64')

        if 'unix' in compiler:
            print("compiling on unix")
            for extension in self.extensions:
                if extension.name == "pylsewave.cynum":
                    extension.extra_compile_args.append('-O2')
                    extension.extra_compile_args.append('-fopenmp')
                    extension.extra_compile_args.append('-std=c++11')
                    extension.extra_link_args = ['-lgomp']

        super().build_extensions()


setup(name="pylsewave",
      version='1.0.2',
      license="GNU GPL v3.0",
      packages=find_packages(),
      description='A python package for pulse wave dynamics and/or any hyperbolic system of PDEs',
      author='Georgios E. Ragkousis',
      author_email='giorgosragos@gmail.com',
      url = "https://giorag.bitbucket.io/pylsewave/pyw_doc.html",
          # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Education',
        'Intended Audience :: Developers',
        'Intended Audience :: Healthcare Industry',
        'Intended Audience :: Information Technology',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Mathematics',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Operating System :: OS Independent'],
      keywords='pdes fdm pulsewave blood-vessels',
      long_description=read('README.md'),
	  long_description_content_type='text/markdown',
      requires=['numpy', "scipy", 'matplotlib'],
      install_requires=['numpy', 'scipy', 'matplotlib'],
      cmdclass={ 'build_ext': build_ext_compiler_check },
      ext_modules=cythonize([ext_cynum], annotate=True))