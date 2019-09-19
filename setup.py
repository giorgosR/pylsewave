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

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

ext_cynum = Extension("pylsewave.cynum",
                      sources=["pylsewave/cynum.pyx",
                               "pylsewave/cypwfdm.cpp",
                               "pylsewave/cypwmesh.cpp",
                               "pylsewave/cypwfuns.cpp"],
                      language='c++',
                      extra_compile_args=['/Ot', '/openmp', '/EHsc', '/GA', '/favor:INTEL64'],
                      # extra_link_args=['-openmp'],
                      # include_path = ["./pylsewave/include/",],
                      include_dirs=["pylsewave/include/", numpy.get_include()]
                      )

setup(name="pylsewave",
      version='1.0.0',
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
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Environment :: Win32 (MS Windows)'],
      keywords='pdes fdm pulsewave blood-vessels',
      long_description=read('README.md'),
	  long_description_content_type='text/markdown',
      requires=['numpy', "scipy", 'matplotlib'],
      install_requires=['numpy', 'scipy', 'matplotlib'],
      ext_modules=cythonize([ext_cynum], annotate=True))