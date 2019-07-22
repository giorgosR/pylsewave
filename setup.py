# setup.py
# build from python setup.py build_ext --build-lib <build directory>
import os
# from setuptools import setup, find_packages
# from distutils.core import Extension
from Cython.Distutils import build_ext
import numpy
from Cython.Build import cythonize
try:
    from setuptools import setup
    from setuptools import find_packages
    from setuptools import Extension
except ImportError:
    from distutils.core import setup
    from distutils.extension import Extension
    from distutils.extension import find_packages

# ext = Extension("mt_random",
#                 sources=["mt_random.pyx", "mt19937ar.c"])
#
# setup(name="mersenne_random",
#       ext_modules = cythonize([ext]))

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
                      include_dirs=["pylsewave/", numpy.get_include()]
                      )

setup(name="pylseWave",
      version='1.0.0',
      packages=find_packages(),
      # package_data={
      #     'cynum': ['./pylsewave/*.pxd']
      # },
      description='A python package for pulse wave dynamics and/or any hyperbolic system of PDEs',
      author='Georgios E. Ragkousis',
      author_email='giorgosragos@gmail.com',
      url = "",
      long_description=read('README.md'),
      ext_modules=cythonize([ext_cynum], annotate=True))
# setup(
#       cmdclass = {'build_ext': build_ext},
#       ext_modules = [Extension("cynum",
#                                ["pulsewavepy/cynum.pyx", "pulsewavepy/cypwfdm.cpp",
#                                 "pulsewavepy/cypwmesh.cpp",
#                                 "pulsewavepy/cypwfuns.cpp"],
#                                language='c++',
#                                include_dirs=["pulsewavepy/", numpy.get_include()]),
#                      Extension("cyfdm",
#                                ["pulsewavepy/cyfdm.pyx"],
#                                 language='c++',
#                                include_dirs=[numpy.get_include()])
#                      # Extension("cymesh",
#                      #           ["pulsewavepy/cymesh.pyx"],
#                      #           language='c++',
#                      #           include_dirs=["pulsewavepy/", numpy.get_include()])
#                      ]
# )


# pysmall = Extension('geom',
# sources = ['pulsewavepy/geom.pyx',
# 'pulsewavepy/Circle.cpp'],
# include_dirs = ['pulsewavepy/include/'],
# language="c++")

# setup(name='geom',
# packages=['pulsewavepy'],
# #      install_requires=['cython==0.17'],
# ext_modules=[pysmall],
# cmdclass = {'build_ext': build_ext})


# from distutils . core import setup
# from distutils . extension import Extension
# from Cython . Distutils import build_ext
# setup (
# cmdclass = {'build_ext ': build_ext },
# ext_modules = [
# Extension("circ", ["circ.pyx", "Circle.cpp"],
# language="c++") ,
# Extension (" trig ", [" trig .pyx"],
# libraries =["m"]),
# ])