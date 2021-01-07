# pylsewave

![pylsewavelogo](./doc/sphinx-rootdir/figures/pylsewave_logo64x64.png)

[![Build Status](https://dev.azure.com/giorgosragos/giorgosr-pylsewave/_apis/build/status/giorgosR.pylsewave?branchName=master)](https://dev.azure.com/giorgosragos/giorgosr-pylsewave/_build/latest?definitionId=4&branchName=master)
![Docker Cloud Build Status](https://img.shields.io/docker/cloud/build/giorgosr/pylsewave)
[![Documentation Status](https://readthedocs.org/projects/pylsewave/badge/?version=latest)](https://pylsewave.readthedocs.io/en/latest/?badge=latest)
[![Anaconda-Server Badge](https://anaconda.org/giorgosr/pylsewave/badges/version.svg)](https://anaconda.org/giorgosr/pylsewave)
[![Anaconda-Server Badge](https://anaconda.org/giorgosr/pylsewave/badges/platforms.svg)](https://anaconda.org/giorgosr/pylsewave)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3426789.svg)](https://doi.org/10.5281/zenodo.3426789)


# Repository structure

```bash
.
+-- doc # documentation directory
+-- JupyterNbs # jupyter notebooks
+-- data # input data for case studies
+-- pylsewave # the package directory (contains all modules)
+-- test # unit test files
+. README.md
+. LICENCE
+. setup.py
+. pyproject.toml
+. Adan_77_example.py # case study for a whole arterial network (linear-elastic arterial walls)
+. vtk_example.py (example with vtk)
```

### A Python package to solve pulse wave dynamics in arterial networks

The following diagram depicts the rationale along with the structure of the package

![pylsewave toolkit](./JupyterNbs/images/pylsewave.png)

A c/python library to solve 1D pulse wave propagation in blood vessels or any hyperbolic PDE system in the form of

![Hyperbolic system](https://latex.codecogs.com/gif.latex?%5Cfrac%7B%5Cpartial%20%5Cbf%7BU%7D%7D%7B%5Cpartial%20t%7D%20&plus;%20%5Cfrac%7B%5Cpartial%20%5Cbf%7BF%7D%7D%7B%5Cpartial%20x%7D%20%3D%20%5Cbf%7BS%7D)

### Purpose

This library is designed to solve efficiently pulse wave dynamics in human arteries. It is written in python, cython and there are some bits of C++.

### Install pylsewave

To install pylsewave toolkit type:

```bash
conda install -c giorgosR pylsewave
```

### Run the example

There are two examples in this directory:

* Adan_77_example.py

To run the case go to the top level directory of the repository and type:

```bash
python Adan_77_example.py -ivesseldatafile data\Arterial_Network_ADAN56.txt -ibcinflowfile data\inflow_Aorta.txt -oresfile arterial_network_77_vessels -language py
```

`NOTE: You can execute the file with python classes translated via cython by changing the -language py to -language cy.`

* vtk_example.py

This is an example for storing the results in vtk multiblock file. The user should install pyvtk first (it can be found on [PYPI](https://pypi.org/project/vtk/)).

```bash
python vtk_example.py -resfile <resfile> -ovisfile <visfile>
```

### Examples in Jupyter notebooks

There are several case studies in jupyter notebooks under `JupyterNbs` directory. To run the cases, the user should install jupyter either with conda or pip (see the latest install instructions`` in [Jupyter documentaion](https://jupyter.readthedocs.io/en/latest/install.html)):

* conda

```bash
conda install jupyter
```

* pip

```bash
pip install jupyter
```

### Build the docs

To build the documentation, under the `docs` directory type (you will need sphinx):

```bash
make html
```

### Cite pylsewave

Prefered citation style for pylsewave:

Georgios E. Ragkousis (2019). PylseWave: A python package to solve pulse wave dynamics in arterial networks. Zenodo. [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3426789.svg)](https://doi.org/10.5281/zenodo.3426789)