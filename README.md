# crysfipy - Tools for analysis of physical properties arising from crystal-field Hamiltonian

For list of modules and their descriptions see documentation at: 
TODO

Installation
1. Create an empty virtual environment with the python version specified in `setup.py: python-version`
2. Install `crysfipy` locally with `python setup.py` or `pip install .`.
3. TODO PyPi availability

Installation FAQ:
- For a lightweight, fast installation miniconda is a great solution. A throwback is, its package manager conda or pip will nominally install the newest version of packages from their channels, which may produce problems. A good solution is to check the list of packages installed by Anaconda for your platform (https://docs.anaconda.com/anaconda/packages/pkg-docs/) and follow the version of packages from there. Otherwise, one might need to carefully test which version of each packages are compatible with each other.
- setup.py file should contain the compatible combination of packages.
- Install with `python setup.py develop` or  `pip install -e .` in development mode, so that all changes in the source code will work for future local scripts.
- For full development utilities, development mode +packages to build the docs, do `pip install -e .[dev]`.
- One can also just download all source files from this Github project to a location `PATH`, and in each python file include lines
  > import sys
  >
  > sys.path.append(PATH)

- numpy==1.23.1 requires python>3.10 I believe. It didnt want to work with python==3.8