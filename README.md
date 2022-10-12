[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4502942.svg)](https://doi.org/10.5281/zenodo.4502942)

# <img src="img/AmolqcLogo.png" alt="Amolqc Logo" height="128"/><img src="img/AmolqcIconTransparent.png" alt="Amolqc Icon" height="128"/>
**Atoms and Molecules with Quantum Monte Carlo**

Amolqc is an open source software package for electron structure calculations with quantum Monte Carlo (QMC).
The code features strong multideterminant and optimization capabilities.

## Building Amolqc

### Build with CMake and Make

CMake version 3.10 or higher is required.

The same settings as in the make/configs/\*.mk files can be set in a cmake/configs/\*.cmake file.
- create your own file (eg. cmake/configs/myconfig.cmake)

If the environment variable `$FC` is set, the compiler from that variable is taken instead of
the one given in configs.

Then in "Amolqc" do:
```
./configure myconfig
mkdir build
cd build
cmake [-DPFUNIT=OFF] ..
make [-j3]
```

To run the serial testsuite, do:
("-R jas" can be given to run only tests matching "jas".
Other examples: "-R ECP", "-R 05")
```
ctest -L serial [-R jas] [--verbose]
```

To run the parallel testsuite, do:
```
ctest -L parallel [-R jas] [--verbose]
```

The executable is:
"Amolqc/build/bin/amolqc"

Add the following line to cmake/config.cmake to get compile warnings:
```
set(WARNINGS ON)
```

Note: on some machines, cmake and ctest are executed with 'cmake3' and 'ctest3'.

## Running Amolqc

Read the manual in **[doc/user_manual.txt](doc/user_manual.txt)** on how to run Amolqc.

## Contributing to Amolqc

Please read **[CONTRIBUTING.md](CONTRIBUTING.md)** on how to contribute to Amolqc.
