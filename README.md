# DwellRegions-open

Implementation for online and offline computation of dwell regions from trajectory data.
Please see the accompanying paper https://ieeexplore.ieee.org/document/6341396/.

# Building code in this repository

```
# Download code
git clone https://github.com/PayasR/DwellRegions-open.git

# Third party dependencies: 
# 1. libgeographic-dev on Debian-based distros, geographiclib-devel on CentOS/RHEL.
# 2. Boost program options.
# 3. GSL (source included in the source tree).
# On Windows, please use vcpkg to install all the above packages, and a recent version of Visual Studio 2019.

# Build code:
mkdir build
cd build
cmake3 --DCMAKE_BUILD_TYPE=Release ../
make -j6

# Run tests
ctest

# Run online dwell degion code
./dwell_region_exists --help

# Run offline dwell degion code
./offline_dwell_region_exists --help
```

