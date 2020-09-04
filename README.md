# Pipe network
> CB-Cities

[![License](https://img.shields.io/badge/license-MIT-blue.svg)](https://raw.githubusercontent.com/cb-cities/pipe-network/develop/license.md)
[![Developer docs](https://img.shields.io/badge/developer-docs-blue.svg)]()
[![User docs](https://img.shields.io/badge/user-docs-blue.svg)]()
[![CircleCI](https://circleci.com/gh/cb-cities/pipe-network.svg?style=svg)](https://circleci.com/gh/cb-cities/pipe-network)
[![codecov](https://codecov.io/gh/cb-cities/pipe-network/branch/develop/graph/badge.svg)](https://codecov.io/gh/cb-cities/pipe-network)
[![](https://img.shields.io/github/issues-raw/cb-cities/pipe-network.svg)](https://github.com/cb-cities/pipe-network/issues)
[![Project management](https://img.shields.io/badge/projects-view-ff69b4.svg)](https://github.com/cb-cities/pipe-network/projects/)

## Install dependencies

* Docker image for CB-Cities pipe-network code [https://hub.docker.com/r/cbcities/pipenetwork](https://hub.docker.com/r/cbcities/pipenetwork)

* Instructions for running pipe-network docker container: [https://github.com/cb-geo/docker-pipenetwork/blob/master/README.md](https://github.com/cb-cities/pipe-network-container/blob/master/README.md).

### Prerequisite packages
> The following prerequisite packages can be found in the docker image:

* [Eigen](http://eigen.tuxfamily.org/)
* [Intel MKL](https://software.intel.com/content/www/us/en/develop/tools/math-kernel-library.html)

## Compile and Run
> See https://pipe-network-doc.cb-cities.com/ for more detailed instructions. 

0. Run `mkdir build && cd build && cmake -DCMAKE_BUILD_TYPE=Release /path/to/CMakeLists.txt`.

1. Run `make clean && make -jN` (where N is the number of cores).

### Run tests

0. Run `./pipe-network-test -s` (for a verbose output) or `ctest -VV`.

### Run hydraulic simulation 
Run ./pipe-network with the following flags 

-f <filename>,  --file <filename>
     .inp file path for the WDN system

-t <save path>,  --to <save path>
     folder to save the simulation results

-d,  --debug
     debug mode: output intermediate results (residuals)
     
-n <meshname>,  --name <meshname>
     Name for the WDN
     
-h,  --help
     Displays usage information and exits.
