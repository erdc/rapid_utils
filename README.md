# rapid_utils

Description
-----------
This repository contains tools for using the Routing Application for Parallel
computatIon of Discharge (RAPID; David et al., 2011) model. Much of the
functionality is adapted from the RAPIDpy codebase written by Alan Snow.

Download
--------
HTTPS:
$ git clone https://github.com/erdc/rapid_utils.git

SSH:
$ git clone git@github.com:erdc/rapid_utils.git

External dependencies
---------------------
 - geopandas
 - netcdf4
 - numpy
 - python
 - pytest-cov

If conda is installed, the following commands will
create and activate an environment (named "rapid_utils" by default)
that contains the required dependencies.

$ conda env create -f environment.yml
$ conda activate rapid_utils

Installation
------------
Execute one of the following commands in the top-level directory of the
repository, i.e. where the file setup.py is located.

Install package for operational use:
$ pip install .

Install package for development:
$ pip install -e .

Testing
-------
Unit tests:
$ pytest

Unit tests with coverage report:
$ pytest --cov --cov-report=html
(creates directory with HTML files reporting test coverage)

RAPID input
-----------
RAPID requires five input files for a forward model simulation. Each of these
files is discussed briefly below.

**vlat (netCDF)**
This contains a variable named "m3_riv", which is an array of lateral inflow
time series associated with each reach in the computational stream network.
The lateral inflow is the primary dynamic forcing for the RAPID mode. "m3_riv"
has dimensions (time, rivid), where the "time" dimension is the number of time
steps represented and the "rivid" dimension is the number of stream reaches in
the network. Optionally, the variables "rivid", "time", "time_bnds", "lat",
"lon", "crs", and "m3_riv_err" and metadata ("title", "institution", "comment")
may also be included. Optional variables and metadata are transferred to the
RAPID output `Qout_file`.

**rapid_connect (CSV)**
This provides connectivity information as a rectangular array for the stream
network used for a RAPID simulation. The first column of the array provides the
(unique) identifier for a reach in the network; the second column holds the
corresponding downstream reach identifier; and the third column holds the
number of reaches directly upstream of the reach indicated in the first column.
Subsequent columns hold the identifiers for all reaches directly upstream. Once
all directly upstream identifiers are listed, all remaining columns are
populated with "0".

**k (CSV)**
This file contains a single column with the Muskingum `k` parameter associated
with each river reach in the stream network. The values in the column should
have the same ordering as the corresponding identifiers in the first column of
the `rapid_connect` file.

**x (CSV)**
This file contains a single column with the Muskingum `x` parameter associated
with each river reach in the stream network. The values in the column should
have the same ordering as the corresponding identifiers in the first column of
the `rapid_connect` file.

**riv_bas_id (CSV)**
This file contains a single column whose values are the reach identifiers to be
used for a simulation. This may be all of the reaches in the first column of
the `rapid_connect` file or a subdomain. An upstream to downstream ordering
results in greater computational efficiency.





