"""Tools to write input input discharge files for RAPID."""

from datetime import datetime
import warnings
import numpy as np
from netCDF4 import Dataset

class RAPIDInputDischarge:

    """A class for writing input discharge netCDF files (initial conditions
    and forcing) for the Routing Application for Parallel Computation of
    Discharge (RAPID) model.
    """

    def __init__(self,
                 output_filename='rapid_input_discharge.nc',
                 input_file_list=None,
                 discharge=None,
                 time=None,
                 rivid=None,
                 latitude=None,
                 longitude=None,
                 input_discharge=None,
                 input_discharge_file=None,
                 input_discharge_rivid=None,
                 rivid_file=None,
                 input_discharge_time_axis=0,
                 input_discharge_time_index=-1,
                 discharge_datatype='f8',
                 rivid_datatype='i4',
                 time_datatype='i8',
                 latlon_datatype='f8',
                 crs_datatype='i4',
                 time_units='seconds since 1970-01-01 00:00:00',
                 discharge_variable_name='Qout',
                 time_variable_name='time',
                 rivid_variable_name='rivid',
                 latitude_variable_name='lat',
                 longitude_variable_name='lon',
                 integration_type='mean'):
        """
        Parameters
        ----------
        output_filename : str (optional)
            Name of output file.
        input_file_list : list (optional)
            List of input discharge (netCDF) filenames.
        discharge : array_like (optional)
            Discharge values.
        time : array_like (optional)
            Integer times corresponding to discharge.
        rivid : array_like (optional)
            Identifiers for river segments.
        latitude : array_like (optional)
            Latitude coordinate associated with each stream reach.
        longitude : array_like (optional)
            Longitude coordinate associated with each stream reach.
        input_discharge_file : str (optional)
            Name of (netCDF) file containing discharge values.
        rivid_file : str (optional)
            Name of (text) file containing ordered list of river identifiers
            for locations where discharge should be extracted/reported.
        input_discharge_time_axis : int (optional)
            Discharge array axis corresponding to the time dimension.
        input_discharge_time_index : int (optional)
            Index along the time axis where values should be extracted from
            a discharge array.
        discharge_datatype : str (optional)
            Data type to be used for writing discharge variable.
        rivid_datatype : str (optional)
            Data type to be used for writing rivid variable.
        time_datatype : str (optional)
            Data type to be used for writing time variable.
        latlon_datatype : str (optional)
            Data type to be used for writing lat/lon variables.
        crs_datatype : str (optional)
            Data type to be used for writing coordinate reference variable.
        time_units : str (optional)
            Units and datum for time variable.
        discharge_variable_name : str (optional)
            Name of discharge variable to be used in output file.
        time_variable_name : str (optional)
            Name of time variable to be used in output file.
        rivid_variable_name : str (optional)
            Name of rivid variable to be used in output file.
        latitude_variable_name : str (optional)
            Name of latitude variable to be used in output file.
        longitude_variable_name : str (optional)
            Name of longitude variable to be used in output file.
        integration_type : str (optional)
            Method to be used for integrating discharge from multiple files.
        """
        self.output_filename = output_filename
        self.input_file_list = input_file_list
        self.discharge = discharge
        self.time = time
        self.rivid = rivid
        self.latitude = latitude
        self.longitude = longitude
        self.input_discharge = input_discharge
        self.input_discharge_rivid = input_discharge_rivid
        self.input_discharge_file = input_discharge_file
        self.rivid_file = rivid_file
        self.input_discharge_time_axis = input_discharge_time_axis
        self.input_discharge_time_index = input_discharge_time_index
        self.discharge_datatype = discharge_datatype
        self.rivid_datatype = rivid_datatype
        self.time_datatype = time_datatype
        self.latlon_datatype = latlon_datatype
        self.crs_datatype = crs_datatype
        self.time_units = time_units
        self.discharge_variable_name = discharge_variable_name
        self.time_variable_name = time_variable_name
        self.rivid_variable_name = rivid_variable_name
        self.latitude_variable_name = latitude_variable_name
        self.longitude_variable_name = longitude_variable_name
        self.integration_type = integration_type

        # Determine slice to be used for input discharge
        if self.input_discharge_time_axis == 0:
            self.input_discharge_slice = (
                slice(self.input_discharge_time_index,None, None),
                slice(None, None, None))
        elif self.input_discharge_time_axis == 1:
            self.input_discharge_slice = (
                slice(None, None, None),
                slice(self.input_discharge_time_index, None, None))

    def write_nc(self, discharge=None, time=None, rivid=None, latitude=None,
                 longitude=None, output_filename=None):
        """Write variables, dimensions, and attributes to output file.

        Parameters
        ----------
        discharge : array_like (optional)
            Discharge values.
        time : array_like (optional)
            Integer times corresponding to discharge.
        rivid : array_like (optional)
            Identifiers for river segments.
        latitude : array_like (optional)
            Latitude coordinate associated with each stream reach.
        longitude : array_like (optional)
            Longitude coordinate associated with each stream reach.
        output_filename : str (optional)
            Name of output file.
        """
        if discharge is not None:
            self.discharge = discharge
        if time is not None:
            self.time = time
        if rivid is not None:
            self.rivid = rivid
        if latitude is not None:
            self.latitude = latitude
        if longitude is not None:
            self.longitude = longitude
        if output_filename is not None:
            self.output_filename = output_filename

        # Input discharge files should only contain a single time step.
        ntime = 1

        if self.rivid is not None:
            nrivid = len(self.rivid)
        else:
            nrivid = None

        data_out_nc = Dataset(self.output_filename, 'w')

        # Create dimensions
        data_out_nc.createDimension(self.time_variable_name, ntime)
        data_out_nc.createDimension(self.rivid_variable_name, nrivid)

        # Create variables
        # Discharge
        discharge_var = data_out_nc.createVariable(
            self.discharge_variable_name, self.discharge_datatype,
            (self.time_variable_name, self.rivid_variable_name),
            fill_value=0)
        discharge_var.long_name = ('instantaneous river water discharge ' +
                                   'downstream of each river reach')
        discharge_var.units = 'm3 s-1'
        discharge_var.coordinates = f'{self.longitude_variable_name} ' + \
                                    f'{self.latitude_variable_name}'
        discharge_var.grid_mapping = 'crs'
        discharge_var.cell_methods = "time: point"
        if self.discharge is not None:
            discharge_var[:] = self.discharge

        # rivid
        rivid_var = data_out_nc.createVariable(
            self.rivid_variable_name, self.rivid_datatype,
            (self.rivid_variable_name,))
        rivid_var.long_name = 'unique identifier for each river reach'
        rivid_var.units = '1'
        rivid_var.cf_role = 'timeseries_id'
        if self.rivid is not None:
            rivid_var[:] = self.rivid

        # Time
        time_var = data_out_nc.createVariable(
            self.time_variable_name, self.time_datatype,
            (self.time_variable_name,), fill_value=-9999.0)
        time_var.long_name = 'time'
        time_var.standard_name = 'time'
        time_var.units = self.time_units
        time_var.axis = 'T'
        time_var.calendar = 'gregorian'
        time_var.bounds = 'time_bnds'
        if self.time is not None:
            time_var[:] = self.time
        else:
            time_var[:] = np.array([-9999.0])

        # Longitude
        lon_var = data_out_nc.createVariable(
            self.longitude_variable_name, self.latlon_datatype,
            (self.rivid_variable_name,), fill_value=-9999.0)
        lon_var.long_name = \
            'longitude of a point related to each river reach'
        lon_var.standard_name = 'longitude'
        lon_var.units = 'degrees_east'
        lon_var.axis = 'X'
        if self.longitude is not None:
            lon_var[:] = self.longitude

        # Latitude
        lat_var = data_out_nc.createVariable(
            self.latitude_variable_name, self.latlon_datatype,
            (self.rivid_variable_name,), fill_value=-9999.0)
        lat_var.long_name = \
            'latitude of a point related to each river reach'
        lat_var.standard_name = 'latitude'
        lat_var.units = 'degrees_north'
        lat_var.axis = 'Y'
        if self.latitude is not None:
            lat_var[:] = self.latitude

        # Coordinate reference system
        crs_var = data_out_nc.createVariable('crs', self.crs_datatype)
        crs_var.grid_mapping_name = 'latitude_longitude'
        crs_var.epsg_code = 'EPSG:4326'  # WGS 84
        crs_var.semi_major_axis = 6378137.0
        crs_var.inverse_flattening = 298.257223563

        # History
        data_out_nc.history = f'date_created: {datetime.utcnow()}'

        # Close file
        data_out_nc.close()

    def parse_input_discharge_file(self, input_discharge_file=None,
                                   discharge_only=False):
        """Extract subset of discharge and (optionally) associated variables
        (time, rivid) from netCDF file.

        Parameters
        ----------
        input_discharge_file : str (optional)
            Name of (netCDF) file containing discharge values.
        discharge_only : bool (optional)
            If True, only extract discharge values. Otherwise, additionally
            extract time and rivid values.
        """
        if input_discharge_file is not None:
            self.input_discharge_file = input_discharge_file

        try:
            dataset = Dataset(self.input_discharge_file)
        except IOError:
            warnings.warn(
                f'Unable to open file {self.input_discharge_file}.')
            dataset = None

        if dataset is not None:
            try:
                self.input_discharge = dataset[
                    self.discharge_variable_name][
                    self.input_discharge_slice].squeeze()
            except (IOError, IndexError):
                warnings.warn(
                    f'Variable "{self.discharge_variable_name}" not found ' +
                    f'in file {self.input_discharge_file}.')

        if not discharge_only:

            try:
                self.input_discharge_rivid = dataset['rivid'][:]
            except (IOError, IndexError):
                warnings.warn(
                    'Variable "rivid" not found ' +
                    f'in file {self.input_discharge_file}.')

            try:
                time_var = dataset['time']
            except IOError:
                warnings.warn(
                    'Variable "time" not found ' +
                    f'in file {self.input_discharge_file}.')

            if time_var is not None:
                try:
                    time_step = time_var[1] - time_var[0]
                except IndexError:
                    time_step = None

                # Increment final simulation time t_n by time step to give
                # value t_(n+1). t_(n+1) is assumed to be the first time step
                # when initializing a new run with discharge from t_n.
                final_time = time_var[-1]

                try:
                    next_time = final_time + time_step
                except TypeError:
                    next_time = np.array([])

                self.time = next_time

    def parse_rivid_file(self, rivid_file=None, delimiter=',', usecols=0):
        """Extract rivid values from a text file. Values in file are assumed
        to be organized in a rectangular array.

        Parameters
        ----------
        rivid_file : str (optional)
            Name of (text) file containing river identifiers in a single column.
        delimiter : str (optional)
            Character indicating how values are separated within `rivid_file`.
        usecols : int (optional)
            Index of
        usecols : int (optiona)
        """
        if rivid_file is not None:
            self.rivid_file = rivid_file

        try:
            self.rivid = np.genfromtxt(
                self.rivid_file, delimiter=delimiter, dtype=int,
                usecols=usecols)
        except (FileNotFoundError, UnicodeDecodeError):
            warnings.warn(f'Unable to read file {self.rivid_file}.')
            self.rivid = None

    def sort_discharge_by_rivid(self, input_discharge_file=None,
                                rivid=None, rivid_file=None):
        """Filter and order the discharge according to the rivid array.

        Parameters
        ----------
        input_discharge_file : str (optional)
            Name of (netCDF) file containing discharge values.
        rivid : array_like (optional)
            River segment identifiers to be used for filtering and ordering.
        rivid_file : str (optional)
            Name of (text) file containing river identifiers in a single column.
        """
        if input_discharge_file is not None:
            self.parse_input_discharge_file(
                input_discharge_file=input_discharge_file)
        if rivid_file is not None:
            self.parse_rivid_file(rivid_file=rivid_file)
        if rivid is not None:
            self.rivid = rivid

        if self.rivid is None:
            warnings.warn('No rivid array specified.')
            has_input = False
        elif self.input_discharge_rivid is None:
            warnings.warn('No input discharge rivid array specified.')
            has_input = False
        elif self.input_discharge is None:
            warnings.warn('No input discharge specified.')
            has_input = False
        else:
            has_input = True

        is_sorted = np.array_equal(self.rivid, self.input_discharge_rivid)

        if has_input and not is_sorted:
            sorted_rivid_indices = np.argsort(self.rivid)
            unsort_rivid_indices = np.argsort(sorted_rivid_indices)

            # Cast attributes `rivid`, `input_discharge`, and
            # `input_discharge_rivid` to numpy arrays to allow
            # array indexing.
            sorted_rivid = np.array(self.rivid)[sorted_rivid_indices]

            sorted_input_rivid_indices = np.argsort(
                self.input_discharge_rivid)
            sorted_input_rivid = np.array(self.input_discharge_rivid)[
                sorted_input_rivid_indices]
            sorted_input_discharge = np.array(self.input_discharge)[
                sorted_input_rivid_indices]

            extract_indices = np.isin(sorted_input_rivid, sorted_rivid)
            extracted_discharge = sorted_input_discharge[extract_indices]
            extracted_discharge = extracted_discharge[unsort_rivid_indices]

            self.discharge = extracted_discharge

    def calculate_mean_input_discharge(self):
        """Calculate mean discharge for a list of netCDF files. It is assumed
        that the discharge variable array in each file has the same
        (time, rivid) coordinates.
        """
        try:
            nfile = len(self.input_file_list)
        except TypeError:
            nfile = None

        if nfile is not None:
            self.parse_input_discharge_file(
                input_discharge_file=self.input_file_list[0])

            input_discharge = self.input_discharge

            for filename in self.input_file_list[1:]:
                self.parse_input_discharge_file(
                    input_discharge_file=filename,
                    discharge_only=True)
                input_discharge += self.input_discharge

            self.input_discharge = input_discharge / nfile

    def integrate_over_files(self, integration_type=None):
        """Integrate discharge values over a list of files.

        Parameters
        ----------
        integration_type : str (optional)
            Identifier for integration method to be used.
        """
        if integration_type is not None:
            self.integration_type = integration_type

        if self.integration_type.lower() == 'mean':
            self.calculate_mean_input_discharge()
        else:
            warnings.warn(f'Integration type {self.integration_type} not ' +
                          'recognized.')

    def main(self, output_filename=None):
        """Main method for the RAPIDInputDischarge class. Writes a netCDF file
        containing discharge and accompanying spatial variables for a single
        time based on discharge and spatial information provided to the class
        constructor.

        Parameters
        ----------
        output_filename : str (optional)
            Name of output file.
        """
        if output_filename is not None:
            self.output_filename = output_filename

        try:
            nfile = len(self.input_file_list)
        except TypeError:
            nfile = None

        if nfile is not None:
            self.integrate_over_files()
        elif self.input_discharge_file is not None:
            self.parse_input_discharge_file()

        if self.rivid_file is not None:
            self.parse_rivid_file()

        if self.rivid is None:
            self.rivid = self.input_discharge_rivid

        self.sort_discharge_by_rivid()

        self.write_nc()
