"""Tools to write input input discharge files for RAPID."""

import numpy as np
from netCDF4 import Dataset
from datetime import datetime

class RAPIDInputDischarge:

    """A class for writing input discharge netCDF files for the 
    Routing Application for Parallel Computation of Discharge (RAPID) model.
    """

    def __init__(self,
                 output_filename='rapid_input_discharge.nc',
                 discharge=None,
                 time=None,
                 rivid=None,
                 latitude=None, longitude=None,
                 discharge_datatype='f8', rivid_datatype='i4',
                 time_datatype='i8', latlon_datatype='f8',
                 crs_datatype='i4',
                 time_units='seconds since 1970-01-01 00:00:00',
                 discharge_variable_name='Qout',
                 time_variable_name='time',
                 rivid_variable_name='rivid',
                 latitude_variable_name='lat',
                 longitude_variable_name='longitude'):
        """
        Parameters
        ----------
        output_filename : str
            Name of output file.
        discharge : array_like
            Discharge values.
        time : array_like
            Integer times corresponding to discharge.
        rivid : array_like
            Identifiers for river segments.
        latitude : array_like
            Latitude coordinates.
        longitude : array_like
            Longitude coordinates.
        discharge_datatype : str
            Data type to be used for writing discharge variable.
        rivid_datatype : str
            Data type to be used for writing rivid variable.
        time_datatype : str
            Data type to be used for writing time variable.
        latlon_datatype : str
            Data type to be used for writing lat/lon variables.
        crs_datatype : str
            Data type to be used for writing coordinate reference variable.
        time_units : str
            Units and datum for time' variable.
        discharge_variable_name : str
            Name of discharge variable to be used in output file.
        time_variable_name : str
            Name of time variable to be used in output file.
        rivid_variable_name : str
            Name of rivid variable to be used in output file.
        latitude_variable_name : str
            Name of latitude variable to be used in output file.
        longitude_variable_name : str
            Name of longitude variable to be used in output file.
        """
        self.output_filename = output_filename
        self.discharge = discharge
        self.time = time
        self.rivid = rivid
        self.latitude = latitude
        self.longitude = longitude
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

    def write_nc(self, discharge=None, time=None, rivid=None, latitude=None,
                 longitude=None):
        """
        Write variables, dimensions, and attributes to output file.
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

        if self.time is not None:
            nt = len(self.time)
        else:
            nt = None

        if self.rivid is not None:
            nr = len(self.rivid)
        else:
            nr = None

        data_out_nc = Dataset(self.output_filename, 'w')

        # create dimensions
        data_out_nc.createDimension(self.time_variable_name, nt)
        data_out_nc.createDimension(self.rivid_variable_name, nr)

        # create variables
        # discharge
        discharge_var = data_out_nc.createVariable(
            self.discharge_variable_name, self.discharge_datatype,
            (self.time_variable_name, self.rivid_variable_name),
            fill_value=0)
        discharge_var.long_name = ('instantaneous river water discharge ' +
                                   'downstream of each river reach')
        discharge_var.units = 'm3 s-1'
        discharge_var.coordinates = 'lon lat'
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

        # time
        time_var = data_out_nc.createVariable(
            self.time_variable_name, self.time_datatype, ('time',))
        time_var.long_name = self.time_variable_name,
        time_var.standard_name = 'time'
        time_var.units = self.time_units
        time_var.axis = 'T'
        time_var.calendar = 'gregorian'
        time_var.bounds = 'time_bnds'
        if self.time is not None:
            time_var[:] = self.time

        # longitude
        lon_var = data_out_nc.createVariable('lon', 'f8', ('rivid',),
                                             fill_value=-9999.0)
        lon_var.long_name = \
            'longitude of a point related to each river reach'
        lon_var.standard_name = 'longitude'
        lon_var.units = 'degrees_east'
        lon_var.axis = 'X'
        if self.longitude is not None:
            lon_var[:] = self.longitude

        # latitude
        lat_var = data_out_nc.createVariable('lat', 'f8', ('rivid',),
                                             fill_value=-9999.0)
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

        # close file
        data_out_nc.close()

class InitialFlows(RAPIDInputDischarge):

    def __init__(self,
                 output_filename='qinit.nc'
                 input_discharge_file=None,
                 connectivity_file=None, 
                 input_discharge_time_axis=0,
                 input_discharge_time_index=-1,
                 discharge_datatype='f8',
                 rivid_datatype='i4',
                 time_datatype='i8',
                 latlon_datatype='f8',
                 crs_datatype='i4',
                 time_units='seconds since 1970-01-01 00:00:00',
                 discharge_variable_name='Qout',
                 rivid_variable_name='rivid',
                 latitude_variable_name='lat',
                 longitude_variable_name='lon',
                 time_variable_name='time'):

        self.input_discharge_file = input_discharge_file
        self.connectivity_file = connectivity_file
        self.input_discharge_time_axis = input_discharge_time_axis
        self.input_discharge_time_index = input_discharge_time_index

        # Attributes to be assigned.
        self.input_discharge = None
        self.connectivity = None
        self.input_discharge_rivid = None
        
        super().__init__(output_filename=output_filename,
                         discharge_datatype=discharge_datatype,
                         rivid_datatype=rivid_datatype,
                         time_datatype=time_datatype,
                         latlon_datatype=latlon_datatype
                         crs_datatype=crs_datatype)

    def parse_input_discharge_file(self, input_discharge_file=None):

        if input_discharge_file is not None:
            self.input_discharge_file = input_discharge_file

        if self.input_discharge_time_axis == 0:
            input_discharge_slice = (
                slice(self.input_discharge_time_index,None, None),
                slice(None, None, None))
        elif self.input_discharge_time_axis == 1:
            input_discharge_slice = (
                slice(None, None, None),
                slice(self.input_discharge_time_index, None, None))

        input_time_var = None
        
        try:
            d = Dataset(self.input_discharge_file)
        except IOError:
            warnings.warn(
                f'Unable to open file {self.input_discharge_file}.')
            d = None
        
        if d is not None:
            try:
                self.input_discharge = d[self.discharge_variable_name][
                    input_discharge_slice]
            except (IOError, IndexError):
                warnings.warn(
                    f'Variable "{self.discharge_variable_name}" not found ' +
                    f'in file {self.input_discharge_file}.')

            try:
                self.input_discharge_rivid = d['rivid'][:]
            except (IOError, IndexError):
                warnings.warn(
                    f'Variable "rivid" not found ' +
                    f'in file {self.input_discharge_file}.')

            try:
                time_var = d['time']
            except IOError:
                warnings.warn(
                    f'Variable "time" not found ' +
                    f'in file {self.input_discharge_file}.')
                
        if time_var is not None:
            try:
                time_step = time_var[1] - time_var[0]
            except IndexError:
                time_step = None

            # Increment final simulation time t_n by time step to give value
            # t_(n+1). t_(n+1) is assumed to be the first time step when
            # initializing a new run with discharge from t_n.
            final_time = time_var[-1]

            try:
                next_time = final_time + time_step
            except:
                next_time = np.array([])

            self.time = next_time

        if self.time is None:
            self.time = np.array([])

    def parse_connectivity_file(self, connectivity_file=None, delimiter=','):

        if connectivity_file is not None:
            self.connectivity_file = connectivity_file

        try:
            self.connectivity = np.genfromtxt(
                self.connectivity_file, delimiter=delimiter, dtype=int)
        except:
            warnings.warn(f'Unable to read file {self.connectivity_file}.')
            self.connectivity = None

        if self.connectivity is not None:
            self.rivid = self.connectivity[:,0]

    def sort_discharge_by_rivid(self, input_discharge_file=None,
                                connectivity_file=None):

        if input_discharge_file is not None:
            self.input_discharge_file = input_discharge_file
        if connectivity_file is not None:
            self.connectivity_file = connectivity_file

        if self.input_discharge is None:
            self.parse_input_discharge_file()

        if self.connectivity is None:
            self.parse_connectivity_file()
        
        if np.array_equal(self.rivid, self.input_discharge_rivid):
            d = Dataset(self.output_filename, 'a')
            d[self.discharge_variable_name][:] = self.input_discharge
            d[self.time_variable_name][:] = self.time

if __name__ == '__main__':
   
