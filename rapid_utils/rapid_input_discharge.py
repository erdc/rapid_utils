"""Tools to write input input discharge files for RAPID."""

import numpy as np
from netCDF4 import Dataset
from datetime import datetime

def read_connectivity(self):
    connectivity = np.genfromtxt(self.connectivity_file, delimiter=',')

    return connectivity

class RAPIDInputDischargFile:

    def __init__(self, output_filename, rivid, time=None,
                 longitude=None, latitude=None,
                 discharge_datatype='f8', rivid_datatype='i4',
                 time_datatype='i8', latlon_datatype='f8',
                 crs_datatype='i4',
                 time_units='seconds since 1970-01-01 00:00:00',
                 discharge_variable_name='Qout'):

        self.output_filename = output_filename
        self.rivid = rivid
        self.time = time
        self.longitude = longitude
        self.latitude = latitude
        self.discharge_datatype = discharge_datatype
        self.rivid_datatype = rivid_datatype
        self.time_datatype = time_datatype
        self.latlon_datatype = latlon_datatype
        self.crs_datatype = crs_datatype
        self.time_units = time_units
        

    def initialize_nc(self):
        """
        Write variables, dimensions, and attributes to output file.
        """

        data_out_nc = Dataset(self.output_filename, 'w')

        # create dimensions
        data_out_nc.createDimension('time', len(self.time))
        data_out_nc.createDimension('rivid', len(self.rivid))

        # create variables
        # discharge
        discharge_var = data_out_nc.createVariable(self.discharge_variable_name,
                                                   self.discharge_datatype,
                                                   ('time', 'rivid'),
                                                   fill_value=0)
        discharge_var.long_name = ('instantaneous river water discharge ' +
                                   'downstream of each river reach')
        discharge_var.units = 'm3 s-1'
        discharge_var.coordinates = 'lon lat'
        discharge_var.grid_mapping = 'crs'
        discharge_var.cell_methods = "time: point"

        # rivid
        rivid_var = data_out_nc.createVariable('rivid', self.rivid_datatype,
                                               ('rivid',))
        rivid_var.long_name = 'unique identifier for each river reach'
        rivid_var.units = '1'
        rivid_var.cf_role = 'timeseries_id'
        rivid_var[:] = self.rivid

        # time
        time_var = data_out_nc.createVariable('time', self.time_datatype,
                                              ('time',))
        time_var.long_name = 'time'
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

        # data_out_nc.featureType = 'timeSeries'
        #data_out_nc.institution = modeling_institution

        # close file
        data_out_nc.close()


class InitialFlowsFile(RAPIDInputDischargFile):

    def __init__(self, input_discharge_file, connectivity_file, 
                 output_filename, rivid, time=None,
                 longitude=None, latitude=None,
                 discharge_datatype='f8', rivid_datatype='i4',
                 time_datatype='i8', latlon_datatype='f8',
                 crs_datatype='i4',
                 time_units='seconds since 1970-01-01 00:00:00',
                 discharge_variable_name='Qout'):

        self.input_discharge_file = input_discharge_file
        self.connectivity_file = connectivity_file
        self.discharge_variable_name = discharge_variable_name

        # Attributes to be assigned.
        self.discharge = None
        self.connectivity = None
        
        super().__init__(output_filename, rivid, time, longitude, latitude,
                         discharge_datatype='f8', rivid_datatype='i4',
                         time_datatype='i8', latlon_datatype='f8',
                         crs_datatype='i4')

    def parse_input_discharge_file(self):
                                   
        d = Dataset(self.input_discharge_file)

        self.input_discharge = d[self.discharge_variable_name][-1,:]
        self.input_discharge_rivid = d['rivid'][:]

    def parse_connectivity_file(self, delimiter=','):

        self.connectivity = np.genfromtxt(self.connectivity_file,
                                          delimiter=delimiter)

    def compute_inital_flows(self):

        if self.input_discharge is None:
            self.parse_input_discharge_file

        if self.connectivity is None:
            self.parse_connectivity_file

        upstream_start_index = 3
        connect_rivid = self.connectivity[:,0]
        nrivid = len(connect_rivid)
        
        initial_flows = np.zeros(nrivid)

        if np.array_equal(connect_rivid, self.input_discharge_rivid):
            for idx, row in enumerate(self.connectivity):
                n_upstream = row[2]
                upstream_rivid = row[
                    upstream_start_index:upstream_start_index + n_upstream]
                discharge_indices = np.where(
                    np.isin(self.input_discharge_rivid, upstream_rivid))
                discharge[idx] = self.discharge[discharge_indices].sum()
        
if __name__ == '__main__':
    fname = 'oonygoogoo.nc'
    rivid = [1,2,3]
    time = [1]
    initial_flow_file = InitialFlowsFile(fname, rivid, time)
    initial_flow_file.initialize_nc()
    print(initial_flow_file.__dict__)
    print('ya bit')
