"""Tools for writing various RAPID input CSV files."""

from csv import writer as csv_writer
import numpy as np
import geopandas as gpd


def write_connectivity_file(flowline_file,
                            flowline_id_field_name,
                            downstream_id_field_name,
                            out_csv_file):

    """
    Write a RAPID connectivity CSV.

    Parameters
    ----------
    flowline_file : str
        Path to vector flowline shapefile.
    flowline_id_field_name : str
        Name of the field containing the identifier of each flowline in the
        flowline file.
    downstream_id_field_name : str
        Name of the field containing the identifier of the reach downstream of
        each flowline.
    out_csv_file : str
        Path to output file to which the connectivity CSV will be written.
    """

    flowline_gdf = gpd.read_file(flowline_file)
    flowline_id_array = flowline_gdf[flowline_id_field_name].values
    downstream_id_array = flowline_gdf[downstream_id_field_name].values

    list_write = []
    max_count_upstream = 0

    for flowline_id in flowline_id_array:
        # For each flowline, get a list of all upstream identifiers
        list_upstream_id = flowline_id_array[
            downstream_id_array == flowline_id]

        # Count the number of the upstream identifiers
        count_upstream = len(list_upstream_id)
        if count_upstream > max_count_upstream:
            max_count_upstream = count_upstream

        # For each flowline, record the downstream identifier
        downstream_id = downstream_id_array[
            flowline_id_array == flowline_id][0]

        # Compile the flowline_id, downstream_id, count_upstream, and
        # list_upstream_id for writing
        list_write.append(
            np.concatenate([
                np.array([flowline_id, downstream_id,
                          count_upstream]), list_upstream_id]).astype(int))

    with open(out_csv_file, 'w', encoding='utf-8') as csv_file:
        connectivity_writer = csv_writer(csv_file)
        for row_list in list_write:
            out = np.concatenate([row_list, np.repeat(0, max_count_upstream -
                                                      row_list[2])])
            connectivity_writer.writerow(out.astype(int))


def write_riv_bas_id_file(flowline_file,
                          flowline_id_field_name,
                          out_csv_file):

    """
    Write a RAPID riv_bas_id_file. The order of the output file will match that
    of the input flowline shapefile.

    Parameters
    ----------
    flowline_file : str
        Path to vector flowline shapefile.
    flowline_id_field_name : str
        Name of the field containing the identifier of each flowline in the
        flowline file.
    out_csv_file : str
        Path to output file to which the riv_bas_id_file will be written.
    """

    flowline_gdf = gpd.read_file(flowline_file)
    flowline_id_array = flowline_gdf[flowline_id_field_name].values

    np.savetxt(out_csv_file, flowline_id_array, fmt='%i')


def write_kfac_file(flowline_file,
                    connectivity_file,
                    out_csv_file,
                    flowline_id_field_name,
                    length_field_name,
                    slope_field_name,
                    formula_type,
                    input_length_units='m',
                    input_slope_percent=False):

    """
    Write a Muskingum kfac file containing first guesses (in seconds) of the
    Muskingum k parameter. Three formula types are available, corresponing to
    equations (5)–(7) in Tavakoly et al. 2016
    (https://doi.org/10.1111/1752-1688.12456).

    formula_type = 1 -> Tavakoly et al. 2016 Eq. (5)
    formula_type = 2 -> Tavakoly et al. 2016 Eq. (6)
    formula_type = 3 -> Tavakoly et al. 2016 Eq. (7)

    Parameters
    ----------
    flowline_file : str
        Path to vector flowline shapefile.
    connectivity_file : str
        Path to a RAPID connectivity CSV.
    out_csv_file : str
        Path to output file to which the kfac_file will be written.
    flowline_id_field_name : str
        Name of the field containing the identifier of each flowline in the
        flowline file.
    length_id_field_name : str
        Name of the field containing the length of each flowline (in m or km)
        in the flowline file.
    slope_id_field_name : str
        Name of the field containing the slope of each flowline in the
        flowline file. The slope can be unitless (i.e., m/m or km/km) or
        represented as a percentage (i.e., m/m * 100 or km/km * 100).
    formula_type : int
        The formula from Tavakoly et al. 2016 to be used.
    input_length_units : str, optional
        The units of the input length. Supported units are 'm' for meters and
        'km' for kilometers. Default is 'm' for meters.
    input_slope_percent : bool, optional
        If True, the slope provided on the flowline file is assumed to be a
        percentage and will be divided by 100. Default is False.
    """

    flowline_gdf = gpd.read_file(flowline_file)
    flowline_id_array = flowline_gdf[flowline_id_field_name].values
    flowline_id_list = list(flowline_id_array)

    connect_rivid_array = np.genfromtxt(connectivity_file, delimiter=',',
                                        usecols=0, dtype=int)

    # The kfac file must be ordered the same as the connectivity file
    # Rewrite the length and slope arrays to ensure this requirement is met
    length_array = flowline_gdf[length_field_name].values
    slope_array = flowline_gdf[slope_field_name].values

    sort_idx = []
    for connect_rivid in connect_rivid_array:
        sort_idx.append(flowline_id_list.index(connect_rivid))

    length_array = length_array[sort_idx]
    slope_array = slope_array[sort_idx]

    if input_length_units == 'km':
        length_array *= 1000.
    if input_slope_percent:
        slope_array /= 100.

    # For Kini1, divide each river length Li by the reference wave celerity
    # C0 of 1 km per hour (1000 m per 3600 s):
    if formula_type == 1:
        # Length in m divided by (1000 m per 3600 s)
        kfac_array = length_array / (1000. / 3600.)

    # For Kini2:
    if formula_type == 2:
        k_ini_1_array = length_array / (1000. / 3600.)
        length_slope_array = length_array / (slope_array**0.5)
        eta_array = np.mean(k_ini_1_array) / np.mean(length_slope_array)
        kfac_array = eta_array*length_slope_array

    # For Kini3:
    if formula_type == 3:
        k_ini_1_array = length_array / (1000. / 3600.)
        length_slope_array = length_array/(slope_array**0.5)
        percentile_5 = np.percentile(length_slope_array, 5)
        percentile_95 = np.percentile(length_slope_array, 95)
        length_slope_array[length_slope_array < percentile_5] = percentile_5
        length_slope_array[length_slope_array > percentile_95] = percentile_95
        eta_array = np.mean(k_ini_1_array)/np.mean(length_slope_array)
        kfac_array = eta_array*length_slope_array

    np.savetxt(out_csv_file, kfac_array, fmt='%.1f')


def write_constant_x_file(connectivity_file,
                          out_csv_file,
                          x_value):

    """
    Write a constant Muskingum x file. The Muskingum x parameter will be the
    same for all river reaches).

    Parameters
    ----------
    connectivity_file : str
        Path to a RAPID connectivity CSV.
    out_csv_file : str
        Path to output file to which the constant x file will be written.
    x_value : float
        The constant x value to be used for all rivier reaches (0 <= x <= 0.5)
    """

    connect_rivid_array = np.genfromtxt(connectivity_file, delimiter=',',
                                        usecols=0, dtype=int)
    np.savetxt(out_csv_file, np.repeat(x_value, len(connect_rivid_array)),
               fmt='%.1f')
