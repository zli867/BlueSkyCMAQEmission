import netCDF4 as nc
import pyproj
import numpy as np
import json
from datetime import datetime, timedelta

"""
The main function should be used are verticalHourlyEmission and dailyVerticalHourlyEmission
verticalHourlyEmission is used for prescribed burning which lasts for one day
dailyVerticalHourlyEmission is used for wildfires which lasts for multiple day or multiple month
dailyVerticalHourlyEmission should be a update of verticalHourlyEmission function since the memory requirement does not
increase as the fire duration increase
"""


def datetime_range(start, end, delta):
    time_array = []
    current = start
    time_array.append(current)
    while current < end:
        current += delta
        time_array.append(current)
    return time_array


def CMAQGrid3D(metcros3d):
    """

    :param metcros3d: METCROS3D filename
    :return: a dictionary which includes CMAQ grid definitions (X, Y, Z) and projection definition
    """
    met_ds = nc.Dataset(metcros3d)
    cmaq_height = met_ds['ZF'][:]
    time_data = met_ds['TFLAG'][:]
    lat_1 = met_ds.getncattr('P_ALP')
    lat_2 = met_ds.getncattr('P_BET')
    lat_0 = met_ds.getncattr('YCENT')
    lon_0 = met_ds.getncattr('XCENT')
    crs = pyproj.Proj("+proj=lcc +a=6370000.0 +b=6370000.0 +lat_1=" + str(lat_1)
                      + " +lat_2=" + str(lat_2) + " +lat_0=" + str(lat_0) +
                      " +lon_0=" + str(lon_0))
    xcell = met_ds.getncattr('XCELL')
    ycell = met_ds.getncattr('YCELL')
    xorig = met_ds.getncattr('XORIG')
    yorig = met_ds.getncattr('YORIG')

    ncols = met_ds.getncattr('NCOLS')
    nrows = met_ds.getncattr('NROWS')

    #> for X, Y cell centers
    x_center_range = np.linspace(xorig+xcell/2, (xorig+xcell/2)+xcell*(ncols-1), ncols)
    y_center_range = np.linspace(yorig+ycell/2, (yorig+ycell/2)+ycell*(nrows-1), nrows)

    Xcenters, Ycenters = np.meshgrid(x_center_range, y_center_range)

    #> for X, Y cell boundaries (i.e., cell corners)
    x_bound_range = np.linspace(xorig, xorig+xcell*ncols, ncols+1)
    y_bound_range = np.linspace(yorig, yorig+ycell*nrows, nrows+1)

    Xbounds, Ybounds = np.meshgrid(x_bound_range, y_bound_range)

    x_max = np.max(Xbounds)
    x_min = np.min(Xbounds)
    y_max = np.max(Ybounds)
    y_min = np.min(Ybounds)

    cmaq_time_array = []
    for i in range(0, time_data.shape[0]):
        time_data_tmp = time_data[i, 0, :]
        time_str = str(time_data_tmp[0]) + str(time_data_tmp[1]).rjust(6, '0')
        parsed = datetime.strptime(time_str, '%Y%j%H%M%S')
        cmaq_time_array.append(parsed)

    res_dict = {"crs": crs, "X_ctr": Xcenters, "Y_ctr": Ycenters, "X_bdry": [x_min, x_max], "Y_bdry": [y_min, y_max],
                "height": cmaq_height, "time": cmaq_time_array}
    return res_dict


def verticalHourlyFraction(timeprofile, plumerise, utc_offset, mcip_time, mcip_height):
    negative = False
    if utc_offset[0] == "-":
        delta = datetime.strptime(utc_offset[1:], "%H:%M") - datetime.strptime("00:00", "%H:%M")
        negative = True
    else:
        delta = datetime.strptime(utc_offset, "%H:%M") - datetime.strptime("00:00", "%H:%M")

    # temporal fraction * height fraction
    # generate hourly fraction
    lay = mcip_height.shape[1]
    flaming_frac = np.zeros((len(mcip_time), lay))
    smoldering_frac = np.zeros((len(mcip_time), lay))
    residual_farc = np.zeros((len(mcip_time), lay))
    for key in timeprofile.keys():
        datetime_object = datetime.strptime(key, '%Y-%m-%dT%H:%M:%S')
        # local time to utc
        if negative:
            utc_date = datetime_object + delta
        else:
            utc_date = datetime_object - delta
        time_idx = mcip_time.index(utc_date)
        # calculate vertical fraction
        fractions = np.zeros(lay)
        vertical_heights_bsp = plumerise[key]["heights"]
        max_plumerise = np.max(vertical_heights_bsp)
        min_plumerise = np.min(vertical_heights_bsp)
        mcip_height_tmp = mcip_height[time_idx, :]
        for i in range(0, len(mcip_height_tmp)):
            if i == 0:
                bottom = 0
                top = mcip_height_tmp[i]
            else:
                bottom = mcip_height_tmp[i - 1]
                top = mcip_height_tmp[i]
            if bottom >= min_plumerise and top <= max_plumerise:
                frac_tmp = (top - bottom) / (max_plumerise - min_plumerise)
            elif bottom <= min_plumerise <= top:
                if top <= max_plumerise:
                    frac_tmp = (top - min_plumerise) / (max_plumerise - min_plumerise)
                else:
                    frac_tmp = 1
            elif bottom <= max_plumerise <= top:
                frac_tmp = (max_plumerise - bottom) / (max_plumerise - min_plumerise)
            else:
                frac_tmp = 0
            fractions[i] = frac_tmp
        flaming_frac[time_idx, :] = timeprofile[key]['flaming'] * fractions
        smoldering_frac[time_idx, :] = timeprofile[key]['smoldering'] * fractions
        residual_farc[time_idx, :] = timeprofile[key]['residual'] * fractions
    return {"flaming_hourly_frac": flaming_frac,
            "smoldering_hourly_frac": smoldering_frac,
            "residual_hourly_farc": residual_farc}


def classifiedEmission(fuelbeds, select_species):
    """

    :param fuelbeds: fuelbeds attribute from bluesky output
    :param select_species: selected chemical emission species
    :return: a dictionary includes three phases the emissions (flaming, smoldering and residual)
            of selected chemical species
    """
    # initialize flaming residual and smoldering emissions
    flaming = {}
    for select_specie in select_species:
        flaming[select_specie] = 0
    smoldering = {}
    for select_specie in select_species:
        smoldering[select_specie] = 0
    residual = {}
    for select_specie in select_species:
        residual[select_specie] = 0

    for fuel_idx in range(0, len(fuelbeds)):
        fuelbed_flaming = fuelbeds[fuel_idx]["emissions"]['flaming']
        fuelbed_smoldering = fuelbeds[fuel_idx]["emissions"]['smoldering']
        fuelbed_residual = fuelbeds[fuel_idx]["emissions"]['residual']
        for specie in select_species:
            if specie in fuelbed_flaming.keys():
                flaming[specie] += fuelbed_flaming[specie][0]
            if specie in fuelbed_smoldering.keys():
                smoldering[specie] += fuelbed_smoldering[specie][0]
            if specie in fuelbed_residual.keys():
                residual[specie] += fuelbed_residual[specie][0]

    # generate flaming smoldering residual emission array
    flaming_arry = []
    smoldering_arry = []
    residual_arry = []
    for select_specie in select_species:
        flaming_arry.append(flaming[select_specie])
        smoldering_arry.append(smoldering[select_specie])
        residual_arry.append(residual[select_specie])
    flaming_arry = np.array(flaming_arry)
    smoldering_arry = np.array(smoldering_arry)
    residual_arry = np.array(residual_arry)
    return {"flaming_emission": flaming_arry, "smoldering_emission": smoldering_arry, "residual_emission": residual_arry}


def findSpatialIndex(fire_x, fire_y, X_ctr, Y_ctr):
    """

    :param fire_x: X of fire location in CMAQ projection
    :param fire_y: Y of fire location in CMAQ projection
    :param X_ctr: CMAQ grid X center
    :param Y_ctr: CMAQ grid Y center
    :return: x_idx, y_idx which are the fire location in CMAQ grid
    """
    dist = np.sqrt((X_ctr - fire_x) ** 2 + (Y_ctr - fire_y) ** 2)
    x_idx, y_idx = np.unravel_index(np.argmin(dist, axis=None), dist.shape)
    return x_idx, y_idx


def verticalHourlyEmission(select_species, metcros3d_files, bsp_filename):
    """
    :param select_species:  a list of string, select species for chemistry mapping usage from BlueSky side
    :param metcros3d_files: MCIP 3D Grid Information
    :param bsp_filename: BlueSKY output
    :return: a emission tensor (species, MCIP time length, LAY, MCIP X length, MCIP Y length)
    """
    # Unit tons/hour
    cmaq_3d_info = CMAQGrid3D(metcros3d_files)
    # CMAQ grid definition
    crs = cmaq_3d_info["crs"]
    x_min = cmaq_3d_info["X_bdry"][0]
    x_max = cmaq_3d_info["X_bdry"][1]
    y_min = cmaq_3d_info["Y_bdry"][0]
    y_max = cmaq_3d_info["Y_bdry"][1]
    Xcenters = cmaq_3d_info["X_ctr"]
    Ycenters = cmaq_3d_info["Y_ctr"]
    mcip_time = cmaq_3d_info["time"]
    mcip_height = cmaq_3d_info["height"]
    LAY = mcip_height.shape[1]
    # read fire data
    with open(bsp_filename) as jsfile:
        bluesky_data = json.load(jsfile)

    bsp_fires = bluesky_data["fires"]
    # species, time, layer, x, y
    emission_tensor = np.zeros((len(select_species), len(mcip_time), LAY, Xcenters.shape[0], Xcenters.shape[1]))
    for bsp_fire in bsp_fires:
        lat = bsp_fire["activity"][0]['active_areas'][0]["specified_points"][0]["lat"]
        lon = bsp_fire["activity"][0]['active_areas'][0]["specified_points"][0]["lng"]
        fire_id = bsp_fire["id"]
        time_profile = bsp_fire["activity"][0]["active_areas"][0]["timeprofile"]
        area_acres = bsp_fire["activity"][0]['active_areas'][0]["specified_points"][0]["area"]
        fuelbeds = bsp_fire["activity"][0]['active_areas'][0]["specified_points"][0]["fuelbeds"]
        emission = bsp_fire["activity"][0]['active_areas'][0]["specified_points"][0]["emissions"]["summary"]
        utc_offset = bsp_fire["activity"][0]['active_areas'][0]["utc_offset"]
        plumerise = bsp_fire["activity"][0]['active_areas'][0]["specified_points"][0]["plumerise"]

        # whether in grid boundary
        fire_x, fire_y = crs(lon, lat)
        if x_min <= fire_x <= x_max and y_min <= fire_y <= y_max:
            # Find the index of fire point in grid
            x_idx, y_idx = findSpatialIndex(fire_x, fire_y, Xcenters, Ycenters)
            mcip_loc_height = mcip_height[:, :, x_idx, y_idx]
            # get fraction of each phase, fraction definition : temporal fraction * height fraction
            emission_frac = verticalHourlyFraction(time_profile, plumerise, utc_offset, mcip_time, mcip_loc_height)
            flaming_hourly_frac = emission_frac["flaming_hourly_frac"]
            smoldering_hourly_frac = emission_frac["smoldering_hourly_frac"]
            residual_hourly_farc = emission_frac["residual_hourly_farc"]
            t_length, h_length = flaming_hourly_frac.shape
            # generate flaming smoldering residual emission array, the array is ordered by species array
            classified_emission = classifiedEmission(fuelbeds, select_species)
            flaming_arry = classified_emission["flaming_emission"]
            smoldering_arry = classified_emission["smoldering_emission"]
            residual_arry = classified_emission["residual_emission"]
            # generate emission matrix (species, time_period, layer)
            emission_matrix = np.zeros((len(select_species), t_length, h_length))
            for i in range(0, len(select_species)):
                emission_matrix[i, :, :] = flaming_hourly_frac * flaming_arry[i] + \
                                        smoldering_hourly_frac * smoldering_arry[i] + \
                                        residual_hourly_farc * residual_arry[i]
            emission_tensor[:, :, :, x_idx, y_idx] += emission_matrix
        else:
            continue
    return emission_tensor


def daily_profile(timeprofile, plumerise, select_utc_date, utc_offset):
    new_timeprofile = {}
    new_plumerise = {}
    daily_mcip_start = datetime(select_utc_date.year, select_utc_date.month, select_utc_date.day, 0)
    daily_mcip_end = datetime(select_utc_date.year, select_utc_date.month, select_utc_date.day, 0) + timedelta(days=1)
    negative = False
    if utc_offset[0] == "-":
        delta = datetime.strptime(utc_offset[1:], "%H:%M") - datetime.strptime("00:00", "%H:%M")
        negative = True
    else:
        delta = datetime.strptime(utc_offset, "%H:%M") - datetime.strptime("00:00", "%H:%M")
    for key in timeprofile.keys():
        datetime_object = datetime.strptime(key, '%Y-%m-%dT%H:%M:%S')
        # local time to utc
        if negative:
            utc_date = datetime_object + delta
        else:
            utc_date = datetime_object - delta
        if utc_date >= daily_mcip_start and utc_date <= daily_mcip_end:
            new_timeprofile[key] = timeprofile[key]
            new_plumerise[key] = plumerise[key]
    return new_timeprofile, new_plumerise


def dailyVerticalHourlyEmission(select_species, metcros3d_file, bsp_filename):
    """
    Get emission tensor fo specific day (metcros3d provided)
    :param select_species: a list of string, select species for chemistry mapping usage from BlueSky side
    :param metcros3d_file: MCIP 3D Grid Information
    :param bsp_filename: BlueSKY output
    :return: a emission tensor (species, MCIP time length, LAY, MCIP X length, MCIP Y length)
    """
    cmaq_3d_info = CMAQGrid3D(metcros3d_file)
    # CMAQ grid definition
    crs = cmaq_3d_info["crs"]
    x_min = cmaq_3d_info["X_bdry"][0]
    x_max = cmaq_3d_info["X_bdry"][1]
    y_min = cmaq_3d_info["Y_bdry"][0]
    y_max = cmaq_3d_info["Y_bdry"][1]
    Xcenters = cmaq_3d_info["X_ctr"]
    Ycenters = cmaq_3d_info["Y_ctr"]
    mcip_time = cmaq_3d_info["time"]
    mcip_height = cmaq_3d_info["height"]
    LAY = mcip_height.shape[1]
    current_select_time = mcip_time[0]
    # read fire data
    with open(bsp_filename) as jsfile:
        bluesky_data = json.load(jsfile)

    bsp_fires = bluesky_data["fires"]
    # species, time, layer, x, y
    emission_tensor = np.zeros((len(select_species), len(mcip_time), LAY, Xcenters.shape[0], Xcenters.shape[1]))
    for bsp_fire in bsp_fires:
        lat = bsp_fire["activity"][0]['active_areas'][0]["specified_points"][0]["lat"]
        lon = bsp_fire["activity"][0]['active_areas'][0]["specified_points"][0]["lng"]
        fire_id = bsp_fire["id"]
        time_profile = bsp_fire["activity"][0]["active_areas"][0]["timeprofile"]
        area_acres = bsp_fire["activity"][0]['active_areas'][0]["specified_points"][0]["area"]
        fuelbeds = bsp_fire["activity"][0]['active_areas'][0]["specified_points"][0]["fuelbeds"]
        emission = bsp_fire["activity"][0]['active_areas'][0]["specified_points"][0]["emissions"]["summary"]
        utc_offset = bsp_fire["activity"][0]['active_areas'][0]["utc_offset"]
        plumerise = bsp_fire["activity"][0]['active_areas'][0]["specified_points"][0]["plumerise"]

        # update the profile to daily profile
        time_profile, plumerise = daily_profile(time_profile, plumerise, current_select_time, utc_offset)

        # whether in grid boundary
        fire_x, fire_y = crs(lon, lat)
        if x_min <= fire_x <= x_max and y_min <= fire_y <= y_max:
            # Find the index of fire point in grid
            x_idx, y_idx = findSpatialIndex(fire_x, fire_y, Xcenters, Ycenters)
            mcip_loc_height = mcip_height[:, :, x_idx, y_idx]
            # get fraction of each phase, fraction definition : temporal fraction * height fraction
            emission_frac = verticalHourlyFraction(time_profile, plumerise, utc_offset, mcip_time, mcip_loc_height)
            flaming_hourly_frac = emission_frac["flaming_hourly_frac"]
            smoldering_hourly_frac = emission_frac["smoldering_hourly_frac"]
            residual_hourly_farc = emission_frac["residual_hourly_farc"]
            t_length, h_length = flaming_hourly_frac.shape
            # generate flaming smoldering residual emission array, the array is ordered by species array
            classified_emission = classifiedEmission(fuelbeds, select_species)
            flaming_arry = classified_emission["flaming_emission"]
            smoldering_arry = classified_emission["smoldering_emission"]
            residual_arry = classified_emission["residual_emission"]
            # generate emission matrix (species, time_period, layer)
            emission_matrix = np.zeros((len(select_species), t_length, h_length))
            for i in range(0, len(select_species)):
                emission_matrix[i, :, :] = flaming_hourly_frac * flaming_arry[i] + \
                                        smoldering_hourly_frac * smoldering_arry[i] + \
                                        residual_hourly_farc * residual_arry[i]
            emission_tensor[:, :, :, x_idx, y_idx] += emission_matrix
        else:
            continue
    return emission_tensor


def get_fire_time(time, utc_offset):
    """

    :param time: a str object of local time
    :param utc_offset: a str object of UTC offset
    :return: a datetime object of UTC time
    """
    negative = False
    if utc_offset[0] == "-":
        delta = datetime.strptime(utc_offset[1:], "%H:%M") - datetime.strptime("00:00", "%H:%M")
        negative = True
    else:
        delta = datetime.strptime(utc_offset, "%H:%M") - datetime.strptime("00:00", "%H:%M")
    datetime_object = datetime.strptime(time, '%Y-%m-%dT%H:%M:%S')
    # local time to utc
    if negative:
        utc_date = datetime_object + delta
    else:
        utc_date = datetime_object - delta
    return utc_date
