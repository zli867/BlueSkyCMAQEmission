from EmissionGenerator import dailyVerticalHourlyEmission, get_fire_time
import netCDF4
import numpy as np
import os
from datetime import datetime, timedelta
import json

bsp_filename = "./input_dir/SE2016_01_30_out.json"
input_dir = "./input_dir/"
work_dir = "./work_dir/"
output_dir = "./results/"
# convert unit
convert_to_g_s = 907184.74 / 3600


def convert_to_mole_s(molecule_weight):
    # convert to g/s -> moles/s
    return convert_to_g_s / molecule_weight

# # SAPRC07
# species_mapping = {
#     'CO': [(1 * convert_to_mole_s(28), 'CO')],
#     # SERA
#     'NO2': [(0.35 * convert_to_mole_s(46), 'NOx')],
#     'NO': [(0.65 * convert_to_mole_s(46), 'NOx')],
#     'NH3': [(1 * convert_to_mole_s(17), 'NH3')],
#     'SO2': [(1 * convert_to_mole_s(64), 'SO2')],
#     'ACYE': [(0.084 * convert_to_mole_s(26.0373), 'VOC')],
#     'ALK1': [(0.1047 * convert_to_mole_s(30.069), 'VOC')],
#     'ALK2': [(0.0035 * convert_to_mole_s(44.0956), 'VOC')],
#     'ALK3': [(0.004108 * convert_to_mole_s(60.8495), 'VOC')],
#     'ALK4': [(0.0255 * convert_to_mole_s(101.5224), 'VOC')],
#     'ALK5': [(0.2443 * convert_to_mole_s(166.6683), 'VOC')],
#     'APIN': [(0.0143 * convert_to_mole_s(208.3353), 'VOC')],
#     'ARO1': [(0.0263 * convert_to_mole_s(166.6883), 'VOC')],
#     'ARO2': [(0.422 * convert_to_mole_s(187.5018), 'VOC')],
#     'B124': [(0.002927 * convert_to_mole_s(187.5018), 'VOC')],
#     'BDE13': [(0.0052 * convert_to_mole_s(54.0904), 'VOC')],
#     'CH4': [(0.0982 * convert_to_mole_s(16.0425), 'VOC')],
#     'CRES': [(0.002397 * convert_to_mole_s(145.8347), 'VOC')],
#     'ETHE': [(0.1911 * convert_to_mole_s(28.0532), 'VOC')],
#     'MVK': [(0.0003175 * convert_to_mole_s(83.3341), 'VOC')],
#     'NROG': [(0.0206 * convert_to_mole_s(41.6671), 'VOC')],
#     'OLE1': [(0.0451 * convert_to_mole_s(74.1464), 'VOC')],
#     'OLE2': [(0.0163 * convert_to_mole_s(83.8603), 'VOC')],
#     'PRD2': [(0.0174 * convert_to_mole_s(125.0012), 'VOC')],
#     'PRPE': [(0.0393 * convert_to_mole_s(42.0797), 'VOC')],
#     'TERP': [(0.0123 * convert_to_mole_s(208.3353), 'VOC')],
#     'PAL': [(4.6 * (10 ** (-4)) * convert_to_g_s, 'PM2.5')],
#     'PCA': [(7.2 * (10 ** (-4)) * convert_to_g_s, 'PM2.5')],
#     'PCL': [(0.00239 * convert_to_g_s, 'PM2.5')],
#     'PEC': [(0.109 * convert_to_g_s, 'PM2.5')],
#     'PFE': [(4.45 * (10 ** (-4)) * convert_to_g_s, 'PM2.5')],
#     'PK': [(0.00135 * convert_to_g_s, 'PM2.5')],
#     'PMFINE': [(0.375 * convert_to_g_s, 'PM2.5')],
#     'PMN': [(1.1 * (10 ** (-4)) * convert_to_g_s, 'PM2.5')],
#     'PMOTHR': [(0.012995 * convert_to_g_s, 'PM2.5')],
#     'PNA': [(0.00135 * convert_to_g_s, 'PM2.5')],
#     'PNCOM': [(0.351 * convert_to_g_s, 'PM2.5')],
#     'PNH4': [(0.00341 * convert_to_g_s, 'PM2.5')],
#     'PNO3': [(0.0107 * convert_to_g_s, 'PM2.5')],
#     'POC': [(0.502 * convert_to_g_s, 'PM2.5')],
#     'PSI': [(1.0 * (10 ** (-4)) * convert_to_g_s, 'PM2.5')],
#     'PSO4': [(0.0033 * convert_to_g_s, 'PM2.5')],
#     'PTI': [(6.7 * (10 ** (-4)) * convert_to_g_s, 'PM2.5')],
#     # SERA
#     'PMC': [(1.2480483 * convert_to_g_s, 'PM2.5'), (-1 * convert_to_g_s, 'PM2.5')]
# }

# CB6
species_mapping = {
    'CO': [(1 * convert_to_mole_s(28), 'CO')],
    'SO2': [(1 * convert_to_mole_s(64), 'SO2')],
    'NH3': [(1 * convert_to_mole_s(17), 'NH3')],
    # SERA
    'PMC': [(1.2480483 * convert_to_g_s, 'PM2.5'), (-1 * convert_to_g_s, 'PM2.5')],
    # SERA
    'NO2': [(0.35 * convert_to_mole_s(46), 'NOx')],
    # SERA
    'NO': [(0.65 * convert_to_mole_s(46), 'NOx')],
    'ALD2_PRIMARY': [(0.026318 * convert_to_mole_s(44.0526), 'VOC')],
    'FORM_PRIMARY': [(0.115995 * convert_to_mole_s(30.026), 'VOC')],
    'SOAALK': [(0.021027 * convert_to_mole_s(80.5347), 'VOC')],
    'ACET': [(0.029243 * convert_to_mole_s(58.0791), 'VOC')],
    'ALD2': [(0.026318 * convert_to_mole_s(44.0526), 'VOC')],
    'ALDX': [(0.081322 * convert_to_mole_s(40.7486), 'VOC')],
    'BENZ': [(0.025761 * convert_to_mole_s(78.1118), 'VOC')],
    'CH4': [(0.247726 * convert_to_mole_s(16.0425), 'VOC')],
    'ETH': [(0.059042 * convert_to_mole_s(28.0532), 'VOC')],
    'ETHA': [(0.03161 * convert_to_mole_s(30.069), 'VOC')],
    'ETHY': [(0.014204 * convert_to_mole_s(26.0373), 'VOC')],
    'ETOH': [(0.000949 * convert_to_mole_s(46.0684), 'VOC')],
    'FORM': [(0.115995 * convert_to_mole_s(30.026), 'VOC')],
    'IOLE': [(0.012975 * convert_to_mole_s(55.6133), 'VOC')],
    'ISOP': [(0.00527 * convert_to_mole_s(68.117), 'VOC')],
    'KET': [(0.007112 * convert_to_mole_s(17.8646), 'VOC')],
    'MEOH': [(0.105412 * convert_to_mole_s(32.0419), 'VOC')],
    'OLE': [(0.092184 * convert_to_mole_s(31.8777), 'VOC')],
    'PAR': [(0.200242 * convert_to_mole_s(21.0484), 'VOC')],
    'PRPA': [(0.013175 * convert_to_mole_s(44.0956), 'VOC')],
    'TERP': [(0.011595 * convert_to_mole_s(136.234), 'VOC')],
    'TOL': [(0.03537 * convert_to_mole_s(93.0188), 'VOC')],
    'UNR': [(0.269171 * convert_to_mole_s(36.5914), 'VOC')],
    'XYLMN': [(0.007905 * convert_to_mole_s(106.165), 'VOC')],
    'PAL': [(4.6 * (10 ** (-4)) * convert_to_g_s, 'PM2.5')],
    'PCA': [(7.2 * (10 ** (-4)) * convert_to_g_s, 'PM2.5')],
    'PCL': [(0.00239 * convert_to_g_s, 'PM2.5')],
    'PEC': [(0.1093 * convert_to_g_s, 'PM2.5')],
    'PFE': [(4.45 * (10 ** (-4)) * convert_to_g_s, 'PM2.5')],
    'PK': [(0.00135 * convert_to_g_s, 'PM2.5')],
    'PMN': [(1.1 * (10 ** (-4)) * convert_to_g_s, 'PM2.5')],
    'PMOTHR': [(0.0125 * convert_to_g_s, 'PM2.5')],
    'PNA': [(0.00135 * convert_to_g_s, 'PM2.5')],
    'PNCOM': [(0.3513 * convert_to_g_s, 'PM2.5')],
    'PNH4': [(0.00341 * convert_to_g_s, 'PM2.5')],
    'PNO3': [(0.0107 * convert_to_g_s, 'PM2.5')],
    'POC': [(0.5019 * convert_to_g_s, 'PM2.5')],
    'PSI': [(1.0 * (10 ** (-4)) * convert_to_g_s, 'PM2.5')],
    'PSO4': [(0.0033 * convert_to_g_s, 'PM2.5')],
    'PTI': [(6.7 * (10 ** (-4)) * convert_to_g_s, 'PM2.5')]
}

# The chemical species of emissions are used for emission mapping BlueSky -> CMAQ
select_species = ["CO", "SO2", "NH3", "NOx", "PM2.5", "VOC"]
with open(bsp_filename) as jsfile:
    bluesky_data = json.load(jsfile)

# To make sure the start time and end time of fires in the bluesky file provided
start_time = None
end_time = None
bsp_fires = bluesky_data["fires"]
for bsp_fire in bsp_fires:
    if start_time is None:
        start_time = get_fire_time(bsp_fire["activity"][0]["active_areas"][0]["start"], bsp_fire["activity"][0]["active_areas"][0]["utc_offset"])
    if end_time is None:
        end_time = get_fire_time(bsp_fire["activity"][0]["active_areas"][0]["end"], bsp_fire["activity"][0]["active_areas"][0]["utc_offset"])
    start_time = min(start_time, get_fire_time(bsp_fire["activity"][0]["active_areas"][0]["start"], bsp_fire["activity"][0]["active_areas"][0]["utc_offset"]))
    end_time = max(end_time, get_fire_time(bsp_fire["activity"][0]["active_areas"][0]["end"], bsp_fire["activity"][0]["active_areas"][0]["utc_offset"]))

current_time = start_time

met_filenames = []
while current_time <= end_time:
    metcros_file = input_dir + "METCRO3D_12US2_" + current_time.strftime("%Y%m%d")
    # generate emission file
    emission_tensor = dailyVerticalHourlyEmission(select_species, metcros_file, bsp_filename)

    current_metfile = metcros_file
    emission_filename = "emis_mole_all_" + current_time.strftime("%Y%m%d") + "_12US2_nobeis_norwc_2016gf_16j.ncf.gz"
    current_emiss_file = input_dir + emission_filename
    os.system("cp " + current_emiss_file + " " + work_dir)
    current_emiss_file = work_dir + emission_filename
    os.system("gunzip " + current_emiss_file)
    current_emiss_file = work_dir + "emis_mole_all_" + current_time.strftime("%Y%m%d") + "_12US2_nobeis_norwc_2016gf_16j.ncf"
    fire_emiss_file = output_dir + "emis_mole_all_" + current_time.strftime("%Y%m%d") + "_12US2_nobeis_norwc_2016gf_16j_fire.ncf"

    met_ds = netCDF4.Dataset(current_metfile)
    species_num, TIMESTEP, LAY, X, Y = emission_tensor.shape
    # Write data to netcdf file
    with netCDF4.Dataset(current_emiss_file) as src, netCDF4.Dataset(fire_emiss_file, "w",
                                                                     format='NETCDF3_64BIT_OFFSET', diskless=True,
                                                                     persist=True) as dst:
        # copy global attributes all at once via dictionary
        globle_dict = src.__dict__
        globle_dict['NLAYS'] = LAY
        globle_dict['VGLVLS'] = met_ds.__dict__['VGLVLS']
        dst.setncatts(globle_dict)
        # copy dimensions
        for name, dimension in src.dimensions.items():
            if name == 'LAY':
                dst.createDimension(name, LAY)
            else:
                dst.createDimension(name, len(dimension))
        # copy all file data except for the excluded
        for name, variable in src.variables.items():
            var_tmp = dst.createVariable(name, variable.datatype, variable.dimensions)
            # REMEMBER DO NOT CHANGE ALL VARIABLES
            if name in species_mapping.keys():
                dst[name][:] = np.zeros(var_tmp.shape)
                # generate emission data
                # Daily emission data, hour is 25
                print("Mapping " + name + " into Emission File")
                emission_specie_name = name
                # Generate emiss based on mapping
                cmaq_emission_data = np.zeros(var_tmp.shape)
                for mapping_tuple in species_mapping[emission_specie_name]:
                    specie_idx = select_species.index(mapping_tuple[1])
                    cmaq_emission_data += mapping_tuple[0] * emission_tensor[specie_idx, :, :, :, :]
                dst[name][:] = cmaq_emission_data
            elif name == "TFLAG":
                # REMEMBER DO NOT CHANGE ALL VARIABLES
                print("Do not revise variable:" + name)
                dst[name][:] = src[name][:]
            else:
                print("Set the emission variable to zeros: " + name)
                dst[name][:] = np.zeros(var_tmp.shape)
            # # copy variable attributes all at once via dictionary
            dst[name].setncatts(src[name].__dict__)

    # zip the output
    os.system("gzip " + fire_emiss_file)
    current_time += timedelta(days=1)
