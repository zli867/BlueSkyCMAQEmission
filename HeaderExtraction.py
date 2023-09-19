import netCDF4
import json


antro_emiss_file = "./input_dir/emis_mole_all_20160130_GA4_nobeis_norwc_2016gf_16j.ncf"
var_dict_info = {}

with netCDF4.Dataset(antro_emiss_file) as src:
    # copy global attributes all at once via dictionary
    globle_dict = src.__dict__
    var_dict_info['VAR-LIST'] = globle_dict['VAR-LIST']
    # copy all file data except for the excluded
    for name, variable in src.variables.items():
        var_dict_info[name] = src[name].__dict__


json_obj = json.dumps(var_dict_info, indent=4)

# Generate BlueSKY input
input_filename = "./HeaderInfo.py"
with open(input_filename, 'w') as fout:
    print(json_obj, file=fout)
