from datetime import datetime, timedelta
import json
from Mechanism import CB6_species_mapping
from EmissionGenerator import write_fire_emissions, get_fire_utc_time

case_type = 'rx'

# CASE 1, prescribed fires which is less than 1 day
if case_type == 'rx':
    bsp_filename = "./input_dir/bsp/SE2016_01_30_out.json"
    metcros_filename = "./input_dir/mcip/METCRO3D_12US2_20160130"
    output_filename = "./results/emis_mole_surface_12US2_20160130_cmaq_cb6_rx_fire_Briggs.ncf"

    write_fire_emissions(bsp_filename, metcros_filename, CB6_species_mapping, output_filename)

# CASE 2, wildfires which lasts longer than one day
if case_type == 'wf':
    bsp_filename = "./input_dir/bsp/WF_NM_sliced_out.json"
    with open(bsp_filename) as jsfile:
        bluesky_data = json.load(jsfile)

    start_time = None
    end_time = None
    bsp_fires = bluesky_data["fires"]
    for bsp_fire in bsp_fires:
        if start_time is None:
            start_time = get_fire_utc_time(bsp_fire["activity"][0]["active_areas"][0]["start"], bsp_fire["activity"][0]["active_areas"][0]["utc_offset"])
        if end_time is None:
            end_time = get_fire_utc_time(bsp_fire["activity"][0]["active_areas"][0]["end"], bsp_fire["activity"][0]["active_areas"][0]["utc_offset"])
        start_time = min(start_time, get_fire_utc_time(bsp_fire["activity"][0]["active_areas"][0]["start"], bsp_fire["activity"][0]["active_areas"][0]["utc_offset"]))
        end_time = max(end_time, get_fire_utc_time(bsp_fire["activity"][0]["active_areas"][0]["end"], bsp_fire["activity"][0]["active_areas"][0]["utc_offset"]))

    current_time = start_time
    while current_time <= end_time:
        current_metcro_3d = "./input_dir/mcip/METCRO3D_12US2_%s" % datetime.strftime(current_time, "%Y%m%d")
        output_filename = "./results/emis_mole_surface_12US2_%s_cmaq_cb6_fire_Briggs.ncf" % datetime.strftime(current_time, "%Y%m%d")
        write_fire_emissions(bsp_filename, current_metcro_3d, CB6_species_mapping, output_filename)
        current_time += timedelta(days=1)