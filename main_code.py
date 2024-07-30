"""
Created on Fri Jun 16 17:16:28 2023
@author: Olivier Chavanne


Modified on Sat Jan 20 11:17:47 2024
@author: Zetong Liu

Modified on Sat Jun 28 14:02:25 2024
@author: Giuliano De Pacsalis

"""


import shapely

import geopandas as gpd
import pandas as pd
from shapely.geometry import box, LinearRing
from shapely.ops import unary_union
import os
import matplotlib.pyplot as plt
import joblib

# Local libraries
from uhiCAD.building import generate_envelope
from uhiCAD.building import generate_buildings
from uhiCAD.building import get_scene_center
from uhiCAD.building import generate_pedestrian
from uhiCAD.building import generate_tree
import uhiCAD.xml as xml
# import enerCAD.result as result
# import enerCAD.network as network
# import uhiCAD.production as prod
# import uhiCAD.KPI as KPI

# URL for RegBL API request
GEOADMIN_BASE_URL = "https://api.geo.admin.ch/rest/services/ech/MapServer/ch.bfs.gebaeude_wohnungs_register/"
    
##################################################
# 
#                  Functions
#
##################################################

def create_xml_root(xml_file_to_copy, climate_file, horizon_file):
    '''
    Parameters                                                          
    ----------
    xml_file_to_copy : TYPE
        DESCRIPTION.
    climate_file : TYPE
        DESCRIPTION.
    horizon_file : TYPE
        DESCRIPTION.

    Returns
    -------
    root : TYPE
        DESCRIPTION.
    district : TYPE
        DESCRIPTION.
    '''
    
    # Write XML file for CitySim :
    print("Writing XML file...")    
    # Add Root 
    root = xml.add_root()
    # Add Simulation days
    xml.add_simulation_days(root)
    # Add Climate
    xml.add_climate(root, climate_file)
    # Add District
    district = xml.add_district(root)
    
    # Horizon
    # read in the tab-separated file as a dataframe
    horizon_df = pd.read_csv(horizon_file, sep='\s+', header=None)
    # assign column names to the dataframe
    horizon_df.columns = ['phi', 'theta']
    # Add Far field obstructions
    xml.add_far_field_obstructions(district, horizon_df)
    
    # Add all the composites and profiles, taken from a source XML
    xml.add_child_from_xml_to_district(district, xml_file_to_copy, 'Composite')
    xml.add_child_from_xml_to_district(district, xml_file_to_copy, 'OccupancyDayProfile')
    xml.add_child_from_xml_to_district(district, xml_file_to_copy, 'OccupancyYearProfile')
    xml.add_child_from_xml_to_district(district, xml_file_to_copy, 'DeviceType')
    xml.add_child_from_xml_to_district(district, xml_file_to_copy, 'ActivityType')
    xml.add_child_from_xml_to_district(district, xml_file_to_copy, 'DHWDayProfile')
    xml.add_child_from_xml_to_district(district, xml_file_to_copy, 'DHWYearProfile')
    
    xml.add_child_from_xml_to_district(district, xml_file_to_copy, 'Building')
    xml.add_child_from_xml_to_district(district, xml_file_to_copy, 'DistrictEnergyCenter')
    
    print("Xml source copied")
    
    return root, district 

def Module_1(gpkg_filepath,XYZfile, GEOADMIN_BASE_URL,
             directory_path, xml_name,
             xml_base_file, climate_file, horizon_file,
             create_geometry_3D=False, calculate_volume_3D=False,
             EGID_column='RegBL_EGID'):
    '''
    Parameters
    ----------
    gpkg_filepath : TYPE
        DESCRIPTION.
    GEOADMIN_BASE_URL : TYPE
        DESCRIPTION.
    directory_path : TYPE
        DESCRIPTION.
    xml_file_to_create : TYPE
        DESCRIPTION.
    xml_base_file : TYPE
        DESCRIPTION.
    climate_file : TYPE
        DESCRIPTION.
    horizon_file : TYPE
        DESCRIPTION.
    create_geometry_3D : TYPE, optional
        DESCRIPTION. The default is False.
    calculate_volume_3D : TYPE, optional
        DESCRIPTION. The default is False.
    EGID_column : TYPE, optional
        DESCRIPTION. The default is 'RegBL_EGID'.

    Returns
    -------
    None.
    '''
    
    ### Exctract geopackage ###
    
    print("Exctracting geopackage layers...")
    
    # MO Cadaster
    MO_all = gpd.read_file(gpkg_filepath, layer = "zone_tout")
    MO_dhn = gpd.read_file(gpkg_filepath, layer = "zone_cad")
    centrale = gpd.read_file(gpkg_filepath, layer = "centrale")
    pedestrian = gpd.read_file(gpkg_filepath, layer = "pedestrian")
    trees = gpd.read_file(gpkg_filepath, layer = "trees")
    new_trees = gpd.read_file(gpkg_filepath, layer = "new_trees")
    EGID_column = 'RegBL_EGID'
    
    # Split Multipolygons into Polygons
    zone_all = MO_all.explode(index_parts=False)
    zone_dhn = MO_dhn.explode(index_parts=False)
    
    # List containing EGID of buildings to simulate
    EGID_list = MO_dhn[EGID_column].tolist()
    
    # Save EGID list of buildings connected to CAD
    df_EGID = pd.DataFrame(EGID_list)
    df_EGID.columns = ['EGID']
    EGID_path = os.path.join(directory_path, 'EGID.csv')     
    df_EGID.to_csv(EGID_path, index=False)
    print("EGID.csv created")
    
    # Swissbuildings3D
    print("Swissbuildings3D processing...")
    try:
        floor_data = gpd.read_file(gpkg_filepath, layer = "floor")
        roof_data = gpd.read_file(gpkg_filepath, layer = "roof")
        wall_data = gpd.read_file(gpkg_filepath, layer = "wall")
        green_data = gpd.read_file(gpkg_filepath, layer = 'green')
        #ori_street_data = gpd.read_file(gpkg_filepath, layer = 'streets')
        #street_data = ori_street_data[~ori_street_data['objektart'].isin( ['Verbindung', 'Platz'])]
        street_data = gpd.read_file(gpkg_filepath, layer = 'streets')
        
        # Filter on the zone with 10m buffer around surrounding square box 
        zone_bounds = MO_all.geometry.buffer(10).values.total_bounds
        zone_box = box(zone_bounds[0], zone_bounds[1], zone_bounds[2], zone_bounds[3])
        
        # Cut swissbuildings3D to zone of concern
        floor_data_intersection = floor_data[floor_data.geometry.intersects(zone_box)]
        roof_data_intersection = roof_data[roof_data.geometry.intersects(zone_box)]
        wall_data_intersection = wall_data[wall_data.geometry.intersects(zone_box)]

        # Split Multipolygons into Polygons
        zone_floor = floor_data_intersection.explode(index_parts=True).reset_index()
        zone_roof = roof_data_intersection.explode(index_parts=True).reset_index()
        zone_wall = wall_data_intersection.explode(index_parts=True).reset_index()
        print('Swissbuildings3D cut to zone of interest \n')
    
    except: print('Error : Swissbuildings3D not provided')
    
    ### Buildings and ground XML processing ###
        
    root, district = create_xml_root(xml_base_file, climate_file, horizon_file)

    ### Footrpints processing ###
    
    try:
        # Get z coordinates of 1st vertex from 1st surface of 1st building's floor polygon as altitude by default for MO footprints
        altitude_default = zone_floor.loc[0].geometry.exterior.coords[0][2]
    except:
        altitude_default = 0
        
    # Create DataFrames containing all necessary information for each building
    print("Creating Buildings GeoDataFrame...")
    footprints, buildings = generate_buildings(zone_all, EGID_list, GEOADMIN_BASE_URL, altitude_default,
                                                create_geometry_3D, calculate_volume_3D, zone_floor, zone_roof, zone_wall)
    print("Buildings GeoDataFrame created \n") 
    
    # Get coordinates of the scene's center to allow more significant digits in further processing
    center_coordinates = get_scene_center(footprints)
        
    ### Ground processing, fist baseline, second scene###
    print("Generating ground...")
    terrain_df = pd.read_table(XYZfile, skiprows=1, delim_whitespace=True, names=['X', 'Y', 'Z'])
    gdf, gdf_intersection, ground_data_baseline, floor_geom3D, zone_geom, zone_geom_qgis = xml.add_ground_from_XYZ(district, terrain_df, zone_box, center_coordinates, buildings)
    
    xml.add_z_to_mo(gdf_intersection, zone_geom, buildings, altitude_default)
    ground_data = xml.close_ground(gdf, gdf_intersection, zone_geom)
    
    #ground_data = xml.fusion_triangles(gdf_full) 
    
    print('Adding Baseline ground in xml file...')
    # Write XML file for baseline
    xml.add_ground_in_xml(ground_data_baseline, district, center_coordinates, groundtype=3, kFactor=0.7, ShortWaveReflectance=0.22)
    #xml.add_ground_in_xml(3, 0.7, 0.22, ground_data_baseline, district, center_coordinates)
    xml.cut(district, ground_data_baseline, MO_dhn, footprints)
    baseline_path = os.path.join(directory_path, f"Baseline")
    os.makedirs(baseline_path, exist_ok=True)
    print('creating xml file \n')
    xml_to_create_path = os.path.join(baseline_path, xml_name+f"_baseline"+".xml")
    xml.write_xml_file(root, xml_to_create_path)
    print(f"{xml_name}.xml files created \n")
    
    # ShortWaveReflectance for Baseline
    b_SWR_df, b_SWA_df = xml.SWR(district)
    
    print('Remove Baseline ground from groundsurface')
    xml.cut_baseline(district, ground_data_baseline, MO_dhn, footprints)
    
    print("Adding actual scene ground in xml file...")
    # id, k-factor, SWR, distric, ground_data,..
    xml.add_ground_in_xml(ground_data, district, center_coordinates, groundtype=31, kFactor=0.1, ShortWaveReflectance=0.4)
    road_index_list, _ = xml.modify_type_road(district, ground_data, road_groundtype=2, road_kfactor=0.1, road_SWR=0.14, modify_data=street_data)
    green_index_list, green_itsctd = xml.modify_type(district, ground_data, groundtype=3, kfactor=0.7, SWR=0.22, modif_data=green_data)
    xml.cut(district, ground_data, MO_dhn, footprints)
    
    # Generate the envelope surfaces
    print("Generating Buildings envelope...")
    envelope, buildings_volume_3D = generate_envelope(footprints, buildings, calculate_volume_3D)
    print("Envelope created \n")
    
    # Merge "volume_3D" and "n_occupants" to main buildings geodataframe according to 'bid'
    merged_buildings = buildings.merge(buildings_volume_3D, left_on='bid', right_on='bid', how='left')    
    if not merged_buildings.empty:
        columns_to_add = ['volume_3D', 'n_occupants']
        for column in columns_to_add:
            buildings[column] = merged_buildings[column]
        print("Buildings 3D volume calculated and merged \n")
    
    print("Adding buildings in xml file...")
    # Add the buildings
    xml.add_all_buildings(district, buildings, envelope, center_coordinates)
    
    print("Creating pedestrian...")
    pedestrian_data, pedestrian_envelope = generate_pedestrian(pedestrian, buildings, ground_data)
    
    print("Adding pedestrians in xml file...")
    xml.add_pedestrians(district, pedestrian_data, pedestrian_envelope, center_coordinates)
    
    print("Creating trees...")
    next_tid = 0
    trees_data, trees_envelope, next_tid = generate_tree(trees, ground_data, next_tid)
    
    print("Adding trees in xml file...")
    xml.add_trees(district, trees_data, trees_envelope, center_coordinates)
    
    print('creating xml file \n')
    # write xml file of default case
    xml_path = os.path.join(directory_path, xml_name+".xml")     
    xml.write_xml_file(root, xml_path) 
    print(f"{xml_name}.xml files created \n")
    
    # ShortWaveReflectance for scene
    SWR_df, SWA_df = xml.SWR(district)

    print("Preparing XML for scenario case...")
    # Write XML file of Scenario 1
    print("Modifiyng ground...")
    sc_id=1
    sidewalk_geneva = gpd.read_file(gpkg_filepath, layer = 'sidewalk_geneva')
    ground_green = gpd.read_file(gpkg_filepath, layer = 'soil_green')
    fosse_tp = gpd.read_file(gpkg_filepath, layer = 'fosse_tp')
    fosse_impluvium = gpd.read_file(gpkg_filepath, layer = 'fosse_impluvium')
    
    sidewalk_geneva_index_list, sidewalk_geneva_itsctd = xml.modify_type(district, ground_data, groundtype=31, kfactor=0.1, SWR=0.4, modif_data=sidewalk_geneva)
    ground_green_index_list, ground_green_itsctd = xml.modify_type(district, ground_data, groundtype=3, kfactor=0.7, SWR=0.22, modif_data=ground_green)
    fosse_tp_index_list, fosse_tp_itsctd = xml.modify_type(district, ground_data, groundtype=7, kfactor=0.6, SWR=0.33, LongWaveEmissivity=0.95, modif_data= fosse_tp)
    fosse_impluvium_index_list, fosse_impluvium_itsctd = xml.modify_type(district, ground_data, groundtype=8, kfactor=0.6, SWR=0.33, LongWaveEmissivity=0.95, modif_data=fosse_impluvium)

    print("Adding extra trees..")
    new_trees_data, new_trees_envelope, next_tid = generate_tree(new_trees, ground_data, next_tid)
    xml.add_trees(district, new_trees_data, new_trees_envelope, center_coordinates)
    
    #approximation comparison
    print('approximation comparison \n') 
    selected_grounds1=ground_data[ground_data['gid'].isin(road_index_list)]
    selected_grounds2=ground_data[ground_data['gid'].isin(green_index_list)]
    selected_grounds3=ground_data[ground_data['gid'].isin(sidewalk_geneva_index_list)]
    selected_grounds4=ground_data[ground_data['gid'].isin(ground_green_index_list)]
    selected_grounds5=ground_data[ground_data['gid'].isin(fosse_tp_index_list)]
    selected_grounds6=ground_data[ground_data['gid'].isin(fosse_impluvium_index_list)]

    # add layers of simulated areas different from base ground
    selected_grounds1.to_file(gpkg_filepath, layer='street_grounds_2m', driver='GPKG')
    selected_grounds2.to_file(gpkg_filepath, layer='green_grounds_2m', driver='GPKG')
    selected_grounds3.to_file(gpkg_filepath, layer='sidewalk_geneva_2m', driver='GPKG')
    selected_grounds4.to_file(gpkg_filepath, layer='ground_green_2m', driver='GPKG')
    selected_grounds5.to_file(gpkg_filepath, layer='fosse_tp_2m', driver='GPKG')
    selected_grounds6.to_file(gpkg_filepath, layer='fosse_impluvium_2m', driver='GPKG')

    print('creating xml file \n')
    scenario_path = os.path.join(directory_path, f"Scenario_{sc_id}")
    os.makedirs(scenario_path, exist_ok=True)
    print(f"_sc_{sc_id}"+".xml" "files created \n")
    
    
    xml_to_create_path = os.path.join(scenario_path, xml_name+f"_sc_{sc_id}"+".xml")
    xml.write_xml_file(root, xml_to_create_path)
    
    # ShortWaveReflectance for scenario
    sc1_SWR_df, sc1_SWA_df = xml.SWR(district)

    return envelope, ground_data, buildings, zone_dhn, centrale, SWA_df, b_SWA_df, sc1_SWA_df, ground_data_baseline, pedestrian_data, zone_geom, zone_geom_qgis
   
def simulate_citysim(directory_path, xml_file, citysim_filepath):
    '''
    Parameters
    ----------
    xml_file : TYPE
        DESCRIPTION.

    Returns
    -------
    None.
    '''
    
    import subprocess
    import time
    start = time.time()
    print('Process started')
    print(f'Simulation of {xml_file}.xml...')

    #run CitySim.exe with xml file
    xml_path = os.path.join(directory_path, xml_file+".xml")
    result = subprocess.run([citysim_filepath, '-q', f"{xml_path}"])
    
    end = time.time()
    duration = end - start
    m, s = divmod(duration, 60)
    print('Simulation ended. Time :', "%.0f" %m,'min', "%.0f" %s,'s \n')
    

#------------------Part 2 iterating----------------------------------------------------------

def Module_2(directory_path, xml_name, gpkg_filepath, root, district,
             ground_data, climate_file, horizon_file,
             scenarios_list):   

    for i in range(len(scenarios_list)):
        sc_id = scenarios_list[i]
        sc1_data = gpd.read_file(gpkg_filepath, layer = f'scenario{sc_id}')
        root_copy, district_copy = root, district
        ### Scenarios XML processing ###
        # xml_to_copy_path = os.path.join(directory_path, xml_name+'.xml' )
        # root, district = create_xml_root(xml_to_copy_path, climate_file, horizon_file)
        road_index_list, green_index_list, _, _  = xml.modify_type(district_copy, ground_data, sc1_data)
        # Write XML file
        scenario_path = os.path.join(directory_path, f"Scenario_{sc_id}")
        os.makedirs(scenario_path, exist_ok=True)
        xml_to_create_path = os.path.join(scenario_path, xml_name+f"_sc_{sc_id}"+".xml")
        # xml.cut(district, ground_data, MO_dhn, footprints)

        xml.write_xml_file(root, xml_to_create_path)
        # print(f'{xml_DHN}_sc_{sc_id}.xml file created \n')

#--------------------- KPI calculation

def Module_KPI(ground_data, buffered_streets, itsctd_greens, road_index_list, green_index_list, scenarios, sc_id):

    real_str_area = buffered_streets['geometry'].area.sum()
    sim_str_area = ground_data.loc[ground_data['gid'].isin(road_index_list), 'geometry'].area.sum()
    str_error = abs(real_str_area-sim_str_area)/real_str_area
    real_gr_area = itsctd_greens['geometry'].area.sum()
    sim_gr_area = ground_data.loc[ground_data['gid'].isin(green_index_list), 'geometry'].area.sum()
    gr_error = abs(real_gr_area-sim_gr_area)/real_gr_area

    print('KPI calculated')
    
    return str_error,  gr_error


##################################################
# 
#         Information to provide
#
##################################################

# Geopackage filepath
gpkg_filepath = r"C:\Users\Giuliano\Documents\EPFL\Master\Projet de master\UHI_CH_sp-main\UHI_CH\Moutier_example.gpkg"

# Create geometry with swissbuildings3D
create_geometry_3D = True                           #TODO

# Calculate volume from swissbuildings3D
calculate_volume_3D = True
# CitySim.exe filepath
citysim_filepath = r"C:\Program Files\CitySim\CitySim.exe" #TODO

# XML name to export 
directory_path = r"C:\Users\Giuliano\Documents\EPFL\Master\Projet de master\UHI_CH_sp-main\test_moutier_example_14.07.2024_14.47"                              #TODO
os.makedirs(directory_path, exist_ok=True)                           
                                      
# XML source files
xml_base_file = r"C:\Users\Giuliano\Documents\EPFL\Master\Projet de master\UHI_CH_sp-main\xml_base.xml"    #TODO                   
horizon_file = r"C:\Users\Giuliano\Documents\EPFL\Master\Projet de master\Geneve\horizon_46.2016875747436_6.153443183166218.hor"   #TODO 
#r"C:\Users\Giuliano\Documents\EPFL\Master\Projet de master\Moutier\horizon_47.27779275971179_7.367307268533389.hor"
#r"C:\Users\Giuliano\Documents\EPFL\Master\Projet de master\UHI_CH_sp-main\UHI_CH\Lausanne.hor"        
XYZfile = r"C:\Users\Giuliano\Documents\EPFL\Master\Projet de master\Moutier\swissalti3d_2019_2594-1236_2_2056_5728.xyz\SWISSALTI3D_2_XYZ_CHLV95_LN02_2594_1236.xyz" #TODO     
#Geneve
#r"C:\Users\Giuliano\Documents\EPFL\Master\Projet de master\Genève\Swissalit3d\SWISSALTI3D_2_XYZ_CHLV95_LN02_2500_1117.xyz"
#Moutier
#r"C:\Users\Giuliano\Documents\EPFL\Master\Projet de master\Moutier\swissalti3d_2019_2594-1236_2_2056_5728.xyz\SWISSALTI3D_2_XYZ_CHLV95_LN02_2594_1236.xyz"

# Scenarios to simulate
scenarios_list = [1]                                #TODO
do_plot = True

def main(): 
    
    # Generate individual buildings XML
    print('***Module 1*** \n')
    Year_of_cli=['Contemporary'] #, 'RCP45_2030', 'RCP45_2040'] 
    for year in Year_of_cli:
        subdirectory_path = os.path.join(directory_path, f"{year}")
        os.makedirs(subdirectory_path, exist_ok=True)
        xml_name = f'_{year}' 
        #xml_name = directory_path+f'_{year}' 
        climate_file = rf"C:\Users\Giuliano\Documents\EPFL\Master\Projet de master\Geneve\climat_geneve\Contemporain\CH_-_Geneve_-_Rive_2-hour.cli"  #TODO
        #rf"C:\Users\Giuliano\Documents\EPFL\Master\Projet de master\Genève\climat_geneve\Contemporain\CH_-_Geneve_-_Rive_2-hour.cli"
        
        # Height of meteo statio
        h_meteo_station = xml.height_meteo(climate_file)
        # dataframe and hottesat day in the year.
        # Hottest day = day with highest mean air temperature. Take tmax of hottest day
        climate_df, Ta_max_values, Ta_max_index = xml.hottest_day(climate_file)
    

        envelope, ground_data, buildings, zone_dhn, centrale, SWA_df, b_SWA_df, sc1_SWA_df, ground_data_baseline, pedestrian_data, zone_geom, zone_geom_qgis = Module_1(gpkg_filepath, XYZfile, GEOADMIN_BASE_URL,  
                                                subdirectory_path, xml_name,
                                                xml_base_file, climate_file, horizon_file,
                                                create_geometry_3D, calculate_volume_3D,
                                                EGID_column='RegBL_EGID')
        
        # Simulation for actual scene
        print('***Actual Scene*** \n')
        simulate_citysim(subdirectory_path, xml_name, citysim_filepath)
        
        #TS
        TS_grounds, TS_buildings = xml.clean_sort_OUT_files(subdirectory_path, xml_name, "_TS.out")
        #LW
        LW_grounds, LW_buildings = xml.clean_sort_OUT_files(subdirectory_path, xml_name, "_LW.out")
        #SW
        SW_grounds, SW_buildings = xml.clean_sort_OUT_files(subdirectory_path, xml_name, "_SW.out")
        columns_to_keep = TS_grounds.columns

        # Filter SW_grounds to keep only columns taht are in TS_grounds
        SW_grounds_filtered = SW_grounds[columns_to_keep]

        VF_file = os.path.join(subdirectory_path, xml_name+"_VF.out")
        VF_df = pd.read_csv(VF_file, delimiter='\t', encoding='latin1')
        VF_df.columns = VF_df.columns.str.replace('Â', '', regex=True)

        VF_df = VF_df.rename(columns={'#': 'id'})
        VF_grounds = VF_df[VF_df['id'].str.contains('NA', na=False)]
        VF_grounds = VF_df[VF_df['Type'] == 'Ground']
        VF_grounds['id'] = VF_grounds['id'].str.split(':').str[-1].str.strip()
        VF_grounds = VF_grounds.rename(columns={'id': 'gid'})
        VF_grounds['gid'] = VF_grounds['gid'].astype(int)

        CM_file = os.path.join(subdirectory_path, xml_name+"_CM.out")
        CM_df = pd.read_csv(CM_file, delimiter='\t', encoding='latin1')
        # Remove the # from the name of the first column
        CM_df.columns = [col.lstrip('#') for col in CM_df.columns]

        Year_RPB_file = os.path.join(subdirectory_path, xml_name+"_YearlyResultsPerBuilding.out")
        Year_RPB_df = pd.read_csv(Year_RPB_file, delimiter='\t', encoding='latin1')
        Year_RPB_df.columns = [col.lstrip('#') for col in Year_RPB_df.columns]
        #Check avec ground data
        ground_data_filtered = ground_data[ground_data['gid'].isin(columns_to_keep)]

        # wb = without buildings
        hc_df = xml.concevtion_coefficient(climate_df, ground_data_filtered, h_meteo_station)
        Tsat_df = xml.T_sol_air_2(TS_grounds, LW_grounds, SW_grounds, SWA_df, hc_df)
        print('Tsolair calculation completed')
        
        ground_SAT = xml.add_hd_data_gpkg_2(Ta_max_index, Ta_max_values, ground_data_filtered, Tsat_df, TS_grounds, SW_grounds_filtered, LW_grounds, hc_df, SWA_df, VF_grounds)
        # Roof buildings, r = roof
        TS_r_df, LW_r_df, SW_r_df, SWA_r_df, hc_r_df, Tsat_r_df = xml.extract_roof_output_2(TS_buildings, LW_buildings, SW_buildings, envelope)
        roof_SAT = xml.add_hd_data_roof_gpkg(Ta_max_index, Ta_max_values, zone_geom_qgis, Tsat_r_df, TS_r_df, SW_r_df, LW_r_df, hc_r_df, SWA_r_df)

        total_SAT = pd.concat([ground_SAT, roof_SAT], ignore_index=True)
        total_SAT.to_file(gpkg_filepath, layer=f'Tsat_@_Ta_max_{year}', driver='GPKG')
        print('Tsol-air added in GPKG')

        pedestrian_GPKG_TMax = xml.add_pedestrian_GPKG_Tmax_2(CM_df, Year_RPB_df, pedestrian_data, Ta_max_index, Ta_max_values)
        pedestrian_GPKG_TMax.to_file(gpkg_filepath, layer=f'Pedestrian_@_Ta_max_{year}', driver='GPKG')
        print('Thermal comfort indices added in GPKG')

        _, _, all_AST = xml.AST_3(envelope, TS_grounds, TS_buildings, ground_data_filtered, zone_geom)
        all_AST.to_file(gpkg_filepath, layer=f'all_AST_{year}', driver='GPKG')
        print('AST added in GPKG')

        # Simulation for baseline
        print('***Baseline*** \n')
        b_path = os.path.join(subdirectory_path, f"Baseline")
        simulate_citysim(b_path, f'{xml_name}_baseline', citysim_filepath)

        #Post process baseline
        b_path = os.path.join(subdirectory_path, f"Baseline")
        xml_name_baseline = f'{xml_name}_baseline'

        #TS
        b_TS_grounds, b_TS_buildings = xml.clean_sort_OUT_files(b_path, xml_name_baseline, "_TS.out")
        # LW
        b_LW_grounds, b_LW_buildings = xml.clean_sort_OUT_files(b_path, xml_name_baseline, "_LW.out")
        # SW
        b_SW_grounds, b_SW_buildings = xml.clean_sort_OUT_files(b_path, xml_name_baseline, "_SW.out")
        columns_to_keep = b_TS_grounds.columns
        # Filter SW_grounds to keep only columns taht are in TS_grounds
        b_SW_grounds_filtered = b_SW_grounds[columns_to_keep]

        b_VF_file = os.path.join(b_path, xml_name_baseline+"_VF.out")
        b_VF_df = pd.read_csv(b_VF_file, delimiter='\t', encoding='latin1')
        b_VF_df.columns = b_VF_df.columns.str.replace('Â', '', regex=True)

        b_VF_df = b_VF_df.rename(columns={'#': 'id'})
        b_VF_grounds = b_VF_df[b_VF_df['id'].str.contains('NA', na=False)]
        b_VF_grounds = b_VF_df[b_VF_df['Type'] == 'Ground']
        b_VF_grounds['id'] = b_VF_grounds['id'].str.split(':').str[-1].str.strip()
        b_VF_grounds = b_VF_grounds.rename(columns={'id': 'gid'})
        b_VF_grounds['gid'] = b_VF_grounds['gid'].astype(int)

        #Check with ground data
        ground_data_baseline_filtered = ground_data_baseline[ground_data_baseline['gid'].isin(columns_to_keep)]

        # wb = without buildings
        b_hc_df = xml.concevtion_coefficient(climate_df, ground_data_baseline_filtered, h_meteo_station)
        b_Tsat_df = xml.T_sol_air_2(b_TS_grounds, b_LW_grounds, b_SW_grounds, b_SWA_df, b_hc_df)
        print('Tsol-air calculation completed for baseline')

        b_ground_SAT = xml.add_hd_data_gpkg_2(Ta_max_index, Ta_max_values, ground_data_baseline_filtered, b_Tsat_df, b_TS_grounds, b_SW_grounds_filtered, b_LW_grounds, b_hc_df, b_SWA_df, b_VF_grounds)
        b_ground_SAT.to_file(gpkg_filepath, layer=f'Tsat_@_Ta_max_{year}_baseline', driver='GPKG')
        print('Tsol-air added in GPKG for baseline')

        b_all_AST = xml.AST_baseline_3(envelope, b_TS_grounds, ground_data_baseline_filtered)
        b_all_AST.to_file(gpkg_filepath, layer=f'all_AST_{year}_baseline', driver='GPKG')
        print('AST added in gpkg for baseline')

        # Scenarii simulations
        for i in range(len(scenarios_list)):
            sc_id = scenarios_list[i]
            print(f'***Scenario {sc_id}*** \n')
            scenario_path = os.path.join(subdirectory_path, f"Scenario_{sc_id}")
            simulate_citysim(scenario_path, f'{xml_name}_sc_{sc_id}', citysim_filepath)
            
            #TS
            sc1_TS_grounds, sc1_TS_buildings = xml.clean_sort_OUT_files(scenario_path, f'{xml_name}_sc_{sc_id}', "_TS.out")
            #LW
            sc1_LW_grounds, sc1_LW_buildings = xml.clean_sort_OUT_files(scenario_path, f'{xml_name}_sc_{sc_id}', "_LW.out")
            #SW
            sc1_SW_grounds, sc1_SW_buildings = xml.clean_sort_OUT_files(scenario_path, f'{xml_name}_sc_{sc_id}', "_SW.out")
            columns_to_keep = sc1_TS_grounds.columns
            # Filter SW_grounds to keep only columns taht are in TS_grounds
            sc1_SW_grounds_filtered = sc1_SW_grounds[columns_to_keep]

            sc1_VF_file = os.path.join(scenario_path, f'{xml_name}_sc_{sc_id}' + "_VF.out")
            sc1_VF_df = pd.read_csv(sc1_VF_file, delimiter='\t', encoding='latin1')
            sc1_VF_df.columns = sc1_VF_df.columns.str.replace('Â', '', regex=True)

            sc1_VF_df = sc1_VF_df.rename(columns={'#': 'id'})
            sc1_VF_grounds = sc1_VF_df[sc1_VF_df['id'].str.contains('NA', na=False)]
            sc1_VF_grounds = sc1_VF_df[sc1_VF_df['Type'] == 'Ground']
            sc1_VF_grounds['id'] = sc1_VF_grounds['id'].str.split(':').str[-1].str.strip()
            sc1_VF_grounds = sc1_VF_grounds.rename(columns={'id': 'gid'})
            sc1_VF_grounds['gid'] = sc1_VF_grounds['gid'].astype(int)
            
            sc1_CM_file = os.path.join(scenario_path, f'{xml_name}_sc_{sc_id}'+"_CM.out")
            sc1_CM_df = pd.read_csv(sc1_CM_file, delimiter='\t', encoding='latin1')
            sc1_CM_df.columns = [col.lstrip('#') for col in sc1_CM_df.columns]
            
            sc1_Year_RPB_file = os.path.join(scenario_path, f'{xml_name}_sc_{sc_id}' + "_YearlyResultsPerBuilding.out")
            sc1_Year_RPB_df = pd.read_csv(sc1_Year_RPB_file, delimiter='\t', encoding='latin1')
            sc1_Year_RPB_df.columns = [col.lstrip('#') for col in sc1_Year_RPB_df.columns]

            # wb = without buildings
            #sc1_hc_df = xml.concevtion_coefficient(climate_df, ground_data_filtered, h_meteo_station) # no need of computing it twice
            sc1_Tsat_df = xml.T_sol_air_2(sc1_TS_grounds, sc1_LW_grounds, sc1_SW_grounds, sc1_SWA_df, hc_df)
            print('Tsol-air calculation completed for scenario')

            sc1_ground_SAT = xml.add_hd_data_gpkg_2(Ta_max_index, Ta_max_values, ground_data_filtered, sc1_Tsat_df, sc1_TS_grounds, sc1_SW_grounds_filtered, sc1_LW_grounds, hc_df, sc1_SWA_df, sc1_VF_grounds)
            # Roof buildings, r = roof
            sc1_TS_r_df, sc1_LW_r_df, sc1_SW_r_df, sc1_SWA_r_df, sc1_hc_r_df, sc1_Tsat_r_df = xml.extract_roof_output_2(sc1_TS_buildings, sc1_LW_buildings, sc1_SW_buildings, envelope)
            sc1_roof_SAT = xml.add_hd_data_roof_gpkg(Ta_max_index, Ta_max_values, zone_geom_qgis, sc1_Tsat_r_df, sc1_TS_r_df, sc1_SW_r_df, sc1_LW_r_df, sc1_hc_r_df, sc1_SWA_r_df)
            sc1_total_SAT = pd.concat([sc1_ground_SAT, sc1_roof_SAT], ignore_index=True)
            sc1_total_SAT.to_file(gpkg_filepath, layer=f'Tsat_@_Ta_max_{year}_sc_{sc_id}', driver='GPKG')
            print('Tsol-air added in GPKG for scenario')

            # Post processing: Pedestrian comfort
            sc1_pedestrian_GPKG_TMax = xml.add_pedestrian_GPKG_Tmax_2(sc1_CM_df, sc1_Year_RPB_df, pedestrian_data, Ta_max_index, Ta_max_values)
            sc1_pedestrian_GPKG_TMax.to_file(gpkg_filepath, layer=f'Pedestrian_@_Ta_max_{year}_sc_{sc_id}', driver='GPKG')
            print('Thermal comfort indices added in GPKG for scenario')

            _, _, scenario_all_AST = xml.AST_3(envelope, sc1_TS_grounds, sc1_TS_buildings, ground_data_filtered, zone_geom)
            scenario_all_AST.to_file(gpkg_filepath, layer=f'all_AST_{year}_sc_{sc_id}', driver='GPKG')
            #merged_df = gpd.GeoDataFrame(pd.merge(all_AST, scenario_all_AST, on='geometry', suffixes=('_dfT', '_s1')))
            #merged_df['T_difference'] = merged_df['T_s1'] - merged_df['T_dfT']
            print('AST added in gpkg scenario') 
            
            print(f"Scenario {sc_id} at {year}processed \n")
        print(f"year {year} processed \n")
    print("***Overall processing finished***")
    print(f"Find all results and graphs in directory : {directory_path}")

if __name__ == "__main__":
    plt.close("all")
    
    import subprocess
    import time
    start_overall = time.time()
    print('Main code started')
    print('-----------------')
    
    main()    

    print('-----------------')
    print('Main code ended')
    print('-----------------')
    end_overall = time.time()
    duration_overall = end_overall - start_overall
    m, s = divmod(duration_overall, 60)
    print('Overall run time :', "%.0f" %m,'min', "%.0f" %s,'s \n')
    
    
    

