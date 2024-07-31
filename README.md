# Urban Heat Island and Pedestrian Comfort in Switzerland (CH)
An open-source framework for the quantification of Urban Heat Islands and Pedestrian Comfort in Switzerland

## Introduction
In this project, Swiss open-source datasets are used to quantify urban heat islands (UHI) and pedestrian comfort, and to demonstrate the application of different scenarios on a case study representative of the Swiss landscape.

## Usage of codes
### Open-source files and import to QGIS
Choose the geographical zone (tile, canton) for the area of study. 
Step 1 : Download the files, Step 2 : Import them to QGIS

Swissbuildings3D : 
1. Geodatabase format from https://www.swisstopo.admin.ch/en/geodata/landscape/buildings3d3.html#download
2. Import *Floor/Roof/Wall/streets/green* layers

MO cadaster
1. GeoPackage format from https://geodienste.ch/services/av
2. Import *lcsf* layer, filter by "Genre"="batiment"

Ground Surface
1. XYZ file format from https://www.swisstopo.admin.ch/en/geodata/height/alti3d.html
2. Choose resolution of 2 meter or o.5 meter and load it with pandas.read_table

Ground Types 
1. Geopackage format from https://www.swisstopo.admin.ch/de/landschaftsmodell-swisstlm3d
2. Import *tlm_bb_bodenbedeckung* layer, filter by 'objektart' IN ('Wald', 'Gehoelzflaeche', 'Gebueschwald', 'Wald offen')
3. Import *tlm_strassen_strasse* layer, filter by "objektart" NOT IN ('Verbindung','Platz')

### QGIS
Create a new GeoPackage file. Add layers with the "export features" option of QGIS (export without "fid" field), with names :
- *zone_cad* : MO features of buildings connected to DHN
- *zone_tout* : MO features of all buildings in the area of study
- *centrale* : point feature of thermal heating station coordinates
- *floor* : swissbuildings3D features of floors of all buildings in the area of study
- *wall* : swissbuildings3D features of walls of all buildings in the area of study
- *roof* : swissbuildings3D features of roofs of all buildings in the area of study
- *streets*: swisstlm 3D features of streets in the area of study
- *green*: features of green areas in the area of study
- *trees*: point features of existing trees in the study area (These data can come from open-source datasets if they exist for the study area)
- *pedestrian*: point features of pedestrian in the study area

Layers for creating scenarios:
- *sidewalk_geneva*: polygons representing the sidewalk in Geneva
- *fosse_tp*: polygons representing soil and stone pits
- *fosse_impluvium*: polygons representing Stockholm pits
- *soil_green*: polygons representing green spaces
- *new_trees*: point features of new trees in the study area

### Code custom modifications
In addition to the newly created GeoPackage, a climatic and a horizon file must be provided and imported on the same directory as where the software code was cloned.

Before launching the "main_code.py" from command line, it must be custom modified for each simulation, either from a text files reading app or a python editor app such as *Spyder* or *Visual Code Studio*.
Each modification is signaled with *#TODO*.
- gpkg_filepath = r"---.gpkg" : File path of the GeoPackage containing all necessary layers
- create_geometry_3D = True/False (default = False) : Activates the simulation with thermal envelope from Swissbuildings3D geometries (much longer simulation)
- calculate_volume_3D = True/False (default = False) : Activates the volume calculation from Swissbuildings3D geometries
- citysim_filepath = r"---/CitySim.exe" : File path of the CitySim solver
- directory_path = r"---" : Name of the new directory to be created by the simulation
- climate_file = r"---.cli" : File path of the climate file
- horizon_file = r"---.hor" : File path of the horizon file

### Python libraries
The required libraries to import are listed in the *requirements.txt* file

### Results
The resulted .xml file is under the directory path defined in "main_code". There are three different result layers with prefix: "All_AST", "Tsat_@_Tmax", "Pedestrian_@_Ta_max" in the .gpkg file used.

### Bakcground
Olivier Chavanne established a framework to process and analyze information of all buildings in the area of study and simulate the surface temperature of those buildings through Citysim. Zetong Liu added the ground in the area of study as well as new pavements: Aspahlt, Concrete, Green and Pavement. To measure UHIs and pedestrian comfort, it is necessary to add the ground around buildings, include trees and pedestrians in the simulation file, and introduce new indicators to quantify UHIs and pedestrian comfort.

### Trees and pedestrian
The addition of trees and pedestrians is based on [Coccolo (2017)](https://infoscience.epfl.ch/entities/publication/7a844472-d333-4867-8364-48772278a7a5). The aim is to convert points from the geopackage file into pedestrian and tree geometries in the simulation file. The pedestrian geometries are created to calculate the Mean Radiant Temperature.

### UHI inidcators
Two new indicators are added to quantify urban heat islands. They are selected for the warmest average annual day at its maximum temperature. These indicators are surface temperature "Ts" and soil-air temperature "Tsol-air", which aim to capture the convective effects of the air.

### Pedestrian comfort indicators
The Universal Thermal Climate Index (UTCI) is implemented because it can assess huma thermal perception in all climates. The annual number of hours of thermal comfort is calculated using this indicator.

## Reference
1. https://github.com/ochavanne/CAD-O
2. https://github.com/ZetongLiu/UHI_CH_sp
3. https://github.com/kaemco/CitySim-Solver/wiki/The-district#groundsurface
4. Silvia Coccolo, Bioclimatic Design of Sustainable Campuses using Advanced Optimisation
Methods, 2017, doi: [https : / / doi . org / 10 . 5075 / epfl - thesis - 7756]([url](https://doi.org/10.5075/epfl-thesis-7756)), URL: [http://infoscience.epfl.ch/record/231147]([url](http://infoscience.epfl.ch/record/231147))
  


