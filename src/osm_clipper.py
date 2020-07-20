"""
Source code for osm_clipper

Copyright (C) 2020 Elco Koks. All versions released under the MIT license.
"""

import os
import numpy
import pandas
import pygeos
import geopandas
import urllib.request
import zipfile
from tqdm import tqdm
from multiprocessing import Pool,cpu_count

from shapely.geometry import MultiPolygon
from shapely.wkb import loads
from geopy.distance import geodesic


def planet_osm(data_path):
    """
    This function will download the planet file from the OSM servers. 
    """
    osm_path_in = os.path.join(data_path,'planet_osm')

    # create directory to save planet osm file if that directory does not exit yet.
    if not os.path.exists(osm_path_in):
        os.makedirs(osm_path_in)
    
    # if planet file is not downloaded yet, download it.
    if 'planet-latest.osm.pbf' not in os.listdir(osm_path_in):
        
        url = 'https://planet.openstreetmap.org/pbf/planet-latest.osm.pbf'
        urllib.request.urlretrieve(url, os.path.join(osm_path_in,'planet-latest.osm.pbf'))
    
    else:
        print('Planet file is already downloaded')

def country_osm(country,data_path):
    """
    This function will download the country file from the GeoFabrik servers.

    Note: some countries are not available no GeoFabrik. This will require you to 
        download the planet file and clip the country from their 
    """
    osm_path_in = os.path.join(data_path,'country_osm')

    # create directory to save planet osm file if that directory does not exit yet.
    if not os.path.exists(osm_path_in):
        os.makedirs(osm_path_in)
    
    # if planet file is not downloaded yet, download it.
    try:
        if '{}.osm.pbf'.format(country) not in os.listdir(osm_path_in):
            
            url = 'http://download.geofabrik.de/{}/{}-latest.osm.pbf'.format(country_list()[country][0],country_list()[country][1])
            urllib.request.urlretrieve(url, os.path.join(osm_path_in,'{}.osm.pbf'.format(country)))
        
        else:
            print('Country file is already downloaded')  
    except:
        print('{} is not available on GeoFabrik, you will have to clip it yourself from the planet file!'.format(country))      

def gadm36_planet(data_path):
    """
    This function will download the GADM36 file. 
    """

    gadm_path_in = os.path.join(data_path,'GADM36')

    # create directory to save planet osm file if that directory does not exit yet.
    if not os.path.exists(gadm_path_in):
        os.makedirs(gadm_path_in)
    
    # if GADM file is not downloaded yet, download it.
    if 'gadm36_levels.gpkg' not in os.listdir(gadm_path_in):
        
        url = 'https://biogeo.ucdavis.edu/data/gadm3.6/gadm36_levels_gpkg.zip'
        urllib.request.urlretrieve(url, os.path.join(gadm_path_in,'gadm36_levels_gpkg.zip'))
        with zipfile.ZipFile(os.path.join(gadm_path_in,'gadm36_levels_gpkg.zip'), 'r') as zip_ref:
            zip_ref.extractall(gadm_path_in)
        os.remove(os.path.join(gadm_path_in,'gadm36_levels_gpkg.zip'))
    
    else:
        print('GADM36 planet file is already downloaded')

def gadm36_country(country_code,data_path):
    """
    This function will download the GADM36 file. 
    """

    gadm_path_in = os.path.join(data_path,'GADM36')

    # create directory to save planet osm file if that directory does not exit yet.
    if not os.path.exists(gadm_path_in):
        os.makedirs(gadm_path_in)
    
    # if GADM file is not downloaded yet, download it.
    if 'gadm36_{}.gpkg'.format(country_code) not in os.listdir(gadm_path_in):
        
        url = 'https://biogeo.ucdavis.edu/data/gadm3.6/gpkg/gadm36_{}_gpkg.zip'.format(country_code)
        urllib.request.urlretrieve(url, os.path.join(gadm_path_in,'gadm36_{}_gpkg.zip').format(country_code))
        with zipfile.ZipFile(os.path.join(gadm_path_in,'gadm36_{}_gpkg.zip').format(country_code), 'r') as zip_ref:
            zip_ref.extractall(gadm_path_in)
        os.remove(os.path.join(gadm_path_in,'gadm36_{}_gpkg.zip').format(country_code))
    
    else:
        print('GADM36 country file is already downloaded')
         
 
def global_shapefiles(data_path,regionalized=False,assigned_level=1):
    """ 
    This function will simplify shapes and add necessary columns, to make further processing more quickly
    
    For now, we will make use of the latest GADM data, split by level: https://gadm.org/download_world.html

    Optional Arguments:
        *regionalized*  : Default is **False**. Set to **True** will also create the global_regions.shp file.
    """

    gadm_path = os.path.join(data_path,'GADM36','gadm36_levels.gpkg')
  
    # path to country GADM file
    if regionalized == False:
        
        # load country file
        gadm_level0 = pandas.DataFrame(geopandas.read_file(gadm_path,layer='level0'))

        #convert to pygeos
        tqdm.pandas(desc='Convert geometries to pygeos')
        gadm_level0['geometry'] = gadm_level0.geometry.progress_apply(lambda x: pygeos.from_shapely(x))

        # remove antarctica, no roads there anyways
        gadm_level0 = gadm_level0.loc[~gadm_level0['NAME_0'].isin(['Antarctica'])]
        
        # remove tiny shapes to reduce size substantially
        tqdm.pandas(desc='Remove tiny shapes')
        gadm_level0['geometry'] = gadm_level0.progress_apply(remove_tiny_shapes,axis=1)

        #simplify geometry
        tqdm.pandas(desc='Simplify geometry')
        gadm_level0.geometry = gadm_level0.geometry.progress_apply(lambda x: pygeos.simplify(
            pygeos.buffer(
                pygeos.simplify(
                    x,tolerance = 0.005, preserve_topology=True),0.01),tolerance = 0.005, preserve_topology=True))  
        
        #save to new country file
        glob_ctry_path = os.path.join(data_path,'cleaned_shapes','global_countries.gpkg')
        tqdm.pandas(desc='Convert geometries back to shapely')
        gadm_level0.geometry = gadm_level0.geometry.progress_apply(lambda x: loads(pygeos.to_wkb(x)))
        geopandas.GeoDataFrame(gadm_level0).to_file(glob_ctry_path,layer='level0', driver="GPKG")
          
    else:
        # this is dependent on the country file, so check whether that one is already created:
        glob_ctry_path = os.path.join(data_path,'cleaned_shapes','global_countries.gpkg')
        if os.path.exists(glob_ctry_path):
            gadm_level0 = geopandas.read_file(os.path.join(glob_ctry_path),layer='level0')
        else:
            print('ERROR: You need to create the country file first')   
            return None
        
        # load region file
        gadm_level_x = pandas.DataFrame(geopandas.read_file(gadm_path,layer='level{}'.format(assigned_level)))

        #convert to pygeos
        tqdm.pandas(desc='Convert geometries to pygeos')
        gadm_level_x['geometry'] = gadm_level_x.geometry.progress_apply(lambda x: pygeos.from_shapely(x))
        
        # remove tiny shapes to reduce size substantially
        tqdm.pandas(desc='Remove tiny shapes')
        gadm_level_x['geometry'] =   gadm_level_x.progress_apply(remove_tiny_shapes,axis=1)
    
         #simplify geometry
        tqdm.pandas(desc='Simplify geometry')
        gadm_level_x.geometry = gadm_level_x.geometry.progress_apply(lambda x: pygeos.simplify(
            pygeos.buffer(
                pygeos.simplify(
                    x,tolerance = 0.005, preserve_topology=True),0.01),tolerance = 0.005, preserve_topology=True))  

        # add some missing geometries from countries with no subregions
        get_missing_countries = list(set(list(gadm_level0.GID_0.unique())).difference(list(gadm_level_x.GID_0.unique())))
        
        #TO DO: GID_2 and lower tiers should first be filled by a tier above, rather then by the country file
        mis_country = gadm_level0.loc[gadm_level0['GID_0'].isin(get_missing_countries)]#
        if assigned_level==1:
            mis_country['GID_1'] = mis_country['GID_0']+'.'+str(0)+'_'+str(1)
        elif assigned_level==2:
            mis_country['GID_2'] = mis_country['GID_0']+'.'+str(0)+'.'+str(0)+'_'+str(1)
        elif assigned_level==3:
            mis_country['GID_3'] = mis_country['GID_0']+'.'+str(0)+'.'+str(0)+'.'+str(0)+'_'+str(1)
        elif assigned_level==4:
            mis_country['GID_4'] = mis_country['GID_0']+'.'+str(0)+'.'+str(0)+'.'+str(0)+'.'+str(0)+'_'+str(1)
        elif assigned_level==5:
            mis_country['GID_5'] = mis_country['GID_0']+'.'+str(0)+'.'+str(0)+'.'+str(0)+'.'+str(0)+'.'+str(0)+'_'+str(1)
        
        tqdm.pandas(desc='Convert geometries back to shapely')
        gadm_level_x.geometry = gadm_level_x.geometry.progress_apply(lambda x: loads(pygeos.to_wkb(x)))       

        # concat missing country to gadm levels 
        gadm_level_x = geopandas.GeoDataFrame( pandas.concat( [gadm_level_x,mis_country] ,ignore_index=True) )
        gadm_level_x.reset_index(drop=True,inplace=True)

        #save to new country file
        gadm_level_x.to_file(os.path.join(data_path,'cleaned_shapes','global_regions.gpkg'),layer='level{}'.format(assigned_level), driver="GPKG")
   
def remove_tiny_shapes(x,regionalized=False):
    """This function will remove the small shapes of multipolygons. Will reduce the size of the file.
    
    Arguments:
        *x* : a geometry feature (Polygon) to simplify. Countries which are very large will see larger (unhabitated) islands being removed.
    
    Optional Arguments:
        *regionalized*  : Default is **False**. Set to **True** will use lower threshold settings (default: **False**).
        
    Returns:
        *MultiPolygon* : a shapely geometry MultiPolygon without tiny shapes.
        
    """
    
    # if its a single polygon, just return the polygon geometry
    if pygeos.geometry.get_type_id(x.geometry) == 3: # 'Polygon':
        return x.geometry
    
    # if its a multipolygon, we start trying to simplify and remove shapes if its too big.
    elif pygeos.geometry.get_type_id(x.geometry) == 6: # 'MultiPolygon':
        
        if regionalized == False:
            area1 = 0.1
            area2 = 250
                
        elif regionalized == True:
            area1 = 0.01
            area2 = 50           

        # dont remove shapes if total area is already very small
        if pygeos.area(x.geometry) < area1:
            return x.geometry
        # remove bigger shapes if country is really big

        if x['GID_0'] in ['CHL','IDN']:
            threshold = 0.01
        elif x['GID_0'] in ['RUS','GRL','CAN','USA']:
            if regionalized == True:
                threshold = 0.01
            else:
                threshold = 0.01

        elif pygeos.area(x.geometry) > area2:
            threshold = 0.1
        else:
            threshold = 0.001

        # save remaining polygons as new multipolygon for the specific country
        new_geom = []
        for index_ in range(pygeos.geometry.get_num_geometries(x.geometry)):
            if pygeos.area(pygeos.geometry.get_geometry(x.geometry,index_)) > threshold:
                new_geom.append(pygeos.geometry.get_geometry(x.geometry,index_))
        
        return pygeos.creation.multipolygons(numpy.array(new_geom))

def poly_files(data_path,global_shape,regionalized=False):

    """
    This function will create the .poly files from the world shapefile. These
    .poly files are used to extract data from the openstreetmap files.
    
    This function is adapted from the OSMPoly function in QGIS.
    
    Arguments:
        *data_path* : base path to location of all files.
        
        *global_shape*: exact path to the global shapefile used to create the poly files.
        
    Optional Arguments:
        *save_shape_file* : Default is **False**. Set to **True** will the new shapefile with the 
        countries that we include in this analysis will be saved.     
        
        *regionalized*  : Default is **False**. Set to **True** will perform the analysis 
        on a regional level.
    
    Returns:
        *.poly file* for each country in a new dir in the working directory.
    """     
    
# =============================================================================
#     """ Create output dir for .poly files if it is doesnt exist yet"""
# =============================================================================
    poly_dir = os.path.join(data_path,'country_poly_files')
    
    if regionalized == True:
        poly_dir = os.path.join(data_path,'regional_poly_files')
    
    if not os.path.exists(poly_dir):
        os.makedirs(poly_dir)

# =============================================================================
#   """Load country shapes and country list and only keep the required countries"""
# =============================================================================
    wb_poly = geopandas.read_file(global_shape)
    
    # filter polygon file
    if regionalized == True:
        wb_poly = wb_poly.loc[wb_poly['GID_0'] != '-']

    else:
        wb_poly = wb_poly.loc[wb_poly['GID_0'] != '-']
   
    wb_poly.crs = {'init' :'epsg:4326'}

# =============================================================================
#   """ The important part of this function: create .poly files to clip the country 
#   data from the openstreetmap file """    
# =============================================================================
    num = 0
    # iterate over the counties (rows) in the world shapefile
    for f in wb_poly.iterrows():
        f = f[1]
        num = num + 1
        geom=f.geometry

#        try:
        # this will create a list of the different subpolygons
        if geom.geom_type == 'MultiPolygon':
            polygons = geom
        
        # the list will be lenght 1 if it is just one polygon
        elif geom.geom_type == 'Polygon':
            polygons = [geom]

        # define the name of the output file, based on the ISO3 code
        ctry = f['GID_0']
        if regionalized == True:
            attr=f['GID_1']
        else:
            attr=f['GID_0']
        
        # start writing the .poly file
        f = open(poly_dir + "/" + attr +'.poly', 'w')
        f.write(attr + "\n")

        i = 0
        
        # loop over the different polygons, get their exterior and write the 
        # coordinates of the ring to the .poly file
        for polygon in polygons:

            if ctry == 'CAN':
                dist = geodesic(reversed(polygon.centroid.coords[:1][0]), (83.24,-79.80), ellipsoid='WGS-84').kilometers
                if dist < 2000:
                    continue

            if ctry == 'RUS':
                dist = geodesic(reversed(polygon.centroid.coords[:1][0]), (82.26,58.89), ellipsoid='WGS-84').kilometers
                if dist < 500:
                    continue
                
            polygon = numpy.array(polygon.exterior)

            j = 0
            f.write(str(i) + "\n")

            for ring in polygon:
                j = j + 1
                f.write("    " + str(ring[0]) + "     " + str(ring[1]) +"\n")

            i = i + 1
            # close the ring of one subpolygon if done
            f.write("END" +"\n")

        # close the file when done
        f.write("END" +"\n")
        f.close()
#        except:
#            print(f['GID_1'])

def clip_osm_osmconvert(data_path,planet_path,area_poly,area_pbf):
    """ Clip the an area osm file from the larger continent (or planet) file and save to a new osm.pbf file. 
    This is much faster compared to clipping the osm.pbf file while extracting through ogr2ogr.
    
    This function uses the osmconvert tool, which can be found at http://wiki.openstreetmap.org/wiki/Osmconvert. 
    
    Either add the directory where this executable is located to your environmental variables or just put it in the 'scripts' directory.
    
    Arguments:
        *continent_osm*: path string to the osm.pbf file of the continent associated with the country.
        
        *area_poly*: path string to the .poly file, made through the 'create_poly_files' function.
        
        *area_pbf*: path string indicating the final output dir and output name of the new .osm.pbf file.
        
    Returns:
        a clipped .osm.pbf file.
    """ 
    print('{} started!'.format(area_pbf))

    osm_convert_path = os.path.join('osmconvert64-0.8.8p')

    try: 
        if (os.path.exists(area_pbf) is not True):
            os.system('{}  {} -B={} --complete-ways -o={}'.format(osm_convert_path,planet_path,area_poly,area_pbf))
        print('{} finished!'.format(area_pbf))

    except:
        print('{} did not finish!'.format(area_pbf))

def clip_osm_osmosis(data_path,planet_path,area_poly,area_pbf):
    """ Clip the an area osm file from the larger continent (or planet) file and save to a new osm.pbf file. 
    This is much faster compared to clipping the osm.pbf file while extracting through ogr2ogr.
    
    This function uses the osmconvert tool, which can be found at http://wiki.openstreetmap.org/wiki/Osmconvert. 
    
    Either add the directory where this executable is located to your environmental variables or just put it in the 'scripts' directory.
    
    Arguments:
        *continent_osm*: path string to the osm.pbf file of the continent associated with the country.
        
        *area_poly*: path string to the .poly file, made through the 'create_poly_files' function.
        
        *area_pbf*: path string indicating the final output dir and output name of the new .osm.pbf file.
        
    Returns:
        a clipped .osm.pbf file.
    """ 
    print('{} started!'.format(area_pbf))

    #osmosis_convert_path = os.path.join("..","osmosis","bin","osmosis.bat")
    osmosis_convert_path = os.path.join("osmosis")

    try: 
        if (os.path.exists(area_pbf) is not True):
            os.system('{} --read-xml file="{}" --bounding-polygon file="{}" --write-xml file="{}"'.format(osmosis_convert_path,planet_osm,area_poly,area_pbf))
        print('{} finished!'.format(area_pbf))

    except:
        print('{} did not finish!'.format(area_pbf))

def single_country(country,data_path,regionalized=False,create_poly_files=False,osm_convert=True,geofabrik=False):
    """    
    Obtain a country osm.pbf file. This can be done by extracting it from the planet file or by downloading it from GeoFabrik.
    
    This function has the option to extract individual regions. 
    
    Arguments:
        *country* : The country for which we want extract the data.
    
    Keyword Arguments:
        *regionalized* : Default is **False**. Set to **True** will parallelize the extraction over all regions within a country.
        
        *create_poly_files* : Default is **False**. Set to **True** will create new .poly files. 

        *osm_convert* : Default is **True**. Set to **False** will use osmosis. 

        *geofabrik* : Default is **False**. Set to **True** will try to obtain country file from geofabrik. 

    """
      
    if geofabrik==False:

        # path to planet file
        planet_path = os.path.join(data_path,'planet_osm','planet-latest.osm.pbf')

        # global shapefile path
        if regionalized == True:
            world_path = os.path.join(data_path,'input_data','global_regions.gpkg')
        else:
            world_path = os.path.join(data_path,'input_data','global_countries.gpkg')

        if not os.path.exists(os.path.join(data_path,'country_poly_files')):
            os.makedirs(os.path.join(data_path,'country_poly_files'))
    else:
        world_path = os.path.join(data_path,'input_data','global_regions.gpkg')

    # create poly files for all countries
    if (create_poly_files == True) & (geofabrik==False):
        poly_files(data_path,world_path,regionalized=regionalized)
    elif (create_poly_files == True) & (geofabrik==True):
        poly_files(data_path,world_path,regionalized=regionalized)

    if not os.path.exists(os.path.join(data_path,'country_osm')):
        os.makedirs(os.path.join(data_path,'country_osm'))
            
    ctry_poly = os.path.join(data_path,'country_poly_files','{}.poly'.format(country))
    ctry_pbf = os.path.join(data_path,'country_osm','{}.osm.pbf'.format(country))


    if (regionalized == False) & (geofabrik==True):
        country_osm(country,data_path)

    elif (regionalized == False) & (geofabrik==False):
        if osm_convert == True:
            try:
                clip_osm_osmconvert(data_path,planet_path,ctry_poly,ctry_pbf)
            except:
                print("NOTE: osmconvert is not correctly installed. Please check your environmental variables settings.")
        else:
            try:
                clip_osm_osmosis(data_path,planet_path,ctry_poly,ctry_pbf)
            except:
                print("NOTE: osmosis is not correctly installed. Please check your environmental variables settings.") 

    elif regionalized == True:
        
        country_pbf_exists = os.path.exists(ctry_pbf) 

        if (country_pbf_exists==False) & (geofabrik == True):
            country_osm(country,data_path)
        elif (country_pbf_exists==False) & (osm_convert == True):
            try:
                clip_osm_osmconvert(data_path,planet_path,ctry_poly,ctry_pbf)
            except:
                print("NOTE: osmconvert is not correctly installed. Please check your environmental variables settings.")
        elif (country_pbf_exists==False) & (osm_convert == True):
            try:
                clip_osm_osmosis(data_path,planet_path,ctry_poly,ctry_pbf)
            except:
                print("NOTE: osmosis is not correctly installed. Please check your environmental variables settings.") 
        else:
            print('Country file is already downloaded')  
           
        if not os.path.exists(os.path.join(data_path,'regional_poly_files')):
            os.makedirs(os.path.join(data_path,'regional_poly_files'))

        if not os.path.exists(os.path.join(data_path,'region_osm_admin1')):
            os.makedirs(os.path.join(data_path,'region_osm_admin1'))
        
        get_poly_files = [x for x in os.listdir(os.path.join(data_path,'regional_poly_files')) if x.startswith(country)]
        polyPaths = [os.path.join(data_path,'regional_poly_files',x) for x in get_poly_files]
        area_pbfs = [os.path.join(data_path,'region_osm_admin1',x.split('.')[0]+'_'+x.split('.')[1]+'.osm.pbf') for x in get_poly_files]
        data_paths = [data_path]*len(polyPaths)
        planet_paths = [ctry_pbf]*len(polyPaths)
      
        # and run all regions parallel to each other
        pool = Pool(cpu_count()-1)
        if osm_convert == True:
            pool.starmap(clip_osm_osmconvert, zip(data_paths,planet_paths,polyPaths,area_pbfs),chunksize=1) 
        else:
            pool.starmap(clip_osm_osmosis, zip(data_paths,planet_paths,polyPaths,area_pbfs),chunksize=1) 


def all_countries(subset = [], regionalized=False,reversed_order=False,osm_convert=True):
    """    
    Clip all countries from the planet osm file and save them to individual osm.pbf files
    
    Optional Arguments:
        *subset* : allow for a pre-defined subset of countries. REquires ISO3 codes. Will run all countries if left empty.
        
        *regionalized* : Default is **False**. Set to **True** if you want to have the regions of a country as well.
        
        *reversed_order* : Default is **False**. Set to **True**  to work backwards for a second process of the same country set to prevent overlapping calculations.
    
    Returns:
        clipped osm.pbf files for the defined set of countries (either the whole world by default or the specified subset)
    
    """
    # set data path
    data_path = os.path.join('..','data')
    
    # path to planet file
    planet_path = os.path.join(data_path,'planet_osm','planet-latest.osm.pbf')
   
    # global shapefile path
    if regionalized == True:
        world_path = os.path.join(data_path,'cleaned_shapes','global_regions.gpkg')
    else:
        world_path = os.path.join(data_path,'cleaned_shapes','global_countries.gpkg')

    # create poly files for all countries
    poly_files(data_path,world_path,regionalized=regionalized)
    
    # prepare lists for multiprocessing
    if not os.path.exists(os.path.join(data_path,'country_poly_files')):
        os.makedirs(os.path.join(data_path,'country_poly_files'))

    if not os.path.exists(os.path.join(data_path,'country_osm')):
        os.makedirs(os.path.join(data_path,'country_osm'))

    if regionalized == False:

        get_poly_files = os.listdir(os.path.join(data_path,'country_poly_files'))
        if len(subset) > 0:
            polyPaths = [os.path.join(data_path,'country_poly_files',x) for x in get_poly_files if x[:3] in subset]
            area_pbfs = [os.path.join(data_path,'region_osm_admin1',x.split('.')[0]+'.osm.pbf') for x in get_poly_files if x[:3] in subset]
        else:
            polyPaths = [os.path.join(data_path,'country_poly_files',x) for x in get_poly_files]
            area_pbfs = [os.path.join(data_path,'region_osm_admin1',x.split('.')[0]+'.osm.pbf') for x in get_poly_files]

        big_osm_paths = [planet_path]*len(polyPaths)
        
    elif regionalized == True:

        if not os.path.exists(os.path.join(data_path,'regional_poly_files')):
            os.makedirs(os.path.join(data_path,'regional_poly_files'))

        if not os.path.exists(os.path.join(data_path,'region_osm')):
            os.makedirs(os.path.join(data_path,'region_osm_admin1'))

        get_poly_files = os.listdir(os.path.join(data_path,'regional_poly_files'))
        if len(subset) > 0:
            polyPaths = [os.path.join(data_path,'regional_poly_files',x) for x in get_poly_files if x[:3] in subset]
            area_pbfs = [os.path.join(data_path,'region_osm_admin1',x.split('.')[0]+'.osm.pbf') for x in get_poly_files if x[:3] in subset]
            big_osm_paths = [os.path.join(data_path,'country_osm',x[:3]+'.osm.pbf') for x in get_poly_files if x[:3] in subset]
        else:
            polyPaths = [os.path.join(data_path,'regional_poly_files',x) for x in get_poly_files]
            area_pbfs = [os.path.join(data_path,'region_osm_admin1',x.split('.')[0]+'.osm.pbf') for x in get_poly_files]
            big_osm_paths = [os.path.join(data_path,'country_osm',x[:3]+'.osm.pbf') for x in get_poly_files]
            
    data_paths = [data_path]*len(polyPaths)

    # allow for reversed order if you want to run two at the same time (convenient to work backwards for the second process, to prevent overlapping calculation)   
    if reversed_order == True:
        polyPaths = polyPaths[::-1]
        area_pbfs = area_pbfs[::-1]
        big_osm_paths = big_osm_paths[::-1]

    # extract all country osm files through multiprocesing
    pool = Pool(cpu_count()-1)
    if osm_convert==True:
        pool.starmap(clip_osm_osmconvert, zip(data_paths,big_osm_paths,polyPaths,area_pbfs),chunksize=1) 
    else:
        pool.starmap(clip_osm_osmosis, zip(data_paths,big_osm_paths,polyPaths,area_pbfs),chunksize=1) 

def country_list():
    """
    Dictioniary to find match ISO-3 with required strings to download country from GeoFabrik.
    """

    return {
       'AFG' : ('asia','afghanistan'),
       'ALB' : ('europe','albania'),
       'DZA' : ('africa','algeria'),
       'AND' : ('europe','andorra'),
       'AGO' : ('africa','angola'),
       'BEN' : ('africa', 'benin'),
       'BWA' : ('africa', 'botswana'),
       'BFA' : ('africa', 'burkina-faso'),       
       'BDI' : ('africa', 'burundi'),
       'CMR' : ('africa', 'cameroon'),
       'BEN' : ('africa', 'canary-islands'),
       'CPV' : ('africa', 'cape-verde'),
       'CAF' : ('africa', 'central-african-republic'),
       'TCD' : ('africa', 'chad'),
       'COM' : ('africa', 'comores'),
       'COG' : ('africa', 'congo-brazzaville'),
       'COD' : ('africa', 'conga-democratic-republic'),
       'DJI' : ('africa', 'djibouti'),      
       'EGY' : ('africa', 'egypt'),
       'GNQ' : ('africa', 'equatorial-guinea'),
       'ERI' : ('africa', 'eritrea'),
       'ETH' : ('africa', 'ethopia'),
       'GAB' : ('africa', 'gabon'),
       'GMB' : ('africa', 'senegal-and-gambia'), #TOGETHER WITH SENEGAL
       'GHA' : ('africa', 'ghana'),
       'GIN' : ('africa', 'guinea'),
       'GNB' : ('africa', 'guinea-bissau'),
       'CIV' : ('africa', 'ivory-coast'),               
       'KEN' : ('africa', 'kenya'),      
       'LSO' : ('africa', 'lesotho'),
       'LBR' : ('africa', 'liberia'),
       'LBY' : ('africa', 'libya'),
       'MDG' : ('africa', 'madagascar'),
       'MWI' : ('africa', 'malawi'),
       'MLI' : ('africa', 'mali'),
       'MRT' : ('africa', 'mauritania'),
       'MAR' : ('africa', 'morocco'),
       'MOZ' : ('africa', 'mozambique'),     
       'NAM' : ('africa', 'namibia'),               
       'NER' : ('africa', 'niger'),      
       'NGA' : ('africa', 'nigeria'),
       'RWA' : ('africa', 'rwanda'),
       'SHN' : ('africa', 'saint-helena-ascension-and-tristan-da-cunha'),
       'STP' : ('africa', 'sao-tome-and-principe'),
       'SEN' : ('africa', 'senegal-and-gambia'), #TOGETHER WITH THE GAMBIA
       'SYC' : ('africa', 'seychelles'),
       'SLE' : ('africa', 'sierra-leone'),
       'SOM' : ('africa', 'somalia'),
       'ZAF' : ('africa', 'south-africa'),         
       'SDN' : ('africa', 'sudan'),    
       'SSD' : ('africa', 'south-sudan'),     
       'SWZ' : ('africa', 'swaziland'),               
       'TZA' : ('africa', 'tanzania'),      
       'TGO' : ('africa', 'togo'),
       'TUN' : ('africa', 'tunisia'),
       'UGA' : ('africa', 'uganda'),
       'ZMB' : ('africa', 'zambia'),
       'ZWE' : ('africa', 'zimbabwe'),
       'ARM' : ('asia', 'armenia'),
       'AZE' : ('asia', 'azerbaijan'),
       'BGD' : ('asia', 'bangladesh'),
       'BTN' : ('asia', 'bhutan'),                
       'KHM' : ('asia', 'cambodia'),
       'CHN' : ('asia', 'china'),
       'SAU' : ('asia', 'gcc-states'), #Together with Kuwait, the United Arab Emirates, Qatar, Bahrain, and Oman
       'KWT' : ('asia', 'gcc-states'), #Together with Saudi Arabia, the United Arab Emirates, Qatar, Bahrain, and Oman
       'ARE' : ('asia', 'gcc-states'), #Together with Saudi Arabia, Kuwait, Qatar, Bahrain, and Oman
       'QAT' : ('asia', 'gcc-states'), #Together with Saudi Arabia, Kuwait, the United Arab Emirates, Bahrain, and Oman
       'OMN' : ('asia', 'gcc-states'), #Together with Saudi Arabia, Kuwait, the United Arab Emirates, Qatar and Oman
       'BHR' : ('asia', 'gcc-states'), #Together with Saudi Arabia, Kuwait, the United Arab Emirates, Qatar and Bahrain
       'IND' : ('asia', 'india'),     
       'IDN' : ('asia', 'indonesia'),
       'IRN' : ('asia', 'iran'),
       'IRQ' : ('asia', 'iraq'),
       'ISR' : ('asia', 'israel-and-palestine'),       # TOGETHER WITH PALESTINE
       'PSE' : ('asia', 'israel-and-palestine'),       # TOGETHER WITH ISRAEL
       'JPN' : ('asia', 'japan'),
       'JOR' : ('asia', 'jordan'),
       'KAZ' : ('asia', 'kazakhstan'),
       'KGZ' : ('asia', 'kyrgyzstan'),             
       'LAO' : ('asia', 'laos'),
       'LBN' : ('asia', 'lebanon'),
       'MYS' : ('asia', 'malaysia-singapore-brunei'), # TOGETHER WITH SINGAPORE AND BRUNEI
       'SGP' : ('asia', 'malaysia-singapore-brunei'), # TOGETHER WITH MALAYSIA AND BRUNEI
       'BRN' : ('asia', 'malaysia-singapore-brunei'), # TOGETHER WITH MALAYSIA AND SINGAPORE
       'MDV' : ('asia', 'maldives'),                
       'MNG' : ('asia', 'mongolia'),
       'MMR' : ('asia', 'myanmar'),
       'NPL' : ('asia', 'nepal'),
       'PRK' : ('asia', 'north-korea'),       
       'PAK' : ('asia', 'pakistan'),
       'PHL' : ('asia', 'philippines'),                
       'RUS' : ('asia', 'russia'),
       'KOR' : ('asia', 'south-korea'),
       'LKA' : ('asia', 'sri-lanka'),
       'SYR' : ('asia', 'syria'),  
       'TWN' : ('asia', 'taiwan'),
       'TJK' : ('asia', 'tajikistan'),       
       'THA' : ('asia', 'thailand'),
       'TKM' : ('asia', 'turkmenistan'),                
       'UZB' : ('asia', 'uzbezkistan'),
       'VNM' : ('asia', 'vietnam'),
       'YEM' : ('asia', 'yemen'),
       'BHS' : ('central-america', 'bahamas'),   
       'BLZ' : ('central-america', 'belize'),                                                        
       'CUB' : ('central-america', 'cuba'),                                                        
       'GTM' : ('central-america', 'guatemala'),                                                        
       'HTI' : ('central-america', 'haiti-and-domrep'),  # TOGETHER WITH DOMINICAN REPUBLIC   
       'DOM' : ('central-america', 'haiti-and-domrep'),  # TOGETHER WITH HAITI                    
       'JAM' : ('central-america', 'jamaica'),                                                        
       'NIC' : ('central-america', 'nicaragua'),                                                        
       'AUT' : ('europe', 'austria'),                                                        
       'BLR' : ('europe', 'belarus'),                                                        
       'BEL' : ('europe', 'belgium'),                                                        
       'BIH' : ('europe', 'bosnia-herzegovina'),                                                        
       'BGR' : ('europe', 'bulgaria'),                                                        
       'HRV' : ('europe', 'croatia'),                                                        
       'CYP' : ('europe', 'cyprus'),                                                        
       'CZE' : ('europe', 'czech-republic'),                                                        
       'DNK' : ('europe', 'denmark'),                                                        
       'EST' : ('europe', 'estonia'),                                                        
       'FRO' : ('europe', 'faroe-islands'),                                                        
       'FIN' : ('europe', 'finland'),                                                        
       'FRA' : ('europe', 'france'),                                                        
       'GEO' : ('europe', 'georgia'),                                                        
       'DEU' : ('europe', 'germany'),                                                        
       'GBR' : ('europe', 'great-britain'),        # DOES NOT INCLUDE NORTHERN ISLAND                                                
       'GRC' : ('europe', 'greece'),                                                        
       'HUN' : ('europe', 'hungary'),                                                        
       'ISL' : ('europe', 'iceland'),                                                        
       'IRL' : ('europe', 'ireland-and-northern-ireland'),                                                        
       'IMN' : ('europe', 'isle-of-man'),                                                        
       'ITA' : ('europe', 'italy'),                                                        
       'LVA' : ('europe', 'latvia'),                                                        
       'LIE' : ('europe', 'liechtenstein'),    
       'LTU' : ('europe', 'lithuania'),                                                        
       'LUX' : ('europe', 'luxembourg'),                                                        
       'MKD' : ('europe', 'macedonia'),    
       'MLT' : ('europe', 'malta'),                                                        
       'MDA' : ('europe', 'moldova'),                                                        
       'MCO' : ('europe', 'monaco'),           
       'MNE' : ('europe', 'montenegro'),           
       'NLD' : ('europe', 'netherlands'),           
       'NOR' : ('europe', 'norway'),           
       'POL' : ('europe', 'poland'),           
       'PRT' : ('europe', 'portugal'),           
       'ROU' : ('europe', 'romania'),           
       'RUS' : ('europe', 'russia'),           
       'SRB' : ('europe', 'serbia'),           
       'SVK' : ('europe', 'slovakia'),           
       'SVN' : ('europe', 'slovenia'),           
       'ESP' : ('europe', 'spain'),           
       'SWE' : ('europe', 'sweden'),           
       'CHE' : ('europe', 'switzerland'),           
       'TUR' : ('europe', 'turkey'),           
       'UKR' : ('europe', 'ukraine'),           
       'CAN' : ('north-america', 'canada'),           
       'GRL' : ('north-america', 'greenland'),           
       'MEX' : ('north-america', 'mexico'),           
       'USA' : ('north-america', 'us'),           
       'AUS' : ('oceania', 'australia'),           
       'COK' : ('oceania', 'cook-islands'),           
       'FJI' : ('oceania', 'fiji'),           
       'KIR' : ('oceania', 'kiribati'),           
       'MHL' : ('oceania', 'marshall-islands'),           
       'FSM' : ('oceania', 'micronesia'),           
       'NRU' : ('oceania', 'nauru'),           
       'NCL' : ('oceania', 'new-caledonia'),           
       'NZL' : ('oceania', 'new-zealand'),           
       'NIU' : ('oceania', 'niue'),           
       'PLW' : ('oceania', 'palau'),           
       'PNG' : ('oceania', 'papua-new-guinea'),           
       'WSM' : ('oceania', 'samoa'),           
       'SLB' : ('oceania', 'solomon-islands'),           
       'TON' : ('oceania', 'tonga'),           
       'TUV' : ('oceania', 'tuvalu'),           
       'VUT' : ('oceania', 'vanuatu'),           
       'ARG' : ('south-america', 'argentina'),           
       'BOL' : ('south-america', 'bolivia'),           
       'BRA' : ('south-america', 'brazil'),           
       'CHL' : ('south-america', 'chile'),           
       'COL' : ('south-america', 'colombia'),           
       'ECU' : ('south-america', 'ecuador'),           
       'PRY' : ('south-america', 'paraguay'),           
       'PER' : ('south-america', 'peru'),
       'SUR' : ('south-america', 'suriname'),           
       'URY' : ('south-america', 'uruguay'),           
       'VEN' : ('south-america', 'venezuela'),           
    }