import os,sys
from tqdm import tqdm
import geopandas as gpd
import pygeos
sys.path.append(os.path.join( '..','src'))
from osm_clipper import remove_tiny_shapes,poly_files,planet_osm,single_country

#from osm_clipper import global_shapefiles

if __name__ == '__main__':

    #planet_osm()
    single_country('VAT')

    # data_path = os.path.join('..','data')
    # country_gadm_path = os.path.join(data_path,'GADM36','gadm36_levels.gpkg')
    # gadm_level0 = gpd.read_file(country_gadm_path,layer='level0')

    # for gadm in gadm_level0.itertuples(index=False):
    #     print(gadm.GID_0)
    #     pygeos.from_shapely(gadm.geometry)

    # remove antarctica, no roads there anyways
    #gadm_level0 = gadm_level0.loc[~gadm_level0['NAME_0'].isin(['Antarctica'])]
    
    # remove tiny shapes to reduce size substantially
    #gadm_level0['geometry'] = gadm_level0.progress_apply(remove_tiny_shapes,axis=1)

    # # simplify geometries
    # gadm_level0['geometry'] = gadm_level0.simplify(tolerance = 0.005, preserve_topology=True).buffer(0.01).simplify(tolerance = 0.005, preserve_topology=True)
    
    # #save to new country file
    # glob_ctry_path = os.path.join(data_path,'cleaned_shapes','global_countries.gpkg')
    # gadm_level0.to_file(glob_ctry_path, layer='countries', driver="GPKG")
     
    # poly_files(data_path,glob_ctry_path) 