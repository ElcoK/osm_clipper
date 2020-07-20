import os,sys
sys.path.append(os.path.join('..','src'))
from osm_clipper import planet_osm,single_country,global_shapefiles,gadm36_planet,poly_files,all_countries,gadm36_country

if __name__ == '__main__':

    data_path = os.path.join('..','data')

    #download planet file
    #planet_osm(data_path)

    #download GADM file
    #gadm36_planet(data_path)

    #get a cleaned global shapefile of the countries
    #global_shapefiles(data_path,regionalized=True)
    single_country('SRB',data_path,regionalized=True,geofabrik=True)

    # create poly files
    #global_shape = os.path.join(data_path,'cleaned_shapes','global_regions.gpkg')
    #poly_files(data_path,global_shape,regionalized=True)

    #run specified country
    #all_countries()