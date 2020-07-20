import os,sys
from osm_clipper import single_country,global_shapefiles,gadm36_download,poly_files

if __name__ == '__main__':

    data_path = os.path.join('..','data')

    #download planet file
    planet_osm(data_path)

    #download GADM file
    gadm36_download(data_path)

    #get a cleaned global shapefile of the countries
    global_shapefiles(data_path)

    # create poly files
    global_shape = os.path.join(data_path,'cleaned_shapes','global_countries.gpkg')
    poly_files(data_path,global_shape,regionalized=False)

    #run specified country
    single_country(sys.argv[1])