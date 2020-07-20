# OpenStreetMap clipper
Small python tool to clip pre-defined areas from large osm.pbf files, using [osmconvert](https://wiki.openstreetmap.org/wiki/Osmconvert) or [osmosis](https://wiki.openstreetmap.org/wiki/Osmosis). 

The tool is at the moment specifically written to split the OpenStreetMap planet file into countries and regions, based on the [GADM36](https://gadm.org/) classification. Future versions will allow for user-defined input as well.  

**Please note:** This package is still in development phase. In case of any problems, or if you have any suggestions for improvements, please raise an *issue*. 

## Installation

1. Clone the repository or download the package on your computer and extract the folder.
2. Go to the osm_clipper folder in your command prompt or bash.
3. Type ``python setup.py install`` and it should install itself into your python environment.
4. Now you can import the package like any other package!

## Dependency on osmcovert or osmosis
In essence, *osm_clipper* is just a wrapper around **osmconvert** and **osmosis**. To make it work, make sure you add either osmconvert or osmosis to your environmental variables. 

### How to set environmental variables on Windows:
1. Search for *Edit the System Environmental Variables* and click on it.
2. Click on the *Environmental Variables* button.
3. Search in your System variables for **Path**.
4. Click on the *New* button and add the location of **osmconvert** and/or **osmosis**.

NOTE: the default is **osmconvert**. If you want to use osmosis, specify **osmconvert==False** in the *single_countries* or *all_countries* functions.

### License
Copyright (C) 2020 Elco Koks. All versions released under the [MIT license](LICENSE).


![IVM](http://ivm.vu.nl/en/Images/IVM_logo_rgb2_tcm234-851594.svg)