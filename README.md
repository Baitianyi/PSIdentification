# PSIdentification
## 1. Introduction
PSIdentification automatically identify the grade, elevation and geographic reach of planation surface using digial elevation data (DEM) and vector profile lines data.

Input data: DEM, vector profile lines data (in shapefile format), buffer width and window size.

Outputdata: the grade and elevation of planation surface (in TXT format), the elevation profile charts and elevation-area statistics charts (in JPG format), the range of tectonic uplift (in shapefile format), the geographic reach of remnant planation surface (in grid-DEM format).

Arcgis Engine 10.1 (or above) and python 2.x is needed when running this project.

Next part aims to intoduce the computing environment.

## 2. Environment
### 2.1 Hardware environment
PSIdentification was run on a computer with 2 cores (1.8 GHz each) and 4 GB RAM, which is a basic running requirement.

### 2.2 Set up computing environment
(1) PSIdentification was interpreted with Python 2.7 and Arcpy (ver 10.1 or above) site-package. Anaconda is recommended for installing Python. Arcgis Engine 10.1 (or above) is needed when running this project.

(2) Install Python GDAL (https://pypi.org/project/GDAL/).

(3) If you want to use your own data, put your DEM file and vector line file into C:\Users\l\Desktop\Identification or change the DEM and line paths in PSIdentification_main.

(4) the buffer width and window size parameters are set in PSIdentification_main.

## 3. Packages
### 3.1 PSIdentification_main
PSIdentification_main is a main project that needs to be run.

### 3.2 peakextraction
Peakextraction is a module which extracts peak points from DEM.

The output peak point data (peak_dem.shp) and attribute csv data (peak.csv) is stored in the peakdata folder.

### 3.3 weightedThiessenPolygons
WeightedThiessenPolygons is a module which constructs Weighted Thiessen polygon.

It was created by Colin Spohn, James Madison University

(Spohn, C. (2015). Creation of weighted Thiessen polygons using Python. Senior Honors Projects, 2010-current. 96. https://commons.lib.jmu.edu/honors201019/96)

and Modified by the authors.

The output weighted Thiessen polygon data (WV.shp) is stored in the Polygons folder.

### 3.4 generation
Generation is a module which generates Buffer areas, constructs original geomorphic digital elevation model (OGDEM) and draws the charts.
The output buffer (buffer.shp) , OGDEM (OGDEM.tif) are stored under the OGDEM folder. The output charts (x.jpg) are stored under the charts folder.

### 3.5 identification
Identification is a module which identifies the grade, elevation and geographic reach of planation surface.
The output grade, height results (level_elevation.txt) and the geographic reach of planation surface data (uplift.shp; remnantPS.tif) are stored in the Results folder.
