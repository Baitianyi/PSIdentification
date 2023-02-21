# -*- coding: utf-8 -*-
# Author: Tianyi Bai; Anbo Li; Xianli Xie; Jiayuan Cai
# Last update: February 21 2023

# Import
import os, sys
#import shutil
from osgeo import gdal

arcpy_path = [r'D:\Program Files (x86)\ArcGIS\Desktop10.2\arcpy',
              r'D:\Program Files (x86)\ArcGIS\Desktop10.2\bin']
sys.path.extend(arcpy_path)

import arcpy
from utility import peakextraction, weightedThiessenPolygons, generation, identification

# Create the file
def mkdir(path):
    isExists = os.path.exists(path)
    if not isExists:
        os.makedirs(path)
    else:
        #shutil.rmtree(path)
        #os.makedirs(path)
        return False

def main():
    # Input the file
    path = r"C:\Users\l\Desktop\Identification"
    windowsize = 11
    width = 1000
    
    # Input the DEM
    DEM = path + r'\DEM.tif'
    
    # Check SRS
    sr = arcpy.Describe(DEM).spatialReference
    print("Spatial Reference System:", sr.name)
    # Check license
    print("Spatial Analyst Extension Available:")
    print(arcpy.CheckOutExtension("spatial"))
    
    # Create files for output
    outpath1 = path + r"\peakdata"
    peak = outpath1 + r"\peak_dem.shp"
    outcsv = outpath1 + r"\peak.csv"
    
    outpath2 = path + r"\OGDEM"
    bufferpath = outpath2 + r'\buffer'
    OGDEM = outpath2 + r"\OGDEM.tif"
    weightedpoly = path + r"\OGDEM\polyvalue.shp"
    imgpath = path + r'\charts'
    
    results = path + r"\Results"
    text = results + r"\level_elevation.txt"
    remnant = results + r"\remnant"
    
    mkdir(outpath1)
    mkdir(outpath2)
    mkdir(imgpath)
    mkdir(results)
    
    # extracte peak points
    peak = peakextraction.extraction(DEM, windowsize, outpath1)
    
    # output attribute table to csv
    peakextraction.SaveShpAsCSV(peak, outcsv)
    
    # Create Weighted Thiessen Polygons
    weightedThiessenPolygons.CreateWeighted(path)
    
    # Generate OGDEM
    generation.OGDEM(DEM, path, outpath2)
    
    # If no profile line data is input, the charts drawing algorithm is not run
    if os.path.exists(path + r'\profile_lines.shp'):
        line = path + r'\profile_lines.shp'
        # Generate buffer
        x1s, x2s, y1s, y2s = generation.create_poly(line, bufferpath, width)
        
        # Elevation profile chart and Elevation-area statistics chart
        h = generation.get_chart(x1s, x2s, y1s, y2s, DEM, OGDEM, imgpath)
        
        # Identify the grade and elevation of planation surface
        identification.clustering(h, results)
        
    else:
        # Identify the grade and elevation of planation surface
        ds_dem = gdal.Open(OGDEM)
        if ds_dem is None:
            print('Cannot find dem')
            sys.exit(1)
        
        geotrans = ds_dem.GetGeoTransform()
        pixelWidth = geotrans[1]
        pixelHeight = geotrans[5]
        i = -1
        
        # Elevation-area statistics chart
        height = generation.get_height(OGDEM)
        generation.get_statistics(height, pixelWidth, pixelHeight, imgpath, i)
        
        # Identify the grade and elevation of planation surface
        identification.clustering(height, results)
    
    # Identify the range of tectonic uplift
    identification.upliftrange(text, weightedpoly, results)
    
    # Identify the geographic reach of remnant planation surface
    identification.remnantPS(text, DEM, remnant)
    
if __name__ == "__main__":
    main()