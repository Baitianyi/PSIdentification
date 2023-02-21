# -*- coding: utf-8 -*-
# Author: Tianyi Bai; Anbo Li; Xianli Xie; Jiayuan Cai
# Last update: February 21 2023

# peakextraction
# Implementation peak points extraction

# import
import os
import arcpy
import pandas as pd
import copy

def extraction(DEM, windowsize, outpath):
    # 1. Fill
    fill = arcpy.sa.Fill(DEM)
    
    # 2. FocalStatistics
    neighborhood = arcpy.sa.NbrRectangle(windowsize, windowsize, "CELL")
    Maxpoint = arcpy.sa.FocalStatistics(fill, neighborhood, "MAXIMUM", "DATA")
    
    # 3. 3D Analyst
    arcpy.CheckOutExtension("3D")
    restercalc_dem = outpath + r"\restercalc_dem.tif"
    arcpy.Minus_3d(Maxpoint,DEM,restercalc_dem)
    
    # 4. Reclass
    RemapValue = arcpy.sa.RemapValue([[0,1]])
    # Reclassification:value 0 to 1
    outReclass  = arcpy.sa.Reclassify(restercalc_dem, "Value", RemapValue, "NODATA")
    
    # 5. Extract values to points
    Polygons = outpath + r"\polys.shp"
    peak11 = outpath + r"\peaks.shp"
    peak = outpath + r"\peak_dem.shp"
    arcpy.RasterToPolygon_conversion(outReclass, Polygons, "NO_SIMPLIFY", "Value")
    arcpy.management.FeatureToPoint(Polygons, peak11, "INSIDE")
    
    arcpy.sa.ExtractValuesToPoints(peak11, DEM, peak, "NONE", "VALUE_ONLY")
    
    arcpy.AddXY_management(peak)
    
    os.remove(restercalc_dem)
    
    return peak

# output attribute table to csv
def SaveShpAsCSV(ShpFile,outcsv):
    fieldList = []
    list1 = []
    data = []
    desc = arcpy.Describe(ShpFile)
    for field in desc.fields:
        fieldList.append(str(field.name))
        cursor = arcpy.SearchCursor(ShpFile)
    for row in cursor:
        for field in fieldList:
            list1.append(row.getValue(str(field)))
        data.append(copy.copy(list1))
        list1[:] = []
    
    data = list(data)
    columns = list(fieldList)
    file_data = pd.DataFrame(data=data, columns=columns)
    file_data.to_csv(outcsv)
    
    