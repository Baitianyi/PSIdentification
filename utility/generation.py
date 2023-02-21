# -*- coding: utf-8 -*-
# Author: Tianyi Bai; Anbo Li; Xianli Xie; Jiayuan Cai
# Last update: February 21 2023

# generation
# Implementation Buffer generation, OGDEM construction and charts drawing

import sys
import arcpy
import copy
import math
import gdal
import numpy as np
import matplotlib.pyplot as plt
from collections import Counter

# Read profile lines
def read_line(feature):
    line = feature.geometry()
    #获取矢量线点数
    #print(line)
    #print(line.GetPointCount())
    n = line.GetPointCount()
    #获取折线上各转折点坐标
    pois = []
    for i in range(n):
        pois.append([line.GetX(i),line.GetY(i)])
    return pois

# Read boundary points of the buffers
def bufferpoint(pois, buffer_dis):
    x1 = []
    x2 = []
    y1 = []
    y2 = []
    for i in range(len(pois)-1):
        # First point
        if i == 0:
            if pois[i][1] != pois[i+1][1]:
                k = (pois[i+1][0]-pois[i][0])/(pois[i][1]-pois[i+1][1])
                # addition of horizontal and vertical coordinates
                addx1 = math.sqrt(buffer_dis*buffer_dis/(1+k*k))
                addy1 = k * addx1
                
                if (i+1) != (len(pois)-1):
                    angle1 = math.atan2(pois[i][1]-pois[i+1][1], pois[i][0]-pois[i+1][0])
                    angle2 = math.atan2(pois[i+2][1]-pois[i+1][1], pois[i+2][0]-pois[i+1][0])
                    angle = (angle1+angle2)/2
                    if angle != (math.pi/2):
                        k = math.tan(angle)
                        buffer1 = abs(buffer_dis/math.sin(angle1-angle))
                        addx2 = math.sqrt(buffer1*buffer1/(1+k*k))
                        addy2 = k * addx2
                    else:
                        addx2 = addx1
                        addy2 = addy1
                elif (i+1) == (len(pois)-1):
                    addx2 = addx1
                    addy2 = addy1
                    
            elif pois[i][1] == pois[i+1][1]:
                addx1 = buffer_dis
                addy1 = 0
                if (i+1) != (len(pois)-1):
                    angle1 = math.atan2(pois[i][1]-pois[i+1][1], pois[i][0]-pois[i+1][0])
                    angle2 = math.atan2(pois[i+2][1]-pois[i+1][1], pois[i+2][0]-pois[i+1][0])
                    angle = (angle1+angle2)/2
                    if angle != (math.pi/2):
                        k = math.tan(angle)
                        buffer1 = abs(buffer_dis/math.sin(angle1-angle))
                        addx2 = math.sqrt(buffer1*buffer1/(1+k*k))
                        addy2 = k * addx2
                    else:
                        addx2 = addx1
                        addy2 = addy1
                elif (i+1) == (len(pois)-1):
                    addx2 = addx1
                    addy2 = addy1
        # middle point           
        elif i != 0:
            if (i+1) != (len(pois)-1):
                angle1 = math.atan2(pois[i][1]-pois[i+1][1], pois[i][0]-pois[i+1][0])
                angle2 = math.atan2(pois[i+2][1]-pois[i+1][1], pois[i+2][0]-pois[i+1][0])
                angle = (angle1+angle2)/2
                if angle != (math.pi/2):
                    k = math.tan(angle)
                    buffer1 = abs(buffer_dis/math.sin(angle1-angle))
                    addx2 = math.sqrt(buffer1*buffer1/(1+k*k))
                    addy2 = k * addx2
                else:
                    k = (pois[i+1][0]-pois[i][0])/(pois[i][1]-pois[i+1][1])
                    addx2 = math.sqrt(buffer_dis*buffer_dis/(1+k*k))
                    addy2 = k * addx2
            elif (i+1) == (len(pois)-1):
                k = (pois[i+1][0]-pois[i][0])/(pois[i][1]-pois[i+1][1])
                addx2 = math.sqrt(buffer_dis*buffer_dis/(1+k*k))
                addy2 = k * addx2
            
        # Output the coordinates
        if i == 0:
            x1.append(pois[i][0]+addx1)
            x2.append(pois[i][0]-addx1)
            y1.append(pois[i][1]+addy1)
            y2.append(pois[i][1]-addy1)
                
        x1.append(pois[i+1][0]+addx2)
        x2.append(pois[i+1][0]-addx2)
        y1.append(pois[i+1][1]+addy2)
        y2.append(pois[i+1][1]-addy2)
    return x1, x2, y1, y2

# Original geomorphic digital elevation model 
def OGDEM(DEM, path, outpath):
    peaks = path + r"\peakdata\peak_dem.shp"
    polys = path + r"\Polygons\WV.shp"
    # Get polygons with attributes
    outfc1 = outpath + r"\polyvalue.shp"
    arcpy.SpatialJoin_analysis(polys, peaks, outfc1, "JOIN_ONE_TO_ONE", "KEEP_ALL")
    
    # Get generate OGDEM
    outfc2 = outpath + r"\OGDEM.tif"
    arcpy.PolygonToRaster_conversion(outfc1, "RASTERVALU", outfc2, "MAXIMUM_AREA", "", DEM)
    
# Buffer
def create_poly(line, outputfile_buffer, buffer_dis):
    # Define boundary points
    x1s = []
    x2s = []
    y1s = []
    y2s = []
    polygons = []
    polygonShapes = []
    # Read shapefile lines
    for line in arcpy.da.SearchCursor(line, ["SHAPE@"]):
        # 获取组成折线的每个折点
        pois = []
        for part in line[0]:
            for p in part:
                pois.append([p.X, p.Y])
        
        # Get the boundary points of each line
        x1, x2, y1, y2 = bufferpoint(pois, buffer_dis)
        
        polygon = []
        for j in range(len(x1)):
            polygon.append([x1[j], y1[j]])
        for k in range(len(x2)-1,-1,-1):
            polygon.append([x2[k], y2[k]])
            
        polygons.append(polygon)
        x1s.append(x1)
        x2s.append(x2)
        y1s.append(y1)
        y2s.append(y2)
    
    # Output the buffer
    for shape in polygons:
        polygonShapes.append(
                arcpy.Polygon(
                        arcpy.Array([arcpy.Point(*eachPair) for eachPair in shape])))
    
    arcpy.CopyFeatures_management(polygonShapes,outputfile_buffer)
    arcpy.management.DefineProjection(outputfile_buffer + r".shp", "3857")
    return x1s, x2s, y1s, y2s
    
# Sample area elevation statistics
def get_h(band, pixelWidth, pixelHeight, xOrigin, yOrigin, dis, x1, y1, x2, y2):
    # Define the elevation and distance
    height_max = []
    height = []
    dis_max = []
    
    for i in range(len(x1)-1):
        length = math.sqrt((x1[i]-x1[i+1])*(x1[i]-x1[i+1])+(y1[i]-y1[i+1])*(y1[i]-y1[i+1]))
        if abs(pixelWidth) == abs(pixelHeight):
            n = int(length/pixelWidth)
        else:
            n = int(length/abs(max(pixelWidth,pixelHeight)))
        
        X1Values = np.linspace(x1[i],x1[i+1],n)
        Y1Values = np.linspace(y1[i],y1[i+1],n)
        X2Values = np.linspace(x2[i],x2[i+1],n)
        Y2Values = np.linspace(y2[i],y2[i+1],n)
        
        # Get the max elevation value of each row
        for j in range(len(X1Values)):
            length1 = math.sqrt((X1Values[j]-X2Values[j])*(X1Values[j]-X2Values[j])+
                                (Y1Values[j]-Y2Values[j])*(Y1Values[j]-Y2Values[j]))
            if abs(pixelWidth) == abs(pixelHeight):
                n1 = int(length1/pixelWidth)
            else:
                n1 = int(length1/abs(max(pixelWidth,pixelHeight)))
                
            # Line interpolation
            poi_x = np.linspace(X1Values[j],X2Values[j],n1)
            poi_y = np.linspace(Y1Values[j],Y2Values[j],n1)
            
            h = []
            for k in range(len(poi_x)):
                x = poi_x[k]
                y = poi_y[k]
                # Get the distance
                if k != 0:
                    dis1 = dis + abs(max(pixelWidth,pixelHeight))/1000
                else:
                    dis1 = 0
                    
                xset = int((x-xOrigin)/pixelWidth)
                yset = int((y-yOrigin)/pixelHeight)
                
                data = band.ReadAsArray(xset,yset,1,1)
                if data is None:
                    continue
                value = data[0,0]
                if value < 0:
                    continue
                height.append(float(value))
                h.append(float(value))
            dis = dis1
            dis_max.append(dis)
            height_max.append(max(h))
    return height, dis_max, height_max

# The whole area elevation statistics
# Get all raster elevation values
def get_height(DEM):
    height = []
    ds_dem = gdal.Open(DEM)
    if ds_dem is None:
        print('Cannot find dem')
        sys.exit(1)
        
    # Get the number of rows, columns and band
    row = ds_dem.RasterYSize
    column = ds_dem.RasterXSize
    
    data = ds_dem.ReadAsArray(0,0,column,row).astype(np.float32)
    h = np.ravel(data)
        
    for i in range(len(h)):
        if h[i] is None:
                    continue
        value = h[i]
        if value < 0:
            continue
        height.append(float(value))
    return height

# Obtain sample area parameters
def get_sample(DEM, x1, y1, x2, y2):
    ds_dem = gdal.Open(DEM)
    if ds_dem is None:
        print('Cannot find dem')
        sys.exit(1)
    band = ds_dem.GetRasterBand(1)
    
    geotrans = ds_dem.GetGeoTransform()
    pixelWidth = geotrans[1]
    pixelHeight = geotrans[5]
    xOrigin = geotrans[0]
    yOrigin = geotrans[3]
    dis = 0
    
    # The sample area is a rectangle
    if len(x1) == 2:
        height, dis_max, height_max = get_h(
                band, pixelWidth, pixelHeight, xOrigin, yOrigin, dis, x1, y1, x2, y2)
    
    # The sample area is made up of section lines
    if len(x1) > 2:
        height_max = []
        height = []
        dis_max = []
        cx1 = copy.deepcopy(x1)
        cx2 = copy.deepcopy(x2)
        cy1 = copy.deepcopy(y1)
        cy2 = copy.deepcopy(y2)
        for j in range(len(x1)-1):
            partx1 = []
            partx2 = []
            party1 = []
            party2 = []
            if (j + 1) != (len(x1)-1):
                partx1.append(x1[j])
                partx2.append(x2[j])
                party1.append(y1[j])
                party2.append(y2[j])
                d1 = (cx1[j]-cx1[j+1])*(cx1[j]-cx1[j+1])+(cy1[j]-cy1[j+1])*(cy1[j]-cy1[j+1])
                d2 = (cx2[j]-cx2[j+1])*(cx2[j]-cx2[j+1])+(cy2[j]-cy2[j+1])*(cy2[j]-cy2[j+1])
                if d1 > d2:
                    if cx1[j+1] > cx1[j]:
                        x = cx1[j] + math.sqrt(d2/(math.pow(((cy2[j]-cy2[j+1])/(cx2[j]-cx2[j+1])),2)+1))
                    else:
                        x = cx1[j] - math.sqrt(d2/(math.pow(((cy2[j]-cy2[j+1])/(cx2[j]-cx2[j+1])),2)+1))
                    y = (cy2[j]-cy2[j+1])/(cx2[j]-cx2[j+1])*(x-cx1[j])+cy1[j]
                    partx1.append(x)
                    partx2.append(cx2[j+1])
                    party1.append(y)
                    party2.append(cy2[j+1])
                    cx2[j+1] = x
                    cy2[j+1] = y
                elif d1 < d2:
                    if cx2[j+1] > cx2[j]:
                        x = cx2[j] + math.sqrt(d2/(math.pow(((cy1[j]-cy1[j+1])/(cx1[j]-cx1[j+1])),2)+1))
                    else:
                        x = cx2[j] - math.sqrt(d2/(math.pow(((cy1[j]-cy1[j+1])/(cx1[j]-cx1[j+1])),2)+1))
                    y = (cy1[j]-cy1[j+1])/(cx1[j]-cx1[j+1])*(x-cx2[j])+cy2[j]
                    partx1.append(cx1[j+1])
                    partx2.append(x)
                    party1.append(cy1[j+1])
                    party2.append(y)
                    cx1[j+1] = x
                    cy1[j+1] = y
            elif (j + 1) == (len(x1)-1):
                partx1.append(x1[j])
                partx1.append(x1[j+1])
                partx2.append(x2[j])
                partx2.append(x2[j+1])
                party1.append(y1[j])
                party1.append(y1[j+1])
                party2.append(y2[j])
                party2.append(y2[j+1])
            h, dm, hm = get_h(
                    band, pixelWidth, pixelHeight, xOrigin, yOrigin, dis, partx1, party1, partx2, party2)
            height = height + h
            dis_max = dis_max + dm
            height_max = height_max + hm
            dis = max(dis_max)
        
    return height, pixelWidth, pixelHeight, dis_max, height_max

# Drawing profile charts
def get_profile(dis0, height0, dis1, height1, imgpath, i):
    config = {
            "font.family": 'serif',
            "font.size": 10,
            "mathtext.fontset": 'stix',
            "font.weight": 'black',
         }
    plt.rcParams.update(config)
    plt.rcParams['axes.unicode_minus'] = False
    plt.figure()#初始化
    plt.plot(dis0,height0,color='royalblue',label = 'DEM elevation line')
    plt.plot(dis1,height1,color='black',linestyle='--',label = 'OGDEM elevation line')
    plt.grid(alpha=0.5,linestyle='--')
    plt.title('Profile Line '+chr(i+65)+': Elevation Profile Chart', fontsize = 15)
    plt.xlabel('Distance(km)', fontsize = 12)
    plt.ylabel('Elevation(m)', fontsize = 12)
    plt.legend(loc=4)
    plt.savefig(imgpath+'\Profile Line '+chr(i+65)+'-Elevation Profile Chart.jpg', dpi=500)
    plt.show()
    
# Drawing statistics chart
def get_statistics(height, pixelWidth, pixelHeight, imgpath, i):
    c = Counter(height)
    h = []
    area = []
    # Sort
    d = sorted(c.items(), key=lambda x: x[0], reverse=False)
    for j in d:
        if j[0] > 0:
            h.append(float(j[0]))
            area.append(abs(int(j[1])*pixelWidth*pixelHeight/1000000))
    
    plt.figure()#初始化
    plt.bar(h, area, 20, color='royalblue', label = 'area')
    plt.xlabel('Elevation(m)', fontsize = 10)
    plt.ylabel('Area(km^2)', fontsize = 10)
    plt.grid(alpha=0.5,linestyle='--')
    if i == -1:
        plt.title('Elevation-Area Statistics Chart', fontsize = 13)
        plt.savefig(imgpath+'\Elevation-Area Statistics Chart.jpg', dpi=500)
    else:
        plt.title('Sample Area '+chr(i+65)+': Elevation-Area Statistics Chart', fontsize = 13)
        plt.savefig(imgpath+'\Sample Area'+chr(i+65)+'-Elevation-Area Statistics Chart.jpg', dpi=500)
    
    plt.show()
    
# Get charts
def get_chart(x1s, x2s, y1s, y2s, DEM, OGDEM, imgpath):
    h = []
    for i in range(len(x1s)):
        x1 = x1s[i]
        x2 = x2s[i]
        y1 = y1s[i]
        y2 = y2s[i]
        
        height0, pixelWidth0, pixelHeight0, dis_max0, height_max0 = get_sample(DEM, x1, y1, x2, y2)
        
        height, pixelWidth, pixelHeight, dis_max, height_max= get_sample(OGDEM, x1, y1, x2, y2)
        
        get_profile(dis_max0, height_max0, dis_max, height_max, imgpath, i)
        
        get_statistics(height, pixelWidth, pixelHeight, imgpath, i)
        
        h = h + height
        
    return h