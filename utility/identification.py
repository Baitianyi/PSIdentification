# -*- coding: utf-8 -*-
# Author: Tianyi Bai; Anbo Li; Xianli Xie; Jiayuan Cai
# Last update: February 21 2023

# identification
# Implementation the grade, elevation and geographic reach identification of planation surface

# Import
import os, sys
import arcpy
import gdal
import math
import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
import re

# Numbers are not displayed in scientific notation
np.set_printoptions(suppress=True)

# Select the optimal cluster number
def K_SSE(loan, outpath):
    # Define SSE list
    sse_list = [ ]

    num = KMeans(n_clusters = 1)
    i = 0
    sse = num.fit(loan).inertia_
    while sse > 10:
        sse = sse / 10
        i = i+1
    for k in range(1,15):
        # set the classification number k
        classify = KMeans(n_clusters = k)
        classify = classify.fit(loan)
        # Obtain SSE
        sse_list.append(classify.inertia_)
        if k > 1:
            j = k - 1
            angle = math.atan2((sse_list[j-1]-sse_list[j])/pow(10, i), 1)/math.pi*180
            if angle < 5:
                classify = KMeans(n_clusters = k+1)
                classify = classify.fit(loan)
                sse_list.append(classify.inertia_)
                break
    config = {"font.size": 10,
            "mathtext.fontset": 'stix',
            "font.weight": 'black'}
    
    plt.rcParams.update(config)
    plt.rcParams['axes.unicode_minus'] = False
    plt.figure()
    plt.plot(np.arange(1, len(sse_list)+1).astype(dtype=np.str), sse_list, color='black', marker = 'v', markersize = '5')
    plt.xlabel('Number of clustering', fontsize = 15)
    plt.ylabel('SSE', fontsize = 15)
    plt.tick_params(labelsize=13)
    for a,b in zip(np.arange(1, len(sse_list)+1).astype(dtype=np.str),sse_list):
        plt.text(a, b+0.001, '%.2f' % (b/pow(10,i)), ha='center', va= 'bottom',fontsize=13)
    plt.grid(alpha=0.5,linestyle='--')
    plt.savefig(outpath + r'\SSE.jpg', dpi=500)
    plt.show()
    
    return j

# Identify the grade and elevation of planation surface
def clustering(height, outpath):
    loan = np.array([height]).T
    
    # Obtain the optimal cluster number
    k = K_SSE(loan, outpath)
    
    # set the cluster number k
    classify = KMeans(n_clusters = k)
    
    classify = classify.fit(loan)
    
    rank = []
    labels = classify.labels_
    for i in range(k):
        index = np.where(labels == i)
        rank.append(max(loan[index])[0])
    rank.sort()
    
    # Print the identification results of the grade and elevations
    txtpath = outpath + r"\level_elevation.txt"
    if os.path.exists(txtpath):
        os.remove(txtpath)
    text = open(outpath + r"\level_elevation.txt", "w")
    text.write('Level,Elevation')
    for i in range(len(rank)):
        h = (int(rank[i]/100)+1)*100
        text.write('\n')
        text.writelines(str(i+1)+','+str(h))
    text.close()
    
# Identify the range of tectonic uplift
def upliftrange(textpath, weightedpoly, outpath):
    upliftpath = outpath + r"\uplift.shp"
    text = open(textpath, "r")
    # Read the first line
    grade = 0
    height = []
    line = text.readline()
    line = text.readline()
    while line:
        parts = re.split(',', line)
        height.append(int(parts[1]))
        grade = grade + 1
        line = text.readline()
    text.close()
    height = sorted(height)
    # generate the tectonic uplift of each grade
    inputs = []
    arcpy.MakeFeatureLayer_management(weightedpoly, "poly_lyr")
    for i in range(grade):
        if i == 0:
            arcpy.SelectLayerByAttribute_management("poly_lyr", 'NEW_SELECTION', '"RASTERVALU" <= '+str(height[i]))
        elif i == (grade-1):
            arcpy.SelectLayerByAttribute_management("poly_lyr", 'NEW_SELECTION', '"RASTERVALU" > '+str(height[i-1]))
        else:
            arcpy.SelectLayerByAttribute_management(
                    "poly_lyr", 'NEW_SELECTION', '"RASTERVALU" > '+str(height[i-1])+' AND "RASTERVALU" <= '+str(height[i]))
        
        # Write the selected features to a new feature class
        dissolve = outpath + r'\uplift'+chr(i+65) + r".shp"
        if os.path.exists(dissolve):
            os.remove(dissolve)
        arcpy.management.Dissolve("poly_lyr", dissolve)
        arcpy.AddField_management(dissolve, "Level", "LONG")
        arcpy.CalculateField_management(dissolve, "Level", str(i+1))
        arcpy.AddField_management(dissolve, "Elevation", "LONG")
        arcpy.CalculateField_management(dissolve, "Elevation", str(height[i]))
        inputs.append(dissolve)
    
    if os.path.exists(upliftpath):
        os.remove(upliftpath)
    arcpy.Merge_management(inputs, upliftpath)
    
# Identify the geographic reach of remnant planation surface
def remnantPS(textpath, DEM, remnantpath):
    text = open(textpath, "r")
    
    # Read the first line
    level = 0
    height = []
    line = text.readline()
    line = text.readline()
    while line:
        parts = re.split(',', line)
        height.append(int(parts[1]))
        level = level + 1
        line = text.readline()
    text.close()
    height = sorted(height)
    
    for i in range(level):
        if i == 0:
            attExtract = arcpy.sa.ExtractByAttributes(DEM, 'VALUE <= '+str(height[i]))
        elif i == (level-1):
            attExtract = arcpy.sa.ExtractByAttributes(DEM, 'VALUE > '+str(height[i-1]))
        else:
            attExtract = arcpy.sa.ExtractByAttributes(DEM, 'VALUE > '+str(height[i-1])+' AND VALUE <= '+str(height[i]))
        
        path = remnantpath + 'PS' + str(i+1) + '.tif'
        attExtract.save(path)