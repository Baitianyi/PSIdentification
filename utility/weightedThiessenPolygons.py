# -*- coding: utf-8 -*-
# Created by Colin Spohn, James Madison University
# Modified by Tianyi Bai, Anbo Li, Xianli Xie and Jiayuan Cai
# Last update: February 21 2023

# weightedThiessenPolygons
# Implementation Weighted Thiessen polygon construction

# Import
import os
from scipy.spatial import Delaunay
import arcpy
import re
import math
from operator import itemgetter

# Get the intersection of lines
def lineIntersection(line1, line2):
    xdiff = (line1[0][0] - line1[1][0], line2[0][0] - line2[1][0])
    ydiff = (line1[0][1] - line1[1][1], line2[0][1] - line2[1][1])
    
    def det(a, b):
        return a[0] * b[1] - a[1] * b[0]
    
    div = det(xdiff, ydiff)
    if div == 0:
        return []
    
    # Define the intersection coordinate x, y
    d = (det(*line1), det(*line2))
    x = (det(d, xdiff) / div)
    y = (det(d, ydiff) / div)
    
    # Check if intersect is on both segments
    # Check line 1 x values
    if round(line1[0][0], 1) == round(line1[1][0], 1):
        x1Check = round(line1[0][0], 1) == round(x, 1)
    elif line1[0][0] <= line1[1][0]:
        x1Check = line1[0][0] <= x <= line1[1][0]
    else :
        x1Check = line1[1][0] <= x <= line1[0][0]
        
    # Check line 2 x values
    if round(line2[0][0], 1) == round(line2[1][0], 1):
        x2Check = round(line2[0][0], 1) == round(x, 1)
    elif line2[0][0] <= line2[1][0]:
        x2Check = line2[0][0] <= x <= line2[1][0]
    else:
        x2Check = line2[1][0] <= x <= line2[0][0]
        
    # Check line 1 y values
    if round(line1[0][1], 1) == round(line1[1][1], 1):
        y1Check = round(line1[0][1], 1) == round(y, 1)
    elif line1[0][1] <= line1[1][1]:
        y1Check = line1[0][1] <= y <= line1[1][1]
    else:
        y1Check = line1[1][1] <= y <= line1[0][1]
        
    # Check line 2 y values
    if round(line2[0][1], 1) == round(line2[1][1], 1):
        y2Check = round(line2[0][1], 1) == round(y, 1)
    elif line2[0][1] <= line2[1][1]:
        y2Check = line2[0][1] <= y <= line2[1][1]
    else:
        y2Check = line2[1][1] <= y <= line2[0][1]
        
    if (x1Check and x2Check and y1Check and y2Check):
        return [x, y]
    else:
        return []

# GeT polygons' boundary
def whittle(point, edges):
    points = []
    
    # Add initial line
    firstPoint = (edges[0][0], edges[0][1])
    secondPoint = (edges[0][4], edges[0][5])
    
    points.append(firstPoint)
    points.append(secondPoint)
    
    # Cycle through remaining lines
    i = 1
    n = len(edges)
    while i < n:
        # Create lines left and right of midpoint
        point1 = (edges[i][0], edges[i][1])
        point2 = (edges[i][2], edges[i][3])
        point3 = (edges[i][2], edges[i][3])
        point4 = (edges[i][4], edges[i][5])
        
        line1 = []
        line2 = []
        line1.append(point1)
        line1.append(point2)
        line2.append(point3)
        line2.append(point4)
        
        leftEndPoint = [edges[i][0], edges[i][1], len(points), 100000]
        rightEndPoint = [edges[i][4], edges[i][5], len(points), 100000]
        
        leftIntersects = []
        rightIntersects = []
        # Edges
        leftIntersects.append(leftEndPoint)
        rightIntersects.append(rightEndPoint)
        
        # Check new line against each existing line segment
        j = 0
        m = len(points) - 1
        while j < m:
            intersect1 = []
            intersect2 = []
            
            point5 = (points[j][0], points[j][1])
            point6 = (points[j+1][0], points[j+1][1])
            
            line3 = []
            line3.append(point5)
            line3.append(point6)
            
            intersect1 = lineIntersection(line1, line3)
            intersect2 = lineIntersection(line2, line3)
            
            if len(intersect1) > 0:
                # Calculate the distance between intersect and midpoint
                x = edges[i][2] - intersect1[0]
                y = edges[i][3] - intersect1[1]
                length = math.sqrt(math.pow(x, 2) + math.pow(y, 2))
                intersect1.append(j+1)
                intersect1.append(length)
                leftIntersects.append(intersect1)
            
            if len(intersect2) > 0:
                x = edges[i][2] - intersect2[0]
                y = edges[i][3] - intersect2[1]
                length = math.sqrt(math.pow(x, 2) + math.pow(y, 2))
                intersect2.append(j+1)
                intersect2.append(length)
                rightIntersects.append(intersect2)  
            j = j + 1
            
        # Repeat for last point to first point
        intersect1 = []
        intersect2 = []
        point5 = (points[len(points)-1][0], points[len(points)-1][1])
        point6 = (points[0][0], points[0][1])
        
        line3 = []
        line3.append(point5)
        line3.append(point6)
        
        intersect1 = lineIntersection(line1, line3)
        intersect2 = lineIntersection(line2, line3)
        
        if len(intersect1) > 0:
            x = edges[i][2] - intersect1[0]
            y = edges[i][3] - intersect1[1]
            length = math.sqrt(math.pow(x, 2) + math.pow(y, 2))
            intersect1.append(len(points))
            intersect1.append(length)
            leftIntersects.append(intersect1)
            
        if len(intersect2) > 0:
            x = edges[i][2] - intersect2[0]
            y = edges[i][3] - intersect2[1]
            length = math.sqrt(math.pow(x, 2) + math.pow(y, 2))
            intersect2.append(len(points))
            intersect2.append(length)
            rightIntersects.append(intersect2)
            
        # Get shortest intersecting segment on top
        leftIntersects = sorted(leftIntersects, key=itemgetter(3))
        rightIntersects = sorted(rightIntersects, key=itemgetter(3))
        
        # Add left intersect
        leftPoint = (leftIntersects[0][0], leftIntersects[0][1])
        points.insert(leftIntersects[0][2], leftPoint)
        
        # Add right intersect
        rightPoint = (rightIntersects[0][0], rightIntersects[0][1])
        if leftIntersects[0][2] <= rightIntersects[0][2]:
            points.insert(rightIntersects[0][2]+1, rightPoint)
        else :
            points.insert(rightIntersects[0][2], rightPoint)
            
        # Add new line if only one intersect
        if len(leftIntersects) == 1 and len(rightIntersects) > 1:
            points.remove(leftPoint)
            points.insert(rightIntersects[0][2]+1, leftPoint)
            
        # Add new line if only one intersect
        if len(rightIntersects) == 1 and len(leftIntersects) > 1:
            points.remove(rightPoint)
            points.insert(leftIntersects[0][2]+1, rightPoint)
            
        # Remove
        j = 0
        m = len(points)
        while j < m:
            if (points[j] != leftPoint) and (points[j] != rightPoint):
                newLine = [leftPoint, rightPoint]
                midLine = [(point[0], point[1]), (points[j][0], points[j][1])]
                
                checkIntersect = []
                checkIntersect = lineIntersection(newLine, midLine)
                if len(checkIntersect) > 0:
                    points.remove(points[j])
                    j = j - 1
                    m = m - 1
            j = j + 1
            
        i = i + 1
        
    pointsCleaned = []
    checked = set()
    for pair in points:
        if pair not in checked:
            pointsCleaned.append([pair[0], pair[1]])
            checked.add(pair)
            
    return pointsCleaned

def sort_angle(center, poi):
    point = []
    polygon = []
    for i in range(len(poi)):
        yd = poi[i][1] - center[1]
        xd = poi[i][0] - center[0]
        if xd == 0:
            if yd >0:
                angle = 90
            else:
                angle = 270
        else:
            if yd > 0:
                angle = math.atan2(yd, xd) / math.pi * 180
            elif yd == 0:
                if xd > 0:
                    angle = 0
                else:
                    angle = 180
            else:
                angle = math.atan2(yd, xd) / math.pi * 180 + 360
        point.append([poi[i][0], poi[i][1], angle])
    point = sorted(point, key=itemgetter(2), reverse=True)
    for x, y, a in point:
        polygon.append([x, y])
    return polygon

def midpoints(coord, point):
    # Create values to be used in loops
    i = 0
    n = len(coord)
    newCoords = []
    
    while i < n:
        x = point[0] - coord[i][0] # center x minus other x
        y = point[1] - coord[i][1] # center y minus other y
        
        midX = coord[i][0] + ((coord[i][2]*x)/(point[2]+coord[i][2]))
        midY = coord[i][1] + ((coord[i][2]*y)/(point[2]+coord[i][2]))
        
        x2 = point[0] - midX
        y2 = point[1] - midY
        
        length = math.sqrt(math.pow(x2, 2) + math.pow(y2, 2))
        
        slope, angle = getSlope(point[0], point[1], coord[i][0], coord[i][1])
        
        if slope == 9999999999999999.0:
            perpendicular = 0
        elif slope == 0:
            perpendicular = 9999999999999999.0
        else :
            perpendicular = -1.0/slope
            
        midPoint = (midX, midY, perpendicular, length, angle)
        newCoords.append(midPoint)
        i = i + 1
    return newCoords

def getSlope(x1, y1, x2, y2):
    yd = y2 - y1
    xd = x2 - x1
    if x2 == x1:
        slope = 9999999999999999.0
        if yd >0:
            angle = 90
        else:
            angle = 270
    else:
        slope = (y2-y1)/(x2-x1)
        if yd > 0:
            angle = math.atan2(yd, xd) / math.pi * 180
        elif yd == 0:
            if xd > 0:
                angle = 0
            else:
                angle = 180
        else:
            angle = math.atan2(yd, xd) / math.pi * 180 + 360
    
    return slope, angle

def createLines(coords, point):
    midPoints = midpoints(coords, point)
    edges = []
    
    # Create values to be used in loops
    i = 0
    n = len(midPoints)
    
    # Take each midpoint and turn it into a possible edge line
    while i < n:
        line = extend(midPoints[i][0], midPoints[i][1], midPoints[i][2], midPoints[i][3], midPoints[i][4])
        edges.append(line)
        i = i + 1
        
    return edges

# Extend a line, from a single point and slope, in both directions
def extend(x, y, slope, distance, angle):
    if slope < 0:
        positiveX = (math.sqrt(10000000000.0/(1+math.pow(slope, 2)))-x)*(-1)
        positiveY = (-math.sqrt(10000000000.0-math.pow(x-positiveX, 2))-y)*(-1)
        negativeX = x + (x - positiveX)
        negativeY = (y + (y - positiveY))
        line = [positiveX, positiveY, x, y, negativeX, negativeY, angle, distance]
    elif slope == 9999999999999999.0:
        positiveX = x
        positiveY = y + 100000
        negativeX = x
        negativeY = y - 100000
        line = [positiveX, positiveY, x, y, negativeX, negativeY, angle, distance]
    elif slope ==0:
        positiveX = x + 100000
        positiveY = y
        negativeX = x - 100000
        negativeY = y
        line = [positiveX, positiveY, x, y, negativeX, negativeY, angle, distance]
    else:
        positiveX = (math.sqrt(10000000000.0/(1+math.pow(slope, 2)))-x)*(-1)
        positiveY = (math.sqrt(10000000000.0-math.pow(x-positiveX, 2))-y)*(-1)
        negativeX = x + (x - positiveX)
        negativeY = (y + (y - positiveY))
        line = [positiveX, positiveY, x, y, negativeX, negativeY, angle, distance]
    return line

def CreateWeighted(path):
    # Input points
    csv = path + r"\peakdata\peak.csv"
    fo = open(csv, "r")
        
    # Create a file geodatabase.
    outpath = path + r"\Polygons"
    if os.path.exists(outpath):
        print('exist')
    else:
        os.mkdir(outpath)

    # Create geodatabase tables to build the points from.
    GDB = path + r"\Polygons"
    arcpy.CreateTable_management(GDB, "Coordinates")
    
    # Define the list to store coordinates
    coord = []
    p = []
    
    # Read the first line
    line = fo.readline()
    line = fo.readline()
    parts = re.split(',', line)
    coord.append([float(parts[7]), float(parts[8]), float(parts[6])])
    p.append([float(parts[7]), float(parts[8])])
    
    # Define the range of the points
    xmin = float(parts[7])
    xmax = float(parts[7])
    ymin = float(parts[8])
    ymax = float(parts[8])
    
    line = fo.readline()
    
    # Loop through the text file adding new sets of coordinates to the list
    while line:
        parts = re.split(',', line)
        coord.append([float(parts[7]), float(parts[8]), float(parts[6])])
        p.append([float(parts[7]), float(parts[8])])
        if xmin > float(parts[7]):
            xmin = float(parts[7])
        elif xmax < float(parts[7]):
            xmax = float(parts[7])
        if ymin > float(parts[8]):
            ymin = float(parts[8])
        elif ymax < float(parts[8]):
            ymax = float(parts[8])
            
        line = fo.readline()
        
    # Create TIN
    tri = Delaunay(p)
    
    # Read all neibour points of a peak point
    neib = []
    l = tri.vertex_neighbor_vertices
    for i in range(len(l[0])-1):
        neib.append(list(l[1][l[0][i]:l[0][i+1]]))
        
    # Create values to be used in loops
    polygons = []
    polygonShapes = []
    
    for i in range(len(coord)):
        extraCoords = []
        for n in neib[i]:
            extraCoords.append(coord[n])
        edges = createLines(extraCoords, coord[i])
        
        edges = sorted(edges, key=itemgetter(7), reverse=True)
        
        # Whittle
        poly_poi = whittle(coord[i], edges)
        
        # Sort by angle
        polygon = sort_angle(coord[i], poly_poi)
        polygons.append(polygon)
    
    for shape in polygons:
        polygonShapes.append(
                arcpy.Polygon(
                        arcpy.Array([arcpy.Point(*eachPair) for eachPair in shape])))
        
    weight = GDB + r"/weight"
    arcpy.CopyFeatures_management(polygonShapes,weight)
    arcpy.management.DefineProjection(weight + r".shp", "3857")
    
    # Clip
    clip = arcpy.Polygon(arcpy.Array(
            [arcpy.Point(xmin, ymin),arcpy.Point(xmax, ymin),
             arcpy.Point(xmax, ymax),arcpy.Point(xmin, ymax)]))
    
    arcpy.Clip_analysis(GDB + r"/weight.shp", clip, GDB + r"/WV.shp")
    
    