# -*- coding: utf-8 -*-
"""
Created on Fri Jul  6 13:35:57 2018

@author: mrwan
"""

import numpy as np
from osgeo import gdal    

lkfile="x:/bziniti/MSB_next_idsmethane/stordalen_2018/modis_landsat/landsat/2014221_LC8_197012_ndvi-cloudmask.tif"
mkfile="x:/bziniti/MSB_next_idsmethane/stordalen_2018/modis_landsat/modis30_landsatgrid_B/2014221_MCD_ndvi.tif"
m0file="x:/bziniti/MSB_next_idsmethane/stordalen_2018/modis_landsat/modis30_landsatgrid_B/2014222_MCD_ndvi.tif"

def read_raster(infile,Band):
  fp = gdal.Open(infile)
  band = fp.GetRasterBand(Band)
  data = band.ReadAsArray()
  
  band.GetScale()
  band.GetOffset()
  band.GetNoDataValue()
  return data

def window(wholeTable,x0=0,y0=0,window_dimention=9):
    window=np.zeros([window_dimention,window_dimention],dtype=int)
    row=0
    while(row<window_dimention):
        col=0
        while(col<window_dimention):
            window[row][col]=wholeTable[x0+row][y0+col]
            col+=1
        row+=1
    return window

def filter_wlk(wlk):
    zero_one_filter=np.zeros([wlk.shape[0],wlk.shape[1]])
    #threshold=(np.max(wlk)+np.min(wlk))/2
    center_point_index=int((wlk.shape[0]-1)/2)
    threshold_up=wlk[center_point_index][center_point_index]+500
    threshold_low=wlk[center_point_index][center_point_index]-500
    row=0
    while(row<wlk.shape[0]):
        col=0
        while(col<wlk.shape[1]):
            if wlk[row][col]>threshold_low and wlk[row][col]<threshold_up:
                zero_one_filter[row][col]=1
            col+=1
        row+=1
    return zero_one_filter

def computeDiff(window1,window2):
    diff_window=np.zeros([window1.shape[0],window2.shape[1]])
    row=0
    while(row<window1.shape[0]):
        col=0
        while(col<window1.shape[1]):
            diff_window[row][col]=abs(window1[row][col]-window2[row][col])
            col+=1
        row+=1
    return diff_window

def computeDistance(wlk):
    dist_window=np.zeros([wlk.shape[0],wlk.shape[1]])
    center_x=(wlk.shape[0]+1)/2
    center_y=(wlk.shape[1]+1)/2
    row=0
    while(row<wlk.shape[0]):
        col=0
        while(col<wlk.shape[1]):
            dist_window[row][col] = np.sqrt( (row - center_x)**2 + (col - center_y)**2)
            col+=1
        row+=1
    return dist_window    

def computeCombinedWeight(spec_diff, temp_diff, dist_pixel):
    #combined_pixel=(spec_diff+1)*(temp_diff+1)*dist_pixel
    #combined_pixel=spec_diff*temp_diff*0.01*dist_pixel
    
    B=10000;combined_pixel=np.log(spec_diff*B+1)*np.log(temp_diff+B+1)*dist_pixel
    combined_sum = np.sum(1/combined_pixel[np.nonzero(combined_pixel)])
    
    weight_window = np.ones(spec_diff.shape)
    row = 0
    while(row < spec_diff.shape[0]):
        col = 0
        while(col < spec_diff.shape[1]):
            if combined_pixel[row][col] != 0:
                weight_window[row][col] = (1/combined_pixel[row][col]) / (combined_sum)
            col += 1
        row += 1

    return weight_window


def generateSinglePredictionValue(weigt_window,wlk_filtered,wmk_filtered,wm0_filtered):
    return np.sum(weigt_window*(wm0_filtered + wlk_filtered - wmk_filtered))


Lk=read_raster(lkfile,1);Lk[Lk<-10000]=0
Mk=read_raster(mkfile,1);Mk[Mk<-10000]=0
M0=read_raster(m0file,1);M0[M0<-10000]=0



window_dimention = 9

# =============================================================================
# row,col=10,10
# wlk = window(Lk,row,col,window_dimention)
# wmk = window(Mk,row,col,window_dimention)
# wm0 = window(M0,row,col,window_dimention)
# 
# zero_one_filter = filter_wlk(wlk)
# wmk_filtered = wmk*zero_one_filter
# wlk_filtered = wlk*zero_one_filter
# wm0_filtered = wm0*zero_one_filter
# 
# spec_diff = computeDiff(wlk_filtered,wmk_filtered)
# temp_diff = computeDiff(wmk_filtered,wm0_filtered)
# dist_pixel = computeDistance(wlk_filtered)
# 
# weight = computeCombinedWeight(spec_diff,temp_diff,dist_pixel)
# 
# center_value = generateSinglePredictionValue(weight,wlk_filtered,wmk_filtered,wm0_filtered)
# print(center_value)
# print(Lk[row+2][col+2])
# =============================================================================



result = np.zeros(Lk.shape)
row = 0
while(row < Lk.shape[0] - window_dimention):
    col = 0
    while(col < Lk.shape[1] - window_dimention):
        wlk = window(Lk,row,col,window_dimention)
        wmk = window(Mk,row,col,window_dimention)
        wm0 = window(M0,row,col,window_dimention)
        
        zero_one_filter = filter_wlk(wlk)
        wmk_filtered = wmk*zero_one_filter
        wlk_filtered = wlk*zero_one_filter
        wm0_filtered = wm0*zero_one_filter
        
        spec_diff = computeDiff(wlk_filtered,wmk_filtered)
        temp_diff = computeDiff(wmk_filtered,wm0_filtered)
        dist_pixel = computeDistance(wlk_filtered)
        
        weight = computeCombinedWeight(spec_diff,temp_diff,dist_pixel)
        center_value = generateSinglePredictionValue(weight,wlk_filtered,wmk_filtered,wm0_filtered)
        result[int(row+(window_dimention+1)/2)-1][int(col+(window_dimention+1)/2)-1] = center_value
        
        col += 1
    row += 1


def write_raster(outfile, data, proj, geo):
   gdal.GDT_UInt8 = gdal.GDT_Byte
   np_dtype = str(data.dtype)
   dtype = eval('gdal.GDT_' + np_dtype.title().replace('Ui','UI'))
   driver = gdal.GetDriverByName('GTiff')
   ny, nx = data.shape
   tfh = driver.Create(outfile, nx, ny, 1, dtype, [])
   tfh.SetProjection(proj)
   tfh.SetGeoTransform(geo)
   tband = tfh.GetRasterBand(1)
   tband.WriteArray(data)
   tfh = None

fp=gdal.Open(lkfile,1)
proj = fp.GetProjection()
geo = fp.GetGeoTransform()
outfile='C:/Users/mrwan/Desktop/AGS_py_code/tifs/11.tif'

#result[result>10000]=0

write_raster(outfile,result,proj,geo)
