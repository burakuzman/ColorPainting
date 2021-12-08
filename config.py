from ColorConversion import *

#INPUT PARAMETERS
#################

pigments=[\
#RGB codes for pigment colors
#red
[255,0,0,255],\
#blue
[0,255,0,255],\
#green
[0,0,255,255],\
#cyan
[0,255,255,255],\
#magenta
[255,0,255,255],\
#yellow
[255,255,0,255],\
#black
[0,0,0,255],\
#white
[255,255,255,255]\
]

otherpigments=[\
#RGB codes for pigment colors
#blue
[0,255,0,255],\
#green
[0,0,255,255],\
#cyan
[0,255,255,255],\
#magenta
[255,0,255,255],\
#yellow
[255,255,0,255],\
]

Pagesize=(8.0,8.0) #in inches

design_zone=(0.75,0.75) #size of the design space in inches (L * W)
eval_zone=(1.0,1.0) #size of the evaluation space in inches (L * W)
eval_offset_x=0.75 #in inches. This is how much the design and eval zones move at each iteration.
eval_offset_y=0.75

error_criteria=70 #percent
#dotradius=0.085 #inches or mm

pos_bit_depth=10 #this it the bit depth used to code/decode position genes
col_bit_depth=4 #this it the bit depth used to code/decode position genes
no_of_shapes=250 #this is the number of shapes crated for each zone
n_color=len(pigments)

import cv2
import numpy as np
import math
import sys
import os
import random
import time

def calcdotraster(dotradius,ppi_x,ppi_y):
    #Calculate the pixels within a circle with respect to the center of the circle, and assign alpha
    #This function converts a disk into rastered pixels with alpha anti-aliasing. The output is in format (xi, xj, alpha) where xi and xj
    #are incremental pixel steps from the center pixel
    r=math.sqrt(1/(ppi_x*ppi_x)+1/(ppi_y*ppi_y))/2 #this is in inches/mm - the radius of a circle that circumscribes a square pixel
    R=dotradius #this is the radius of the circle for dot made by the pen
    a_x=(int(dotradius*ppi_x)+1) #this is in pixels
    a_y=(int(dotradius*ppi_y)+1)
    rasterboxcorner=(-int(a_x),-int(a_y)) #this is in pixels
    pixlist=[]
    for j in range(int(2*a_y)):
        for i in range(int(2*a_x)):
            evaluatedpixel_center=((rasterboxcorner[0]+i+0.5)/ppi_x,(rasterboxcorner[1]+j+0.5)/ppi_y) # this is in inches/mm
            d=math.sqrt( evaluatedpixel_center[1]**2+evaluatedpixel_center[0]**2)
            squeval=(4*d*d*R*R-(d*d-r*r+R*R)**2)
            if d<=R-r:
                alfa = 1 #alfa is a float in range [0 to 1]
                pixlist.append((int(evaluatedpixel_center[0]*ppi_x), int(evaluatedpixel_center[1]*ppi_y),alfa))
            elif d>=r+R:
                alfa=0
            elif squeval>0:
                h=1/(2*d)*math.sqrt(squeval)
                teta = math.asin(h / R)
                if d-R*math.cos(teta)>0:
                    beta=math.asin(h/r)
                else:
                    beta=math.pi-math.asin(h/r)
                alfa=(r*r*(beta-math.sin(2*beta)/2)+R*R*(teta-math.sin(2*teta)/2))/(math.pi*r*r) #alfa is a float in range [0 to 1]
                pixlist.append((int(evaluatedpixel_center[0]*ppi_x), int(evaluatedpixel_center[1]*ppi_y),alfa))
    return pixlist

def alpha_over(fg,bg):
    #this does the alpha blending with fg over bg. requires fg and bg to be both in form RGBA.
    #returns res as RGBA
    fg_alfa=fg[3]/255
    #bg_alfa=bg[3]/255
    bg_alfa=1.0
    res_alfa=fg_alfa+bg_alfa*(1-fg_alfa)
    fg_RGB=fg[0:3]
    bg_RGB=bg[0:3]
    fg_XYZ=np.array(RGB2XYZ(fg_RGB))
    bg_XYZ=np.array(RGB2XYZ(bg_RGB))
    res_XYZ=(fg_XYZ*fg_alfa+bg_XYZ*bg_alfa*(1-fg_alfa))/res_alfa
    res_RGBA=list(XYZ2RGB(res_XYZ.tolist()))+[int(res_alfa*255)]
    return tuple(res_RGBA)

class brushstroke:
    id_count=0
    dotradius=0.085
    br_ppi_x=None
    br_ppi_y=None
    br_pixelmap=None
    def __init__(self,centerpoint,zz=1,bcolor=(0,0,0,0)):
        #brushstroke.id_count=brushstroke.id_count+1
        self.centerpoint=centerpoint #centerpoinst is in (x,y) format where x and y are in inches or mm in global coord
        #self.z=brushstroke.id_count #this will be an integer indicating depth 1 is the background,
        self.z=zz
        self.color=bcolor  #this will be in RGBA mode
        #self._pixelmap=pixelmap()

    @property
    def pixelmap(self):
        if brushstroke.br_pixelmap==None:
            br_pixelmap = calcdotraster(brushstroke.dotradius, brushstroke.br_ppi_x, brushstroke.br_ppi_y)
            brushstroke.br_pixelmap=br_pixelmap
        else:
            br_pixelmap=brushstroke.br_pixelmap
        return br_pixelmap

    def exportpixels(self):
        pixlist=[]
        ppi_x=brushstroke.br_ppi_x
        ppi_y=brushstroke.br_ppi_y
        cx = int(self.centerpoint[0] * ppi_x) #in pixels
        cy = int(self.centerpoint[1] * ppi_y) #in pixels
        for increments in self.pixelmap:
            delx=increments[0]
            dely=increments[1]
            alphachannel=increments[2]*self.color[3] #this adjusts the color alpha for coverage percentage. result is
            # in range [0 to 255]
            tempix=(cx+delx,cy+dely,(self.color[0],self.color[1],self.color[2],int(alphachannel)))
            pixlist.append(tempix) #tempix is in format (cx,cy, RGBA) cx and cy are in pixels, cx and cy are in global
            # coordinates. RGBA is in (255,255,255,255) format
        return pixlist

    def delete(self):
        del self

def roi_errorcalc(source,canvas,nregions):
    (height,width)=source.shape[0:2]
    if (height,width) != canvas.shape[0:2]:
        print ('****roi_error_calc -- REGIONS ARE DIFFERENT SIZE***')
    cell_h=int(height/nregions) #pixels
    cell_w=int(width/nregions) #pixels
    cum_error=0
    y_pos = 0
    max_error=0
    while y_pos < height:
        x_pos = 0
        while x_pos < width:
            roi_x1 = x_pos #pixels
            roi_y1 = y_pos #pixels
            roi_x2 = x_pos + cell_w #pixels
            roi_y2 = y_pos + cell_h #pixels
            #print (roi_x1,',',roi_y1,'**',roi_x2,',',roi_y2)
            source_ROI = source[roi_y1:roi_y2, roi_x1:roi_x2]
            k1=source_ROI.shape[0]*source_ROI.shape[1]
            sourcepixels=source_ROI.reshape(k1,1,3)
            sourceXYZ=[RGB2XYZ(pixel.tolist()[0]) for pixel in sourcepixels]
            sourcexyz2=np.array(sourceXYZ)
            source_val=np.mean(sourcexyz2, axis=0)
            source_rgb=XYZ2RGB(source_val)

            canvas_ROI = canvas[roi_y1:roi_y2, roi_x1:roi_x2]
            k2 = canvas_ROI.shape[0] * canvas_ROI.shape[1]
            canvaspixels=canvas_ROI.reshape(k2,1,3)
            canvasXYZ=[RGB2XYZ(pixel.tolist()[0]) for pixel in canvaspixels]
            canvasxyz2=np.array(canvasXYZ)
            canvas_val=np.mean(canvasxyz2, axis=0)
            canvas_rgb=XYZ2RGB(canvas_val)

            cell_error=math.sqrt( (canvas_rgb[0]-source_rgb[0])**2+(canvas_rgb[1]-source_rgb[1])**2 +(canvas_rgb[2]-source_rgb[2])**2)
            if cell_error>max_error:
                max_error=cell_error
            cum_error=cum_error+cell_error
            x_pos = x_pos + cell_w
        y_pos = y_pos + cell_h
    cum_error=cum_error*100/(nregions*nregions*255) #on a scale of 0 to 100
    max_error=max_error*100/(nregions*nregions*255) #on a scale of 0 to 100
    return (cum_error,max_error)
