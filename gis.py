'''
Created on 23 feb. 2018

@author: thomasgumbricht
'''


from __future__ import division
import os
import sys

from gdalconst import *
from osgeo import gdal, ogr, osr
from itertools import chain
import array as arr
import numpy as np
#import geoimagine.support.karttur_dt as mj_dt
'''
import mj_datetime_v50 as mj_dt

from shapely import geometry
from itertools import chain
import array as arr
import csv
import poplib
#from dbfpy import dbf
from shutil import move
'''
class Composition:
    def __init__(self,celltype,cellnull):
        """The constructoris just an empty container.""" 
        self.celltype = celltype
        self.cellnull = cellnull

class RasterPalette():
    def __init__(self):
        """The constructoris just an empty container.""" 
            
    def SetTuplePalette(self,paletteT):
        self.paletteT = paletteT
        self.FixGDALPalette()
        ct = gdal.ColorTable() 
        for row in self.PcR:
            pass

        for c in range(1,len(self.PcR)):
            ct.CreateColorRamp(self.PcR[c-1][0],self.PcR[c-1][1],self.PcR[c][0],self.PcR[c][1])
        self.colortable = ct
       
    def FixGDALPalette(self):
        #paletteT = mj_psycopg2.SelectPalettePub(palette)
        if len(self.paletteT) == 0:
            NOPALETTE
            return False
        PcR = []
        AT = []
        ATvalue = []
        ATattr = []
        for c in range(len(self.paletteT)):
            cr = (self.paletteT[c][0],(self.paletteT[c][1],self.paletteT[c][2],self.paletteT[c][3]))
            PcR.append(cr)
            ATvalue.append(self.paletteT[c][0])
            ATattr.append(self.paletteT[c][5])
        for a in range(self.paletteT[c][0]): #Up to the highest attibute
            if a in ATvalue:
                x = ATvalue.index(a)
                AT.append(ATattr[x])
            else:
                AT.append('')
        self.PcR = PcR
        self.AT = AT
        self.maxAT = self.paletteT[c][0]
        #return PcR,AT,paletteT[c][0]
        
class FieldDef():
    def __init__(self, name,fieldDefD):
        '''{'type','transfer','source','width','precision','keyfield/field'}'''
        self.name = name
        if 'type' in fieldDefD:
            self.type = fieldDefD['type']
        if 'transfer' in fieldDefD:
            self.transfer = fieldDefD['transfer']
        if 'source' in fieldDefD:
            self.source = fieldDefD['source']
        if 'width' in fieldDefD:
            self.width = fieldDefD['width']
        else:
            self.width = 8
        if 'precision' in fieldDefD:
            self.precision = fieldDefD['precision']
        else:
            self.precision = 0
        if 'keyfield' in fieldDefD:
            self.keyfield = fieldDefD['keyfield']
        if 'field' in fieldDefD:
            self.keyfield = fieldDefD['field']
        else:
            self.keyfield = False
                       
class MjProj:
    def __init__(self):
        """The constructoris just an empty container.""" 
         
    def SetProj(self,spatialRef):
        self.proj_cs = spatialRef
        
    def SetGeoTrans(self,gt):
        self.gt = gt
        
    def SetFromProj4(self,proj4):
        self.proj4 = proj4
        srs = osr.SpatialReference()  
        srs.ImportFromProj4(self.proj4)
        self.proj_cs = osr.SpatialReference()
        self.proj_cs.ImportFromWkt(srs.ExportToWkt())
        
    def SetProj4(self): 
        srs = osr.SpatialReference()
        self.WKT = self.proj_cs.ExportToWkt()
        srs.ImportFromWkt(self.WKT)
        self.proj4 = srs.ExportToProj4()
        
    def SetFromWKT(self,WKT):
        self.WKT = WKT
        self.proj_cs = osr.SpatialReference()
        self.proj_cs.ImportFromWkt(WKT)
        
    def SetFromEPSG(self,epsg):
        self.epsg = epsg      
        self.proj_cs = osr.SpatialReference() 
        if epsg == 6842:
            WKT = 'PROJCS["Sinusoidal_Sanson_Flamsteed",GEOGCS["GCS_Unknown",DATUM["D_unknown",SPHEROID["Unknown",6371007.181,"inf"]],PRIMEM["Greenwich",0],UNIT["Degree",0.017453292519943295]],PROJECTION["Sinusoidal"],PARAMETER["central_meridian",0],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["Meter",1]]'
            self.SetFromWKT(WKT)
        else:
            self.proj_cs.ImportFromEPSG(epsg)
 
    def ReadSpatialRef(self):
        self.SetProj4()
        print (self.proj4)
        print (self.proj_cs)
        self.SetFromProj4(self.proj4)
        print (self.proj_cs)
        self.epsg = int( self.proj_cs.GetAttrValue("AUTHORITY", 1) )
        self.spatialunit = self.proj_cs.GetAttrValue("UNIT", 0)
        self.datum = self.proj_cs.GetAttrValue("DATUM", 0)
        self.spheroid = self.proj_cs.GetAttrValue("SPHEROID", 0)
        self.projection = self.proj_cs.GetAttrValue("PROJECTION", 0)
        self.setEpsg = False
        if self.epsg == None and self.datum in ['WGS_1984','WGS_84'] and self.spheroid in ['WGS_1984','WGS_84'] and self.spatialunit.lower() == 'degree':
            self.epsg = 4326
            self.setEpsg = True
        elif self.epsg == 4326:
            pass
        elif self.proj4.replace(' ','') == '+proj=sinu+lon_0=0+x_0=0+y_0=0+a=6371007.181+b=6371007.181+units=m+no_defs':
            self.epsg = 6842
            self.setEpsg = True
        elif self.projection == None:
            self.epsg = False

            
class RasterDataSource: 
    def __init__(self):
        """The constructoris just an empty container.""" 
        
    def OpenGDALRead(self,FPN): 
        self.rasterFPN = FPN
        if os.path.exists(self.rasterFPN):
            self.datasource = gdal.Open(FPN, GA_ReadOnly)
            if self.datasource is None:
                exitstr = 'Exiting - Failed to open raster file %s' %(self.rasterFPN)
                sys.exit(exitstr)
        else:
            exitstr = 'Raster datasource %s does not exist' %(self.rasterFPN)
            sys.exit(exitstr)
                                                      
    def OpenGDALEdit(self,FPN): 
        self.rasterFPN = FPN
        if os.path.exists(self.rasterFPN):
            self.datasource = gdal.Open(FPN)
            if self.datasource is None:
                exitstr = 'Exiting - Failed to open raster file %s' %(self.rasterFPN)
                sys.exit(exitstr)
        else:
            exitstr = 'Raster datasource %s does not exist' %(self.rasterFPN)
            sys.exit(exitstr)
       
    def CreateGDALraster(self,FPN,layer,of = 'GTiff'):
        #self.datasource = gdal.GetDriverByName('GTiff').Create(FPN, x_res, y_res, 1, gdal.GDT_Byte)
        driver = gdal.GetDriverByName( of )
        #metadata = driver.GetMetadata()
        #if metadata.has_key(gdal.DCAP_CREATE) \
        #    and metadata[gdal.DCAP_CREATE] == 'YES':
        if layer.comp.celltype.lower() in ['byte','uint8']:
            self.datasource = driver.Create( FPN, layer.cols, layer.lins, 1, gdal.GDT_Byte )   
            layer.comp.cellnull = int(layer.comp.cellnull) 
        elif layer.comp.celltype.lower() == 'int16':  
            self.datasource = driver.Create( FPN, layer.cols, layer.lins, 1, gdal.GDT_Int16 ) 
            layer.comp.cellnull = int(layer.comp.cellnull)
        elif layer.comp.celltype.lower() == 'uint16':
            self.datasource = driver.Create( FPN, layer.cols, layer.lins, 1, gdal.GDT_UInt16 )
            layer.comp.cellnull = int(layer.comp.cellnull)
        elif layer.comp.celltype.lower() == 'float32':
            self.datasource = driver.Create( FPN, layer.cols, layer.lins, 1, gdal.GDT_Float32 )
            layer.comp.cellnull = float(layer.comp.cellnull)
        else:
            exitstr ='numpy type not defined',layer.comp.celltype
            sys.exit(exitstr) 
        self.datasource.SetGeoTransform( layer.geotrans )
        self.datasource.SetProjection( layer.projection  ) 
        if hasattr(layer.comp,'palette') and layer.comp.palette:
            palette = RasterPalette()
            palette.SetTuplePalette(layer.comp.palette)
            self.datasource.GetRasterBand(1).SetColorTable(palette.colortable) 
        self.datasource.GetRasterBand(1).WriteArray( layer.BAND )
        self.datasource.GetRasterBand(1).SetNoDataValue(layer.comp.cellnull)
        '''
            #PcR,AT,maxAT = self.FixGDALPalette(self.palette)
            #ct = gdal.ColorTable() 
            
            For discrete colors
            #ct.CreateColorRamp(0,(178,223,138),5,(255,127,0))
            #ct.CreateColorRamp(Pcr)
            for c in PcR:
        
                ct.SetColorEntry(c[0],c[1])
        
            #for color ramps
            for c in range(1,len(PcR)):
                ct.CreateColorRamp(PcR[c-1][0],PcR[c-1][1],PcR[c][0],PcR[c][1])
            self.datasource.GetRasterBand(1).SetColorTable(ct) 
        
            rat = gdal.RasterAttributeTable()
            rat.CreateColumn("Value", GFT_String, GFT_String)
            for i in range(maxAT): 
                rat.SetValueAsString(i, 0, AT[i])
            self.datasource.GetRasterBand(1).SetDefaultRAT(rat)
        '''
        self.datasource.FlushCache()
        self.datasource = None
        #else:
        #    print ( 'Driver %s does not supports CreateCopy() method.' % of )
 
    def CreateGDALDS(self,FPN,layer,of = 'GTiff'):
        #self.datasource = gdal.GetDriverByName('GTiff').Create(FPN, x_res, y_res, 1, gdal.GDT_Byte)
        driver = gdal.GetDriverByName( of )
        metadata = driver.GetMetadata()
        if metadata.has_key(gdal.DCAP_CREATE) \
            and metadata[gdal.DCAP_CREATE] == 'YES':
            if layer.comp.celltype.lower() in ['byte','uint8']:
                self.datasource = driver.Create( FPN, layer.cols, layer.lins, 1, gdal.GDT_Byte )   
                layer.comp.cellnull = int(layer.comp.cellnull) 
            elif layer.comp.celltype.lower() == 'int16':  
                self.datasource = driver.Create( FPN, layer.cols, layer.lins, 1, gdal.GDT_Int16 ) 
                layer.comp.cellnull = int(layer.comp.cellnull)
            elif layer.comp.celltype.lower() == 'uint16':
                self.datasource = driver.Create( FPN, layer.cols, layer.lins, 1, gdal.GDT_UInt16 )
                layer.comp.cellnull = int(layer.comp.cellnull)
            elif layer.comp.celltype.lower() == 'float32':
                self.datasource = driver.Create( FPN, layer.cols, layer.lins, 1, gdal.GDT_Float32 )
                layer.comp.cellnull = float(layer.comp.cellnull)
            else:
                exitstr ='numpy type not defined',layer.comp.celltype
                sys.exit(exitstr) 
            self.datasource.SetGeoTransform( layer.geotrans )
            self.datasource.SetProjection( layer.projection  )
            self.datasource.GetRasterBand(1).SetNoDataValue(layer.comp.cellnull)
            if hasattr(layer.comp,'palette') and layer.comp.palette:
                PcR,AT,maxAT = self.FixGDALPalette(self.comp.palette)
                ct = gdal.ColorTable() 
                for c in range(1,len(PcR)):
                    ct.CreateColorRamp(PcR[c-1][0],PcR[c-1][1],PcR[c][0],PcR[c][1])
                self.datasource.GetRasterBand(1).SetColorTable(ct)
                    
    def WriteBlock(self,col,row,arr):
        self.datasource.GetRasterBand(1).WriteArray(arr,col, row)
        #self.datasource.GetRasterBand(1).WriteArray
 
    def arrTo2Dnp(self,layer):
        if layer.comp.celltype.lower() in ['byte','uint8']:
            band2DArray = np.asarray(layer.BAND, dtype=np.int8)
        elif layer.comp.celltype.lower() == 'int16':  
            band2DArray = np.asarray(layer.BAND, dtype=np.int16)    
        elif layer.comp.celltype.lower() == 'uint16':
            band2DArray = np.asarray(layer.BAND, dtype=np.uint16)
        else:
            ( 'print numpy type not defined',layer.comp.celltype )
            sys.exit()   
        #reshape the 1D array to a 2D image
        band2DArray.shape = (-1, layer.cols)
        layer.BAND = band2DArray
    
    def SetGeoTransform(self,gt): 
        self.datasource.SetGeoTransform(gt) 
    
    def SetProjection(self,proj): 
        self.datasource.SetProjection(proj)  
            
    def CloseDS(self):
        #close the datasource
        self.datasource = None

class RasterLayer:
    def __init__(self):
        """The constructoris just an empty container.""" 
        
    def SetDS(self,DS): 
        self.datasource =  DS.datasource
        
    def GetLayer(self,bandnr):
        self.layer = self.datasource.GetRasterBand(bandnr)
        
    def GetSpatialRef(self): 
        self.spatialRef = self.proj_cs = self.datasource.GetProjection()
        self.gt = self.datasource.GetGeoTransform()
        
        #self.spatialRef = self.proj_cs = self.layer.GetSpatialRef() 
    def SetSpatialRef(self, tarProj):
        self.spatialRef = tarProj.proj_cs
        self.gt = tarProj.gt
        
    def GetGeometry(self): 
        #Get the necessary band information
        self.lins = self.datasource.RasterYSize
        self.cols = self.datasource.RasterXSize
        self.cellnull = self.layer.GetNoDataValue()
        self.celltype = gdal.GetDataTypeName(self.layer.DataType)
        self.projection = self.datasource.GetProjection()
        self.geotrans = self.datasource.GetGeoTransform()

        #Get the extent of the image
        self.ext=[]
        xarr=[0,self.cols]
        yarr=[0,self.lins]
        for px in xarr:
            for py in yarr:
                x=self.gt[0]+(px*self.gt[1])+(py*self.gt[2])
                y=self.gt[3]+(px*self.gt[4])+(py*self.gt[5])
                self.ext.append([x,y])
            yarr.reverse()
        self.bounds = (self.ext[0][0], self.ext[2][1], self.ext[2][0],self.ext[0][1])
        #Get the spatial resolution
        cellsize = [(self.ext[2][0]-self.ext[0][0])/self.cols, (self.ext[0][1]-self.ext[2][1])/self.lins] 
        if cellsize[0] != cellsize[1]:
            pass
        self.cellsize = cellsize[0]
 
    def ReadBand(self):
        self.NPBAND = self.layer.ReadAsArray()
        
        
    def ReadBlock(self, col,row,ncols,nrows):
        
        NPBLOCK = self.layer.ReadAsArray(col, row, ncols, nrows)
        self.NPBLOCK = NPBLOCK.ravel()


    def Flatten2DtoArr(self,structcodeD):
        self.GetGeometry()
        structcode  = structcodeD[self.cellType]

        bandList = list(chain.from_iterable(self.NPBAND))   
        self.bandArray = arr.array(structcode)
        self.bandArray.fromlist(bandList)
  
    def ColorTable(self,colortable):
        print ('    Setting color table (0 if success):'),
        print ( self.layer.SetColorTable(colortable) )
        
    def CloseLayer(self):
        #close the layer
        self.layer = None

class VectorDataSource:  
    def __init__(self):
        """The constructoris just an empty container.""" 
        
    def OpenESRIshapeRead(self,FPN): 
        self.ESRIshapeFPN = FPN
        if os.path.exists(self.ESRIshapeFPN):
            driver = ogr.GetDriverByName("ESRI Shapefile")
            self.datasource = driver.Open(self.ESRIshapeFPN, 0)
        else:
            exitstr = 'ESRI datasource %s does not exist' %(self.ESRIshapeFPN)
            sys.exit(exitstr)
            
    def OpenESRIshapeEdit(self,FPN): 
        self.ESRIshapeFPN = FPN
        if os.path.exists(self.ESRIshapeFPN):
            driver = ogr.GetDriverByName("ESRI Shapefile")
            self.datasource = driver.Open(self.ESRIshapeFPN, 1)
        else:
            exitstr = 'ESRI datasource %s does not exist' %(self.ESRIshapeFPN)
            sys.exit(exitstr)
            
    def CreateESRIdatasource(self,FPN):
        self.FPN = FPN
        driver = ogr.GetDriverByName("ESRI Shapefile")
        if os.path.exists(self.FPN):
            driver.DeleteDataSource(self.FPN)
        self.datasource = driver.CreateDataSource(self.FPN)
        if self.datasource is None:
            print ( 'Could not create file' )
            sys.exit(1)
                   
    def CloseDS(self):
        self.datasource = None 
  
class VectorLayer:
    def __init__(self):
        """The constructoris just an empty container.""" 
        
    def SetDS(self,DS): 
        self.datasource =  DS.datasource
             
    def GetLayer(self):
        geomTypeD = {1:'point',2:'line',3:'polygon'} 
        self.layer = self.datasource.GetLayer()
        self.layerid = self.layer.GetName()

        print ('self.layerid',self.layerid)
        print ('geomType', self.layer.GetGeomType() )
        self.geomtype = geomTypeD[self.layer.GetGeomType()]
        '''
        self.layerDef = self.layer.GetLayerDefn()

        print ( self.layerDef.GetGeomType() )
        
        first_feat = self.layer.GetFeature(1)

        self.geomtype = first_feat.geometry().GetGeometryName()
        '''
        return self.GetLayerDef()
        
    def GetLayerDef(self):
        layerDefn = self.layer.GetLayerDefn()
        fieldDefD = {}  
        fieldDefL = []
        for i in range( layerDefn.GetFieldCount() ):            
            n = layerDefn.GetFieldDefn(i).GetName()
            t = layerDefn.GetFieldDefn(i).GetType()
            #fieldDefD['typename'] = layerDefn.GetFieldDefn(i).GetFieldTypeName()
            t = layerDefn.GetFieldDefn(i).GetFieldTypeName(t)
            w = layerDefn.GetFieldDefn(i).GetWidth()
            p = layerDefn.GetFieldDefn(i).GetPrecision()
            fieldDefD[n] = {'type':t.lower(), 'width':w, 'precision': p, 'transfer': 'constant'}      
            #fieldDefD['source'] = self.feature.GetField(fieldDefD['name'])
        for key in fieldDefD:
            fieldDefL.append(FieldDef(key,fieldDefD[key]))
        return fieldDefL
                        
    def CreateOgrLayer(self,geomtype,tarLayerId):
        self.geomtype = geomtype
        if self.geomtype.lower() == 'point':
            self.layer = self.datasource.CreateLayer(tarLayerId, self.spatialRef, geom_type=ogr.wkbPoint) 
        elif self.geomtype.lower() == 'line':
            self.layer = self.datasource.CreateLayer(tarLayerId, self.spatialRef, geom_type=ogr.wkbLineString) 
        elif self.geomtype.lower() == 'polygon':
            self.layer = self.datasource.CreateLayer(tarLayerId, self.spatialRef, geom_type=ogr.wkbPolygon) 
        elif self.geomtype.lower() == 'multipolygon':
            self.layer = self.datasource.CreateLayer(tarLayerId, self.spatialRef, geom_type=ogr.wkbMultiPolygon)     
   
        else:
            exitstr = 'Can not recoginze geometry type %s for creating layer' %(self.geomtype) 
            print (exitstr)
            BALLE
            sys.exit(exitstr)
        if self.layer is None:
            BALLE
            print ( 'Could not create layer' )
            sys.exit() 
            
    def AddFieldDef(self,fieldDef):
        #Get all field
        self.layerDefn = self.layer.GetLayerDefn()
        if fieldDef.name not in self.schemaD:
            #Check if field exists:
            if fieldDef.type == 'string':
                new_field = ogr.FieldDefn(fieldDef.name, ogr.OFTString)
                new_field.SetWidth(fieldDef.width)
            elif fieldDef.type == 'integer':
                new_field = ogr.FieldDefn(fieldDef.name, ogr.OFTInteger)
            elif fieldDef.type == 'real':
                new_field = ogr.FieldDefn(fieldDef.name, ogr.OFTReal)
            elif fieldDef.type == 'date':
                new_field = ogr.FieldDefn(fieldDef.name, ogr.OFTDate)
            elif fieldDef.type == 'time':
                new_field = ogr.FieldDefn(fieldDef.name, ogr.OFTTime)
            else:
                printstr =  'field type not recognized %s' %(fieldDef.type) 
                print (printstr)
                BALLE           
            self.layer.CreateField(new_field)
        else:
            for i in range( self.layerDefn.GetFieldCount() ):
                if self.layerDefn.GetFieldDefn(i).GetName() == fieldDef.name:

                    fieldTypeCode = self.layerDefn.GetFieldDefn(i).GetType()
                    fieldType = self.layerDefn.GetFieldDefn(i).GetFieldTypeName(fieldTypeCode)
                    fieldWidth = self.layerDefn.GetFieldDefn(i).GetWidth()
                    fieldPrecision = self.layerDefn.GetFieldDefn(i).GetPrecision()

                    if fieldType.lower() != fieldDef.type.lower():
                        sys.exit('Existing field of different type already exists')
         
    def AddField(self,f):
        FITTA
              
    def FillAllFields(self,fieldDefL,keyfield):
        for f in fieldDefL:
            if f.transfer == 'constant':
                for feature in self.layer:  
                    feature.SetField(f.name, f.source)
                    self.layer.SetFeature(feature)
            elif f.transfer == 'copy':
                #all features have the same value
                for feature in self.layer:  
                    copy = feature.GetField(f.source)
                    feature.SetField(f.name, copy)
                    self.layer.SetFeature(feature)
            elif f.transfer == 'dict':
                #all features have the same value
                for feature in self.layer:  
                    featkey = feature.GetField(keyfield)
                    fielddata = f.source[featkey]
                    feature.SetField(f.name, fielddata)
                    self.layer.SetFeature(feature)
            else:
                exitstr = 'missing transfer type %s' %(f.transfer)
                sys.exit(exitstr)
       
    def AddPtDataFromDict(self,xkey,ykey,dictL):
        from datetime import datetime
        featDefn = self.layer.GetLayerDefn()
        for D in dictL:
            tarfeat = ogr.Feature(featDefn)
            #set the geometry
            geom = ShapelyPointGeom( (D[xkey], D[ykey]) )
            geom.ShapelyToOgrGeom()
            tarfeat.SetGeometry(geom.ogrGeom)
                    
            for key in D:  
                print ('feature', key, D[key])
                tarfeat.SetField(key, D[key])         
            self.layer.CreateFeature(tarfeat) 
       
    def GetSpatialRef(self): 
        self.spatialRef = self.proj_cs = self.layer.GetSpatialRef() 

    def SetSpatialRef(self, tarProj):
        self.spatialRef = tarProj.proj_cs
        
    def SetCSSpatialRef(self, proj_cs):
        self.spatialRef = proj_cs
        
    def GetLayerFields(self):  
        self.schemaD = {}
        ldefn = self.layer.GetLayerDefn()
        for n in range(ldefn.GetFieldCount()):
            fdefn = ldefn.GetFieldDefn(n)
            self.schemaD[fdefn.name]= [fdefn.GetType(),fdefn.GetWidth(),fdefn.GetPrecision()]
        
    def CompareFields(self,srcFieldId,srcFieldName,fieldDefL):
        res = True
        for feature in self.layer:  
            item = feature.GetField(srcFieldId)
            for f in fieldDefL:
                if f.transfer == 'dict':
                    if item in f.source:
                        pass
                    else:
                        print ('missing', item,feature.GetField(srcFieldName) )
                        res = False
        return res  

class ogrFeature:
    def __init__(self,layer):
        """The constructoris just an empty container.""" 
        self.layer = layer.layer
         
    def SetFeature(self,feature):
        self.feature = feature 
        
    def CreateOgrFeature(self, geom, fieldDefL, dictL = False, key = False):
        featDefn = self.layer.GetLayerDefn()
        tarfeat = ogr.Feature(featDefn)
        tarfeat.SetGeometry(geom.ogrGeom)
        
        for f in fieldDefL:
            if f.transfer in ['constant','fixed','field']:
                tarfeat.SetField(f.name,f.source)       
            elif f.transfer == 'copy':
                pass
                STOP
            elif f.transfer == 'db':
                tarfeat.SetField(f.name,f.source)
            elif f.transfer == 'dict':
                tarfeat.SetField(f.name,f.source[key])
            else:
                exitstr = 'Unrecognized transfer command: %s' %(f.transfer)
                sys.exit(exitstr)
        self.layer.CreateFeature(tarfeat)
        
    def CopyOgrFeature(self,srcfeat,fieldDefL):
        featDefn = self.layer.GetLayerDefn()
        tarfeat = ogr.Feature(featDefn)
        tarfeat.SetGeometry(srcfeat.GetGeometryRef())
        for f in fieldDefL:        
            tarfeat.SetField(f.name, f.source)           
        self.layer.CreateFeature(tarfeat)      
        
    def AddSourceValtoFieldDef(self,fieldDefL):
        for f in fieldDefL:
            if f.transfer == 'fixed':
                pass
            elif f.transfer == 'field':

                f.source = self.feature.GetField(f.field)
            else:
                f.source = self.feature.GetField(f.name)
    
class Geometry:
    def __init__(self):
        """The constructoris just an empty container.""" 
        
    def GeomFromFeature(self,feature): 
        self.ogrGeom = feature.GetGeometryRef() 
        self.OgrGeomToShapely()     
        
    def SetShapelyGeom(self,geom):   
        self.shapelyGeom = geom
    
    def CascadeUnion(self,shapelyGeomL):
        from shapely.ops import cascaded_union
        self.SetShapelyGeom(cascaded_union(shapelyGeomL))
        #return cascaded_union(shapelyGeomL)
        
    def OgrGeomToShapely(self):
        from shapely.wkb import loads
        self.shapelyGeom = loads(self.ogrGeom.ExportToWkb())
        
    def ShapelyToOgrGeom(self):
        self.ogrGeom = ogr.CreateGeometryFromWkb(self.shapelyGeom.wkb)
        
    def ShapelyBuffer(self,buff):
        if not self.shapelyGeom.is_valid:
            print ( '    mending invalid geometry by buffer = 0' )
            self.shapelyGeom = self.shapelyGeom.buffer(0)
            if not self.shapelyGeom.is_valid:
                exitstr = 'Can not fix invalid geometry in ShapelyBuffer'
                sys.exit(exitstr)
                
    def ShapelyPolysToMultiPoly(self,multiPolyL):
        from shapely.geometry.multipolygon import MultiPolygon
        self.shapelyGeom = MultiPolygon(multiPolyL)  
        
    def MultiPolyGeomFromGeomL(self,geomL):
        singlegeomL = []
        for testgeom in geomL:
            if testgeom.geom_type == 'MultiPolygon':
                for polygeom in [polygeom for polygeom in testgeom]:
                    if len(list(polygeom.exterior.coords)) > 2:
                        singlegeomL.append(polygeom)
            else:
                if len(list(testgeom.exterior.coords)) > 2:
                    singlegeomL.append(testgeom)
        self.ShapelyPolysToMultiPoly(self,singlegeomL)
        
    def SplitPolygonsToSingleGeomL(self):
        singlegeomL = []
        if self.shapelyGeom.geom_type == 'MultiPolygon':
            for polygeom in [polygeom for polygeom in self.shapelyGeom]:
                if len(list(polygeom.exterior.coords)) > 2:
                    singlegeomL.append(polygeom)
        else:
            if len(list(self.shapelyGeom.exterior.coords)) > 2:
                    singlegeomL.append(self.shapelyGeom)
        return singlegeomL
             
    def ShapelyRemoveSlivers(self):
        eps = 0.000001
        from shapely.geometry import JOIN_STYLE
        self.ShapelyBuffer(0)
        #remove slivers
        self.shapelyGeom = self.shapelyGeom.buffer(eps, 1, join_style=JOIN_STYLE.mitre).buffer(-eps, 1, join_style=JOIN_STYLE.mitre)
                
    def ShapelyContainCut(self,cut):
        if cut.shapelyGeom.contains(self.shapelyGeom):
            pass
        #geom remains intact
        elif cut.shapelyGeom.intersects(self.shapelyGeom):
            self.shapelyGeom = cut.shapelyGeom.intersection(self.shapelyGeom)
        else:
            self.shapelyGeom = False
                
    def ShapelyOutsideCut(self,cut):
        if cut.shapelyGeom.contains(self.shapelyGeom):
            self.shapelyGeom = False
        elif cut.shapelyGeom.intersects(self.shapelyGeom):
            self.shapelyGeom = self.shapelyGeom.difference(cut.shapelyGeom)
        else:
            pass
        #geom remains intact
        
    def GeoTransform(self,srcProj,tarProj):
        self.geoTransform = osr.CoordinateTransformation(srcProj.proj_cs,tarProj.proj_cs)
        res = self.ogrGeom.Transform(self.geoTransform)
        if res != 0:
            sys.exit('Geotransformation failed')
        #set shapley to reproejcted
        self.OgrGeomToShapely()
            
    def GeoTransformCoordsFIX(self,srcProj,tarProj,ptL): 
        #Special version that allows geotransformation outised lat lon boundaries (i.e. for SIN outised erait sphere
        transform = osr.CoordinateTransformation(srcProj.proj_cs,tarProj.proj_cs)
        xyL = []
        for pt in ptL:
            xyL.append(transform.TransformPoint(pt[0],pt[1]))
        return xyL
    
    def PointsToPolygonGeom(self,ptL):
        from shapely.geometry import Polygon
        #ptgeomL = []
        #for pt in ptL:
        #    ptgeomL.append(Point(pt))
        self.shapelyGeom = Polygon(ptL)
        #self.shapelyGeom = Polygon(ptgeomL)
        self.ShapelyToOgrGeom()
        
    def ShapelyIntersection(self, otherGeom): 
        return self.shapelyGeom.intersection(otherGeom.shapelyGeom)
        
class ShapelyPointGeom(Geometry):
    def __init__(self,pt):
        from shapely.geometry import Point
        Geometry.__init__(self)
        self.shapelyGeom = Point(pt)
        
class ShapelyMultiPointGeom(Geometry):
    def __init__(self,ptL):
        from shapely.geometry import MultiPoint
        Geometry.__init__(self)
        self.shapelyGeom = MultiPoint(ptL)

class ShapelyLineGeom(Geometry):
    def __init__(self,ptLT):
        from shapely.geometry import LineString
        Geometry.__init__(self)
        self.shapelyGeom = LineString(ptLT)

class ShapelyPolyGeom(Geometry):
    def __init__(self,ptLT):
        from shapely.geometry import Polygon
        Geometry.__init__(self)
        self.shapelyGeom = Polygon(ptLT)
        
def ESRICreateDSLayer(tarShpFPN,spatialRef,geomtype,layerid,fieldDefL):
    '''Creates a standard ESRI datasource and sets the layer and layerDef
    '''
    #Create an instance of datasource for the target
    tarDS = VectorDataSource()
    #Create the datasource
    tarDS.CreateESRIdatasource(tarShpFPN)
    #create an instance of layer
    tarLayer = VectorLayer()
    #link the datasource to the layer
    tarLayer.SetDS(tarDS)
    tarLayer.SetCSSpatialRef(spatialRef)   
    tarLayer.CreateOgrLayer(geomtype, str(layerid))
    tarLayer.GetLayerFields()
    for fieldDef in fieldDefL:
        tarLayer.AddFieldDef(fieldDef)
    return tarDS, tarLayer  

def ESRIOpenGetLayer(srcShpFPN,mode):
    '''Opens a standard ESRI datasource and gets the layer and layerDef
    '''
    #Create an instance of datasource for the source data
    srcDS = VectorDataSource()
    #open the source data for reading
    if mode == 'read':
        srcDS.OpenESRIshapeRead(srcShpFPN)
    elif mode == 'edit':
        srcDS.OpenESRIshapeEdit(srcShpFPN)
    #create a layer instance
    srcLayer = VectorLayer()
    #set the datasource of the layer
    srcLayer.SetDS(srcDS)
    #get the layer and the fieldDefintions 
    fieldDefL = srcLayer.GetLayer()
    #get the spatialref
    srcLayer.GetSpatialRef()
    

    return srcDS,srcLayer,fieldDefL 

def RasterCreateWithFirstLayer(dstRastFPN,layer):
    dstDS = RasterDataSource()
    dstDS.CreateGDALDS(dstRastFPN,layer)
    return dstDS 

def RasterOpenGetFirstLayer(srcRastFPN,modeD):
    '''Opens a standard GDAL Raster datasource and gets the first layer
    '''
    #Create an instance of datasource for the source data
    srcDS = RasterDataSource()
    #open the source data for reading    
    if modeD['mode'] == 'read':
        srcDS.OpenGDALRead(srcRastFPN)
        #create a layer instance
        srcLayer = RasterLayer()
        #set the datasource of the layer
        srcLayer.SetDS(srcDS)
    elif modeD['mode'] == 'edit':
        print (BALLE)
        srcDS.OpenGDALEdit(srcRastFPN)
        srcLayer = RasterLayer()
        srcLayer.SetDS(srcDS)
    else:
        print ( 'Can not understand edit mode', modeD['mode'] )
    srcLayer.GetLayer(1)
    #get the spatialref
    srcLayer.GetSpatialRef()
    srcLayer.GetGeometry()
    return srcDS,srcLayer 

def ReadRasterArray(srcRastFPN,readD):
    srcDS,srcLayer = RasterOpenGetFirstLayer(srcRastFPN,readD['mode'])
    srcLayer.ReadBand()
    srcDS.CloseDS()
    return srcLayer.NPBAND

def CreateDSWriteRasterArray(layer,writeD):
    tarDS = RasterDataSource()
    tarDS.CreateGDALraster(layer.FPN,layer,writeD['of'])   
     
def CreateESRIPolygonPtL(tarShpFPN, fieldDefL, ptL, proj_cs, layerid):
    #create the datasource
    tarDS,tarLayer = ESRICreateDSLayer(tarShpFPN, proj_cs, 'polygon', layerid, fieldDefL)
    #create the geometry
    geom = ShapelyPolyGeom(ptL)
    geom.ShapelyToOgrGeom()
    tarFeat = ogrFeature(tarLayer)
    tarFeat.CreateOgrFeature(geom, fieldDefL)
    tarDS.CloseDS()
    
def CreateESRIPolygonGeom(tarShpFPN, fieldDefL, geom, proj_cs, layerid):
    #create the datasource
    tarDS,tarLayer = ESRICreateDSLayer(tarShpFPN, proj_cs, 'polygon', layerid, fieldDefL)
    #create the geometry
    geom.ShapelyToOgrGeom()
    tarFeat = ogrFeature(tarLayer)
    tarFeat.CreateOgrFeature(geom, fieldDefL)
    tarDS.CloseDS()
    
def GetFeatureBounds(srcShpFPN,fieldId):
    srcDS,srcLayer,srcFieldDefL = ESRIOpenGetLayer(srcShpFPN,'read')
    boundsD = {}
    for feature in srcLayer.layer:
        srcFeat = ogrFeature(srcLayer)
        srcFeat.SetFeature(feature)
        fid = feature.GetField( str(fieldId) )
        geom = Geometry()
        geom.GeomFromFeature(feature)      
        boundsD[fid] = geom.shapelyGeom.bounds
    srcDS.CloseDS()
    return boundsD

def GetVectorProjection(srcShpFPN):
    filetype = os.path.splitext(srcShpFPN)[1]
    if filetype.lower() == '.shp':
        srcDS,srcLayer,srcFieldDefL = ESRIOpenGetLayer(srcShpFPN,'read')
        srcLayer.GetSpatialRef()
        spatialRef = MjProj()
        if srcLayer.spatialRef != None:
            spatialRef.SetProj(srcLayer.spatialRef)
            spatialRef.ReadSpatialRef()
        else:
            spatialRef.epsg = False
        srcDS.CloseDS()
        return spatialRef
    else:
        errorstr = 'urecognized format for retrieveing spatial reference: %s' %(filetype)
        sys.exit(errorstr)
    
def CopyRasterSetMeta(srcRasterFPN, tarRasterFPN, gdalformat, **kwargs): 
    # Opens source dataset
    srcDS,srcLayer = RasterOpenGetFirstLayer(srcRasterFPN,'read')
    #spatialRef = MjProj()
    #TGTODO HERE I ASSUME THAT RASTER PROJ IS IN WKT FORMAT
    srcLayer.GetSpatialRef()
    driver = gdal.GetDriverByName(gdalformat) 
    # Open destination dataset
    tardatasouce = driver.CreateCopy(tarRasterFPN, srcDS.datasource, 0)
    srcDS.CloseDS()
    tarDS = RasterDataSource()
    tarDS.datasource = tardatasouce
    # Set location
    tarDS.SetGeoTransform(srcLayer.gt)
    # Set projection
    tarDS.SetProjection(srcLayer.proj_cs)
    #Get the tar
    tarLayer = RasterLayer()
    tarLayer.SetDS(tarDS)
    tarLayer.GetLayer(1)
    if kwargs is not None:
        for key, value in kwargs.items():
            if key == 'cellnull' and value:
                tarLayer.SetNoDataValue(value)
            elif key == 'palette' and value:
                palette = RasterPalette()
                palette.SetTuplePalette(value)
                tarLayer.ColorTable(palette.colortable)
    # Close files
    tarDS.CloseDS()
 
def ReplaceRasterDS(srcRasterFPN, *args, **kwargs):
    FP = os.path.split(srcRasterFPN)[0]
    tmpFPN = os.path.join(FP,'temp.tmp')
    #set the gdal format
    if args is not None and len(args) > 0:
        if args[0].lower() in ['tif','.tif']:
            gdalformat = 'GTiff'
        else:     
            gdalformat = args[0]
    else:
        gdalformat = 'GTiff'

    if kwargs is None:
        printstr = 'No metadata to change for %s' %(srcRasterFPN)  
        print (printstr)
        return
    #Get the celltype of the existing layer
    spatialRef,srcLayer = GetRasterMetaData(srcRasterFPN)
    print ('cellType',srcLayer.cellType )


    changecelltype = anychange = False
    for key, value in kwargs.items():
        if key == 'celltype' and value:
            if value.lower() != srcLayer.cellType.lower():
                changecelltype = True
            if value:
                anychange = True
    if not anychange:
        printstr = 'No metadata to change for %s' %(srcRasterFPN) 
        print ( printstr )
        return
    if changecelltype:
        NOTYET
    else:   
        CopyRasterSetMeta(srcRasterFPN, tmpFPN, gdalformat, **kwargs)
    os.remove(srcRasterFPN)
    os.rename(tmpFPN, srcRasterFPN)
     
def GetRasterMetaData(srcRasterFPN): 
    modeD = {'mode':'read'}
    srcDS,srcLayer = RasterOpenGetFirstLayer(srcRasterFPN,modeD)                      
    #srcRast = RasterDataSource()
    spatialRef = MjProj()
    #TGTODO HERE I ASSUME THAT RASTER PROJ IS IN WKT FORMAT
    spatialRef.SetFromWKT(srcLayer.spatialRef)
    spatialRef.SetProj4()
    #spatialRef.SetProj(srcLayer.spatialRef)
    spatialRef.ReadSpatialRef()
    srcLayer.GetGeometry()
    #close the layer
    srcLayer.CloseLayer()
    srcDS.CloseDS()
    return spatialRef, srcLayer

def SetRasterMeta(srcRasterFPN,*args,**kwargs):
    srcDS,srcLayer = RasterOpenGetFirstLayer(srcRasterFPN,'edit')
    print ( 'xsize', srcLayer.layer.XSize )
    if kwargs is not None:
        for key, value in kwargs.items():
            if key == 'cellnull' and value:
                srcLayer.SetNoDataValue(value)
            elif key == 'palette' and value:
                palette = RasterPalette()
                palette.SetTuplePalette(value)
                srcLayer.ColorTable(palette.colortable)
    #srcLayer.CloseLayer()
    srcDS.CloseDS()
    BALLE
        
def GetFeatureAttributeList(srcShpFPN, fieldL, idfield):
    #Get all fields
    srcDS,srcLayer,srcFieldDefL = ESRIOpenGetLayer(srcShpFPN,'edit')
    fieldNameL = []
    fieldD = {}
    for srcFieldDef in srcFieldDefL:
        fieldNameL.append(srcFieldDef.name)
    for f in fieldL:
        if not f in fieldNameL:
            printstr = 'the requested field %s does not exists in the file %s' %(f, srcShpFPN)
            print ( printstr )
            return
    for feature in srcLayer.layer:
        values_list = [feature.GetField(str(j)) for j in fieldL]
        fid = feature.GetField(str(idfield))
        fieldD[fid] = dict(zip(fieldL,values_list))
    srcDS.CloseDS()
    return fieldD

def ExtractFeaturesToNewDS(srcShpFPN,tarShpFPN,fieldname,valueLL,fieldDefL = False, dictL = False, key = False):
    #Open ESRI DS for reading
    srcDS,srcLayer,srcFieldDefL = ESRIOpenGetLayer(srcShpFPN,'read')
    if not fieldDefL:
        tarFieldDefL = srcFieldDefL
    else:
        tarFieldDefL = fieldDefL
    #Create an instance of datasource for the target


    tarDS,tarLayer = ESRICreateDSLayer(tarShpFPN, srcLayer.proj_cs, srcLayer.geomtype, srcLayer.layerid, fieldDefL)
    #loop over all the features in the source file        
    for feature in srcLayer.layer:
        #create an instance of feature from the srcLayer
        srcFeat = ogrFeature(srcLayer)
        #set the feature
        srcFeat.SetFeature(feature)
        #Set the fieldDef source from the source feature
        srcFeat.AddSourceValtoFieldDef(srcFieldDefL)
        extract = False
        for valueL in valueLL:
            for f in srcFieldDefL:
                if f.name == fieldname and f.source in valueL:
                    extract = True
                if extract:
                    #create an instance of feature for the target
                    tarFeat = ogrFeature(tarLayer) 
                    tarFeat.CopyOgrFeature(feature,tarFieldDefL)  
                    break #as the correct field was found      
    srcDS.CloseDS()
    tarDS.CloseDS()
    
def FeatureFieldsToNewDS(srcShpFPN,tarShpFPN,fieldDefL):
    #Open ESRI DS for reading
    srcDS,srcLayer,srcFieldDefL = ESRIOpenGetLayer(srcShpFPN,'read')
    #tarFieldDefL = fieldDefL
    #Create an instance of datasource for the target
    tarDS,tarLayer = ESRICreateDSLayer(tarShpFPN, srcLayer.proj_cs, srcLayer.geomtype, srcLayer.layerid, fieldDefL)
    #loop over all the features in the source file        
    for feature in srcLayer.layer:
        #create an instance of feature from the srcLayer
        srcFeat = ogrFeature(srcLayer)
        #set the feature
        srcFeat.SetFeature(feature)
        #Set the fieldDef source from the source feature
        srcFeat.AddSourceValtoFieldDef(srcFieldDefL)
        tarFeat = ogrFeature(tarLayer)
        tarFeat.CopyOgrFeature(feature,fieldDefL)    
    srcDS.CloseDS()
    tarDS.CloseDS()

def LinkPtToPoly(ptL,wrsShpFPN):
    import shapely.wkt
    from shapely.geometry import Point
    #Open ESRI DS for reading
    wrsDS,wrsLayer,wrsFieldDefL = ESRIOpenGetLayer(wrsShpFPN,'read')
    polygons = []
    # For each feature in the layer, read the polygon and store it
    for i in range(wrsLayer.layer.GetFeatureCount()):
        # Get the feature, and its path and row attributes
        feature = wrsLayer.layer.GetFeature(i)
        path = feature['PATH']
        row = feature['ROW']
        geom = feature.GetGeometryRef()
        shape = shapely.wkt.loads(geom.ExportToWkt())
        polygons.append((shape, path, row))
    wrsDS.CloseDS()
    #linkD = {}
    #Loop over the points and get the wrs
    for idpt in ptL:
        pt = Point(idpt[1], idpt[2])
        res = []
        # Iterate through every polgon
        for poly in polygons:
            # If the point is within the polygon then
            # append the current path/row to the results
            # list
            if pt.within(poly[0]):
                res.append({'path': poly[1], 'row': poly[2]})
        #linkD[ idpt[0] ] = res
        return res

class LinkRegions:
    def __init__(self,wrsShpFPN,pathField,rowField,wrsrecs = False):
        """The constructoris just an empty container.""" 
        import shapely.wkt
        #from shapely.geometry import Polygon
        #Open ESRI DS for reading
        wrsDS,wrsLayer = ESRIOpenGetLayer(wrsShpFPN,'read')[0:2]
        self.wrspolygons = []
        incl = True
        # For each feature in the layer, read the polygon and store it
        for i in range(wrsLayer.layer.GetFeatureCount()):
            # Get the feature, and its path and row attributes
            feature = wrsLayer.layer.GetFeature(i)
            path = feature[pathField]
            row = feature[rowField]
            if wrsrecs:
                incl = False
                for wrs in wrsrecs:  
                    if wrs[0] == path and wrs[1] == row:
                        incl = True
                        break
            if incl:                
                geom = feature.GetGeometryRef()
                shape = shapely.wkt.loads(geom.ExportToWkt())
                self.wrspolygons.append((shape, path, row))      
        wrsDS.CloseDS()
 
    def OverlayRegion(self,srcShpFPN,simp = 0.005):
        import shapely.wkt
        srcDS,srcLayer= ESRIOpenGetLayer(srcShpFPN,'read')[0:2]
        # For each feature in the layer, read 
        if srcLayer.layer.GetFeatureCount() != 1:
            print ( 'error in featurecount' )
            return False
        res = []

        # Get the feature, and its path and row attributes
        feature = srcLayer.layer.GetFeature(0)
        geom = feature.GetGeometryRef()
        shape = shapely.wkt.loads(geom.ExportToWkt()) 
        #simpShape = shape.simplify(0.05, preserve_topology=False)
        if simp:
            simpShape = shape.simplify(simp, preserve_topology=False)
        else:
            simpShape = shape
        for poly in self.wrspolygons:
            # If the point is within the polygon then
            # append the current path/row to the results
            if simpShape.intersects(poly[0]):
                res.append({'path': poly[1], 'row': poly[2]})          
        srcDS.CloseDS()
        return res
    
    def LinkPtToPoly(self, ptL):
        from shapely.geometry import Point
        #linkD = {}
        res = []
        for idpt in ptL:
            pt = Point(idpt[1], idpt[2])
            
            # Iterate through every polgon
            for poly in self.wrspolygons:
                # If the point is within the polygon then
                # append the current path/row to the results
                # list
                if pt.within(poly[0]):
                    res.append({'path': poly[1], 'row': poly[2], 'ascdesc':'D', 'wrs':'2', 'plotid':idpt[0]})

            if len(res) == 0:
                exitstr = 'No wrs locations found for %s, the default region is not corrext' %()
                sys.exit(exitstr)
            #linkD[ idpt[0] ] = res
        return res
         
def LinkSingleRegionsToWRS(srcShpFPN, wrsShpFPN, simp = 0.005):
    linkwrs = LinkRegionsToWRS(wrsShpFPN)
    return linkwrs.OverlayRegion(srcShpFPN, simp)
    '''
    #from shapely.geometry import Polygon
    #Open ESRI DS for reading
    wrsDS,wrsLayer = ESRIOpenGetLayer(wrsShpFPN,'read')[0:2]
    polygons = []
    # For each feature in the layer, read the polygon and store it
    for i in range(wrsLayer.layer.GetFeatureCount()):
        # Get the feature, and its path and row attributes
        feature = wrsLayer.layer.GetFeature(i)
        path = feature['PATH']
        row = feature['ROW']
        geom = feature.GetGeometryRef()
        shape = shapely.wkt.loads(geom.ExportToWkt())
        polygons.append((shape, path, row))
    wrsDS.CloseDS()
    
    srcDS,srcLayer= ESRIOpenGetLayer(srcShpFPN,'read')[0:2]

    # For each feature in the layer, read 
    res = []
    for i in range(srcLayer.layer.GetFeatureCount()):
        # Get the feature, and its path and row attributes
        feature = srcLayer.layer.GetFeature(i)
        geom = feature.GetGeometryRef()

        shape = shapely.wkt.loads(geom.ExportToWkt()) 
        #simpShape = shape.simplify(0.05, preserve_topology=False)
        if simp:
            simpShape = shape.simplify(simp, preserve_topology=False)
        else:
            simpShape = shape
        for poly in polygons:
            # If the point is within the polygon then
            # append the current path/row to the results
            if simpShape.intersects(poly[0]):
                res.append({'path': poly[1], 'row': poly[2]})          
    srcDS.CloseDS()
    return res
    '''
 
def ReprojectPtOld(pt,cs_src,cs_tar):
    
    transform = osr.CoordinateTransformation(cs_src,cs_tar)
    pt.ogrGeom.Transform(transform) 
    
def ReprojectPt(pt,cs_src,cs_tar):
    ptgeom = ShapelyPointGeom(pt)
    ptgeom.ShapelyToOgrGeom()
    transform = osr.CoordinateTransformation(cs_src,cs_tar)
    ptgeom.ogrGeom.Transform(transform) 
    ptgeom.OgrGeomToShapely()
    return list(ptgeom.shapelyGeom.coords)[0]
    
def LatLonPtToRasterLinCol(pt,cs_tar,gt):
    spatialRef = MjProj()
    spatialRef.SetFromEPSG(4326)
    pt = ShapelyPointGeom(pt)
    pt.ShapelyToOgrGeom()
    ReprojectPt(pt,spatialRef.proj_cs,cs_tar)
    pt.OgrGeomToShapely()
    #gt = (minX, pixelSize, 0, maxY, 0, -pixelSize)
    #lin = int(round( (gt[3]-pt.shapelyGeom.coords[0][1]) / -gt[5]))
    #col = int(round( (pt.shapelyGeom.coords[0][0]-gt[0]) /gt[1]))
    lin = (gt[3]-pt.shapelyGeom.coords[0][1]) / -gt[5]
    col = (pt.shapelyGeom.coords[0][0]-gt[0]) /gt[1]
    return lin,col

def LatLonPtProject(pt,cs_tar,gt):
    spatialRef = MjProj()
    spatialRef.SetFromEPSG(4326)
    pt = ShapelyPointGeom(pt)
    pt.ShapelyToOgrGeom()
    ReprojectPt(pt,spatialRef.proj_cs,cs_tar)
    pt.OgrGeomToShapely()
    return (pt.shapelyGeom.coords[0][0],pt.shapelyGeom.coords[0][1])

def LatLonPtToRasterCoord(pt,cs_tar,gt):
    spatialRef = MjProj()
    spatialRef.SetFromEPSG(4326)
    pt = ShapelyPointGeom(pt)
    pt.ShapelyToOgrGeom()
    ReprojectPt(pt,spatialRef.proj_cs,cs_tar)
    pt.OgrGeomToShapely()
    #gt = (minX, pixelSize, 0, maxY, 0, -pixelSize)
    lin = int(round( (gt[3]-pt.shapelyGeom.coords[0][1]) / -gt[5]))
    col = int(round( (pt.shapelyGeom.coords[0][0]-gt[0]) /gt[1]))
    return lin,col
   
def ReprojectBounds(bounds,cs_src,cs_tar):
    transform = osr.CoordinateTransformation(cs_src,cs_tar)
    ptL = []
    for pt in bounds:
        ptgeom = ShapelyPointGeom(pt)
        ptgeom.ShapelyToOgrGeom()
        #ptgeom.
        #point = ogr.CreateGeometryFromWkt("POINT (1120351.57 741921.42)")
        ptgeom.ogrGeom.Transform(transform)
        ptgeom.OgrGeomToShapely()
        ptcoord = list(ptgeom.shapelyGeom.coords)[0]
        ptL.extend([ptcoord[0],ptcoord[1]])     
    paramL = ['ullon','ullat','urlon','urlat','lrlon','lrlat','lllon','lllat']
    return dict(zip(paramL,ptL))
         
def ConvertPtLHdrToDict(ptL, headerL):
    import datetime
    dictL = []
    for item in ptL:
        litem = list(item)
        for x,i in enumerate(litem):
            if isinstance(i, datetime.date):
                litem[x] = i.strftime("%Y/%m/%d")
            if isinstance(i, datetime.datetime):
                pass
                #should look like this
                #feat['datetimefield'] = '2014/06/07 14:20:00'
        dictL.append( dict(zip(headerL,litem)))
    return dictL

def SetLayerDefFromPGTab(headerL,tabFormat):
    typD = {'text':'string','character':'string','variable character':'string','time without time zone':'time', 'date':'date','double precision':'real','real':'real','integer':'integer'}

    fieldDefL = []
    for tabform in tabFormat:
        #tabformat: name, typ, stringlen, precision
        fieldDefD = {}
        if tabform[0] in headerL:
            #skip if time, not dbf format unless numeric
            #if tabform[1][0:4] == 'time':
                
            #set name
            name = tabform[0]
            #set type
            fieldDefD['type'] = typD[tabform[1]]
            #set the string len
            if tabform[2] == None and fieldDefD['type'] == 'string':
                fieldDefD['width'] = 254
            elif tabform[2] != None:
                fieldDefD['width'] = tabform[2]
            #set the precision:
            if tabform[3] != None:
                if tabform[1] == 'double precision':
                    fieldDefD['precision'] = 16
                elif tabform[1] == 'real':
                    fieldDefD['precision'] = 16
            fieldDefL.append(FieldDef(name,fieldDefD))
    return fieldDefL
           
def SpecimenToLayer(tarShpFPN, ptL, headerL, tabFormat, append, inregion):
    if append and os.path.exists(tarShpFPN):
        tarDS,tarLayer,tarFieldDefL = ESRIOpenGetLayer(tarShpFPN,'edit')
    else:   
        fieldDefL = SetLayerDefFromPGTab(headerL,tabFormat)
        layerid = 'test'
        spatialRef = MjProj()
        spatialRef.SetFromEPSG(4326)
        #Create the shapefile
        tarDS,tarLayer = ESRICreateDSLayer(tarShpFPN, spatialRef.proj_cs, 'point', layerid, fieldDefL)       
    ptDL = ConvertPtLHdrToDict(ptL, headerL)
    tarLayer.AddPtDataFromDict('lon','lat',ptDL)
    tarDS.CloseDS()
        
def CreateEmptyVectorDS(tarShpFPN, fieldDefL, proj_cs, layerid, ogrtype):
    #create the datasource
    tarDS,tarLayer = ESRICreateDSLayer(tarShpFPN, proj_cs, ogrtype, layerid, fieldDefL)
    return tarDS, tarLayer
    '''
    #create the geometry
    #geom = ShapelyPolyGeom(ptL)
    #geom.ShapelyToOgrGeom()
    tarFeat = ogrFeature(tarLayer)
    tarFeat.CreateOgrFeature(geom, fieldDefL)
    return tarDS()

    tarDS.CloseDS()
    fieldDefL = SetLayerDefFromPGTab(headerL,tabFormat)
        layerid = 'test'
        spatialRef = MjProj()
        spatialRef.SetFromEPSG(4326)
        #Create the shapefile
        tarDS,tarLayer = ESRICreateDSLayer(tarShpFPN, spatialRef.proj_cs, 'point', layerid, fieldDefL)
    fittan h       
    ptDL = ConvertPtLHdrToDict(ptL, headerL)
    tarLayer.AddPtDataFromDict('lon','lat',ptDL)
    tarDS.CloseDS()
    '''


if __name__ == '__main__':
    paa
    #TRMM test
    #srcRasterFPN = '/Volumes/WETLAND/climate_download/3B43/3B43.19980101.7.HDF'
