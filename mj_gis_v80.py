
from __future__ import division
import os
import sys

from gdalconst import *
from osgeo import gdal, ogr, osr
from itertools import chain
import array as arr
import numpy as np
#import mj_datetime_v70 as mj_dt
import glob
from os import listdir
#import fiona
import json

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
        #print self.proj_cs
        #print 'epsg',self.proj_cs.GetAttrValue("AUTHORITY", 1)
        self.SetProj4()
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

        '''
        elif self.projection.lower() == 'sinusoidal' and self.datum == 'Not_specified_based_on_custom_spheroid':
            self.epsg = 6842
        elif self.projection.lower() == 'sinusoidal' and self.epsg in ['9901','9001']:
            self.epsg = 6842
        '''
            
    def ReprojectGeom(self, geom, cs_tar):
        #ptgeom = ShapelyPointGeom(pt)
        #ptgeom.ShapelyToOgrGeom()
        #print ('self.proj_cs',self.proj_cs)
        #print ('')
        #print ('cs_tar',cs_tar.proj_cs)

        transform = osr.CoordinateTransformation(self.proj_cs,cs_tar.proj_cs)
        geom.ogrGeom.Transform(transform)
        geom.OgrGeomToShapely()
         
        return geom

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
            print ('Could not create file')
            sys.exit(1)
                   
    def CloseDS(self):
        self.datasource = None
        
class FieldDef:
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
        
class VectorLayer:
    def __init__(self):
        """The constructoris just an empty container.""" 
        
    #VectorDataSource.__init(self)   
    def SetDS(self,DS): 
        self.datasource =  DS.datasource
             
    def GetLayer(self):
        geomTypeD = {1:'point',2:'line',3:'polygon'} 
        self.layer = self.datasource.GetLayer()
        self.geomtype = geomTypeD[self.layer.GetGeomType()]
        self.layerid = self.layer.GetName()
        self.GetLayerDef()
        
    def GetLayerDef(self):
        layerDefn = self.layer.GetLayerDefn()
        fieldDefD = {}  
        self.fieldDefL = []
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
            self.fieldDefL.append(FieldDef(key,fieldDefD[key]))
                   
    def CreateOgrLayer(self,geomtype,tarLayerId):
        #print geomtype
        self.geomtype = geomtype
        if self.geomtype.lower() == 'point':
            self.layer = self.datasource.CreateLayer(tarLayerId, self.spatialRef, geom_type=ogr.wkbPoint) 
        elif self.geomtype.lower() == 'line':
            self.layer = self.datasource.CreateLayer(tarLayerId, self.spatialRef, geom_type=ogr.wkbLineString) 
        elif self.geomtype.lower() == 'polygon':
            self.layer = self.datasource.CreateLayer(tarLayerId, self.spatialRef, geom_type=ogr.wkbPolygon) 
        else:
            exitstr = 'Can not recoginze geometry type %s for creating layer' %(self.geomtype) 
            sys.exit(exitstr)
        if self.layer is None:
            sys.exit('Could not create layer') 
            
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
                printstr = ('field type not recognized %s' %(fieldDef.type) )
                sys.exit(printstr)
            self.layer.CreateField(new_field)
        else:
            for i in range( self.layerDefn.GetFieldCount() ):
                STOPPAUNKNOWN
                if self.layerDefn.GetFieldDefn(i).GetName() == fieldDef.name:
                    fieldTypeCode = self.layerDefn.GetFieldDefn(i).GetType()
                    fieldType = self.layerDefn.GetFieldDefn(i).GetFieldTypeName(fieldTypeCode)
                    fieldWidth = self.layerDefn.GetFieldDefn(i).GetWidth()
                    fieldPrecision = self.layerDefn.GetFieldDefn(i).GetPrecision()

                    if fieldType.lower() != fieldDef.type.lower():
                        sys.exit('Existing field of different type already exists')
         
    def AddField(self,f):
        pass
              
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
                tarfeat.SetField(key, D[key])         
            self.layer.CreateFeature(tarfeat) 
       
    def GetSpatialRef(self): 
        self.spatialRef = self.layer.GetSpatialRef() 

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
                        print ('missing', item,feature.GetField(srcFieldName))
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
        tarfeat.SetGeometry(geom)
        
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
        self.ShapelyToOgrGeom()
        
    def LoadWKTGeom(self,geom):
        import shapely.wkt
        self.shapelyGeom = shapely.wkt.loads(geom)
        self.ShapelyToOgrGeom()
        
    def BoundsToPoly(self):
        from shapely.geometry import Polygon
        minx, miny, maxx, maxy = self.shapelyGeom.bounds
        ptLT = ( (minx,maxy), (maxx,maxy), (maxx,miny), (minx,miny), (minx,maxy) )
        return Polygon(ptLT)
    
    def PointsToMultiPointGeom(self,ptL):
        from shapely.geometry import Point,MultiPoint
        ptgeomL = []
        for pt in ptL:
            ptgeomL.append(Point(pt))
        self.shapelyGeom = MultiPoint(ptgeomL)
        self.ShapelyToOgrGeom()
        
    def PointsToPolygonGeom(self,ptL):
        from shapely.geometry import Polygon
        #ptgeomL = []
        #for pt in ptL:
        #    ptgeomL.append(Point(pt))
        self.shapelyGeom = Polygon(ptL)
        #self.shapelyGeom = Polygon(ptgeomL)
        self.ShapelyToOgrGeom()
        
    def CascadeUnion(self,shapelyGeomL):
        from shapely.ops import cascaded_union
        self.SetShapelyGeom(cascaded_union(shapelyGeomL))
        #return cascaded_union(shapelyGeomL)
        
    def MultPolyToSinglePoly(self,multipoly):
        from shapely.ops import cascaded_union
        from shapely.geometry import Polygon
        '''
        self.shapelyGeom = cascaded_union([
            Polygon(component.exterior) for component in multipoly
        ])
        '''
        polyL = [item for item in multipoly]
        #print ( len(polyL) )
        self.shapelyGeom = polyL[0]
        for item in polyL:
            self.shapelyGeom.union(item)
        self.ShapelyToOgrGeom()
        
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
        
    def ShapelyIntersection(self, otherGeom): 
        return self.shapelyGeom.intersection(otherGeom.shapelyGeom)
        
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
        
    def ShapelyWithin(self,cut):
        if self.shapelyGeom.within(cut.shapelyGeom):
            return True
        else:
            return False
        
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

def ESRIOpenGetLayer(srcShpFPN,mode='read'):
    '''Opens a standard ESRI datasource and gets the layer and layerDef
    '''
    #Create an instance of datasource for the source data
    srcDS = VectorDataSource()
    #open the source data for reading
    if mode[0] == 'r':

        srcDS.OpenESRIshapeRead(srcShpFPN)
    else:
        srcDS.OpenESRIshapeEdit(srcShpFPN)
    #create a layer instance
    srcLayer = VectorLayer()
    #set the datasource of the layer
    srcLayer.SetDS(srcDS)
    #get the layer and the fieldDefintions 
    srcLayer.GetLayer()
    #get the spatialref
    srcLayer.GetSpatialRef()
    return srcDS,srcLayer

def CreateESRIPolygonGeom(tarShpFPN, fieldDefL, geom, proj_cs, layerid):
    #create the datasource
    tarDS,tarLayer = ESRICreateDSLayer(tarShpFPN, proj_cs, 'polygon', layerid, fieldDefL)
    #create the geometry
    geom.ShapelyToOgrGeom()
    tarFeat = ogrFeature(tarLayer)
    tarFeat.CreateOgrFeature(geom.ogrGeom, fieldDefL)
    tarDS.CloseDS()

def CreateVectorAttributeDef(fieldDD): 
        fieldDefD = {}
        fieldDefL =[]
        for key in fieldDD:

            fieldD = fieldDD[key]
            if 'width' in fieldD:
                width = fieldD['width']
            else:
                width = 8
            if 'precision' in fieldD:
                precision = fieldD['precision']
            else:
                precision = 0
            if 'keyfield' in fieldD:
                keyfield = fieldD['keyfield']
            elif 'field' in fieldD:
                keyfield = fieldD['field']
            else:
                keyfield = False
            fieldDefD[key] = {'type':fieldD['type'].lower(), 'width':width, 'precision': precision, 'transfer': fieldD['transfer'].lower(), 'source':fieldD['source'], 'keyfield':keyfield}      

        for key in fieldDefD:
            fieldDefL.append(FieldDef(key,fieldDefD[key]))
        return fieldDefL
     
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
    
def ImportKMLtoShape(kmlFPN,shpFPN):
    oscmd = '/Library/Frameworks/GDAL.framework/Versions/2.1/Programs/ogr2ogr -skipfailures %(dst)s %(src)s' %{'dst':shpFPN, 'src':kmlFPN}

    ERRORCHECK
    
def ExportToGeoJson(srcFPN,dstFPN):
    features = []
    crs = None
    with fiona.collection(srcFPN, "r") as source:
        for feat in source:
        #    feat['properties'].update(...) # with your attributes
            features.append(feat)
        crs = " ".join("+%s=%s" % (k,v) for k,v in source.crs.items())
    
    my_layer = {
        "type": "FeatureCollection",
        "features": features,
        "crs": {
            "type": "link", 
            "properties": {"href": "my_layer.crs", "type": "proj4"} }}
    
    with open(dstFPN, "w") as f:
        f.write(json.dumps(my_layer))
    #with open("my_layer.crs", "w") as f:
    #    f.write(crs)
     
if __name__ == '__main__':

    #Create the target file
    tarVectorFPN = '/Volumes/karttur2tb/sites/All_sites3_4326.shp'
    tarProj = MjProj()
    #Set target projection to LatLon
    tarProj.SetFromProj4('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')
    #print ('tarproj',tarProj.proj_cs)

    
    #spatialRef = MjProj()
    #spatialRef.SetFromEPSG(4326)
    #pt = ShapelyPointGeom(pt)
    #pt.ShapelyToOgrGeom()
    #ReprojectPt(pt,spatialRef.proj_cs,cs_tar)

    fieldD = {}
    #for key in self.process.xparamTagD['fielddef']:
    #    f = self.process.xparamTagD['fielddef'][key]
    fieldD['id'] = {'type':'integer','width':8,'precision':0,'transfer':'constant','source':0 }
    fieldD['origin'] = {'type':'string','width':32,'precision':0,'transfer':'constant','source':'ABext_LAEA' }
    fieldD['searchid'] = {'type':'string','width':128,'precision':0,'transfer':'constant','source':'any' }
    #fieldD['area(ha)'] = {'type':'real','width':32,'precision':2,'transfer':'constant','source':'ABext_LAEA' }
        
    fieldDefL = CreateVectorAttributeDef(fieldD)
    
    

    #tarDS,tarLayer = ESRICreateDSLayer(tarVectorFPN, tarProj.proj_cs, srcLayer.geomtype, srcLayer.layerid, fieldDefL)
    tarDS,tarLayer = ESRICreateDSLayer(tarVectorFPN, tarProj.proj_cs, 'polygon', 'wetland', fieldDefL)
    
    
    #srcRasterFPN = '/Volumes/WETLAND/climate_download/3B43/3B43.19980101.7.HDF'
    #srcVectorFPN = '/Volumes/karttur2tb/data_raw/shp_arctic_wetland_sites/ABext_LAEA.shp'
    srcVectorFP = '/Volumes/karttur2tb/data_raw/shp_arctic_wetland_sites'
    
    searhPath = '%(s)s.*.shp' %{'s':srcVectorFP}
    shpFPNL = glob.glob("searhPath")

    shpFNL = [f for f in os.listdir(srcVectorFP) if os.path.isfile(os.path.join(srcVectorFP, f)) and os.path.splitext(f)[1] == '.shp']
    count = -1
    for shpFN in shpFNL:
        count += 1
        #if shpFN == 'ABext_LAEA.shp':
        #    continue
        print ('opening',shpFN)
        shpFPN = os.path.join(srcVectorFP,shpFN)
        srcDS,srcLayer,fieldDefL = ESRIOpenGetLayer(shpFPN)
        

        fieldD = {}
        #for key in self.process.xparamTagD['fielddef']:
        #    f = self.process.xparamTagD['fielddef'][key]
        fieldD['id'] = {'type':'integer','width':8,'precision':0,'transfer':'constant','source':count }
        fieldD['origin'] = {'type':'string','width':32,'precision':0,'transfer':'constant','source':shpFN }
        ktid = '%(origin)s-%(id)d' %{'origin':shpFN,'id':count}
        fieldD['searchid'] = {'type':'string','width':128,'precision':0,'transfer':'constant','source':ktid }
        #fieldD['area(ha)'] = {'type':'real','width':32,'precision':2,'transfer':'constant','source':'ABext_LAEA' }
        
        fieldDefL = CreateVectorAttributeDef(fieldD)

        srcDS.CloseDS()

        srcProj = MjProj()
        srcProj.SetProj(srcLayer.spatialRef)
        #srcProj.SetProj(srcLayer.spatialRef)
        #print ('srcLayer.spatialRef',srcLayer.spatialRef)
        #srcProj.SetFromWKT(srcLayer.spatialRef)
        #srcProj = srcLayer.spatialRef
        #srcProj.SetFromProj4('+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs')
        #print srcProj.proj4
        #print ('srcProj.proj_cs',srcProj.proj_cs)
        
        for feature in srcLayer.layer:
            #feature = srcLayer.layer.GetFeature(0)
            #create an instance of feature from the srcLayer
            #srcFeat = ogrFeature(srcLayer)
            #get the src feature
            #srcFeat.SetFeature(feature)
            geom = feature.GetGeometryRef()
            print ('geom',geom)
            #print 'srcLayer.spatialRef.ReprojectGeom',srcLayer.spatialRef.ReprojectGeom
            geom = srcProj.ReprojectGeom(geom, tarProj.proj_cs)
            #print ('srcLayer.spatialRef',srcLayer.spatialRef)
            
            #print ('geom',geom)

            #Set the fieldDef 
            tarFeat = ogrFeature(tarLayer)
            tarFeat.CreateOgrFeature(geom, fieldDefL)
            #tarFeat.CreateOgrFeature(geom, fieldDefL)
        
            #srcFeat.AddSourceValtoFieldDef(fieldDefL)
            '''
            for valueL in valueLL:
                for f in fieldDefL:
            '''
            #create an instance of feature for the target
            #tarFeat = ogrFeature(tarLayer) 
            #tarFeat.CopyOgrFeature(feature,fieldDefL)  
    

    tarDS.CloseDS()
