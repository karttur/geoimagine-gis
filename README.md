# geoimagine-gis

Karttur Geoimagine GIS python project

## Introduction

Karttur's GeoImagine Framework is an attempt for semi-automated processing of spatial data. 
The GIS package contains the processing for handling translations between spatial data formats, resolutions and 
projections. 


The package contains 4 modules:


- \_\_init.py\_\_
- gis.py
- mj\_gis\_v80.py
- version.py

### gis.py

gis.py is an older version, but that is still called by some other GeoImagine Framework packages.

#### GeoImagine dependencies

None

#### Classes

- Composition
- RasterPalette
- FieldDef
- MjProj
- RasterDataSource
- RasterLayer
- VectorDataSource
- VectorLayer
- ogrFeature
- Geometry
- ShapelyPointGeom(Geometry)
- ShapelyMultiPointGeom(Geometry)
- ShapelyLineGeom(Geometry)
- ShapelyPolyGeom(Geometry)

### mj\_gis\_v80.py

#### GeoImagine dependencies

None

#### Classes

- VectorDataSource
- VectorLayer
- ogrFeature
- Geometry