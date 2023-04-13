###!/usr/bin/env python
# -*- coding: utf-8 -*-

####################################################################################################
####################################################################################################
#Este script realiza las siguientes funciones:
#1. Descargar datos MODIS de vegetación (NDVI) y temperatura superficial (LST)
#2. Reproyectar los datos a EPSG:4326
#3. Generar Mosaicos de los datos de NDVI y LST para Argentina, Uruguay y Paraguay
#4. Calcular TVDI
#############################################01-04-2019#############################################
####################################################################################################


import os
import glob
import sys
from osgeo import gdal
import numpy as np
from datetime import datetime, timedelta

from funciones import descarga, mosaico, filldata, recorte, recortexTIFF, getIntMaskImg, getNodata
from funciones import calc_indice, subirtiff, filtro_media, lead_to_tvdi
from funciones import hdfToTiff

# path inputs
###path = "/mnt/datos/Repositorio/sancor/" #Directorio principal de procesamiento
path = "/home/jrubio/Documentos/VENVS/SANCOR5/" #Directorio principal de procesamiento
path_estadisticos = path + "estadisticos/"
path_anomalias = path + "anomalias/"
path_mask = path + "datos/mascara/"

shp = path_mask + "mask.shp" #Archivo shapefile con máscara de zona de interés
uw_mask = path_mask + "mask_urbano_agua_i5.tif" #Raster máscara zona de zonas urbanas y cuerpos de agua
uw_mask2 = path_mask + "mask_urbano_agua_5.tif" #Raster máscara zona de zonas urbanas y cuerpos de agua
#uw_mask = path_mask + "mask_urbano_agua_i3.tif" #Raster máscara zona de zonas urbanas y cuerpos de agua
#uw_mask2 = path_mask + "mask_urbano_agua_3.tif" #Raster máscara zona de zonas urbanas y cuerpos de agua
dem = path + "datos/dem500m/dem_500m.tif" #Modelo de Elevación Digital para corrección de temperatura superficial

#path outputs
path_descarga = path + "datos/crudos/" #Directorio para descarga de información
path_proc = path + "procesamiento/" #Directorio para almacenar datos temporales de procesamiento
path_xml = path + "datos/xml/" #Directorio con archivos xml para reproyección de las imágenes
path_indices = path + "indices/" #Directorio donde se almacenan los índices calculados


ndvi_nodata = -3000
lst_nodata = 0
#red_nodata = -28672
#mir_nodata = -28672

factor_ndvi = 0.0001 #Factor de escala para el producto MOD11A2
factor_lst = 0.02 #Factor de escala para el producto MOD13A1
resx,resy = 0.005086568914507000154,0.005086568914507000154 #Resolución para remuestreo de datos de temperatura
resxLST,resyLST = 0.008626647658857677925,0.008626647658857677925

dataset = ['LST_Day_1km', '500m 16 days NDVI']

#Crear rutas para procesamiento de datos
if not os.path.exists(path+"datos/"):os.mkdir(path+"datos/")
if not os.path.exists(path+"datos/crudos/"):os.mkdir(path+"datos/crudos/")
if not os.path.exists(path+"datos/xml/"):os.mkdir(path+"datos/xml/")
if not os.path.exists(path+"indices/"):os.mkdir(path+"indices/")

if not os.path.exists(path+"procesamiento/"):os.mkdir(path+"procesamiento/")
if not os.path.exists(path_proc+"lst_temp/"):os.mkdir(path_proc+"lst_temp/")
if not os.path.exists(path_proc+"lst8_temp/"):os.mkdir(path_proc+"lst8_temp/")
if not os.path.exists(path_proc+"ndvi_temp/"):os.mkdir(path_proc+"ndvi_temp/")
if not os.path.exists(path_proc+"mosaicos/"):os.mkdir(path_proc+"mosaicos/")
if not os.path.exists(path_proc+"recortes/"):os.mkdir(path_proc+"recortes/")

#Crear rutas para almacenamiento de información
#if not os.path.exists(path_indices+"ndwi/"):os.mkdir(path_indices+"ndwi/")
if not os.path.exists(path_indices+"tvdi/"):os.mkdir(path_indices+"tvdi/")
if not os.path.exists(path_indices+"ndvi/"):os.mkdir(path_indices+"ndvi/")
if not os.path.exists(path_indices+"lst/"):os.mkdir(path_indices+"lst/")

for aa in range(2023,2024):
	for d in range(65,66,8):
		fecha = str(aa)+"%03d"%d
		fecha_lst8 = (datetime.strptime(fecha, "%Y%j") + timedelta(days=8)).strftime("%Y%j")

		###########################################################################################################
		######################################### Cálculo del índice TVDI #########################################
		###########################################################################################################

		############################
		#Escalar datos
		lst_escal = path_indices+"lst/lst_Celsius_" + fecha + "_500m_fill.tif"
		lst8_escal = path_indices+"lst/lst_Celsius_" + fecha_lst8 + "_500m_fill.tif"
		ndvi_escal = path_indices+"ndvi/ndvi_" + fecha + "_500m.tif"

		print("""
--------------------
Recortando ndvi escalado
--------------------
		""")
		recortexTIFF(ndvi_escal, uw_mask, ndvi_nodata, ndvi_escal, "Float32")
		
		lead_to_tvdi(ndvi_escal, lst_escal, fecha)
		lead_to_tvdi(ndvi_escal, lst8_escal, fecha_lst8)

