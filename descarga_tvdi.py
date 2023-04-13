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

from funciones import descarga, mosaico, filldata, recorte, recortexTIFF
from funciones import subirtiff, filtro_media
from funciones import hdfToTiff, getNodata

tiles = ["h09v09","h10v07","h10v08","h10v09","h10v10","h11v07","h11v08","h11v09","h11v10","h11v11","h11v12","h12v08","h12v09","h12v10","h12v11","h12v12","h12v13","h13v09","h13v10","h13v11","h13v12","h13v13","h13v14","h14v09","h14v10","h14v11","h14v14"]

url = "https://ladsweb.modaps.eosdis.nasa.gov/archive/allData/61/" #Dirección de descarga de los datos
prods = ["MOD11A2","MOD13A1"] #Productos requeridos
#appkey = "51B6A2FE-4FC6-11E9-84A2-8E39E106C194" #Clave de autenticación al servidor de descarga
appkey = "ai5lLnJ1YmlvOmFuSjFZbWx2UUdOdmJtRmxMbWR2ZGk1aGNnPT06MTYyODI1MjU2MTpiYTI4ZmMyZDMyNjc0Y2MzYjc0NmFkODUzMWJjOWI1YmNkMmNhMTA5"

# path inputs
###path = "/mnt/datos/Repositorio/sancor/" #Directorio principal de procesamiento
path = "/home/jrubio/Documentos/VENVS/SANCOR5/" #Directorio principal de procesamiento
path_estadisticos = path + "estadisticos/"
path_anomalias = path + "anomalias/"
path_mask = path + "datos/mascara/"

shp = path_mask + "mask.shp" #Archivo shapefile con máscara de zona de interés
uw_mask = path_mask + "mask_urbano_agua_i5.tif" #Raster máscara zona de zonas urbanas y cuerpos de agua
uw_mask2 = path_mask + "mask_urbano_agua_5.tif" #Raster máscara zona de zonas urbanas y cuerpos de agua
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
#if not os.path.exists(path_indices+"lst8/"):os.mkdir(path_indices+"lst/")

for aa in range(2023,2024):
	for d in range(65,66,16):
		fecha = str(aa)+"%03d"%d
		fecha_lst8 = (datetime.strptime(fecha, "%Y%j") + timedelta(days=8)).strftime("%Y%j")

		###########################################################################################################
		##########################################Cálculo del índice TVDI##########################################
		###########################################################################################################
		#Descargar datos
#3#		descarga(url,appkey,tiles,prods,fecha,path_descarga,path_xml)
#3#		descarga(url,appkey,tiles,[prods[0]],fecha_lst8,path_descarga,path_xml)

		#Reproyectar capas de vegetación y temperatura superficial
		lista_lst = glob.glob("%s/%s/%s/%s/*%s*.hdf"%(path_descarga,prods[0],fecha[0:4],fecha[4:7],fecha))
		lista_lst8 = glob.glob("%s/%s/%s/%s/*%s*.hdf"%(path_descarga,prods[0],fecha_lst8[0:4],fecha_lst8[4:7],fecha_lst8))
		lista_veg = glob.glob("%s/%s/%s/%s/*%s*.hdf"%(path_descarga,prods[1],fecha[0:4],fecha[4:7],fecha))

		for vnp in lista_lst:
			resLST = hdfToTiff(vnp, path_proc+'lst_temp/', dataset[0], lst_nodata)

		for vnp in lista_lst8:
			res = hdfToTiff(vnp, path_proc+'lst_temp/', dataset[0], lst_nodata)

		for vnp in lista_veg:
			res = hdfToTiff(vnp, path_proc+'ndvi_temp/', dataset[1], ndvi_nodata)

		#Realizar el Mosaico de las imágenes
		mosaico(glob.glob(path_proc+"lst_temp/*"+fecha+".*.tif"),path_proc+"mosaicos/mosaico_"+fecha+"_lst.tif",lst_nodata,resxLST,resyLST)
		mosaico(glob.glob(path_proc+"lst_temp/*"+fecha_lst8+".*.tif"),path_proc+"mosaicos/mosaico_"+fecha_lst8+"_lst.tif",lst_nodata,resxLST,resyLST)
		mosaico(glob.glob(path_proc+"ndvi_temp/*"+fecha+".*.tif"),path_proc+"mosaicos/mosaico_"+fecha+"_ndvi.tif",ndvi_nodata,resx,resy)

		#Interpolación de líneas faltantes
		lst_fill = path_proc+"mosaicos/mosaico_"+fecha+"_lst_fill.tif"
		lst8_fill = path_proc+"mosaicos/mosaico_"+fecha_lst8+"_lst_fill.tif"
		ndvi_fill = path_proc+"mosaicos/mosaico_"+fecha+"_ndvi_fill.tif"

		filldata(path_proc+"mosaicos/mosaico_"+fecha+"_lst.tif", lst_fill)
		filldata(path_proc+"mosaicos/mosaico_"+fecha_lst8+"_lst.tif", lst8_fill)
		filldata(path_proc+"mosaicos/mosaico_"+fecha+"_ndvi.tif", ndvi_fill)

		#Escalar datos
		lst_escal = path_indices+"lst/lst_Celsius_"+fecha+"_500m.tif"
		lst8_escal = path_indices+"lst/lst_Celsius_"+fecha_lst8+"_500m.tif"
		ndvi_escal = path_indices+"ndvi/ndvi_"+fecha+"_500m.tif"

		os.system('gdal_calc.py -A %s --outfile=%s --calc="A*%.2f-273.15" --NoDataValue %d --type="Float32" --overwrite'%(lst_fill, lst_escal, factor_lst, lst_nodata))
		os.system('gdal_calc.py -A %s --outfile=%s --calc="A*%.2f-273.15" --NoDataValue %d --type="Float32" --overwrite'%(lst8_fill, lst8_escal, factor_lst, lst_nodata))
		os.system('gdal_calc.py -A %s --outfile=%s --calc="A*%.4f" --NoDataValue %d --type="Float32" --overwrite'%(ndvi_fill, ndvi_escal, factor_ndvi, ndvi_nodata))

		#Rellenar datos
		temp1 = lst_escal[:-4] + "_temp1.tif"
		temp2 = lst_escal[:-4] + "_temp2.tif"

		#Llevar de 1000 a 500mt y al tamaño adecuado
		os.system('gdalwarp -te -99.9999999747252275 -60.0015558891651537 -29.9986385732798908 19.9999999982009413 -tr %s %s %s %s '%(resx,resy,lst_escal,temp1))
		os.system('gdal_translate -outsize 13762 15728 %s %s'%(temp1,lst_escal))
		os.system("rm " + temp1)
		os.system('gdalwarp -te -99.9999999747252275 -60.0015558891651537 -29.9986385732798908 19.9999999982009413 -tr %s %s %s %s '%(resx,resy,lst8_escal,temp1))
		os.system('gdal_translate -outsize 13762 15728 %s %s'%(temp1,lst8_escal))
		os.system("rm " + temp1)

		#Corregir temperatura
		lst_corr = path_indices+"lst/lst_Celsius_"+fecha+"_500m_corr.tif"
		lst8_corr = path_indices+"lst/lst_Celsius_"+fecha_lst8+"_500m_corr.tif"
		os.system('gdal_calc.py -A %s -B %s --outfile=%s --calc="A+B*0.006" --NoDataValue %d --type="Float32" --overwrite'%(lst_escal, dem, lst_corr, lst_nodata))#La corrección de lst 6grados/km vertical
		os.system('gdal_calc.py -A %s -B %s --outfile=%s --calc="A+B*0.006" --NoDataValue %d --type="Float32" --overwrite'%(lst8_escal, dem, lst8_corr, lst_nodata))#La corrección de lst 6grados/km vertical

		#Enmascarado
		recortexTIFF(lst_escal, uw_mask, lst_nodata, lst_escal, "Float32")
		recortexTIFF(lst_corr, uw_mask, lst_nodata, lst_corr, "Float32")
		recortexTIFF(lst8_escal, uw_mask, lst_nodata, lst8_escal, "Float32")
		recortexTIFF(lst8_corr, uw_mask, lst_nodata, lst8_corr, "Float32")
		recortexTIFF(ndvi_escal, uw_mask, ndvi_nodata, ndvi_escal, "Float32")

		#Comprimir archivos TIFF
#		os.system('gdal_calc.py -A %s --outfile=%s --calc="A*1" --NoDataValue %d --type="Float32" --co="COMPRESS=DEFLATE" --overwrite'%(lst_escal, lst_escal, lst_nodata))
#		os.system('gdal_calc.py -A %s --outfile=%s --calc="A*1" --NoDataValue %d --type="Float32" --co="COMPRESS=DEFLATE" --overwrite'%(lst_corr, lst_corr, lst_nodata))
#		os.system('gdal_calc.py -A %s --outfile=%s --calc="A*1" --NoDataValue %d --type="Float32" --co="COMPRESS=DEFLATE" --overwrite'%(lst8_escal, lst8_escal, lst_nodata))
#		os.system('gdal_calc.py -A %s --outfile=%s --calc="A*1" --NoDataValue %d --type="Float32" --co="COMPRESS=DEFLATE" --overwrite'%(lst8_corr, lst8_corr, lst_nodata))
#		os.system('gdal_calc.py -A %s --outfile=%s --calc="A*1" --NoDataValue %d --type="Float32" --co="COMPRESS=DEFLATE" --overwrite'%(ndvi_escal, ndvi_escal, ndvi_nodata))
		
		#RECORTAR LST Y LLEVAR A ENTEROS
		lst_mask = path_indices+"lst/lst_"+fecha+"_500m_mask.tif"
		lst_int = path_indices+"lst/lst_"+fecha+"_500m_int.tif"
		lst_nodata = getNodata(lst_corr).astype(int)
		os.system('gdal_calc.py -A %s --outfile=%s --calc="A*1" --NoDataValue=%s --type="Int16" --overwrite'%(lst_corr, lst_int, lst_nodata))
		recortexTIFF(lst_int, uw_mask, lst_nodata, lst_mask, "Int16")
		os.system('gdal_calc.py -A %s --outfile=%s --calc="A*1" --NoDataValue=%s --type="Int16" --overwrite'%(lst_corr, lst_int, lst_nodata))
		recortexTIFF(lst_int, uw_mask, lst_nodata, lst_mask, "Int16")
		os.system("rm " + lst_int)

		#RECORTAR NDVI Y LLEVAR A ENTEROS
		ndvi_mask = path_indices+"ndvi/ndvi_"+fecha+"_500m_mask.tif"
		ndvi_int = path_indices+"ndvi/ndvi_"+fecha+"_500m_int.tif"
		ndvi_nodata = getNodata(ndvi_escal).astype('int')
		os.system('gdal_calc.py -A %s --outfile=%s --calc="A*1000" --NoDataValue=%s --type="Int16" --overwrite'%(ndvi_escal, ndvi_int, ndvi_nodata))
		recortexTIFF(ndvi_int, uw_mask, ndvi_nodata, ndvi_mask, "Int16")
		os.system("rm " + ndvi_int)
		

		#Borrar archivos temporales
		print("Borrando archivos temporales")
		os.system("rm " + path_indices + "lst/*temp*.tif")
		os.system("rm " + path_proc + "lst_temp.tif")
		os.system("rm " + path_proc + "ndvi_temp/*.tif")
		os.system("rm " + path_proc + "mosaicos/*.tif")
