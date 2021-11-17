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
import gdal
import numpy as np
from datetime import datetime, timedelta

from funciones import descarga, mosaico, filldata, recorte
from funciones import fn_tvdi, subirtiff, filtro_media
from funciones import hdfToTiff, fillgaps, anomalia

tiles = ["h11v11","h11v12","h12v10","h12v11","h12v12","h12v13","h13v11","h13v12","h13v13","h13v14","h14v14"]

url = "https://ladsweb.modaps.eosdis.nasa.gov/archive/allData/6/" #Dirección de descarga de los datos
prods = ["MOD11A2","MOD13A1","MOD09A1"] #Productos requeridos
#appkey = "51B6A2FE-4FC6-11E9-84A2-8E39E106C194" #Clave de autenticación al servidor de descarga
appkey = "ai5lLnJ1YmlvOmFuSjFZbWx2UUdOdmJtRmxMbWR2ZGk1aGNnPT06MTYyODI1MjU2MTpiYTI4ZmMyZDMyNjc0Y2MzYjc0NmFkODUzMWJjOWI1YmNkMmNhMTA5"

# path inputs
path = "/mnt/datos/Repositorio/sancor/" #Directorio principal de procesamiento
path_estadisticos = path + "estadisticos/"
path_anomalias = path + "anomalias/"
shp = path + "datos/mascara/mask.shp" #Archivo shapefile con máscara de zona de interés
uw_mask = path + "datos/mascara/mask_urbano_agua_i3.tif" #Raster máscara zona de zonas urbanas y cuerpos de agua
uw_mask2 = path + "datos/mascara/mask_urbano_agua_3.tif" #Raster máscara zona de zonas urbanas y cuerpos de agua
dem = path + "datos/dem500m/dem_500m.tif" #Modelo de Elevación Digital para corrección de temperatura superficial

#path outputs
path_descarga = path + "datos/crudos/" #Directorio para descarga de información
path_proc = path + "procesamiento/" #Directorio para almacenar datos temporales de procesamiento
path_xml = path + "datos/xml/" #Directorio con archivos xml para reproyección de las imágenes
path_indices = path + "indices/" #Directorio donde se almacenan los índices calculados


ndvi_nodata = -3000
lst_nodata = 0
red_nodata = -28672
mir_nodata = -28672

factor_ndvi = 0.0001 #Factor de escala para el producto MOD11A2
factor_lst = 0.02 #Factor de escala para el producto MOD13A1
#factor_red = 0.0001 #Factor de escala para el producto MOD11A2
#factor_mir = 0.0001 #Factor de escala para el producto MOD11A2
resx,resy = 0.005086568914507,0.005086568914507 #Resolución para remuestreo de datos de temperatura
resxLST,resyLST = 0.008626647658857677925,0.008626647658857677925


# create session with the user credentials that will be used to authenticate access to the data
username = "xxx"
password= "xxx"

dataset = ['LST_Day_1km', '500m 16 days NDVI']

#Estilos
estilo_lst = "lst"
estilo_ndvi = "ndvi"
estilo_tvdi = "tvdi"
estilo_anom = "anomaliasTVDI"
estilo_ndwi = "ndwi"

#Workspace
ws_lst = "lst"
ws_ndvi = "ndvi"
ws_tvdi = "tvdi"
ws_anom = "anomaliasTVDI"
ws_ndwi = "ndwi"

#Datos Conexión
user = "xxx" #Usuario
passw = "xxx" #Contraseña
ipgeo = "xx.xx.xx.xx" #IP del servidor
port = "8080" #Puerto

#Crear rutas para procesamiento de datos
if not os.path.exists(path+"datos/"):os.mkdir(path+"datos/")
if not os.path.exists(path+"datos/crudos/"):os.mkdir(path+"datos/crudos/")
if not os.path.exists(path+"datos/xml/"):os.mkdir(path+"datos/xml/")
if not os.path.exists(path+"indices/"):os.mkdir(path+"indices/")

if not os.path.exists(path+"procesamiento/"):os.mkdir(path+"procesamiento/")
if not os.path.exists(path_proc+"lst_temp/"):os.mkdir(path_proc+"lst_temp/")
if not os.path.exists(path_proc+"lst8_temp/"):os.mkdir(path_proc+"lst8_temp/")
if not os.path.exists(path_proc+"ndvi_temp/"):os.mkdir(path_proc+"ndvi_temp/")
if not os.path.exists(path_proc+"red_temp/"):os.mkdir(path_proc+"red_temp/")
if not os.path.exists(path_proc+"mir_temp/"):os.mkdir(path_proc+"mir_temp/")
if not os.path.exists(path_proc+"red8_temp/"):os.mkdir(path_proc+"red8_temp/")
if not os.path.exists(path_proc+"mir8_temp/"):os.mkdir(path_proc+"mir8_temp/")
if not os.path.exists(path_proc+"mosaicos/"):os.mkdir(path_proc+"mosaicos/")
if not os.path.exists(path_proc+"recortes/"):os.mkdir(path_proc+"recortes/")

#Crear rutas para almacenamiento de información
if not os.path.exists(path_indices+"ndwi/"):os.mkdir(path_indices+"ndwi/")
if not os.path.exists(path_indices+"tvdi/"):os.mkdir(path_indices+"tvdi/")
if not os.path.exists(path_indices+"ndvi/"):os.mkdir(path_indices+"ndvi/")
if not os.path.exists(path_indices+"lst/"):os.mkdir(path_indices+"lst/")

fecha = sys.argv[1] #Fecha de procesamiento
fecha_lst8 = (datetime.strptime(fecha, "%Y%j") + timedelta(days=8)).strftime("%Y%j")


###########################################################################################################
##########################################Cálculo del índice TVDI##########################################
###########################################################################################################
#Descargar datos
descarga2(url,username,password,tiles,prods,fecha,path_descarga,path_xml)
descarga2(url,username,password,tiles,[prods[0],prods[2]],fecha_lst8,path_descarga,path_xml)

#Reproyectar capas de vegetación y temperatura superficial

lista_lst8 = glob.glob("%s/%s/%s/%s/*%s*.hdf"%(path_descarga,prods[0],fecha_lst8[0:4],fecha_lst8[4:7],fecha_lst8))
lista_lst = glob.glob("%s/%s/%s/%s/*%s*.hdf"%(path_descarga,prods[0],fecha[0:4],fecha[4:7],fecha))
lista_veg = glob.glob("%s/%s/%s/%s/*%s*.hdf"%(path_descarga,prods[1],fecha[0:4],fecha[4:7],fecha))
lista_ndwi = glob.glob("%s/%s/%s/%s/*%s*.hdf"%(path_descarga,prods[2],fecha[0:4],fecha[4:7],fecha))
lista_ndwi8 = glob.glob("%s/%s/%s/%s/*%s*.hdf"%(path_descarga,prods[2],fecha_lst8[0:4],fecha_lst8[4:7],fecha_lst8))

for vnp in lista_lst:
	resLST = hdfToTiff(vnp, path_proc+'lst_temp/', dataset[0], lst_nodata,1)

for vnp in lista_lst8:
	res = hdfToTiff(vnp, path_proc+'lst8_temp/', dataset[0], lst_nodata,1)

for vnp in lista_veg:
	res = hdfToTiff(vnp, path_proc+'ndvi_temp/', dataset[1], ndvi_nodata,1)

for vnp in lista_ndwi:
	res = hdfToTiff(vnp, path_proc+'red_temp/', dataset[2], red_nodata,1)

for vnp in lista_ndwi8:
	res = hdfToTiff(vnp, path_proc+'red8_temp/', dataset[2], red_nodata,1)

for vnp in lista_ndwi:
	res = hdfToTiff(vnp, path_proc+'mir_temp/', dataset[2], mir_nodata,6)

for vnp in lista_ndwi8:
	res = hdfToTiff(vnp, path_proc+'mir8_temp/', dataset[2], mir_nodata,6)


#Realizar el Mosaico de las imágenes
mosaico(glob.glob(path_proc+"lst8_temp/*"+fecha_lst8+".*.tif"),path_proc+"mosaicos/mosaico_"+fecha_lst8+"_lst8.tif",lst_nodata,resxLST,resyLST)
mosaico(glob.glob(path_proc+"lst_temp/*"+fecha+".*.tif"),path_proc+"mosaicos/mosaico_"+fecha+"_lst.tif",lst_nodata,resxLST,resyLST)
mosaico(glob.glob(path_proc+"ndvi_temp/*"+fecha+".*.tif"),path_proc+"mosaicos/mosaico_"+fecha+"_ndvi.tif",ndvi_nodata,resx,resy)
mosaico(glob.glob(path_proc+"red_temp/*"+fecha+".*.tif"),path_proc+"mosaicos/mosaico_"+fecha+"_red.tif",red_nodata,resx,resy)
mosaico(glob.glob(path_proc+"red8_temp/*"+fecha+".*.tif"),path_proc+"mosaicos/mosaico_"+fecha+"_red8.tif",red_nodata,resx,resy)
mosaico(glob.glob(path_proc+"mir_temp/*"+fecha+".*.tif"),path_proc+"mosaicos/mosaico_"+fecha+"_mir.tif",red_nodata,resx,resy)
mosaico(glob.glob(path_proc+"mir8_temp/*"+fecha+".*.tif"),path_proc+"mosaicos/mosaico_"+fecha+"_mir8.tif",red_nodata,resx,resy)

#Interpolación de líneas faltantes
lst8_fill = path_proc+"mosaicos/mosaico_"+fecha_lst8+"_lst8_fill.tif"
lst_fill = path_proc+"mosaicos/mosaico_"+fecha+"_lst_fill.tif"
ndvi_fill = path_proc+"mosaicos/mosaico_"+fecha+"_ndvi_fill.tif"
red_fill = path_proc+"mosaicos/mosaico_"+fecha+"_red_fill.tif"
mir_fill = path_proc+"mosaicos/mosaico_"+fecha+"_mir_fill.tif"
red8_fill = path_proc+"mosaicos/mosaico_"+fecha+"_red8_fill.tif"
mir8_fill = path_proc+"mosaicos/mosaico_"+fecha+"_mir8_fill.tif"

filldata(path_proc+"mosaicos/mosaico_"+fecha_lst8+"_lst8.tif", lst8_fill)
filldata(path_proc+"mosaicos/mosaico_"+fecha+"_lst.tif", lst_fill)
filldata(path_proc+"mosaicos/mosaico_"+fecha+"_ndvi.tif", ndvi_fill)
filldata(path_proc+"mosaicos/mosaico_"+fecha+"_red.tif",red_fill)
filldata(path_proc+"mosaicos/mosaico_"+fecha+"_mir.tif",mir_fill)
filldata(path_proc+"mosaicos/mosaico_"+fecha+"_red8.tif",red_fill)
filldata(path_proc+"mosaicos/mosaico_"+fecha+"_mir8.tif",mir_fill)

#Recorte de zona de interés
lst8_rec500 = path_proc+"recortes/lst8_"+fecha_lst8+"_500m.tif"
lst_rec500 = path_proc+"recortes/lst_"+fecha+"_500m.tif"
ndvi_rec500 = path_proc+"recortes/ndvi_"+fecha+"_500m.tif"
red_rec500 = path_proc+"recortes/red_"+fecha+"_500m.tif"
mir_rec500 = path_proc+"recortes/mir_"+fecha+"_500m.tif"
red8_rec500 = path_proc+"recortes/red8_"+fecha+"_500m.tif"
mir8_rec500 = path_proc+"recortes/mir8_"+fecha+"_500m.tif"

recorte(shp, resx, resy, lst8_fill, lst8_rec500, lst_nodata)
recorte(shp, resx, resy, lst_fill, lst_rec500, lst_nodata)
recorte(shp, resx, resy, ndvi_fill, ndvi_rec500, ndvi_nodata)
recorte(shp, resx, resy,red_fill,red_rec500,red_nodata)
recorte(shp, resx, resy,mir_fill,mir_rec500,mir_nodata)
recorte(shp, resx, resy,red8_fill,red8_rec500,red_nodata)
recorte(shp, resx, resy,mir8_fill,mir8_rec500,mir_nodata)

#Escalar datos
lst8_escal = path_indices+"lst8/lst8_Celsius_"+fecha_lst8+"_500m.tif"
lst_escal = path_indices+"lst/lst_Celsius_"+fecha+"_500m.tif"
ndvi_escal = path_indices+"ndvi/ndvi_"+fecha+"_500m.tif"

os.system('gdal_calc.py -A %s --outfile=%s --calc="A*%.2f-273.15" --NoDataValue %d --type="Float32" --overwrite'%(lst8_rec500, lst8_escal, factor_lst, lst_nodata))
os.system('gdal_calc.py -A %s --outfile=%s --calc="A*%.2f-273.15" --NoDataValue %d --type="Float32" --overwrite'%(lst_rec500, lst_escal, factor_lst, lst_nodata))
os.system('gdal_calc.py -A %s --outfile=%s --calc="A*%.4f" --NoDataValue %d --type="Float32" --overwrite'%(ndvi_rec500, ndvi_escal, factor_ndvi, ndvi_nodata))

#Rellenar datos
lst_fill = lst_escal[:-4] + "_fill.tif"
fillgaps(lst_escal, lst8_escal, lst_fill)

#Corregir temperatura
lst_corr = path_indices+"lst/lst_Celsius_"+fecha+"_500m_corr.tif"
os.system('gdal_calc.py -A %s -B %s --outfile=%s --calc="A+B*0.006" --NoDataValue %d --type="Float32" --overwrite'%(lst_fill, dem, lst_corr, lst_nodata))#La corrección de lst 6grados/km vertical

#TVDI
salida_tvdi = path_indices + 'tvdi/tvdi_%s_500m.tif'%fecha
res = fn_tvdi(ndvi_escal, lst_corr, salida_tvdi)

tvdi_mask = salida_tvdi[:-4] + '_mask.tif'
os.system('gdal_calc.py -A %s -B %s --outfile=%s --calc="A*B" --NoDataValue=-9999 --type="Float32" --overwrite'%(salida_tvdi, uw_mask, tvdi_mask))


#ANOMALIA
print("Calculando la anomalía")
path_stats_mean = path_estadisticos + "tvdi/mean/"
path_stats_sd = path_estadisticos + "tvdi/std/"
path_salida = path_anomalias + "tvdi/anomalias_m/"

media = path_stats_mean + "tvdi_%s_fill_mean_2000_2019_500m.tif"%fecha[-3:]
sd = path_stats_sd + "tvdi_%s_fill_std_2000_2019_500m.tif"%fecha[-3:]
tvdi_anom = path_salida + "anom_%s_500fill_m.tif"%fecha

anomalia(tvdi_mask, media, sd, tvdi_anom)

#ENMASCARAMIENTO Y REMUESTREO A 10 KM USANDO LA MEDIA
##TVDI
print("Aplicando filtro de media a salidas")
smooth_tvdi_path = path_indices + "tvdi/tvdi_fill_m_10km/" #Directorio resultado de procesamiento
out_tvdi_10km_filt = filtro_media(tvdi_mask, smooth_tvdi_path, uw_mask, uw_mask2)
##ANOMALIA
smooth_anom_path = path_anomalias + "tvdi/anomalias_m_10km/" #Directorio resultado de procesamiento
out_anom_10km_filt = filtro_media(tvdi_anom, smooth_anom_path, uw_mask, uw_mask2)


#Subir datos al GeoServer
subirtiff(out_tvdi_10km_filt, estilo_tvdi, ws_tvdi, user, passw, ipgeo, port)
subirtiff(out_anom_10km_filt, estilo_anom, ws_anom, user, passw, ipgeo, port)
subirtiff(lst_escal, estilo_lst, ws_lst, user, passw, ipgeo, port)
subirtiff(ndvi_escal, estilo_ndvi, ws_ndvi, user, passw, ipgeo, port)


#Borrar archivos temporales
print("Borrando archivos temporales")
os.system("rm "+path_proc+"lst8_temp/*%s.*"%fecha_lst8)
os.system("rm "+path_proc+"ndvi_temp/*%s.*"%fecha)
os.system("rm "+path_proc+"lst_temp/*%s.*"%fecha)
os.system("rm "+path_proc+"mosaicos/*")
os.system("rm "+path_proc+"recortes/*")
os.remove(lst_corr)
os.remove(salida_tvdi)
os.remove(tvdi_mask)
