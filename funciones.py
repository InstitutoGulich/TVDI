#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import glob
from osgeo import gdal
#import gdal
import numpy as np
from osgeo.gdalconst import *
from scipy import stats
from datetime import datetime, timedelta
from astropy.convolution import convolve, Box2DKernel

from multiprocessing import Pool


#Funciones de descarga y procesamiento de datos para el proyecto SANCOR

def extraccion3(paths, prods, dem, mask, fecha, sufix):
	path = paths[0]
	path_descarga = paths[1]
	path_proc = paths[2]
	path_indices = paths[3]

	for product in prods:
		prod = product["name"]
		bands = product["bands"]

		for band in bands:
			band_name = band["name"]
			band_dir = band["dir"]
			band_nodata = band["nodata"]
			band_factor = band["factor"]
			band_res = band["res"]

			# Llevar de hdf a tif
			lista_hdf = glob.glob("%s/%s/%s/%s/*%s*.hdf"%(path_descarga,prod,fecha[0:4],fecha[4:7],fecha))
	
			for vnp in lista_hdf:
				hdfToTiff(vnp, path_proc + band_dir + "_temp/", band_dir, band_nodata)

			# Realizar el Mosaico de las imágenes
			mosaico(glob.glob(path_proc + band_dir + "_temp*/" + fecha+".*.tif"),path_proc + "mosaicos/mosaico_" + fecha + "_" + band_dir + ".tif", band_nodata, band_res, band_res)
	
			# Interpolación de líneas faltantes
			band_fill = path_proc +" mosaicos/mosaico_" + fecha + "_" + band_dir + "_fill.tif"
			filldata(path_proc + "mosaicos/mosaico_" + fecha + "_" + band_dir + ".tif", band_fill)
			
			#Recorte de zona de interés
			recortexTIFF_Float32(band_fill, mask, band_nodata, band_fill)
	
			# Escalar datos
			band_escal = path_indices + band_dir + "/" + band_dir + "_" + fecha + sufix, + ".tif"
			os.system('gdal_calc.py -A %s --outfile=%s ---calc="A*%.4f" --NoDataValue %d --type="Float32" --overwrite'%(band_fill, band_escal, band_factor, band_nodata))

			# Llevar al tamaño adecuado
			temp1 = band_escal[:-4] + "_temp1.tif"
	
			# Para SA: -99.9999999747252275 -60.0015558891651537 -29.9986385732798908 19.9999999982009413
			# Para ArUyPy: -73.5777800000000042 -55.0601126068132203 -53.0941669812803099 -19.2913600000000010
			os.system('gdalwarp -te -99.9999999747252275 -60.0015558891651537 -29.9986385732798908 19.9999999982009413 -tr %s %s %s %s '%(band_res,band_res,band_escal,temp1))
			os.system('gdal_translate -outsize 13762 15728 %s %s'%(temp1,band_escal))

	#Borrar archivos temporales
	print("Borrando archivos temporales")
	os.system("rm " + path_indices + "lst/*temp*.tif")


def extraccion2(paths, prods, dirs, dataset, nodata, factor, res, dem, mask, fecha, sufix):
	path = paths[0]
	path_descarga = paths[1]
	path_proc = paths[2]
	path_indices = paths[3]

	for i_p, prod in enumerate(prods):
		bands = dataset[i_p]
		dir_dwn_bands = dirs[i_p]
		nodata_bands = nodata[i_p]
		factor_bands = factor[i_p]
		resol_bands = res[i_p]
		
		for i_b, band in enumerate(bands):
			band_dir = dir_dwn_band[i_b]
			band_nodata = nodata_bands[i_b]
			band_factor = factor_bands[i_b]
			band_resol = resol_bands[i_b]
			
			# Llevar de hdf a tif
			lista_hdf = glob.glob("%s/%s/%s/%s/*%s*.hdf"%(path_descarga,prod,fecha[0:4],fecha[4:7],fecha))
	
			for vnp in lista_hdf:
				hdfToTiff(vnp, path_proc + band_dir + "_temp/", band_dir, band_nodata)

			# Realizar el Mosaico de las imágenes
			mosaico(glob.glob(path_proc + band_dir + "_temp*/" + fecha+".*.tif"),path_proc + "mosaicos/mosaico_" + fecha + "_" + band_dir + ".tif", band_nodata, band_resol, band_resol)
	
			# Interpolación de líneas faltantes
			band_fill = path_proc +" mosaicos/mosaico_" + fecha + "_" + band_dir + "_fill.tif"
			filldata(path_proc + "mosaicos/mosaico_" + fecha + "_" + band_dir + ".tif", band_fill)
			
			#Recorte de zona de interés
			recortexTIFF_Float32(band_fill, mask, band_nodata, band_fill)
	
			# Escalar datos
			band_escal = path_indices + band_dir + "/" + band_dir + "_" + fecha + sufix, + ".tif"
			os.system('gdal_calc.py -A %s --outfile=%s ---calc="A*%.4f" --NoDataValue %d --type="Float32" --overwrite'%(band_fill, band_escal, band_factor, band_nodata))

			# Llevar al tamaño adecuado
			temp1 = band_escal[:-4] + "_temp1.tif"
	
			# Para SA: -99.9999999747252275 -60.0015558891651537 -29.9986385732798908 19.9999999982009413
			# Para ArUyPy: -73.5777800000000042 -55.0601126068132203 -53.0941669812803099 -19.2913600000000010
			os.system('gdalwarp -te -99.9999999747252275 -60.0015558891651537 -29.9986385732798908 19.9999999982009413 -tr %s %s %s %s '%(band_resol,band_resol,band_escal,temp1))
			os.system('gdal_translate -outsize 13762 15728 %s %s'%(temp1,band_escal))

	#Borrar archivos temporales
	print("Borrando archivos temporales")
	os.system("rm " + path_indices + "lst/*temp*.tif")
	

def extraccion(paths, prods, fecha, dataset, nodata, factor, dem, mask):
	# dataset = [LST, Red, NIR, SWIR1]
	
	path = paths[0]
	path_descarga = paths[1]
	path_proc = paths[2]
	path_indices = paths[3]
	
	resx,resy = 0.005086568914507000154,0.005086568914507000154 #Resolución MOD09A1
	resxLST,resyLST = 0.008626647658857677925,0.008626647658857677925 #Resolución MOD11A2

	lst_nodata = nodata[0]
	red_nodata = nodata[1]
	nir_nodata = nodata[2]
	swir_nodata = nodata[3]
	
	lst_factor = factor[0]
	red_factor = factor[1]
	nir_factor = factor[2]
	swir_factor = factor[3]

	# Llevar de hdf a tif
	lista_lst = glob.glob("%s/%s/%s/%s/*%s*.hdf"%(path_descarga,prods[0],fecha[0:4],fecha[4:7],fecha))
	lista_veg = glob.glob("%s/%s/%s/%s/*%s*.hdf"%(path_descarga,prods[1],fecha[0:4],fecha[4:7],fecha))

	for vnp in lista_lst:
		resLST = hdfToTiff(vnp, path_proc+'lst_temp/', dataset[0], lst_nodata)

	for vnp in lista_veg:
		res = hdfToTiff(vnp, path_proc+'red_temp/', dataset[1], red_nodata)
		res = hdfToTiff(vnp, path_proc+'nir_temp/', dataset[2], nir_nodata)
		res = hdfToTiff(vnp, path_proc+'swir_temp/', dataset[3], swir_nodata)

	# Realizar el Mosaico de las imágenes
	mosaico(glob.glob(path_proc+"lst_temp/*"+fecha+".*.tif"),path_proc+"mosaicos/mosaico_"+fecha+"_lst.tif",lst_nodata,resxLST,resyLST)
	mosaico(glob.glob(path_proc+"red_temp/*"+fecha+".*.tif"),path_proc+"mosaicos/mosaico_"+fecha+"_red.tif",red_nodata,resx,resy)
	mosaico(glob.glob(path_proc+"nir_temp/*"+fecha+".*.tif"),path_proc+"mosaicos/mosaico_"+fecha+"_nir.tif",nir_nodata,resx,resy)
	mosaico(glob.glob(path_proc+"swir_temp/*"+fecha+".*.tif"),path_proc+"mosaicos/mosaico_"+fecha+"_swir.tif",swir_nodata,resx,resy)

	# Interpolación de líneas faltantes
	lst_fill = path_proc+"mosaicos/mosaico_"+fecha+"_lst_fill.tif"
	red_fill = path_proc+"mosaicos/mosaico_"+fecha+"_red_fill.tif"
	nir_fill = path_proc+"mosaicos/mosaico_"+fecha+"_nir_fill.tif"
	swir_fill = path_proc+"mosaicos/mosaico_"+fecha+"_swir_fill.tif"

	filldata(path_proc+"mosaicos/mosaico_"+fecha+"_lst.tif", lst_fill)
	filldata(path_proc+"mosaicos/mosaico_"+fecha+"_red.tif", red_fill)
	filldata(path_proc+"mosaicos/mosaico_"+fecha+"_nir.tif", nir_fill)
	filldata(path_proc+"mosaicos/mosaico_"+fecha+"_swir.tif", swir_fill)

	#Recorte de zona de interés
#	lst8_rec500 = path_proc+"recortes/lst8_"+fecha_lst8+"_500m.tif"
#	lst_rec500 = path_proc+"recortes/lst_"+fecha+"_500m.tif"
#	ndvi_rec500 = path_proc+"recortes/ndvi_"+fecha+"_500m.tif"

#	recorte(shp, resx, resy, lst8_fill, lst8_rec500, lst_nodata)
#	recorte(shp, resx, resy, lst_fill, lst_rec500, lst_nodata)
#	recorte(shp, resx, resy, ndvi_fill, ndvi_rec500, ndvi_nodata)
	recortexTIFF_Float32(lst_fill, mask, lst_nodata, lst_fill)
	recortexTIFF_Float32(red_fill, mask, red_nodata, red_fill)
	recortexTIFF_Float32(nir_fill, mask, nir_nodata, nir_fill)
	recortexTIFF_Float32(swir_fill, mask, swir_nodata, swir_fill)
	
	# Escalar datos
	lst_escal = path_indices + "lst/lst_Celsius_"+fecha+"_500m.tif"
	red_escal = path_indices + "red/red_"+fecha+"_500m.tif"
	nir_escal = path_indices + "nir/nir_"+fecha+"_500m.tif"
	swir_escal = path_indices + "swir/swir_"+fecha+"_500m.tif"

	os.system('gdal_calc.py -A %s --outfile=%s --calc="A*%.2f-273.15" --NoDataValue %d --type="Float32" --overwrite'%(lst_fill, lst_escal, lst_factor, lst_nodata))
	os.system('gdal_calc.py -A %s --outfile=%s --calc="A*%.4f" --NoDataValue %d --type="Float32" --overwrite'%(red_fill, red_escal, red_factor, red_nodata))
	os.system('gdal_calc.py -A %s --outfile=%s --calc="A*%.4f" --NoDataValue %d --type="Float32" --overwrite'%(nir_fill, nir_escal, nir_factor, nir_nodata))
	os.system('gdal_calc.py -A %s --outfile=%s --calc="A*%.4f" --NoDataValue %d --type="Float32" --overwrite'%(swir_fill, swir_escal, swir_factor, swir_nodata))

	# Llevar al tamaño adecuado
	temp1 = lst_escal[:-4] + "_temp1.tif"
	temp2 = lst_escal[:-4] + "_temp2.tif"
	temp3 = lst_escal[:-4] + "_temp3.tif"
	temp4 = lst_escal[:-4] + "_temp4.tif"
	
	# Para SA: -99.9999999747252275 -60.0015558891651537 -29.9986385732798908 19.9999999982009413
	# Para ArUyPy: -73.5777800000000042 -55.0601126068132203 -53.0941669812803099 -19.2913600000000010
	os.system('gdalwarp -te -99.9999999747252275 -60.0015558891651537 -29.9986385732798908 19.9999999982009413 -tr %s %s %s %s '%(resx,resy,lst_escal,temp1))
	os.system('gdal_translate -outsize 13762 15728 %s %s'%(temp1,lst_escal))
	os.system('gdalwarp -te --99.9999999747252275 -60.0015558891651537 -29.9986385732798908 19.9999999982009413 -tr %s %s %s %s '%(resx,resy,red_escal,temp2))
	os.system('gdal_translate -outsize 13762 15728 %s %s'%(temp2,red_escal))
	os.system('gdalwarp -te -99.9999999747252275 -60.0015558891651537 -29.9986385732798908 19.9999999982009413 -tr %s %s %s %s '%(resx,resy,nir_escal,temp3))
	os.system('gdal_translate -outsize 13762 15728 %s %s'%(temp3,nir_escal))
	os.system('gdalwarp -te -99.9999999747252275 -60.0015558891651537 -29.9986385732798908 19.9999999982009413 -tr %s %s %s %s '%(resx,resy,swir_escal,temp4))
	os.system('gdal_translate -outsize 13762 15728 %s %s'%(temp4,swir_escal))

	# Corregir temperatura por altura
	lst_corr = path_indices+"lst/lst_Celsius_"+fecha+"_500m_corr.tif"
	os.system('gdal_calc.py -A %s -B %s --outfile=%s --calc="A+B*0.006" --NoDataValue %d --type="Float32" --overwrite'%(lst_escal, dem, lst_corr, lst_nodata))#La corrección de lst 6grados/km vertical

	#Borrar archivos temporales
	print("Borrando archivos temporales")
	os.system("rm " + path_indices + "lst/*temp*.tif")
	
#def descarga(url,appkey,tiles,indices,fecha,path_descarga,path_xml):
#	for ind in indices:
#		for t in tiles:
#			os.system('wget -e robots=off -m -np -R .html,.tmp -nH --cut-dirs=3 -A "%s.A%s.%s.061.*.hdf" "%s%s/%s/%s/" --header "Authorization: Bearer %s" -P %s'%(ind,fecha,t,url,ind,fecha[0:4],fecha[4:7],appkey,path_descarga))
#			if len(glob.glob("%s.A%s.%s.061.*.hdf"%(path_descarga,fecha,ind)))>0:
#				arch = glob.glob("%s/%s/%s/%s/*%s*.hdf"%(path_descarga,ind,fecha[0:4],fecha[4:7],t))[0]
#				os.system("cp %s*%s* %s.xml"%(path_xml,t,arch))


def descarga(url,appkey,tiles,prods,fecha,path_descarga,path_xml):
	for prod in prods:
		for t in tiles:
			os.system('wget -e robots=off -m -np -R .html,.tmp -nH --cut-dirs=3 -A "%s.A%s.%s.061.*.hdf" "%s%s/%s/%s/" --header "Authorization: Bearer %s" -P %s'%(prod,fecha,t,url,prod,fecha[0:4],fecha[4:7],appkey,path_descarga))
		list_files = glob.glob("%s/%s/%s/%s/*.hdf"%(path_descarga,prod,fecha[0:4],fecha[4:7]))
		if len(list_files)<len(tiles):
			print("list_files: ", list_files)
			print("tiles ", tiles)
			#descarga(url,appkey,tiles,prods,fecha,path_descarga,path_xml)
		else:
			arch = glob.glob("%s/%s/%s/%s/*%s*.hdf"%(path_descarga,prod,fecha[0:4],fecha[4:7],t))[0]
			os.system("cp %s*%s* %s.xml"%(path_xml,t,arch))


def descarga3(url,appkey,tiles,prods,fecha,path_descarga,path_xml):
	for prod["name"] in prods:
		for t in tiles:
			os.system('wget -e robots=off -m -np -R .html,.tmp -nH --cut-dirs=3 -A "%s.A%s.%s.061.*.hdf" "%s%s/%s/%s/" --header "Authorization: Bearer %s" -P %s'%(prod,fecha,t,url,prod,fecha[0:4],fecha[4:7],appkey,path_descarga))
		list_files = glob.glob("%s/%s/%s/%s/*.hdf"%(path_descarga,prod,fecha[0:4],fecha[4:7]))
		if len(list_files)<len(tiles):
			print("list_files: ", list_files)
			print("tiles ", tiles)
			#descarga(url,appkey,tiles,prods,fecha,path_descarga,path_xml)
		else:
			arch = glob.glob("%s/%s/%s/%s/*%s*.hdf"%(path_descarga,prod,fecha[0:4],fecha[4:7],t))[0]
			os.system("cp %s*%s* %s.xml"%(path_xml,t,arch))


#######################################################################################################################
########################################Función para exportar los hdf a geotiff########################################
#######################################################################################################################
def hdfToTiff(inputHDF, outputPath, product, nodata):
	ds = gdal.Open(inputHDF) #Abre la imagen
	datasets = ds.GetSubDatasets() #Obtiene los datasets
#	print(datasets)
	for d in datasets: #Se queda con el dataset de interés
		if product in d[0]:
			prod = d[0]

	ds2=gdal.Open(prod) #Abre el dataset
	sds = ds2.GetRasterBand(1).ReadAsArray() #Obtiene la matriz de datos
	geoTs = ds2.GetGeoTransform() #Copia los parámetros de la imagen
	prj = ds2.GetProjection() #Copia la proyección
	dataType = ds2.GetRasterBand(1).DataType #Obtiene el tipo de dato original
	outputTif = outputPath+os.path.basename(inputHDF).replace(".hdf","_sin.tif") #Nombre de imagen de salida

	if not os.path.exists(outputPath): os.makedirs(outputPath) #Crea la ruta si no existe

	# Output raster array to GeoTIFF file            
	driver = gdal.GetDriverByName('GTiff')
	export = driver.Create(outputTif, sds.shape[1], sds.shape[0], 1, dataType)
	band = export.GetRasterBand(1)
	band.WriteArray(sds)
	export.SetGeoTransform(geoTs)
	export.SetProjection(prj) # Set output coordinate referense system information
	band.FlushCache() #Flush band from memory
	export.FlushCache() #Flush memory
	
	del export
	
	os.system("gdalwarp -srcnodata '%s' -dstnodata '%s' '%s' '%s' -t_srs '%s'"%(nodata,nodata, outputTif,outputPath+os.path.basename(inputHDF).replace(".hdf",".tif"),'EPSG:4326'))

	os.system("rm "+outputPath+os.path.basename(inputHDF).replace(".hdf","_sin.tif"))
	return(geoTs[1],geoTs[5])


def mosaico(lista, salida, nodata, resx, resy):
	os.system("gdal_merge.py -n %d -a_nodata %d -ps %.15f %.15f -of GTiff -o '%s' %s"%(nodata,nodata,resx,resy,salida,str(lista).replace(","," ")[1:-1]))


def recorte(shp,resx,resy,entrada,salida,nodata):
	os.system("gdalwarp -q -srcnodata %d -dstnodata %d -cutline %s -crop_to_cutline -tr %.15f %.15f -of GTiff '%s' '%s'"%(nodata, nodata,shp,resx,resy,entrada,salida))


def recortexTIFF(tiff,mask,nodata,salida,type_out):
	if type_out == 'int':
		nodata = getNodata(mask).astype('float')
	os.system("gdal_calc.py -A %s -B %s --NoDataValue=%s --outfile=%s --type=%s --calc='A*B' --co='COMPRESS=DEFLATE' --overwrite"%(tiff,mask,nodata,salida,type_out))

	
def filldata(entrada,salida):
	ds = gdal.Open(entrada)
	ds1 = ds.GetRasterBand(1)
	sds = ds1.ReadAsArray().astype('float')
	nodata = ds1.GetNoDataValue()
	
	for i in range(1,sds.shape[0]-1):
		if (np.nansum(sds[i,:])==0 and np.nansum(sds[i-1,:])!=0 and np.nansum(sds[i+1,:])!=0):
			sds[i,:]=(sds[i-1,:]+sds[i+1,:])/2
		elif (np.nansum(sds[i,:])==nodata*sds.shape[1] and np.nansum(sds[i-1,:])!=nodata*sds.shape[1] and np.nansum(sds[i+1,:])!=nodata*sds.shape[1]):
			sds[i,:]=(sds[i-1,:]+sds[i+1,:])/2
		else:
			pass
	
	geoTs = ds.GetGeoTransform() #Parámetros de la imagen (coordenadas origen y dimensiones)
	driver = gdal.GetDriverByName("GTiff") #Tipo de imagen (geotiff)
	prj = ds.GetProjection() #Sistema de referencia de la imagen (aquí se lee el mismo sistema de referencia de la imagen de entrada)

	#Crear el espacio
	export=driver.Create(salida,sds.shape[1],sds.shape[0],1,GDT_Int32)
	banda=export.GetRasterBand(1) #Cargar la banda creada en el paso anterior
	banda.WriteArray(sds) #Escribir los valores de la matriz calculada
	banda.SetNoDataValue(nodata) #Asignar los parametros de la transformacion a la salida
	export.SetGeoTransform(geoTs) #Asignar los parametros de la transformacion a la salida
	export.SetProjection(prj) #definir la proyección
	banda.FlushCache()#descargar de la memoria virtual al disco
	export.FlushCache()#descargar de la memoria virtual al disco


def subirtiff(img, estilo, workspace, user, passw, ipgeo, port):
	namelayer = os.path.basename(img) #Nombre de la capa a subir con extensión
	coveragestore = namelayer[:-4] #Nombre de la capa a subir sin extensión

	#Línea para crear el almacén de datos
	os.system('curl -u '+user+':'+passw+' -XPOST -H "Content-type: application/xml" -d "<coverageStore><name>"'+coveragestore+'"</name><workspace>"'+workspace+'"</workspace><enabled>true</enabled><type>GeoTIFF</type></coverageStore>" http://'+ipgeo+':'+port+'/geoserver/rest/workspaces/'+workspace+'/coveragestores')
	print("coverage creado")

	#Línea para cargar la capa
	os.system('curl -u '+user+':'+passw+' -XPUT -H "Content-type:image/tiff" --data-binary @'+img+' http://'+ipgeo+':'+port+'/geoserver/rest/workspaces/'+workspace+'/coveragestores/'+coveragestore+'/file.geotiff')
	print("tiff subido")

	#Línea para asignar el estilo a la capa
	os.system('curl -u '+user+':'+passw+' -XPUT -H "Content-type:text/xml" -d "<layer><defaultStyle><name>'+estilo+'</name><workspace>tvdi</workspace></defaultStyle></layer>" http://'+ipgeo+':'+port+'/geoserver/rest/layers/'+workspace+':'+namelayer)
	print("estilo asignado")


def subirtiff2(img, estilo, workspace, user, passw, ipgeo, port):
	basename = os.path.basename(img) #Nombre de la capa a subir con extensión
	namelayer = "_".join(basename.split("_")[:2])
	coveragestore = namelayer #Nombre de la capa a subir sin extensión

	#Línea para crear el almacén de datos
	os.system('curl -u '+user+':'+passw+' -XPOST -H "Content-type: application/xml" -d "<coverageStore><name>"'+coveragestore+'"</name><workspace>"'+workspace+'"</workspace><enabled>true</enabled><type>GeoTIFF</type></coverageStore>" http://'+ipgeo+':'+port+'/geoserver/rest/workspaces/'+workspace+'/coveragestores')
	print("coverage creado")

	#Línea para cargar la capa
	os.system('curl -u '+user+':'+passw+' -XPUT -H "Content-type:image/tiff" --data-binary @'+img+' http://'+ipgeo+':'+port+'/geoserver/rest/workspaces/'+workspace+'/coveragestores/'+coveragestore+'/file.geotiff')
	print("tiff subido")

	#Línea para asignar el estilo a la capa
	os.system('curl -u '+user+':'+passw+' -XPUT -H "Content-type:text/xml" -d "<layer><defaultStyle><name>'+estilo+'</name><workspace>'+workspace+'</workspace></defaultStyle></layer>" http://'+ipgeo+':'+port+'/geoserver/rest/layers/'+workspace+':'+namelayer)
	print('curl -u '+user+':'+passw+' -XPUT -H "Content-type:text/xml" -d "<layer><defaultStyle><name>'+estilo+'</name><workspace>'+workspace+'</workspace></defaultStyle></layer>" http://'+ipgeo+':'+port+'/geoserver/rest/layers/'+workspace+':'+namelayer)
	print("estilo asignado")
	
	
def creartif(file_in,matriz,tipoDato,nodata,salida):
    in_ds = gdal.Open(file_in)
    in_geoTs = in_ds.GetGeoTransform() #Parámetros de la imagen (coordenadas origen y dimensiones)
    in_prj = in_ds.GetProjection() #Sistema de referencia de la imagen (aquí se lee el mismo sistema de referencia de la imagen de entrada)
#    print("creando imagen")

    driver = gdal.GetDriverByName("GTiff") #Tipo de imagen (geotiff)

    #Crear el espacio
    export = driver.Create(salida,matriz.shape[1],matriz.shape[0],1,tipoDato)
    banda = export.GetRasterBand(1) #Cargo la banda creada en el paso anterior
    banda.WriteArray(matriz) #Escribe los valores de NDVI en la imagen
    banda.SetNoDataValue(nodata)

    export.SetGeoTransform(in_geoTs) #Asigna los parametros de la transformacion a la salida
    export.SetProjection(in_prj) #define la proyección

    banda.FlushCache()#descargar de la memoria virtual al disco
    export.FlushCache()#descargar de la memoria virtual al disco

    del in_ds


def filtro_media(
				img,
				path_out,
				mask1,
				mask2,
				scale=0.005086568914507,
				block_xsize=11,
				block_ysize=11			
				):

	start_time = datetime.now()
	print("Aplicando filtro de medias")
	
	in_ds = gdal.Open(mask2)
	in_band = in_ds.GetRasterBand(1)
	mask_data = in_band.ReadAsArray()
	del in_ds

	out_fn = path_out + '_temp.tif'
	
	in_ds = gdal.Open(img)
	in_band = in_ds.GetRasterBand(1)
	in_data = in_band.ReadAsArray()
	nodata = in_data[0,0]
	in_data = np.where(in_data == nodata, np.nan, in_data)
	del in_ds

	smooth_data = convolve(in_data, kernel=Box2DKernel(block_xsize), mask=mask_data, nan_treatment='interpolate', preserve_nan=True)
	creartif(img, smooth_data, gdal.GDT_Float32, np.nan, out_fn)

	basename = os.path.basename(img)
	basename = "_".join(basename.split("_")[:2])+"_filter.tif"   #"_10km.tif"
	out_img = path_out + basename

	# ~ os.system('gdal_translate -tr %s %s -projwin -73.5777816800000011 -19.2913703900000009 -53.0941686612803068 -55.0601229968132202 %s %s'%(scale,scale,out_fn,path_out+'_temp2.tif'))
	# ~ os.system('gdal_translate -outsize 4027 7032 %s %s'%(path_out+'_temp2.tif',path_out+'_temp3.tif'))
	os.system('gdal_calc.py -A %s -B %s --outfile=%s --calc="A*B" --type="Int16" --NoDataValue %d --co="COMPRESS=DEFLATE" --overwrite'%(out_fn, mask1, out_img, nodata))

	end_time = datetime.now()
	print('\nTotal run time:', end_time - start_time, '\n')
		
	os.system("rm "+path_out+"_temp.tif")

	return out_img


def getIntMaskImg(img,img_sal):
	ds1 = gdal.Open(img)
	sds1 = ds1.GetRasterBand(1).ReadAsArray()
	nan_value = sds1[0,0]
	mask_nan = np.where(sds1==nan_value, -99, sds1)
	del ds1
	mask_nan = mask_nan*1000
	creartif(img, mask_nan, gdal.GDT_Int16, -99000, img_sal)


def getNodata(img):
	ds1 = gdal.Open(img)
	sds1 = ds1.GetRasterBand(1).ReadAsArray()
	nan_value = sds1[0,0]
	return nan_value
	

def calc_indice(indice, img1, img2, salida, mask):
	if indice == "tvdi":
		# LST ALtitude correction

		#cargar imágenes
		ds_ndvi = gdal.Open(img1)
		sds_ndvi = ds_ndvi.GetRasterBand(1)
		matriz_ndvi = sds_ndvi.ReadAsArray()
		matriz_ndvi[matriz_ndvi==sds_ndvi.GetNoDataValue()]=np.nan
		matriz_ndvi[matriz_ndvi>100]=np.nan
		del ds_ndvi
		del sds_ndvi

		ds_lst = gdal.Open(img2)
		sds_lst = ds_lst.GetRasterBand(1)
		matriz_lst = sds_lst.ReadAsArray()
		matriz_lst[matriz_lst==sds_lst.GetNoDataValue()]=np.nan
		matriz_lst[matriz_lst<=0]=np.nan

		#Salida de la imagen georreferenciada
		geoTs = ds_lst.GetGeoTransform() #Parámetros de la imagen (coordenadas origen y dimensiones)
		driver = gdal.GetDriverByName("GTiff") #Tipo de imagen (geotiff)
		prj = ds_lst.GetProjection() #Sistema de referencia de la imagen (aquí se lee el mismo sistema de referencia de la imagen de entrada)

		del ds_lst
		del sds_lst

		################################
		filas = matriz_lst.shape[0]
		cols = matriz_lst.shape[1]
		matriz_tvdi = np.zeros((filas,cols))
		ancho_filas = 2000

		#Parámetros
		min_ndvi = 0. #límite inferior del corte del histograma para calcular la línea de temperatura
		cant_px = 1 #Cantidad de píxeles a tomar de cada delta (10%)
		max_ndvi = np.nanmax(matriz_ndvi)
		delta = 0.02

		for i in range(0,filas-ancho_filas+1,2):
			matriz_ndvi_sub = matriz_ndvi[i:i+ancho_filas,:]
			matriz_lst_sub = matriz_lst[i:i+ancho_filas,:]

			#Eliminar valores nulos
			nan_data = ~np.isnan(matriz_lst_sub)
			nan_data2 = ~np.isnan(matriz_ndvi_sub)
			nan_data = nan_data*nan_data2
			lst_C_reshape = matriz_lst_sub[nan_data]
			ndvi_reshape = matriz_ndvi_sub[nan_data]

			lst_regr = []	#Lista para almacenar los valores de temperatura para la regresión
			ndvi_regr = []	#Lista para almacenar los valores de vegetación para la regresión
			tmin = []
			print("fila %s a %s: calculando deltas"%(str(i),str(i+ancho_filas)))

			for v in np.arange(min_ndvi,max_ndvi,delta):
				#Valores que están en el delta definido
				lst_arr = lst_C_reshape[np.where((ndvi_reshape>v) & (ndvi_reshape<=(v+delta)))]
				ndvi_arr = ndvi_reshape[np.where((ndvi_reshape>v) & (ndvi_reshape<=(v+delta)))]

				#Ordenar los valores según mayor temperatura
				indices = lst_arr.argsort()
				lst_arr = lst_arr[indices]
				ndvi_arr = ndvi_arr[indices]

				if len(ndvi_arr)>0:

					tmin.append(lst_arr[0]) #Valores bajos de la dispersión (límite húmedo)
					lst_regr.append(lst_arr[-cant_px:][0]) #Valores altos de la dispersión (límite seco)
					ndvi_regr.append(ndvi_arr[-cant_px:][0]) #Valor equivalente NDVI (límite seco)
				else:
					pass

			tmin = np.array(tmin)
			lst_regr = np.array(lst_regr)
			ndvi_regr = np.array(ndvi_regr)
			print("calculando regresión")
			#Regresión lineal
			slope, intercept, r_value, p_value, std_err = stats.linregress(ndvi_regr,lst_regr)

			print("calculando tvdi")
			#Cálculo del TVDI
			print(i)
			print(ancho_filas)
			tvdi = (matriz_lst_sub-np.mean(tmin))/(intercept+slope*matriz_ndvi_sub-np.mean(tmin))
			if i==0:
				matriz_tvdi[i:i+ancho_filas,:] = tvdi
			elif (i==filas-ancho_filas):
				matriz_tvdi[i+ancho_filas//2-1:,:] = tvdi[ancho_filas//2-1:,:]
			else:
				matriz_tvdi[i+ancho_filas//2-1,:] = tvdi[ancho_filas//2-1,:]
				matriz_tvdi[i+1+ancho_filas//2-1,:] = tvdi[1+ancho_filas//2-1,:]
		print("creando imagen")
		#Crear el espacio
		export=driver.Create(salida,matriz_tvdi.shape[1],matriz_tvdi.shape[0],1,GDT_Float32)
		banda=export.GetRasterBand(1) #Cargo la banda creada en el paso anterior
		banda.WriteArray(matriz_tvdi) #Escribe los valores de NDVI en la imagen
		export.SetGeoTransform(geoTs) #Asigna los parametros de la transformacion a la salida
		export.SetProjection(prj) #define la proyección
		banda.FlushCache()#descargar de la memoria virtual al disco
		export.FlushCache()#descargar de la memoria virtual al disco
	
#		nodata = getNodata(salida)
#		recortexTIFF(salida,mask,nodata,salida,"Float32")
		return 1
	else:
		nodata = getNodata(img2)
		os.system('gdal_calc.py -A %s -B %s --outfile=%s --calc="(A-B)*1.0/(A+B)" --NoDataValue %d --type="Float32" --co="COMPRESS=DEFLATE" --overwrite'%(img1,img2, salida, nodata))
	
	
def calc_anom(indice, path_estadisticos, path_anomalias, fecha, uw_mask, img):
	print("Calculando la anomalía de " + indice)

	path_stats_mean = path_estadisticos + indice + "/mean/"
	path_stats_sd = path_estadisticos + indice + "/std/"
	path_salida = path_anomalias + indice + "/anomalias_m/"
#	rango_anyos = str(fecha[:4]-11) + "_" + str(fecha[:4]-1)
	rango_anyos = "2000_2019"

	mean = path_stats_mean + "%s_%s_fill_mean_%s_500m.tif"%(indice, fecha[-3:], rango_anyos)
	sd = path_stats_sd + "%s_%s_fill_std_%s_500m.tif"%(indice, fecha[-3:], rango_anyos)
		
	anom_out = path_salida + "anom_%s_%s_500_m.tif"%(indice, fecha)
	nodata = getNodata(uw_mask)
		
	os.system('gdal_calc.py '\
			'-A %s '\
			'-B %s '\
			'-C %s '\
			'--outfile=%s '\
                '--calc="(A-B)/(C)" '\
                '--NoDataValue=-9999 '\
                '--type="Float32"'\
					%(img, mean, sd, anom_out))

	recortexTIFF(anom_out, uw_mask, nodata, anom_out, "Float32")


def img_Int16(img_in, img_out, escal, uw_mask):
	# RECORTAR LST Y LLEVAR A ENTEROS
	img_int = img_out[:-8] + "int.tif"
	img_nodata = getNodata(img_in).astype(int)
	os.system('gdal_calc.py -A %s --outfile=%s --calc="A*%s" --NoDataValue=%s --type="Int16" --overwrite'%(img_in, img_int, escal, img_nodata))
	recortexTIFF_Int16(img_int, uw_mask, img_nodata, img_out)
	os.system("rm " + img_int)


def lead_to_tvdi(ndvi_escal, lst_escal, fecha):
	# path inputs
	path = "/home/jrubio/Documentos/VENVS/SANCOR5/" #Directorio principal de procesamiento
	path_estadisticos = path + "estadisticos/"
	path_anomalias = path + "anomalias/"
	path_mask = path + "datos/mascara/"
	path_stats_mean = path_estadisticos + "tvdi/mean/"
	path_stats_sd = path_estadisticos + "tvdi/std/"
	path_salida = path_anomalias + "tvdi/anomalias_m/"

	uw_mask = path_mask + "mask_urbano_agua_i5.tif" #Raster máscara zona de zonas urbanas y cuerpos de agua
	uw_mask2 = path_mask + "mask_urbano_agua_5.tif" #Raster máscara zona de zonas urbanas y cuerpos de agua
	uw_mask3 = path_mask + "mask_urbano_agua_i3_250m.tif" #Raster máscara zona de zonas urbanas y cuerpos de agua
	uw_mask4 = path_mask + "mask_urbano_agua_3_250m.tif" #Raster máscara zona de zonas urbanas y cuerpos de agua
	uw_mask5 = path_mask + "mask_urbano_agua_i3.tif" #Raster máscara zona de zonas urbanas y cuerpos de agua
	uw_mask6 = path_mask + "mask_urbano_agua_3.tif" #Raster máscara zona de zonas urbanas y cuerpos de agua
	pas_mask = path_mask + "mask_zonas_pas.tif"
	crop_mask = path_mask + "mask_cultivos_i3_250m.tif"
	dem = path + "datos/dem500m/dem_500m.tif" #Modelo de Elevación Digital para corrección de temperatura superficial

	#path outputs
	path_indices = path + "indices/" #Directorio donde se almacenan los índices calculados

	ndvi_nodata = -3000
	lst_nodata = 0

	print("""
--------------------
Recortando lst escalado
--------------------
	""")
	recortexTIFF(lst_escal, uw_mask, lst_nodata, lst_escal, "Float32")

	############################
	#Corregir temperatura
	lst_corr = path_indices+"lst/lst_Celsius_"+fecha+"_500m_corr.tif"
	print("""
--------------------
Corregiendo lst x altura
--------------------
	""")
	os.system('gdal_calc.py -A %s -B %s --outfile=%s --calc="A+B*0.006" --NoDataValue %d --type="Float32" --co="COMPRESS=DEFLATE" --overwrite'%(lst_escal, dem, lst_corr, lst_nodata))#La corrección de lst 6grados/km vertical

	############################
	#TVDI
	print("""
--------------------
Calculando tvdi
--------------------
	""")
	###salida: tvdi_fecha_500m_f.tif
	###estado: original
	###tipo: float
	###region: toda
	salida_tvdi = path_indices + 'tvdi/tvdi_%s_500m.tif'%fecha
	res = calc_indice("tvdi", ndvi_escal, lst_corr, salida_tvdi, uw_mask)
	salida_tvdi_f = path_indices + 'tvdi/tvdi_%s_500m_f.tif'%fecha
	os.system("cp " + salida_tvdi + " " + salida_tvdi_f)
	
	print("""
--------------------
Recortando tvdi
--------------------
	""")
	###salida: tvdi_fecha_500m_f.tif
	###estado: original recortado a máscara
	###tipo: float
	###region: toda
	recortexTIFF(salida_tvdi_f, uw_mask, -33, salida_tvdi_f, "Float32")
		
	#TVDI A ENTEROS
	print("""
--------------------
Recortando tvdi entero
--------------------
	""")
	###salida: tvdi_fecha_500m_i.tif
	###estado: original recortado a máscara
	###tipo: int
	###region: toda
	tvdi_int = salida_tvdi[:-4] + '_i.tif'
	nodata = getNodata(salida_tvdi_f).astype('int')*1000
	os.system('gdal_calc.py -A %s --outfile=%s --calc="A*1000" --NoDataValue=%s --type="Int16" --co="COMPRESS=DEFLATE" --overwrite'%(salida_tvdi_f, tvdi_int, nodata))
		
	#FILTRANDO A 5 KM USANDO LA MEDIA
	print("""
--------------------
Filtrando tvdi entero
--------------------
	""")
	###salida: tvdi_fecha_filter.tif
	###estado: filtrado recortado a máscara
	###tipo: int
	###region: toda
	smooth_tvdi_path = path_indices + "tvdi/tvdi_fill_m_5km/" #Directorio resultado de procesamiento
	out_tvdi_5km_filt = filtro_media(tvdi_int, smooth_tvdi_path, uw_mask, uw_mask2)

	print("""
--------------------
Enmascarando tvdi a zonas pas
--------------------
	""")
	###salida: tvdicrop_fecha_500m_f.tif 
	###estado: original recortado a máscara
	###tipo: float
	###region: pas
	tvdi_crop = path_indices + 'tvdi/tvdicrop_%s_500m_f.tif'%fecha
	recortexTIFF(salida_tvdi_f, pas_mask, -33, tvdi_crop, "Float32")

	print("""
--------------------
Limitando y Remuestreando tvdi en zonas pas a 250m
--------------------
	""")
	###salida: tvdicrop_fecha_250m_f.tif 
	###estado: original recortado a máscara
	###tipo: float
	###region: pas
	scale = 0.002543284457288304268
	temp1 = salida_tvdi_f[:-4] + "_temp1.tif"
	
	tvdi_crop_salida = path_indices + 'tvdi/tvdicrop_%s_250m_f.tif'%fecha
	# Para ArUyPy: -73.5777800000000042 -55.0601126068132203 -53.0941669812803099 -19.2913600000000010
	os.system('gdalwarp -te -73.5777800000000042 -55.0601126068132203 -53.0941669812803099 -19.2913600000000010 -tr %s %s %s %s '%(scale,scale,tvdi_crop,temp1))
	os.system('gdal_translate -outsize 8054 14064 %s %s'%(temp1,tvdi_crop_salida))
	os.system('rm '+ temp1)

	print("""
--------------------
Enmascarando tvdi a zonas de cultivos 250m
--------------------
	""")
	###salida: tvdicrop_fecha_250m_f.tif 
	###estado: original recortado a máscara cultivos
	###tipo: float
	###region: pas
	recortexTIFF(tvdi_crop_salida, crop_mask, -99, tvdi_crop_salida, "Float32")

	print("""
--------------------
Recortando tvdi crop entero con filtro
--------------------
	""")
	###salida: tvdicrop_fecha_5km_i.tif 
	###estado: original recortado a máscara cultivos
	###tipo: float
	###region: pas
	tvdicrop_int = tvdi_crop_salida[:-6] + '_i.tif'
	nodata = getNodata(tvdi_crop_salida).astype('int')*1000
	os.system('gdal_calc.py -A %s --outfile=%s --calc="A*1000" --NoDataValue=%s --type="Int16" --co="COMPRESS=DEFLATE" --overwrite'%(tvdi_crop_salida, tvdicrop_int, nodata))
	recortexTIFF(tvdicrop_int, crop_mask, nodata, tvdicrop_int, "Int16")

	print("""
--------------------
Filtrando tvdi crop
--------------------
	""")
	###salida: tvdicrop_fecha_filter.tif 
	###estado: original recortado a máscara cultivos
	###tipo: float
	###region: pas
	smooth_tvdi_path = path_indices + "tvdi/tvdi_fill_m_5km/" #Directorio resultado de procesamiento
	out_tvdi_5km_filt = filtro_media(tvdicrop_int, smooth_tvdi_path, uw_mask3, uw_mask4)

	print("""
--------------------
Limpiando archivos temporales
--------------------
	""")
	
	#3# OJO QUE SE NECESITAN PARA LAS ANOMALIAS
#3#	os.system("rm " + path_indices + "tvdi/tvdi_%s_500m_f.tif"%fecha)
#3#	os.system("rm " + path_indices + "tvdi/tvdi_%s_500m_i.tif"%fecha)
#3#	os.system("rm " + path_indices + "tvdi/tvdi_fill_m_5km/tvdi_%s_5km_i.tif"%fecha)
#3#	os.system("rm " + path_indices + "tvdi/tvdicrop_%s_500m_f.tif"%fecha)
#3#	os.system("rm " + path_indices + "tvdi/tvdicrop_%s_500m_i.tif"%fecha)
#3#	os.system("rm " + path_indices + "tvdi/tvdicrop_%s_250m_f.tif"%fecha)
#3#	os.system("rm " + path_indices + "tvdi/tvdicrop_%s_250m_i.tif"%fecha)








	############################
	#Anomalías TVDI
	print("""
--------------------
Calculando anomalías de tvdi
--------------------
	""")
#	tvdi = path_indices + 'tvdi/tvdi_%s_500m_f.tif'%fecha

#	media = path_stats_mean + "tvdi_%s_fill_mean_2000_2019_500m.tif"%fecha[-3:]
#	sd = path_stats_sd + "tvdi_%s_fill_std_2000_2019_500m.tif"%fecha[-3:]
	
	media = path_stats_mean + "tvdi_%s_fill_mean_2000_2019_250m.tif"%fecha[-3:]
	sd = path_stats_sd + "tvdi_%s_fill_std_2000_2019_250m.tif"%fecha[-3:]
	
#	tvdi_mask = path_indices + "tvdi/tvdi_%s_500m.tif"%fecha
#	tvdi_anom = path_salida + "anom_%s_500m.tif"%fecha
	
	tvdi_mask = path_indices + "tvdi/tvdicrop_%s_250m_f.tif"%fecha
#	tvdi_anom = path_salida + "anom_%s_250m.tif"%fecha
	
	anom_out = path_salida + "anom_%s_250m_crop.tif"%fecha
	nodata = getNodata(uw_mask)
	tvdi_anom = path_salida + "anom_%s_250fill_m.tif"%fecha

#3#	anomalia(tvdi_mask, media, sd, tvdi_anom)
#3#	recortexTIFF(tvdi_anom, uw_mask, nodata, anom_out)

#ENMASCARAMIENTO Y REMUESTREO A 5 KM USANDO LA MEDIA
##ANOMALIA
#3#	smooth_anom_path = path_anomalias + "tvdi/anomalias_m_5km/" #Directorio resultado de procesamiento
#3#	out_anom_5km_filt = filtro_media(anom_out, smooth_anom_path, uw_mask, uw_mask2, scale, 21, 21)
#3#	recortexTIFF(out_anom_5km_filt, uw_mask, nodata, out_anom_5km_filt)

	os.system("rm " + tvdi_anom)




#3#	media = path_stats_mean + "tvdi_%s_fill_mean_2000_2019_500m.tif"%fecha[-3:]
#3#	sd = path_stats_sd + "tvdi_%s_fill_std_2000_2019_500m.tif"%fecha[-3:]
	
#3#	tvdi_mask = path_indices + "tvdi/tvdi_%s_500m.tif"%fecha
#3#	tvdi_anom = path_salida + "anom_%s_500m.tif"%fecha
	
#3#	anom_out = path_salida + "anom_%s_500m.tif"%fecha
#3#	nodata = getNodata(uw_mask)
#3#	tvdi_anom = path_salida + "anom_%s_250fill_m.tif"%fecha

#3#	calc_anom("tvdi", path_estadisticos, path_anomalias, fecha, uw_mask5, tvdi)
#3#	recortexTIFF(tvdi_anom, uw_mask, nodata, anom_out, "Float32")
#3#	os.system('rm '+ tvdi_anom)

#ENMASCARAMIENTO Y REMUESTREO A 5 KM USANDO LA MEDIA
##ANOMALIA
#3#	smooth_anom_path = path_anomalias + "tvdi/anomalias_m_5km/" #Directorio resultado de procesamiento
#3#	out_anom_5km_filt = filtro_media(anom_out, smooth_anom_path, uw_mask, uw_mask2, scale, 21, 21)
#3#	recortexTIFF(out_anom_5km_filt, uw_mask, nodata, out_anom_5km_filt)


#3#	media = path_stats_mean + "tvdi_%s_fill_mean_2000_2019_250m.tif"%fecha[-3:]
#3#	sd = path_stats_sd + "tvdi_%s_fill_std_2000_2019_250m.tif"%fecha[-3:]
	
#3#	tvdi_mask = path_indices + "tvdi/tvdi_%s_250m.tif"%fecha
#3#	tvdi_anom = path_salida + "anom_%s_250m.tif"%fecha
	
#3#	anom_out = path_salida + "anom_%s_250m_crop.tif"%fecha
#3#	nodata = getNodata(uw_mask)
#3#	tvdi_anom = path_salida + "anom_%s_250fill_m.tif"%fecha

#3#	anomalia(tvdi_mask, media, sd, tvdi_anom)
#3#	recortexTIFF(tvdi_anom, uw_mask, nodata, anom_out)

#ENMASCARAMIENTO Y REMUESTREO A 5 KM USANDO LA MEDIA
##ANOMALIA
#3#	smooth_anom_path = path_anomalias + "tvdi/anomalias_m_5km/" #Directorio resultado de procesamiento
#3#	out_anom_5km_filt = filtro_media(anom_out, smooth_anom_path, uw_mask, uw_mask2, scale, 21, 21)
#3#	recortexTIFF(out_anom_5km_filt, uw_mask, nodata, out_anom_5km_filt)





