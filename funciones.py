#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import glob
import gdal
import numpy as np
from gdalconst import *
from scipy import stats
import requests

#Funciones de descarga y procesamiento de datos para el proyecto SANCOR

def descarga(url,appkey,tiles,indices,fecha,path_descarga,path_xml):
    for ind in indices:
        for t in tiles:
            os.system('wget -e robots=off -m -np -R .html,.tmp -nH --cut-dirs=3 -A "%s.A%s.%s.006.*.hdf" "%s%s/%s/%s/" --header "Authorization: Bearer %s" -P %s'%(ind,fecha,t,url,ind,fecha[0:4],fecha[4:7],appkey,path_descarga))
            arch = glob.glob("%s/%s/%s/%s/*%s*.hdf"%(path_descarga,ind,fecha[0:4],fecha[4:7],t))[0]
            os.system("cp %s*%s* %s.xml"%(path_xml,t,arch))
            

def descarga2(url,username,password,tiles,indices,fecha,path_descarga,path_xml):
    for ind in indices:
        for t in tiles:
            os.system('wget -e robots=off -m -np -R .html,.tmp -nH --cut-dirs=3 -A "%s.A%s.%s.006.*.hdf" "%s%s/%s/%s/" --user=%s --password=%s -P %s'%(ind,fecha,t,url,ind,fecha[0:4],fecha[4:7],username,password,path_descarga))
            arch = glob.glob("%s/%s/%s/%s/*%s*.hdf"%(path_descarga,ind,fecha[0:4],fecha[4:7],t))[0]
            os.system("cp %s*%s* %s.xml"%(path_xml,t,arch))
            

#######################################################################################################################
########################################Función para exportar los hdf a geotiff########################################
#######################################################################################################################
def hdfToTiff(inputHDF, outputPath, product, nodata, band):
    ds = gdal.Open(inputHDF) #Abre la imagen
    datasets = ds.GetSubDatasets() #Obtiene los datasets
    for d in datasets: #Se queda con el dataset de interés
        if product in d[0]:
            prod = d[0]

    ds2=gdal.Open(prod) #Abre el dataset
    sds = ds2.GetRasterBand(band).ReadAsArray() #Obtiene la matriz de datos
    geoTs = ds2.GetGeoTransform() #Copia los parámetros de la imagen
    prj = ds2.GetProjection() #Copia la proyección
    dataType = ds2.GetRasterBand(1).DataType #Obtiene el tipo de dato original
    outputTif = outputPath+os.path.basename(inputHDF).replace(".hdf","_sin.tif") #Nombre de imagen de salida

    if not os.path.exists(outputPath): os.makedirs(otputPath) #Crea la ruta si no existe

    # Output raster array to GeoTIFF file            
    driver = gdal.GetDriverByName('GTiff')
    export = driver.Create(outputTif, sds.shape[1], sds.shape[0], 1, dataType)
    band = export.GetRasterBand(band)
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

def fn_tvdi(img_ndvi,img_lst,salida_tvdi):
    #cargar imágenes
    ds_ndvi = gdal.Open(img_ndvi)
    sds_ndvi = ds_ndvi.GetRasterBand(1)
    matriz_ndvi = sds_ndvi.ReadAsArray()
    matriz_ndvi[matriz_ndvi==sds_ndvi.GetNoDataValue()]=np.nan
    
    ds_lst = gdal.Open(img_lst)
    sds_lst = ds_lst.GetRasterBand(1)
    matriz_lst = sds_lst.ReadAsArray()
    matriz_lst[matriz_lst==sds_lst.GetNoDataValue()]=np.nan
    matriz_lst[matriz_lst<=0]=np.nan
    
    ################################
    filas = matriz_lst.shape[0]
    cols = matriz_lst.shape[1]
    matriz_tvdi = np.zeros((filas,cols))
    ancho_filas = 2000
    
    #Parámetros
    min_ndvi = 0. #límite inferior del corte del histograma para calcular la línea de temperatura
    cant_px = 1 #Cantidad de píxeles a tomar de cada delta (10%)
    max_ndvi = np.nanmax(matriz_ndvi)
    delta = 0.01
    
    for i in range(0,filas-ancho_filas+1,1):
        matriz_ndvi_sub = matriz_ndvi[i:i+ancho_filas,:]
        matriz_lst_sub = matriz_lst[i:i+ancho_filas,:]
        
        #Eliminar valores nulos
        nan_data = ~np.isnan(matriz_lst_sub)
        nan_data2 = ~np.isnan(matriz_ndvi_sub)
        nan_data = nan_data*nan_data2
        lst_C_reshape = matriz_lst_sub[nan_data]
        ndvi_reshape = matriz_ndvi_sub[nan_data]

        lst_regr = []    #Lista para almacenar los valores de temperatura para la regresión
        ndvi_regr = []    #Lista para almacenar los valores de vegetación para la regresión
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
    #Salida de la imagen georreferenciada28001139
    geoTs = ds_lst.GetGeoTransform() #Parámetros de la imagen (coordenadas origen y dimensiones)
    driver = gdal.GetDriverByName("GTiff") #Tipo de imagen (geotiff)
    prj = ds_lst.GetProjection() #Sistema de referencia de la imagen (aquí se lee el mismo sistema de referencia de la imagen de entrada)
    print("creando imagen")
    #Crear el espacio
    export=driver.Create(salida_tvdi,matriz_tvdi.shape[1],matriz_tvdi.shape[0],1,GDT_Float32)
    banda=export.GetRasterBand(1) #Cargo la banda creada en el paso anterior
    banda.WriteArray(matriz_tvdi) #Escribe los valores de NDVI en la imagen
    export.SetGeoTransform(geoTs) #Asigna los parametros de la transformacion a la salida
    export.SetProjection(prj) #define la proyección
    banda.FlushCache()#descargar de la memoria virtual al disco
    export.FlushCache()#descargar de la memoria virtual al disco

def subirtiff(img,estilo,user,passw,ipgeo,port):
    namelayer=os.path.basename(img) #Nombre de la capa a subir con extensión
    coveragestore=namelayer[:-4] #Nombre de la capa a subir sin extensión
    workspace=coveragestore.split("_")[0]
        
    #Línea para crear el almacén de datos
    os.system('curl -u '+user+':'+passw+' -XPOST -H "Content-type: application/xml" -d "<coverageStore><name>"'+coveragestore+'"</name><workspace>"'+workspace+'"</workspace><enabled>true</enabled><type>GeoTIFF</type></coverageStore>" http://'+ipgeo+':'+port+'/geoserver/rest/workspaces/'+workspace+'/coveragestores')
    print("coverage creado")
    
    #Línea para cargar la capa
    os.system('curl -u '+user+':'+passw+' -XPUT -H "Content-type:image/tiff" --data-binary @'+img+' http://'+ipgeo+':'+port+'/geoserver/rest/workspaces/'+workspace+'/coveragestores/'+coveragestore+'/file.geotiff')
    print("tiff subido")
    
    #Línea para asignar el estilo a la capa
    os.system('curl -u '+user+':'+passw+' -XPUT -H "Content-type:text/xml" -d "<layer><defaultStyle><name>'+estilo+'</name><workspace>tvdi</workspace></defaultStyle></layer>" http://'+ipgeo+':'+port+'/geoserver/rest/layers/'+workspace+':'+namelayer)
    print("estilo asignado")


def fillgaps(path_tvdi,path_salida,img1,img2):
    img_sal = path_salida+img1[:-4]+"fill.tif"

    ds1 = gdal.Open(path_tvdi+img1)
    sds1 = ds1.GetRasterBand(1).ReadAsArray()

    ds2 = gdal.Open(path_tvdi+img2)
    sds2 = ds2.GetRasterBand(1).ReadAsArray()

    nan_value = sds1[0,0].astype(int)

    mask_nan = np.where(sds1.astype(int)==nan_value)

    sds1[mask_nan] = sds2[mask_nan]

    geoTs = ds1.GetGeoTransform() #Parámetros de la imagen (coordenadas origen y dimensiones)
    driver = gdal.GetDriverByName("GTiff") #Tipo de imagen (geotiff)
    prj = ds1.GetProjection() #Sistema de referencia de la imagen (aquí se lee el mismo sistema de referencia de la imagen de entrada)

    #Crear el espacio
    export=driver.Create(img_sal,sds1.shape[1],sds1.shape[0],1,GDT_Float32)
    banda=export.GetRasterBand(1) #Cargar la banda creada en el paso anterior
    banda.WriteArray(sds1) #Escribir los valores de la matriz calculada
    export.SetGeoTransform(geoTs) #Asignar los parametros de la transformacion a la salida
    export.SetProjection(prj) #definir la proyección
    banda.FlushCache()#descargar de la memoria virtual al disco
    export.FlushCache()#descargar de la memoria virtual al disco

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
