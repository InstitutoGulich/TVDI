from funciones import filtro_media, getNodata
import os

def anomalia(tvdi, mean, sd, output):

        os.system('gdal_calc.py '\
                '-A %s '\
                '-B %s '\
                '-C %s '\
                '--outfile=%s '\
###                '--calc="(A-B)/(C)" '\
###                '--NoDataValue=-9999 '\
###                '--type="Float32"'\
                '--calc="(A-B)/(C)*1000" '\
                '--NoDataValue=-9999 '\
                '--type="Int16"'\
                %(tvdi, mean, sd, output))

def recortexTIFF(tiff,mask,nodata,salida):
	os.system('gdal_calc.py -A %s -B %s --outfile=%s --calc="A*B" --NoDataValue=%s --type="Int16" --overwrite'%(tiff,mask,salida,nodata))
	
path = "/home/jrubio/Documentos/VENVS/SANCOR5/" #Directorio principal de procesamiento
path_estadisticos = path + "estadisticos/"
path_anomalias = path + "anomalias/"
path_indices = path + "indices/"
scale = 0.002543284457288304268 #250m

uw_mask = path + "datos/mascara/mask_urbano_agua_i5_250m.tif" #Raster máscara zona de zonas urbanas y cuerpos de agua
uw_mask2 = path + "datos/mascara/mask_urbano_agua_5_250m.tif" #Raster máscara zona de zonas urbanas y cuerpos de agua
uw_mask = path + "datos/mascara/mask_urbano_agua_i3_250m.tif" #Raster máscara zona de zonas urbanas y cuerpos de agua
uw_mask2 = path + "datos/mascara/mask_urbano_agua_3_250m.tif" #Raster máscara zona de zonas urbanas y cuerpos de agua

uw_mask = path + "datos/mascara/mask_cultivos_i3_250m.tif" #Raster máscara zona de zonas urbanas y cuerpos de agua
uw_mask2 = path + "datos/mascara/mask_cultivos_3_250m.tif" #Raster máscara zona de zonas urbanas y cuerpos de agua

#ANOMALIA
print("Calculando la anomalía")
path_stats_mean = path_estadisticos + "tvdi/mean/"
path_stats_sd = path_estadisticos + "tvdi/std/"
path_salida = path_anomalias + "tvdi/anomalias_m/"

for d in range(65,82,16):
	fecha = "2023%03d"%d
#3#	media = path_stats_mean + "tvdi_%s_fill_mean_2000_2019_500m.tif"%fecha[-3:]
#3#	sd = path_stats_sd + "tvdi_%s_fill_std_2000_2019_500m.tif"%fecha[-3:]
	
	media = path_stats_mean + "tvdi_%s_fill_mean_2000_2019_250m.tif"%fecha[-3:]
	sd = path_stats_sd + "tvdi_%s_fill_std_2000_2019_250m.tif"%fecha[-3:]
	
#3#	tvdi_mask = path_indices + "tvdi/tvdi_%s_500m.tif"%fecha
#3#	tvdi_anom = path_salida + "anom_%s_500m.tif"%fecha
	
	tvdi_mask = path_indices + "tvdi/tvdicrop_%s_250m_f.tif"%fecha
	tvdi_anom = path_salida + "anom_%s_250m.tif"%fecha
	
	nodata = getNodata(uw_mask)

#3#	anom_out = path_salida + "anom_%s_500m.tif"%fecha
#3#	tvdi_anom = path_salida + "anom_%s_500mfill_m.tif"%fecha

	anom_out = path_salida + "anomcrop_%s_250m.tif"%fecha
	tvdi_anom = path_salida + "anomcrop_%s_250mfill_m.tif"%fecha

	anomalia(tvdi_mask, media, sd, tvdi_anom)
	recortexTIFF(tvdi_anom, uw_mask, nodata, anom_out)

#ENMASCARAMIENTO Y REMUESTREO A 5 KM USANDO LA MEDIA
##ANOMALIA
	smooth_anom_path = path_anomalias + "tvdi/anomalias_m_5km/" #Directorio resultado de procesamiento
#3#	out_anom_5km_filt = filtro_media(anom_out, smooth_anom_path, uw_mask, uw_mask2, scale)
	out_anom_5km_filt = filtro_media(anom_out, smooth_anom_path, uw_mask, uw_mask2, scale, 21, 21)
	recortexTIFF(out_anom_5km_filt, uw_mask, nodata, out_anom_5km_filt)
	os.system('gdal_calc.py -A %s -B %s --outfile=%s --calc="A*B" --type="Int16" --NoDataValue %d --co="COMPRESS=DEFLATE" --overwrite'%(out_anom_5km_filt, uw_mask, out_anom_5km_filt, nodata))

	os.system("rm " + tvdi_anom)
