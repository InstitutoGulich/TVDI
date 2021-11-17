#!/usr/bin/env python
# -*- coding: utf-8 -*-

############################################################################################
################# Este script verifica a diario si se encuentra disponible #################
######################## un nuevo dato para procesar el índice TVDI ########################
############################################################################################

import os
from datetime import date, datetime, timedelta
import glob
import numpy as np

path = "/mnt/datos/Repositorio/sancor/indices/tvdi/" #Ruta de los archivos TVDI

fecha_act = datetime.today() #Fecha actual
fecha_final = fecha_act - timedelta(days = 21) #Fecha fin búsqueda 21 = 16 + 5 (días índice + delay procesamiento)
yy_act = fecha_final.year #Año actual
dia_j_act = fecha_final.strftime("%j") #Día juliano actual

dias_interes = np.arange(1,365,16) #Días del año donde hay disponibilidad de imágenes

fecha_ant = fecha_act+timedelta(days=-95) #Fecha desde la cual quiero revisar los datos procesados
dia_j_ant = fecha_ant.strftime("%j") #Día juliano desde el cual quiero revisar los datos procesados

dia_inic = dias_interes[np.argmin(abs(dias_interes-int(dia_j_ant)))] #Día más cercano a la fecha desde la que quiero verificar


if int(dia_j_ant)>int(dia_j_act): #Si el día inicial corresponde a un año anterior
	fecha = datetime.strptime(str(int(yy_act)-1)+"%03d"%dia_inic,"%Y%j")
	fechas_ref = [(fecha+timedelta(days=f)).strftime("%Y%j") for f in range(0,365-int(dia_j_ant),16)]
	fechas_ref.extend([(datetime(yy_act,1,1)+timedelta(days=f)).strftime("%Y%j") for f in range(0,int(dia_j_act),16)])
	print(fechas_ref)

else: #Si el día inicial se encuentra en el año en curso
	fecha = datetime.strptime(str(yy_act)+"%03d"%dia_inic,"%Y%j")
	fechas_ref = [(fecha+timedelta(days=f)).strftime("%Y%j") for f in range(0,int(dia_j_act)-int(dia_inic),16)]

fechas_faltantes = [] #Vector donde voy a almacenar las fechas faltantes

for f_r in fechas_ref: #Verifico por cada fecha la existencia de imágenes
	lista = glob.glob(path+"*%s*"%f_r)
	if lista: #Si ya está procesada no hago nada
		pass
	else:  #En caso de no existir, almaceno la fecha faltante
		fechas_faltantes.append(f_r)  

print(fechas_faltantes)

if len(fechas_faltantes)>0:
	for ff in fechas_faltantes:
		try: #Intento procesar la fecha faltante
			os.system("python3 /mnt/datos/Repositorio/sancor/scripts/01_main_tvdi.py "+ff)
		except: #Si no es posible puede ser que no se han cargado los datos
			pass
else:
	sys.exit(0)

