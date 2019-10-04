#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  4 10:56:15 2019

@author: macbookair
"""
import numpy as np
from matplotlib import pyplot as plt
import cartopy.feature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cartopy.crs as ccrs
import os
#%%


##- abrir .nc con xarray
import xarray as xr

ruta="/Users/macbookair/Desktop/practicasCG/atmo/Practica 1/simulacion 1/" #ruta CARO
os.chdir(ruta)
dir = 'shallow.nc'
dS = xr.open_dataset(dir, decode_times=False)
print(dS)       # visualizo la info del .nc

lat = dS['lat'].values
lon = dS['lon'].values
time = dS['time'].values
fr = dS['fr'].values
ucomp = dS['ucomp'].values
vor = dS['vor'].values

#Pasamos las latitudes/longitudes del dataset a una reticula para graficar

lons, lats = np.meshgrid(lon, lat)

#%% Graficamos el forzante

#plt.imshow(fr[0,:,:]) plt.colorbar() #hago un grafico pra visualizar min y max

cmin = 40000
cmax = 50000
ncont = 11 #me queda cada 1000
clevs = np.linspace(cmin, cmax, ncont)

#Creamos figura
fig=plt.figure(figsize=(6,4),dpi=200)
LONMIN= 0
LONMAX= 359.9
LATMIN = -90
LATMAX = 90
#Definimos proyección
ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=180))
crs_latlon = ccrs.PlateCarree()
ax.set_extent([LONMIN, LONMAX, LATMIN, LATMAX], crs=crs_latlon)

#Graficamos
im=ax.contourf(lons, lats, fr[0,:,:], clevs, cmap=plt.get_cmap("jet"), extend='both', transform=crs_latlon)

#Agregamos barra de colores
plt.colorbar(im, fraction=0.052, pad=0.04, shrink=0.8, aspect=8)

#Características del mapa
ax.add_feature(cartopy.feature.LAND, facecolor='#d9d9d9')
ax.add_feature(cartopy.feature.COASTLINE)
ax.add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=.5)
ax.gridlines(crs=crs_latlon, linewidth=0.3, linestyle='-')
ax.set_xticks(np.arange(0, 360, 60), crs=crs_latlon)
ax.set_yticks(np.arange(-90, 90, 30), crs=crs_latlon)
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)

#Características del mapa
#Titulo
plt.title("Forzante", fontsize=8, y=0.98, loc="center")

#Guardar figura
plt.savefig('Forzante.jpg')


#%%

#%% Graficamos U

plt.imshow(ucomp[49,:,:]) , plt.colorbar() #hago un grafico pra visualizar min y max

cmin = 0
cmax = 1.2
ncont = 7 #me queda cada 1000
clevs = np.linspace(cmin, cmax, ncont)

#Creamos figura
fig=plt.figure(figsize=(6,4),dpi=200)
LONMIN= 0
LONMAX= 359.9
LATMIN = -90
LATMAX = 90
#Definimos proyección
ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=180))
crs_latlon = ccrs.PlateCarree()
ax.set_extent([LONMIN, LONMAX, LATMIN, LATMAX], crs=crs_latlon)

#Graficamos
im=ax.contourf(lons, lats, ucomp[49,:,:], clevs, cmap=plt.get_cmap("jet"), extend='both', transform=crs_latlon)

#Agregamos barra de colores
plt.colorbar(im, fraction=0.052, pad=0.04, shrink=0.8, aspect=8)

#Características del mapa
ax.add_feature(cartopy.feature.LAND, facecolor='#d9d9d9')
ax.add_feature(cartopy.feature.COASTLINE)
ax.add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=.5)
ax.gridlines(crs=crs_latlon, linewidth=0.3, linestyle='-')
ax.set_xticks(np.arange(0, 360, 60), crs=crs_latlon)
ax.set_yticks(np.arange(-90, 90, 30), crs=crs_latlon)
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)

#Características del mapa
#Titulo
plt.title("Viento Zonal", fontsize=8, y=0.98, loc="center")

#Guardar figura
plt.savefig('Viento Zonal.jpg')

#%% Gradiente de la vorticidad

#primero hago el calculo de gradiente de cortividad con rutina dada

from DerY import derivy

# dy=156543.75 porque pase la cantidad de puntos que tengo a grados luego converti grados en metros y obtuve la disancia en y en vez de en puntos en metros180/128*111320

gvor = derivy(vor[49,:,:],156543.75)

#Graficamos gradiente de vorticidad

plt.imshow(gvor) , plt.colorbar() #hago un grafico pra visualizar min y max

cmin = -6e-13
cmax = 6e-13
ncont = 7 #me queda cada 2
clevs = np.linspace(cmin, cmax, ncont)

#Creamos figura
fig=plt.figure(figsize=(6,4),dpi=200)
LONMIN= 0
LONMAX= 359.9
LATMIN = -90
LATMAX = 90
#Definimos proyección
ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=180))
crs_latlon = ccrs.PlateCarree()
ax.set_extent([LONMIN, LONMAX, LATMIN, LATMAX], crs=crs_latlon)

#Graficamos
im=ax.contourf(lons, lats, gvor, clevs, cmap=plt.get_cmap("jet"), extend='both', transform=crs_latlon)

#Agregamos barra de colores
plt.colorbar(im, fraction=0.052, pad=0.04, shrink=0.8, aspect=8)

#Características del mapa
ax.add_feature(cartopy.feature.LAND, facecolor='#d9d9d9')
ax.add_feature(cartopy.feature.COASTLINE)
ax.add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=.5)
ax.gridlines(crs=crs_latlon, linewidth=0.3, linestyle='-')
ax.set_xticks(np.arange(0, 360, 60), crs=crs_latlon)
ax.set_yticks(np.arange(-90, 90, 30), crs=crs_latlon)
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)

#Características del mapa
#Titulo
plt.title("Gradiente de vorticidad", fontsize=8, y=0.98, loc="center")

#Guardar figura
plt.savefig('Gradiente de vorticidad.jpg')

#%% Numero de onda estacionaria
lats_para_k=lats*np.pi/180
lat_para_k=lat*np.pi/180

beta=1e-11
Rt= 6371000
betha=2*7.29e-5*np.cos(lats_para_k)/Rt # hay que calcularla no constante sino variando con la lat
Ks= np.sqrt((betha-gvor)/ucomp[49,:,:])
Ks_betacte=np.sqrt((beta-gvor)/ucomp[49,:,:])

i=1
Ksp=np.empty_like(Ks)
for i in np.arange(0,np.size(lat_para_k)):
    for j in np.arange(0,np.size(lon)):
        Ksp[i,j]=Rt*np.cos(lat_para_k[i])*Ks[i,j]
        
i=1
Ksp_betacte=np.empty_like(Ks_betacte)
for i in np.arange(0,np.size(lat_para_k)):
    for j in np.arange(0,np.size(lon)):
        Ksp_betacte[i,j]=Rt*np.cos(lat_para_k[i])*Ks_betacte[i,j]       
        
Lean=Ksp - Ksp_betacte # la diferencia entre beta cte y no cte
#Graficamos Ks

plt.imshow(Ksp) , plt.colorbar() #hago un grafico pra visualizar min y max
plt.imshow(Ksp_betacte) , plt.colorbar() #hago un grafico pra visualizar min y max
plt.imshow(Lean) , plt.colorbar() #hago un grafico pra visualizar min y max

cmin = 0 
cmax = 10
ncont = 10 #me queda cada 1
clevs = np.linspace(cmin, cmax, ncont)

#Creamos figura
fig=plt.figure(figsize=(6,4),dpi=200)
LONMIN= 0
LONMAX= 359.9
LATMIN = -90
LATMAX = 90
#Definimos proyección
ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=180))
crs_latlon = ccrs.PlateCarree()
ax.set_extent([LONMIN, LONMAX, LATMIN, LATMAX], crs=crs_latlon)

#Graficamos
im=ax.contourf(lons, lats, Ksp , clevs, cmap=plt.get_cmap("jet"), extend='both', transform=crs_latlon)

#Agregamos barra de colores
plt.colorbar(im, fraction=0.052, pad=0.04, shrink=0.8, aspect=8)

#Características del mapa
ax.add_feature(cartopy.feature.LAND, facecolor='#d9d9d9')
ax.add_feature(cartopy.feature.COASTLINE)
ax.add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=.5)
ax.gridlines(crs=crs_latlon, linewidth=0.3, linestyle='-')
ax.set_xticks(np.arange(0, 360, 60), crs=crs_latlon)
ax.set_yticks(np.arange(-90, 90, 30), crs=crs_latlon)
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
#Características del mapa
#Titulo
plt.title("Ks por circulo de latitud", fontsize=8, y=0.98, loc="center")

#Guardar figura
plt.savefig('Ks por circulo de latitud.jpg')

#%% Vamos con la Simulacion 2

ruta2="/Users/macbookair/Desktop/practicasCG/atmo/Practica 1/simulacion 2/" #ruta CARO
os.chdir(ruta2)
dir = 'shallow.nc'
dS = xr.open_dataset(dir, decode_times=False)
print(dS)       # visualizo la info del .nc

lat = dS['lat'].values
lon = dS['lon'].values
time = dS['time'].values
fr = dS['fr'].values
ucomp = dS['ucomp'].values
vor = dS['vor'].values

#Pasamos las latitudes/longitudes del dataset a una reticula para graficar

lons, lats = np.meshgrid(lon, lat)
#%% Graficamos el forzante

plt.imshow(fr[49,:,:]) 
plt.colorbar() #hago un grafico pra visualizar min y max

cmin = 40000
cmax = 110000
ncont = 7 #me queda cada 1000
clevs = np.linspace(cmin, cmax, ncont)

#Creamos figura
fig=plt.figure(figsize=(6,4),dpi=200)
LONMIN= 0
LONMAX= 359.9
LATMIN = -90
LATMAX = 90
#Definimos proyección
ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=180))
crs_latlon = ccrs.PlateCarree()
ax.set_extent([LONMIN, LONMAX, LATMIN, LATMAX], crs=crs_latlon)

#Graficamos
im=ax.contourf(lons, lats, fr[49,:,:], clevs, cmap=plt.get_cmap("jet"), extend='both', transform=crs_latlon)

#Agregamos barra de colores
plt.colorbar(im, fraction=0.052, pad=0.04, shrink=0.8, aspect=8)

#Características del mapa
ax.add_feature(cartopy.feature.LAND, facecolor='#d9d9d9')
ax.add_feature(cartopy.feature.COASTLINE)
ax.add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=.5)
ax.gridlines(crs=crs_latlon, linewidth=0.3, linestyle='-')
ax.set_xticks(np.arange(0, 360, 60), crs=crs_latlon)
ax.set_yticks(np.arange(-90, 90, 30), crs=crs_latlon)
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)

#Características del mapa
#Titulo
plt.title("Forzante", fontsize=8, y=0.98, loc="center")

#Guardar figura
plt.savefig('Forzante.jpg')



#%% Graficamos U

plt.imshow(ucomp[49,:,:]) , plt.colorbar() #hago un grafico pra visualizar min y max

cmin = -10
cmax = 25
ncont = 7 #me queda cada 1000
clevs = np.linspace(cmin, cmax, ncont)

#Creamos figura
fig=plt.figure(figsize=(6,4),dpi=200)
LONMIN= 0
LONMAX= 359.9
LATMIN = -90
LATMAX = 90
#Definimos proyección
ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=180))
crs_latlon = ccrs.PlateCarree()
ax.set_extent([LONMIN, LONMAX, LATMIN, LATMAX], crs=crs_latlon)

#Graficamos
im=ax.contourf(lons, lats, ucomp[49,:,:], clevs, cmap=plt.get_cmap("jet"), extend='both', transform=crs_latlon)

#Agregamos barra de colores
plt.colorbar(im, fraction=0.052, pad=0.04, shrink=0.8, aspect=8)

#Características del mapa
ax.add_feature(cartopy.feature.LAND, facecolor='#d9d9d9')
ax.add_feature(cartopy.feature.COASTLINE)
ax.add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=.5)
ax.gridlines(crs=crs_latlon, linewidth=0.3, linestyle='-')
ax.set_xticks(np.arange(0, 360, 60), crs=crs_latlon)
ax.set_yticks(np.arange(-90, 90, 30), crs=crs_latlon)
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)

#Características del mapa
#Titulo
plt.title("Viento Zonal", fontsize=8, y=0.98, loc="center")

#Guardar figura
plt.savefig('Viento Zonal.jpg')

#%% Gradiente de la vorticidad

#primero hago el calculo de gradiente de cortividad con rutina dada

from DerY import derivy

# dy=156543.75 porque pase la cantidad de puntos que tengo a grados luego converti grados en metros y obtuve la disancia en y en vez de en puntos en metros180/128*111320

gvor = derivy(vor[49,:,:],156543.75)

#Graficamos gradiente de vorticidad

plt.imshow(gvor) , plt.colorbar() #hago un grafico pra visualizar min y max

cmin = -6e-12
cmax = 6e-12
ncont = 7 #me queda cada 2
clevs = np.linspace(cmin, cmax, ncont)

#Creamos figura
fig=plt.figure(figsize=(6,4),dpi=200)
LONMIN= 0
LONMAX= 359.9
LATMIN = -90
LATMAX = 90
#Definimos proyección
ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=180))
crs_latlon = ccrs.PlateCarree()
ax.set_extent([LONMIN, LONMAX, LATMIN, LATMAX], crs=crs_latlon)

#Graficamos
im=ax.contourf(lons, lats, gvor, clevs, cmap=plt.get_cmap("jet"), extend='both', transform=crs_latlon)

#Agregamos barra de colores
plt.colorbar(im, fraction=0.052, pad=0.04, shrink=0.8, aspect=8)

#Características del mapa
ax.add_feature(cartopy.feature.LAND, facecolor='#d9d9d9')
ax.add_feature(cartopy.feature.COASTLINE)
ax.add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=.5)
ax.gridlines(crs=crs_latlon, linewidth=0.3, linestyle='-')
ax.set_xticks(np.arange(0, 360, 60), crs=crs_latlon)
ax.set_yticks(np.arange(-90, 90, 30), crs=crs_latlon)
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)

#Características del mapa
#Titulo
plt.title("Gradiente de vorticidad", fontsize=8, y=0.98, loc="center")

#Guardar figura
plt.savefig('Gradiente de vorticidad.jpg')

#%% Numero de onda estacionaria
lats_para_k=lats*np.pi/180

Rt= 6371000
betha=2*7.29e-5*np.cos(lats_para_k)/Rt # hay que calcularla no constante sino variando con la lat
Ks= np.sqrt((betha-gvor)/ucomp[49,:,:])

i=1
Ksp=np.empty_like(Ks)
for i in np.arange(0,np.size(lat_para_k)):
    for j in np.arange(0,np.size(lon)):
        Ksp[i,j]=Rt*np.cos(lat_para_k[i])*Ks[i,j]
        
        
#Graficamos Ks

plt.imshow(Ksp) , plt.colorbar() #hago un grafico pra visualizar min y max

cmin = 0 
cmax = 20
ncont = 10 #me queda cada 1
clevs = np.linspace(cmin, cmax, ncont)

#Creamos figura
fig=plt.figure(figsize=(6,4),dpi=200)
LONMIN= 0
LONMAX= 359.9
LATMIN = -90
LATMAX = 90
#Definimos proyección
ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=180))
crs_latlon = ccrs.PlateCarree()
ax.set_extent([LONMIN, LONMAX, LATMIN, LATMAX], crs=crs_latlon)

#Graficamos
im=ax.contourf(lons, lats, Ksp , clevs, cmap=plt.get_cmap("jet"), extend='both', transform=crs_latlon)

#Agregamos barra de colores
plt.colorbar(im, fraction=0.052, pad=0.04, shrink=0.8, aspect=8)

#Características del mapa
ax.add_feature(cartopy.feature.LAND, facecolor='#d9d9d9')
ax.add_feature(cartopy.feature.COASTLINE)
ax.add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=.5)
ax.gridlines(crs=crs_latlon, linewidth=0.3, linestyle='-')
ax.set_xticks(np.arange(0, 360, 60), crs=crs_latlon)
ax.set_yticks(np.arange(-90, 90, 30), crs=crs_latlon)
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
#Características del mapa
#Titulo
plt.title("Ks por circulo de latitud", fontsize=8, y=0.98, loc="center")

#Guardar figura
plt.savefig('Ks por circulo de latitud.jpg')