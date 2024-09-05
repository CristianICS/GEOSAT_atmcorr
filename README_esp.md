---
autor: Cristian Iranzo
afiliación: Universidad de Zaragoza
---
# Tratamiento de imágenes GEOSAT

[![DOI](https://zenodo.org/badge/852767833.svg)](https://zenodo.org/doi/10.5281/zenodo.13694180)

Aplicar correcciones atmosféricas a las imágenes GEOSAT. El código es adaptable a todo tipo de imágenes, pero ha sido desarrollado sobre imágenes GEOSAT con las siguietnes características:

- PM4 o *Bundle*: imagen pancromática (1m GSD), y las cuatro bandas de la imagen multiespectral (3m GSD)
- Nivel de tratamiento L1C: Producto calibrado, corregido radiométricamente y ortorrectificado.

| Banda | Min $\lambda$ (nm) | Max $\lambda$ (nm) | GSD (m) |
| ----- | ------------------ | ------------------ | ------- |
| PAN   | 560                | 900                | 0.75    |
| BLUE  | 466                | 525                | 3       |
| GREEN | 532                | 599                | 3       |
| RED   | 640                | 697                | 3       |
| NIR   | 770                | 892                | 3       |

**Importante**: El directorio que contiene la imagen debe mantener los nombres y la estructura original de los archivos `tif` y `dim`:

```text
DE2_PM4_L1C_000000_20210629T102606_20210629T102608_DE2_38083_838B
|_ DE2_MS4*.dim
|_ DE2_MS4*.tif
|_ DE2_PAN*.dim
|_ DE2_PAN*.tif
```

## Instalación

El código se ha ejecutado en [miniconda](https://docs.anaconda.com/free/miniconda/index.html), dentro del siguiente ambiente:

```text
conda create -n geosat_atmcorr python py6s earthengine-api gdal -c conda-forge
```

Los módulos `py6s` y `earthengine-api` son específicos para la corrección 6S.

Una vez instalado miniconda (o anaconda), descargar los dos códigos de python en la misma carpeta:

- `atmcorr_utils.py`: Clases y funciones para aplicar la corrección.
- `geosat_atmcorr.py`: Código que genera las correcciones

Se ha testado con las siguientes versiones:

```text
earthengine-api=0.1.401 (pyhd8ed1ab_0)
gdal=3.9.0 (py312hea5013e_2)
py6s=1.9.2 (pyhd8ed1ab_0)
python=3.12.3 (h2628c8c_0_cpython)
```

## Aplicar la corrección

Se realiza a través de los siguientes comandos:

- `folder`

  Directorio de las imágenes GEOSAT a corregir.

- `folders`

  Directorio principal con las carpetas de las imágenes GEOSAT. Si se incluye este parámetro en lugar de `folder`, se ejecuta una corrección en bucle de todas las imágenes.

  - Permite que haya otros archivos que no son directorios de imágenes, no los tendrá en cuenta.

- `atm_key`

  Código de la corrección atmosférica a utilizar. Elegir entre `DOS`, `COST` y `6S`.

  **Importante**: Todas las correcciones requieren la presencia de una imagen en niveles de radiancia en el directorio, i.e., el producto resultante de aplicar la transformación `ARC`. Es decir, que para aplicar las correcciones deberá ejecutarse primero la transformación anterior, que creará la imagen ARC con el prefijo `<original_imgname>_ARC.tif`.

  ```text
  python D:\geosat\geosat_atmcorrection.py -folder=D:\geosat\DE2_* -atm_key=ARC
  ```

- `geecloud_id`

  Si se aplica la corrección `6S`, deberá utilizarse este comando para incluir el id del proyecto de Google Cloud con el que ejecutar las funciones de la API de GEE.

  *Nota*: Previamente se deberá incluir el comando `earthengine authenticate` en la consola de conda y registrarse con la cuenta de google asociada al id del proyecto.

Al finalizar la corrección se genera un archivo `.log` con los valores aplicados y las fórmulas correspondientes a las correcciones de cada banda.

**Importante**: Todas las rutas deben ser absolutas (comenzando por la unidad de disco), excepto si los códigos de python con las funciones se localizan en la carpeta del ambiente de conda.

**Importante**: Si la imagen a corregir ya existe, no se sobrescribirá y se saltará a la siguiente.

*Nota*: El código debe aplicarse dentro del ambiente que contiene los paquetes requeridos. En el caso del ejemplo anterior:

```text
conda activate geosat_atmcorr
```

*Nota*: El ID de Google Cloud de una cuenta puede consultarse en la [consola de Google Cloud](https://console.cloud.google.com/?hl=es).

## Ejemplos

```text
python D:\geosat\geosat_atmcorrection.py -folders=D:\geosat\data -atm_key=6S -geecloud_id=id
```

## Tareas pendientes

- Funcionalidad para ejecutar 6S con parámetros extraidos de estaciones meteorológicas terrestres (sin necesidad de Google Cloud).
- Englobar todas las funcionalidades en un ejecutable. Utilizar herramientas como [PyInstaller](https://pyinstaller.org/en/stable/index.html)

## Referencias

Chavez, P. S. & others. (1996). Image-based atmospheric corrections-revisited and improved. Photogrammetric Engineering and Remote Sensing, 62(9), 1025–1035.

Fernández, C., de Castro, C., Garcı́a, L., Calleja, M. E., Niño, R., Fraile, S., & Sousa, R. (2023). Evaluación del impacto de la superresolución sobre imágenes multiespectrales GEOSAT-2. Revista de Teledetección, 61, 83–96.

Mahiny, A. S., & Turner, B. J. (2007). A comparison of four common atmospheric correction methods. Photogrammetric Engineering & Remote Sensing, 73(4), 361–368.

[Murphy, S. & Hård, J. `gee-atmcorr-S2` repository](https://github.com/samsammurphy/gee-atmcorr-S2/tree/master)
 
Thuillier, G., Hersé, M., Labs, D., Foujols, T., Peetermans, W., Gillotay, D., Simon, P., & Mandel, H. (2003). The solar spectral irradiance from 200 to 2400 nm as measured by the SOLSPEC spectrometer from the ATLAS and EURECA missions. Solar Physics, 214, 1–22.

Vermote, E. F., Tanré, D., Deuze, J. L., Herman, M., & Morcette, J.-J. (1997). Second simulation of the satellite signal in the solar spectrum, 6S: An overview. IEEE Transactions on Geoscience and Remote Sensing, 35(3), 675–686.

Wilson, R. T. (2013). Py6S: A Python interface to the 6S radiative transfer model. Computers & Geosciences, 51, 166–171.
