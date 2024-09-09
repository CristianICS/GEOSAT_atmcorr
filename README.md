---
author: Cristian Iranzo
address: Universidad de Zaragoza
---
# GEOSAT Image Processing

[![DOI](https://zenodo.org/badge/852767833.svg)](https://zenodo.org/doi/10.5281/zenodo.13694180)

Apply atmospheric corrections to GEOSAT images. The code is adaptable to all types of images, but it has been developed to work with GEOSAT ones and with the following characteristics:

- PM4 or *Bundle*: panchromatic image (1m GSD), and the four bands of the multispectral image (3m GSD).
- Processing level L1C (Ortho): Calibrated and radiometrically corrected scaled ToA radiance, orthorectified and resampled.

| Band | Min $\lambda$ (nm) | Max $\lambda$ (nm) | GSD (m) |
| ---- | ------------------ | ------------------ | ------- |
| PAN  | 560                | 900                | 0.75    |
| BLUE | 466                | 525                | 3       |
| GREEN| 532                | 599                | 3       |
| RED  | 640                | 697                | 3       |
| NIR  | 770                | 892                | 3       |

**Important**: The directory containing the image(s) must maintain the original structure and names of the `tif` and `dim` files:

```text
DE2_PM4_L1C_000000_20210629T102606_20210629T102608_DE2_38083_838B
|_ DE2_MS4*.dim
|_ DE2_MS4*.tif
|_ DE2_PAN*.dim
|_ DE2_PAN*.tif
```

## Installation

The code was executed using [miniconda](https://docs.anaconda.com/free/miniconda/index.html) within the following environment:

```text
conda create -n geosat_atmcorr python py6s earthengine-api gdal -c conda-forge
```

The `py6s` and `earthengine-api` modules are specific to the 6S correction.

Once miniconda (or anaconda) is installed, download the two python scripts into the same folder:

- `atmcorr_utils.py`: Classes and functions to apply the correction.
- `geosat_atmcorr.py`: Code that generates the corrections.

It has been tested with the following versions:

```text
earthengine-api=0.1.401 (pyhd8ed1ab_0)
gdal=3.9.0 (py312hea5013e_2)
py6s=1.9.2 (pyhd8ed1ab_0)
python=3.12.3 (h2628c8c_0_cpython)
```

## Applying the correction

Use the following commands:

- `folder`

  Directory of the GEOSAT images to correct.

- `folders`

  Main directory containing the folders of the GEOSAT images. If this parameter is included instead of `folder`, a loop correction of all the images is executed.

  - It allows for other files that are not image directories, which will be ignored.

- `atm_key`

  Atmospheric correction code to use. Choose from `DOS`, `COST`, and `6S`.

  **Important**: All corrections will require an image in radiance levels, i.e., the resulting product from applying the `ARC` transformation. To apply these, the previous transformation must be executed first, which will create a file with the prefix `<original_imgname>_ARC.tif` in the folder containing the images.

  ```text
  python D:\geosat\geosat_atmcorrection.py -folder=D:\geosat\DE2_* -atm_key=ARC
  ```

- `geecloud_id`

  If the `6S` correction is applied, this command should be used to include the Google Cloud project ID for executing GEE API functions.

  *Note*: You must previously include the `earthengine authenticate` command in the conda console and sign in with the Google account associated with the project ID.

After the correction, a `.log` file is generated with the applied values and the formulas corresponding to the corrections for each band.

**Important**: All paths must be absolute (starting from the disk drive), unless the python scripts with the functions are located in the conda environment folder.

**Important**: If the image to be corrected already exists, it will not be overwritten, and the next one will be processed.

*Note*: The code must be run within the environment containing the required packages. In the example provided:

```text
conda activate geosat_atmcorr
```

*Note*: The Google Cloud account ID can be found in the [Google Cloud Console](https://console.cloud.google.com/?hl=en).

## Unzipping files

The GEOSAT images are downloaded in zip folders. To unzip them, execute the following code in python:

```python
import os
import zipfile

main_folder = r"D:\GEOSAT"

def unzip(folder):
    """Unzip folders"""
    if folder.endswith('.zip'):
        # Select the location to unzip the folder content
        unzip_folder_name = os.path.basename(folder).split('.')[0]
        unzip_folder = os.path.join(os.path.dirname(folder), unzip_folder_name)
        if not os.path.exists(unzip_folder):
            # Unzip
            with zipfile.ZipFile(folder, 'r') as zip_ref:
                zip_ref.extractall(unzip_folder)

for folder in os.listdir(main_folder):
    unzip(os.path.join(main_folder, folder))
```

## Examples

```text
python D:\geosat\geosat_atmcorrection.py -folders=D:\geosat\data -atm_key=6S -geecloud_id=id
```

## Pending Tasks

- Functionality to execute 6S with parameters extracted from ground weather stations (without the need for EearthEngine API).
- Consolidate all functionalities into an executable. Use tools like [PyInstaller](https://pyinstaller.org/en/stable/index.html).

## References

Chavez, P. S. & others. (1996). Image-based atmospheric corrections-revisited and improved. Photogrammetric Engineering and Remote Sensing, 62(9), 1025–1035.

Fernández, C., de Castro, C., Garcı́a, L., Calleja, M. E., Niño, R., Fraile, S., & Sousa, R. (2023). Evaluación del impacto de la superresolución sobre imágenes multiespectrales GEOSAT-2. Revista de Teledetección, 61, 83–96.

Fernández, C., De Castro, C., Calleja, M. E., Sousa, R., Niño, R., García, L., Fraile, S., & Molina, I. (2023). GEOSAT 2 Atmospherically Corrected Images: Algorithm Validation. ECRS 2023, 64. https://doi.org/10.3390/ECRS2023-16296

Mahiny, A. S., & Turner, B. J. (2007). A comparison of four common atmospheric correction methods. Photogrammetric Engineering & Remote Sensing, 73(4), 361–368.

[Murphy, S. & Hård, J. `gee-atmcorr-S2` repository](https://github.com/samsammurphy/gee-atmcorr-S2/tree/master)
 
Thuillier, G., Hersé, M., Labs, D., Foujols, T., Peetermans, W., Gillotay, D., Simon, P., & Mandel, H. (2003). The solar spectral irradiance from 200 to 2400 nm as measured by the SOLSPEC spectrometer from the ATLAS and EURECA missions. Solar Physics, 214, 1–22.

Vermote, E. F., Tanré, D., Deuze, J. L., Herman, M., & Morcette, J.-J. (1997). Second simulation of the satellite signal in the solar spectrum, 6S: An overview. IEEE Transactions on Geoscience and Remote Sensing, 35(3), 675–686.

Wilson, R. T. (2013). Py6S: A Python interface to the 6S radiative transfer model. Computers & Geosciences, 51, 166–171.

