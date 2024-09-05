---
author: Cristian Iranzo
address: Universidad de Zaragoza
---
# GEOSAT Image Processing

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

## Examples

```text
python D:\geosat\geosat_atmcorrection.py -folders=D:\geosat\data -atm_key=6S -geecloud_id=id
```

## Pending Tasks

- Functionality to execute 6S with parameters extracted from ground weather stations (without the need for EearthEngine API).
- Consolidate all functionalities into an executable. Use tools like [PyInstaller](https://pyinstaller.org/en/stable/index.html).

## References



