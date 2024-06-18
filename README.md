# plotting_maps
The aim of this repo is to plot standardised maps of ACS climate hazard data. We will develop and provide examples of mapping climate hazards for Australia so that data can be consistently presented.

Examples will include maps for Australia and for selected states or regions. 

Intended uses include taking netcdf or xarray dataarrays of hazards and indices such as Rx1day, TXx, FFDI and plotting the data on a map of Australia. 

Developed by Gen Tolhurst (gentolhurst@gmail.com), Supervised by Mitchell Black (Mitchell.Black@bom.gov.au).\
Funded by ACS

## Python environment
This code is designed to work with hh5 analysis3-24.04 virtual environment.

In your terminal, this may look like:

```
$ module use /g/data/hh5/public/modules
$ module load conda_concept/analysis3-24.04
```

When starting a new ARE JupyterLab session (https://are.nci.org.au/pun/sys/dashboard/batch_connect/sys/jupyter/ncigadi/session_contexts/new, requires NCI login), selecting the hh5 analysis3-24.04 virtual environment might look like this:

![image](https://github.com/AusClimateService/plotting_maps/assets/45543810/6607e78a-8599-4eeb-8cee-4e910e067d5a)

## Cloning this repo
Before you can ```import acs_plotting_maps``` to use the plotting function ```plot_aus_shapefiles```, you will need to clone a copy of this repository to your own working directory.

If you are working in your home directory, navigate there:
```
$ cd ~/
```

Else, if you are working elsewhere (eg. scratch or project), specify the path:
```
$ cd /path/to/dir/

$ cd /scratch/PROJECT/USER/

$ cd /g/data/PROJECT/USER/
```

Then, you can clone this repository to access the python code and notebooks. \
If you want the new directory to be called anything other than "plotting_maps" please replace the final argument to your choice of directory name:
```
$ git clone https://github.com/AusClimateService/plotting_maps.git plotting_maps
```
You will now be able to access the functions, python scripts, and Jupyter notebooks from your user.

## Usage in Jupyter notebook:
1. **Navigate to the directory you cloned to:**
```
cd ~/plotting_maps
```

2. **Import the ACS plotting maps function.** This will also import dependencies including xarray and pandas
```python 
from acs_plotting_maps import *
```

3. **Load some data.** For example, this will load extratropical storm rx5day rainfall
```python
ds = xr.open_dataset("/g/data/ia39/ncra/extratropical_storms/RX5D_AGCD-05i_ACCESS-CM2_historical_r4i1p1f1_BOM_BARPA-R_v1-r1_annual.nc")
```
This data has three dimensions (time, lon, lat). There is a value for every year from 1960 to 2015. We can only plot 2D, so next we will calculate a statistic to summarise the data

4. **Summarise data into a 2D xr.DataArray.** For example, calculate the median:
```python
da = ds.median(dim="time")
```

5. **Finally, use the plotting function**.\
You will need to specify:
     * the data (and select the variable eg "pr");
     * suitable arguments for the colorbar including cmap, ticks, cbar_label, and cbar_extend;
     * annotations including title, dataset_name, date_range; and
     * where you want the image outfile saved.
   
```python
plot_aus_shapefiles(data = da["pr"],
                    regions = regions_dict['ncra_regions'],
                    cmap = cmap_dict["pr"],
                    ticks = tick_dict['pr_mon'],
                    cbar_label = "rainfall [mm]",
                    cbar_extend = "max",
                    title = "Extratropical storms Rx5day median",
                    dataset_name = "BARPA-R ACCESS-CM2",
                    date_range = '01 January 1960 to 01 January 2015',
                    outfile = "figures/outfile.png",
                   );
```
![Extratropical_storms_Rx5day_median](https://github.com/AusClimateService/plotting_maps/assets/45543810/b5735647-c886-4d35-b230-aee7c8012a0c)

## Suggested regions, colormaps and scales for plotting
Using suggested colormaps and scales will improve the consistency across teams producing similar variables. This will support comparison across different plots.

We have provided dictionaries with suggested region shapefiles, cmap colormaps, and tick intervals. Using the suggested items may help making plotting between users more consistent, but if they are not fit for your purpose, you may specify whatever you like.

### regions_dict
Region shape files are stored here: /g/data/ia39/aus-ref-clim-data-nci/shapefiles/data/

More information here: [https://aus-ref-clim-data-nci.github.io/aus-ref-clim-data-nci/datasets/shapefiles.html]

Regions include Australian: 
- Local Government Areas (LGAs),
- State and Territories,
- land boundary,
- Natural Resource Management (NRM) regions,
- river regions,
- broadacre regions, and
- National Climate Risk Assessment (NCRA) regions.

```
dict_keys(['aus_local_gov', 'aus_states_territories', 'australia', 'nrm_regions', 'river_regions', 'broadacre_regions', 'ncra_regions'])
```

### cmap_dict
These are suggested colormaps matched with possible variables to plot.  This includes color maps for total amount and for anomalies.
![colormaps_aus_maps](https://github.com/AusClimateService/plotting_maps/blob/main/colormaps_aus_maps.png)


### tick_dict
```python
# This dictionary gives some suggestions on the scale of the colour map to use for some variables. Some scales are taken from climate maps on bom.gov.au/climate
tick_dict = {"pr_annual":  [0, 50, 100, 200, 300, 400, 600, 1000, 1500, 2000, 3000, 6000],
             "pr_6mon":    [0, 50, 100, 200, 300, 400, 600,  900, 1200, 1800, 2400, 6000],
             "pr_3mon":  [0, 10,  25,  50, 100, 200, 300,  400,  600,  800, 1200, 2500],
             "pr_mon" :  [0,  1,   5,  10,  25,  50, 100,  200,  300,  400,  600, 1200],
             "pr_hour" :  [0, 1,   2,   5,  10,  15,  20,   30,   50,   75,  100,  200,],
             "pr_days": [0, 2, 3, 5, 10, 20, 30, 40, 50, 75, 100, 125, 150, 175],
             "pr_anom_mon": [-1000, -400, -200, -100, -50, -25, -10, 0, 10, 25, 50, 100, 200, 400, 1000],
             "pr_anom_3mon": [-2000, -600, -400, -200, -100, -50, -25, 0, 25, 50, 100, 200, 400, 600, 2000],
             "pr_anom_6mon": [-3000, -1200, -800, -400, -200, -100, -50, 0, 50, 100, 200, 400, 800, 1200, 3000],
             "pr_anom_ann": [-4000, -1800, -1200, -800, -400, -200, -100, 0, 100, 200, 400, 800, 1200, 1800, 4000],
             "pr_diff_mon": [-1000, -400, -200, -100, -50, -25, -10, 10, 25, 50, 100, 200, 400, 1000],
             "pr_diff_ann": [-3000, -1800, -1200, -800, -400, -200, -100, 100, 200, 400, 800, 1200, 1800, 3000],
             "frost_days": [0, 10, 20, 30, 40, 50, 75, 100, 150, 300],
             "frost_days_mon":[0, 2, 5, 10, 15, 20, 25, 31],
             "tas": np.arange(-9,52, 3),
             "tas_anom_day": np.arange(-14, 14.1, 2),
             "tas_anom_mon": np.arange(-7, 7.1, 1),
             "tas_anom_ann": np.arange(-3.5, 3.6, 0.5),
             "apparent_tas": np.arange(-6, 42, 3),
             "percent": np.arange(0,101,10),
             "xts_freq":[0.00, 0.005, 0.01, 0.02, 0.03, 0.05, 0.07, 0.10, 0.12, 0.15]
            }
```

