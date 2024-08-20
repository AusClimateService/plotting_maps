# plotting_maps
This repo aims to plot standardised maps of ACS climate hazard data. We will develop and provide examples of mapping climate hazards for Australia so that data can be consistently presented.

Examples will include maps for Australia and selected states or regions. 

Intended uses include taking netcdf or xarray dataarrays of hazards and indices such as Rx1day, TXx, FFDI and plotting the data on a map of Australia. 

Developed by Gen Tolhurst (gentolhurst@gmail.com), Supervised by Mitchell Black (Mitchell.Black@bom.gov.au).\
Funded by ACS

## Python environment
This code is designed to work with hh5 analysis3-24.04 virtual environment.

In your terminal, this may look like:

```
$ module use /g/data/hh5/public/modules
$ module load conda/analysis3-24.04
```

When starting a new ARE JupyterLab session (https://are.nci.org.au/pun/sys/dashboard/batch_connect/sys/jupyter/ncigadi/session_contexts/new, requires NCI login), selecting the hh5 analysis3-24.04 virtual environment might look like this:

![image](https://github.com/AusClimateService/plotting_maps/assets/45543810/e0d93235-c0a7-4a24-adb5-8bf99f3febe0)

## Access shapefiles
This code references shapefiles stored in ```/g/data/ia39/```. You will need to be a member of this project to access the data. Request membership https://my.nci.org.au/mancini/project/ia39

See https://github.com/aus-ref-clim-data-nci/shapefiles for more information on the shapefiles.

## Cloning this repo
Before you can ```import acs_plotting_maps``` to use the plotting function ```plot_acs_hazard```, you will need to clone a copy of this repository to your own working directory.

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

Then, you can clone this repository to access the Python code and notebooks. \
If you want the new directory to be called anything other than "plotting_maps" please replace the final argument with your choice of directory name:
```
$ git clone https://github.com/AusClimateService/plotting_maps.git plotting_maps
```
You will now be able to access the functions, python scripts, and Jupyter notebooks from your user.

## Usage in Jupyter Notebook:

See a small, easy-to-follow example here: 
- [https://github.com/AusClimateService/plotting_maps/blob/main/minimal_plotting_example_pr.ipynb]
- [https://github.com/AusClimateService/plotting_maps/blob/main/area_statistics_example.ipynb]

Other examples:
- [https://github.com/AusClimateService/plotting_maps/blob/main/plotting_and_stats_examples.ipynb]
- [https://github.com/AusClimateService/plotting_maps/blob/main/acs_plotting_maps_examples.ipynb]

1. **Navigate to the directory you cloned to:**
```
cd ~/plotting_maps
```

2. **Import the ACS plotting maps function and dictionaries and Xarray.** 
```python 
from acs_plotting_maps import plot_acs_hazard, regions_dict, cmap_dict, tick_dict
import xarray as xr
from acs_area_statistics import get_regions # this line has been updated 19 August 2024
regions = get_regions(["ncra_regions", "australia"]) # this line has been updated 19 August 2024
```

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
plot_acs_hazard(data = da["pr"],
                regions = regions, # this line has been updated 19 August 2024
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

6. **Calculate summary statistics for the range of models.**

```python 
# Import needed packages
from acs_area_statistics import acs_regional_stats, get_regions
regions = get_regions(["ncra_regions", "australia"])
```

**this has changed. Previously**

"# import needed packages

from acs_area_statistics import acs_regional_stats, regions"

For Calculating the NCRA region stats, we want to compare the regional averages based on different models, eg what is the regional mean value from the coolest/driest model relisation, what is the mean, what is the regional mean from the hottest/wettest model for this, we want ds to have the 10th, median and 90th percentile values from each model, then we can find the range of the models and the MMM.

```python
# calculate the stats using the acs_region_fractional_stats function
# Find the min, mean, max value for each region

ds = xr.open_dataset(filename)
mask_frac = regions.mask_3D_frac_approx(ds)
dims = ("lat", "lon",)
how = ["min", "mean", "max"]

da_summary = acs_regional_stats(ds=ds, infile = filename, mask=mask_frac, dims = dims, how = how,)
da_summary

```

The dataframe will be saved to: ```infile.replace(".nc", f"_summary-{'-'.join(how)}_ncra-regions.csv"```

For example only, this would make a dataframe in this format:

|   region | abbrevs   | names                   |   pr_min |   pr_mean |   pr_max |
|---------:|:----------|:------------------------|---------:|----------:|---------:|
|        0 | VIC       | Victoria                |  415.729 |   909.313 |  3005.45 |
|        1 | NT        | Northern Territory      |  397.385 |   941.405 |  3934.81 |
|        2 | TAS       | Tasmania                |  555.644 |  1760.66  |  4631.81 |
|        3 | SA        | South Australia         |  284.455 |   575.952 |  1413.98 |
|        4 | NSW       | New South Wales & ACT   |  294.329 |   768.1   |  3440.04 |
|        5 | WAN       | Western Australia North |  123.651 |   921.906 |  3470.24 |
|        6 | WAS       | Western Australia South |  249.566 |   545.317 |  1819.89 |
|        7 | SQ        | Queensland South        |  287.613 |   584.155 |  1654.74 |
|        8 | NQ        | Queensland North        |  264.447 |   766.444 |  7146.55 |
|        9 | AUS       | Australia               |  123.614 |   742.735 |  7146.55 |



## Suggested regions, colormaps and scales for plotting
Using suggested colormaps and scales will improve the consistency across teams producing similar variables. This will support comparison across different plots.

We have provided dictionaries with suggested region shapefiles, cmap colormaps, and tick intervals. Using the recommended items may help make plotting between users more consistent, but if they are not fit for your purpose, you may specify whatever you like.

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

