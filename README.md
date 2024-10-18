# plotting_maps
This repo has been developed to enable standardised statistics and plotting of ACS climate hazard data to support the Australian Climate Service (ACS) Hazard teams and the National Climate Risk Assessment (NCRA). We have developed Python functions and provided examples of mapping climate hazards for Australia so that data can be consistently and clearly presented.

Examples include maps and stats for Australia with a range of regions including states/territories and NCRA regions. Plotting is also possible for station data and ocean data. The functions are flexible to plot any given lat/lon but are optimised for Australia.

Intended uses include taking netcdf or xarray dataarrays of hazards and indices such as Rx1day, TXx, FFDI and plotting the data on a map of Australia. 

This work has enabled consistent mapping and summary analyses of Hazard Metrics for present and future scenarios at Global Warming Levels (GWLs). Subsequent figures and tables have been presented to national departments and ministers. The figures and tables contribute to the timely delivery of presentations and reports on Australia’s current and future climate hazards.

The code has been developed to be flexible and functional for different hazards; plotting land, ocean and station-based data; creating single- and multi-panel plots; applying different regions (eg states and NCRA regions); and masking areas (eg AGCD mask). The goal was to create code that is easy to apply, well-supported by documentation and notebook tutorials, and effective at producing aesthetic and clear plots. This has enabled collaboration between the ACS team of scientists and science communicators to deliver high-quality products to stakeholders.

Figures have been developed to align with ACS design guidelines and IPCC standards, where possible and applicable.

This work has been developed with advice from ACS, CSIRO, and BOM scientists, ACS communication teams, and stakeholders.

This repo has been developed by Gen Tolhurst (gentolhurst@gmail.com or Gen.Tolhurst@bom.gov.au) and supervised by Mitch Black (Mitchell.Black@bom.gov.au). Work has been undertaken from May 2024 to October 2024.\
Funded by ACS.

## What's possible?
### [acs_plotting_maps.py](https://github.com/AusClimateService/plotting_maps/blob/main/acs_plotting_maps.py) for plots
There's many possibilities built into this function. ```plot_acs_hazard``` is the single plot function. Multiple plots can be made in the same figure using ```plot_acs_hazard_2pp```, ```plot_acs_hazard_3pp```, ```plot_acs_hazard_4pp```, and ```plot_acs_hazard_1plus3```; these multi-panel plots have the same functionalities as the single plot function.

To access docstrings and learn about input arguments, use ```plot_acs_hazard?```. This will describe each parameter you can give to the function to customise your plot.

 - Basic usage: Single plot of Australia eg temperature  [Climate_and_hazards_report](https://github.com/AusClimateService/plotting_maps/blob/main/reports/Climate_and_hazards_report.ipynb) and [acs_plotting_maps_examples.ipynb](https://github.com/AusClimateService/plotting_maps/blob/main/example_notebooks/acs_plotting_maps_examples.ipynb)
<img src="https://github.com/AusClimateService/plotting_maps/blob/main/figures/ch_report/Annual-maximum-daily%0Amaximum-temperature.png" width="300">

 - Plot ocean data: Plots of ocean data eg marine heat waves [acs_plotting_maps_examples.ipynb](https://github.com/AusClimateService/plotting_maps/blob/main/example_notebooks/acs_plotting_maps_examples.ipynb), and  [Climate_and_hazards_report](https://github.com/AusClimateService/plotting_maps/blob/main/reports/Climate_and_hazards_report.ipynb)
<img src="https://github.com/AusClimateService/plotting_maps/blob/main/figures/Ocean-acidification.png" width="300">
 - Plot stations data: Single plot of station data eg coastal flooding [acs_plotting_maps_examples.ipynb](https://github.com/AusClimateService/plotting_maps/blob/main/example_notebooks/acs_plotting_maps_examples.ipynb) , [multi_plots](https://github.com/AusClimateService/plotting_maps/blob/main/example_notebooks/multi_plots.ipynb) and [Climate_and_hazards_report](https://github.com/AusClimateService/plotting_maps/blob/main/reports/Climate_and_hazards_report.ipynb)

<img src="https://github.com/AusClimateService/plotting_maps/blob/main/figures/ch_report/Change-in-frequency-of-flood-days.png" width="300">
  
 - Plot multiple data types in one figure. Gridded data and station data can be plotted on the same plot: eg ocean data and station data (station and gridded data on the same plot)
 - Plot categorical data: Single plot of categorical data eg [aridity](https://github.com/AusClimateService/plotting_maps/blob/main/example_notebooks/minimal_plotting_example_ai.ipynb) and [aridity or fire climate classes](https://github.com/AusClimateService/plotting_maps/blob/main/example_notebooks/plotting_and_stats_examples.ipynb)
 - Plot categorical data with stippling: Single plot of hazard data with stippling eg [multi_plots](https://github.com/AusClimateService/plotting_maps/blob/main/example_notebooks/multi_plots.ipynb), [Fire climate classes](https://github.com/AusClimateService/plotting_maps/blob/main/reports/fire_climate_classes_projections.ipynb)
![Fire-climate-classes-and-shift](https://github.com/user-attachments/assets/606762e1-d9f5-41a7-8df4-1f089f2c8596)

 - Mask remote data sparse regions of Australia using agcd_mask. Particularly for maps that depend on *in situ* rainfall observations. See issue https://github.com/AusClimateService/plotting_maps/issues/28 and [Climate_and_hazards_report](https://github.com/AusClimateService/plotting_maps/blob/main/reports/Climate_and_hazards_report.ipynb).
 - Plot a selected region: Single plot of single state/region [acs_plotting_maps_examples.ipynb](https://github.com/AusClimateService/plotting_maps/blob/main/example_notebooks/acs_plotting_maps_examples.ipynb)
 - Plot a selection of regions: Single plot of multiple selected regions [acs_plotting_maps_examples.ipynb](https://github.com/AusClimateService/plotting_maps/blob/main/example_notebooks/acs_plotting_maps_examples.ipynb)
 - Plot only data below particular latitude [Climate_and_hazards_report](https://github.com/AusClimateService/plotting_maps/blob/main/reports/Climate_and_hazards_report.ipynb)
 - Plot multi-panelled plots with shared colour bars for multiple future global warming levels (GWLs). For example, `plot_acs_hazard_4pp` and `plot_acs_hazard_1plus3` are both four-panel plots for gwl1.2, gwl1.5, gwl2.0, and gwl3.0.  `plot_acs_hazard_1plus3` plots the first (gwl1.2) panel as the baseline and the subsequent 3 gwls as anomalies from this baseline. [Multiplots examples](https://github.com/AusClimateService/plotting_maps/blob/main/example_notebooks/multi_plots.ipynb), [Climate_and_hazards_report](https://github.com/AusClimateService/plotting_maps/blob/main/reports/Climate_and_hazards_report.ipynb), [ncra_briefing_plots](https://github.com/AusClimateService/plotting_maps/blob/main/reports/ncra_briefing_plots.ipynb), [fire_climate_classes_projections](https://github.com/AusClimateService/plotting_maps/blob/main/reports/fire_climate_classes_projections.ipynb)
  ![Average-daily-maximum-temperature](https://github.com/user-attachments/assets/05e67c12-a6d5-478d-a088-7af5c6eaaadd)

 - All the above functionality is available in multi-panelled plots. Functions exist for 1, 2, 3, and 4-panelled plots in vertical or horizontal orientations eg [multi_plots](https://github.com/AusClimateService/plotting_maps/blob/main/example_notebooks/multi_plots.ipynb). (also 2-by-2 “square” for 4pp) The hazard plotting function eg plot_acs_hazard_4pp for four-panelled-plots is constructed using helper functions 

**Limitations**
-	Region shapefiles with many regions (eg LGAs) are very slow to load (big regions, like states, are ok)
-	Stippling can be weird when the mask has fuzzy edges (ie data is noisy), the stippling can get confused about what should be stippled and what shouldn’t be and may put hatches where there shouldn’t be hatches. (problem with contour). This is a problem with stippling for fire climate classes, to overcome this, I coarsened the mask to a larger grid see [multi_plots](https://github.com/AusClimateService/plotting_maps/blob/main/example_notebooks/multi_plots.ipynb), [Fire climate classes](https://github.com/AusClimateService/plotting_maps/blob/main/reports/fire_climate_classes_projections.ipynb)
-	Setting contourf=True to smoothly plot your gridded data can cause errors for particular data and particular projections. This is a known issue with contourf, be careful if you use it (check with contour=False plots). Contour and contourf are quite slow to calculate for noisy high-resolution data. (see issue https://github.com/AusClimateService/plotting_maps/issues/10)
-	Specifying tick_labels for non-categorical data produces unexpected results. The tick_labels argument is designed to label categorical data. It might be misunderstood to allow for labelling only major ticks or for labelling data with the units on the tick labels. Be aware of this. Possibly could change the functionality, if desired. (see issue https://github.com/AusClimateService/plotting_maps/issues/7)

### Colours and design
Using suggested colormaps and scales will improve the consistency across teams producing similar variables. This will support comparison across different plots.

Most colours have been tested for common red-green colorblindness eg Deuteranopia. [Coblis](https://www.color-blindness.com/coblis-color-blindness-simulator/) is a handy tool to understand what your plots look like with a range of colorblind types.

Colorscales follow [IPCC design principles](https://www.ipcc.ch/site/assets/uploads/2022/09/IPCC_AR6_WGI_VisualStyleGuide_2022.pdf) and [ACS design guide (internal BOM document)](https://bom365-my.sharepoint.com/:w:/g/personal/amy_walsh_bom_gov_au/EU0i7YY8nlNHrFo3shk35nwBbl-0A4gFqG9QyxKajo2l1A). Subject matter experts gave guidance on common colourscales used in their field.
ACS has specific guidelines on figure layout and text label sizes etc.

We have provided dictionaries with suggested region shapefiles, cmap colormaps, and tick intervals. Using the recommended items may help make plotting between users more consistent, but if they are not fit for your purpose, you may specify whatever you like.

Below are suggested colormaps matched with possible variables to plot.  This includes color maps for the total amount and anomalies. They are stored as `cmap_dict` in the `acs_plotting_maps` module.

<img src="https://github.com/AusClimateService/plotting_maps/blob/main/figures/colormaps_aus_maps.png" width="300">

### [acs_area_stats.py](https://github.com/AusClimateService/plotting_maps/blob/main/acs_area_stats.py) for area statistics
This module enables calculating a range of statistics for areas defined by shapefiles, including area averages. It is best used for reducing 2D maps into a table of summary statistics for each region or state. The function can be used for more dimensions (eg, lat, lon, time, model) but may be slow and memory intensive depending on the size of the data and the number of regions.  
 - Here's a [verbose example](https://github.com/AusClimateService/plotting_maps/blob/main/example_notebooks/area_statistics_example.ipynb) of using the function.
 - The function works for continuous and numerical variables eg [rainfall](https://github.com/AusClimateService/plotting_maps/blob/main/example_notebooks/minimal_plotting_example_pr.ipynb), [marine heatwaves](https://github.com/AusClimateService/plotting_maps/blob/main/example_notebooks/minimal_plotting_example_ocean.ipynb), [temperature](https://github.com/AusClimateService/plotting_maps/blob/main/example_notebooks/minimal_plotting_example_tx.ipynb)
 - The function also works for calculating stats for [categorical data](https://github.com/AusClimateService/plotting_maps/blob/main/example_notebooks/plotting_and_stats_examples.ipynb), including calculating mode, median (if ordinal), and each category's proportions.
 - The function can calculate the area averages for many models individually or across the multi-member ensemble. eg [ensemble-table](https://github.com/AusClimateService/plotting_maps/blob/main/example_notebooks/ensemble-table.ipynb)
 - The function can [work with any custom shapefile](https://github.com/AusClimateService/plotting_maps/blob/main/example_notebooks/area_statistics_example_basin_gpkg.ipynb)
 - The function can be used for time series extraction for regions, but it can be very memory intensive (TODO set up a workflow to cope with large data input)
 
**Limitations**
-	Stats function cannot work with NaNs
-	Region shapefiles with many regions is very slow (big regions are ok)

### Time Series extraction
For time series extraction of point locations see https://github.com/AusClimateService/TimeSeriesExtraction

### Masks
Shapefiles and masks that define regions can be at /g/data/ia39/shapefiles/data and /g/data/ia39/aus-ref-clim-data-nci/shapefiles/masks/.

These shapefiles and masks can be used to outline some selected regions, calculate area statistics, or any other use you like. 

[More information on the shapefiles](https://github.com/aus-ref-clim-data-nci/shapefiles) is in the readme and example notebooks.

You may apply your own shapefiles or masks. You may need to rename some columns so that functions work as intended.

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

### other

See the github “issues” https://github.com/AusClimateService/plotting_maps/issues?q=is%3Aissue for some history of added functionality etc.


## Getting started:

### Python environment
This code is designed to work with hh5 analysis3-24.04 virtual environment.

In your terminal, this may look like:

```
$ module use /g/data/hh5/public/modules
$ module load conda/analysis3-24.04
```

When starting a new ARE JupyterLab session (https://are.nci.org.au/pun/sys/dashboard/batch_connect/sys/jupyter/ncigadi/session_contexts/new, requires NCI login), selecting the hh5 analysis3-24.04 virtual environment might look like this:

![image](https://github.com/AusClimateService/plotting_maps/assets/45543810/e0d93235-c0a7-4a24-adb5-8bf99f3febe0)

### Access shapefiles
This code references shapefiles stored in ```/g/data/ia39/```. You will need to be a member of this project to access the data. Request membership https://my.nci.org.au/mancini/project/ia39

See https://github.com/aus-ref-clim-data-nci/shapefiles for more information on the shapefiles.

Include the projects you need when you start an ARE session. Eg, storage: "gdata/ia39+gdata/hh5+gdata/mn51"

![image](https://github.com/user-attachments/assets/97b5b23d-4d21-45ab-bbc0-feeff5d74388)


### Cloning this repo
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

### Update to the lastest version of the repo (pull)
Navigate to your existing version of the plotting maps repository (if you don't have an existing version, follow the above directions for cloning).

```
$ cd /path/to/dir/plotting_maps
````

Then pull the latest version using git

```
$ git pull
```

### Usage in Jupyter Notebook:

See small, easy-to-follow examples here: 
- [https://github.com/AusClimateService/plotting_maps/blob/main/example_notebooks/minimal_plotting_example_pr.ipynb]
- [https://github.com/AusClimateService/plotting_maps/blob/main/example_notebooks/area_statistics_example.ipynb]

Other examples:
- [https://github.com/AusClimateService/plotting_maps/blob/main/example_notebooks/plotting_and_stats_examples.ipynb]
- [https://github.com/AusClimateService/plotting_maps/blob/main/example_notebooks/acs_plotting_maps_examples.ipynb]

1. **Navigate to the directory you cloned to:**
```
cd ~/plotting_maps
```

2. **Import the ACS plotting maps function and dictionaries and Xarray.** 
```python 
from acs_plotting_maps import plot_acs_hazard, regions_dict, cmap_dict, tick_dict, plot_acs_hazard_3pp
import xarray as xr
```

3. **Load some data.** For example, this will load extratropical storm rx5day rainfall
```python
ds = xr.open_dataset("/g/data/ia39/australian-climate-service/test-data/CORDEX-CMIP6/bias-adjustment-input/AGCD-05i/BOM/ACCESS-CM2/historical/r4i1p1f1/BARPA-R/v1-r1/day/pr/pr_AGCD-05i_ACCESS-CM2_historical_r4i1p1f1_BOM_BARPA-R_v1-r1_day_19600101-19601231.nc")
```
This data has three dimensions (time, lon, lat). There is a value for every day from 01-01-1960 to 31-12-1960. We can only plot 2D, so next we will calculate a statistic to summarise the data

4. **Summarise data into a 2D xr.DataArray.** For example, calculate the annual sum:
```python
var="pr"
da = ds.sum(dim="time")[var]
```

5. **Finally, use the plotting function**.\
You will need to specify:
     * the data (and select the variable eg "pr");
     * suitable arguments for the colorbar including cmap, ticks, cbar_label, and cbar_extend;
     * annotations including title, dataset_name, date_range; and
     * where you want the image outfile saved.
   
```python
regions = regions_dict['ncra_regions']
plot_acs_hazard(data = da,
                regions = regions,
                ticks=tick_dict['pr_annual'],
                cbar_label="annual rainfall [mm]",
                cbar_extend="max",
                title = "Rainfall",
                dataset_name = ds_pr.source_id,
                date_range=f"{start} to {end}",
                agcd_mask=True,
                cmap_bad="lightgrey",
                watermark="",
                outfile = "~/figures/out.png");
```
![rainfall_plot](https://github.com/user-attachments/assets/112b3911-1807-4d70-b035-acad148eb96c)

**Plot a three-panel plot**
```python
%%time
from plotting_maps.acs_plotting_maps import plot_acs_hazard_3pp

var = "HWAtx"

ds_gwl12 =xr.open_dataset("/g/data/ia39/ncra/heat/data/HWAtx/bias-corrected/ensemble/GWL-average/HWAtx_AGCD-05i_MME50_ssp370_v1-r1-ACS-QME-AGCD-1960-2022_GWL12.nc")
ds_gwl15 = xr.open_dataset("/g/data/ia39/ncra/heat/data/HWAtx/bias-corrected/ensemble/GWL-average/HWAtx_AGCD-05i_MME50_ssp370_v1-r1-ACS-QME-AGCD-1960-2022_GWL15.nc")
ds_gwl20 = xr.open_dataset("/g/data/ia39/ncra/heat/data/HWAtx/bias-corrected/ensemble/GWL-average/HWAtx_AGCD-05i_MME50_ssp370_v1-r1-ACS-QME-AGCD-1960-2022_GWL20.nc")
ds_gwl30 = xr.open_dataset("/g/data/ia39/ncra/heat/data/HWAtx/bias-corrected/ensemble/GWL-average/HWAtx_AGCD-05i_MME50_ssp370_v1-r1-ACS-QME-AGCD-1960-2022_GWL30.nc")

plot_acs_hazard_3pp(ds_gwl15 = ds_gwl15[var], 
                    ds_gwl20 = ds_gwl20[var],
                    ds_gwl30 = ds_gwl30[var],
                    regions = regions_dict['ncra_regions'],
                    cbar_label=f"Temperature [degC]",
                    title=f"Maximum Temperature of Hottest Heatwave for future warming scenarios", 
                    date_range = "Insert subtitle - should include the date range of the data \nand then the dataset below that", 
                    # baseline = "GWL1.2", 
                    dataset_name= "MME50_ssp370",
                    issued_date=None,
                    watermark="EXPERIMENTAL IMAGE ONLY", 
                    watermark_color="k",
                    cmap = cmap_dict["tasmax"],
                    ticks = np.arange(18,53,2),)
```
![Maximum-Temperature-of-Hottest-Heatwave-for-future-warming-scenarios](https://github.com/user-attachments/assets/6741a325-aac5-44bd-b878-d1aab1dab2ba)


6. **Calculate summary statistics for the range of models.**

```python 
# Import needed packages
from acs_area_statistics import acs_regional_stats, get_regions
regions = get_regions(["ncra_regions", "australia"])
```

**this has changed. Previously**

"# import needed packages

from acs_area_statistics import acs_regional_stats, regions"

For Calculating the NCRA region stats, we want to compare the regional averages based on different models, eg what is the regional mean value from the coolest/driest model realisation, what is the mean, what is the regional mean from the hottest/wettest model for this, we want ds to have the 10th, median and 90th percentile values from each model, then we can find the range of the models and the MMM.

```python
# calculate the stats using the acs_region_fractional_stats function
# Find the min, mean, max value for each region

ds = xr.open_dataset(filename)
mask_frac = regions.mask_3D_frac_approx(ds)
dims = ("lat", "lon",)
how = ["min", "mean", "max"]

da_summary = acs_regional_stats(ds=ds, infile = filename, mask=mask_frac, dims = dims, how = how,)
da_summary.to_DateFrame()

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

## FAQs
### Something is not working and I don't know why!
Here are some common suggestions for troubleshooting: 
 -	see “getting started” above and make sure you have followed all the instructions
 -	Check you are using the right venv. This code is designed to work with hh5 analysis3-24.04 virtual environment.
 -	Restart the kernel and rerun all cells from start. Especially if you have made a variety of modifications, you may have renamed a function/variable.
 -	If python can't find the module, check you have the .py module in your working directory. If not cd to the directory with the module.
 -	Make sure you have requested access to all the right gdata projects (eg gdata/ia39)

### An argument I have used before using this code no longer works. What's happening?
During development, priorities and requests have changed what the functions needed to do. As a result, there are a few deprecated features and functionalities. Some things that were needed that are now not required:
 -	 “show_logo”, it was initially requested to have an ACS logo in the figures. The comms team now prefers only the copywrite in the bottom
 -	Contour and contourf are generally not recommended now due to errors in plotting and long computational time. They are left in the function because they can be useful for lower resolution data, eg ocean data.
 -	“infile” is not used. The idea was to use this for well-organised data with a consistent DRS to enable a good plot to be made without lots of keyword inputs. The data we have is not organised consistently enough for this.
 -	“regions_dict” in acs_plotting_maps.py made to module very slow to load Shapefiles can take many seconds to load. It is inefficient to load all these regions even when you don’t use them all. This was replaced with a class
 -	“regions” in acs_area_stats had preloaded shapefiles. Shapefiles can take many seconds to load. It is inefficient to load all these regions even when you don’t use them all. This was replaced with “get_regions”
```python
 from acs_area_statistics import acs_regional_stats, get_regions
regions = get_regions(["ncra_regions", "australia"])
```

### How can I add stippling (hatching) to plots to indicate model agreement?
The plotting scripts can add stippling to the plots using the stippling keyword(s).
You will need to calculate the mask and provide this as a dataarray with "lat" and "lon". The mask must be a True/False boolean mask. It does not have to be the same resolution as the underlying data (you may wish to coarsen the mask if the underlying data is high-resolution and noisy).
See [this link](https://github.com/AusClimateService/plotting_maps/issues/2) for a brief example of applying stippling.

For the multi-panel plots, you can give a mask for each of the plots eg see [fire_climate_classes_projections.ipynb](https://github.com/AusClimateService/plotting_maps/blob/main/reports/fire_climate_classes_projections.ipynb) (you may ignore the "coarsen..." this is needed so smooth out the fuzzy edges of the fire climate classes). In this

Your function will look something like this:
```python
plot_acs_hazard_4pp(ds_gwl12=ds_gwl12[var],
                    ds_gwl15=ds_gwl15[var],
                    ds_gwl20=ds_gwl20[var],
                    ds_gwl30=ds_gwl30[var],
                    stippling_gwl12=stippling_gwl12,
                    stippling_gwl15=stippling_gwl15,
                    stippling_gwl20=stippling_gwl20,
                    stippling_gwl30=stippling_gwl30,
                    regions = regions,
                    title = "Fire Climate Classes",
                    # figsize=(7,2),
                    # baseline="GWL1.2",
                    cmap = cmap_dict["fire_climate"],
                    ticks = tick_dict["fire_climate_ticks"],
                    tick_labels = ["Tropical\nSavanna","Arid grass \n& woodland","Wet Forest","Dry Forest","Grassland",],
                    cbar_label = "classes",
                    dataset_name = "BARPA MRNBC-ACGD",
                    watermark="",
                    orientation="horizontal",
                    issued_date="",
                    );
```

### Why is the stippling weird?
You may need to check that the stippling is in the areas you expect it to be. There is a bug in contourf that causes the stippling to get confused when plotting noisy high-resolution mask. If that is the case, I recommend coarsening the stippling mask 
E.g. 
new_stippling_mask =  stippling_mask.coarsen(lat=2, boundary="pad").mean().coarsen(lon=2, boundary="pad").mean()>0.4

(full example here https://github.com/AusClimateService/plotting_maps/blob/main/reports/fire_climate_classes_projections.ipynb)


### Is there a way to use the 4pp plot with the average conditions for GWL1.2 and the change % for GWL1.5 to GWL3? Or does it only work for plots that use a consistent colourbar?
plot_acs_hazard_1plus3

### How can I change the orientation (eg from vertical to horizontal) of the figures in a multipaneled plot?
orientation

### Can I use my own shapefiles to define regions?
Yes, you can provide any shapefiles you like. We have provided some helpful Australian regions from /g/data/ia39, but the functions are flexible to take custom regions. 

### I want to use a divergent colormap, but the neutral color isn't in the middle of my ticks. What can I do to align the centre of the colormap to zero?
vcentre

### What does gwl mean?
GWL describe global warming levels. These are 20 year periods centred on the year when a climate model is projected to reach a specified global surface temperature above the pre-industrial era. Global climate models reach these temperature thresholds at different years.

For example, the Paris Agreement (2012) refers to global warming levels in its aims:

“…to strengthen the global response to the threat of climate change by keeping a global temperature rise this century well below 2 degrees Celsius above pre-industrial levels and to pursue efforts to limit the temperature increase even further to 1.5 degrees Celsius.”

Find more information here https://github.com/AusClimateService/gwls

The plotting functions have been designed to accommodate present and future global warming levels. This is indicated by argument names containing "gwl12", "gwl15", "gwl20", "gwl30". If you want to use the function for other time periods or scenarios, you can still use these functions. The functions will work for any data in the right format (eg 2D xarray data array with lat and lon).

### I am not using GWLs but I want to use these functions. How can I can the subtitles?


### I only want to plot data below 30S latitude, is there a mask for this?
There is no specific mask for this, but it is easy to adjust your input to achieve this.

If you just want to plot the data below 30S, you can use ```plot_acs_hazard(data=  ds.where(ds["lat"]<-30)[var] , ...)```
![image](https://github.com/user-attachments/assets/549ccef8-3aed-4dfa-9ee0-c08e225ab386)


You may also like to apply a custom mask to the stats function using "clipped" to only select by a lat lon box:
```
import geopandas as gpd
from glob import glob
from shapely.geometry import box
import regionmask
import xarray as xr
from acs_area_statistics import acs_regional_stats, get_regions

# get the shapefile for australia
PATH = "/g/data/ia39/aus-ref-clim-data-nci/shapefiles/data"
shapefile = "australia"
gdf = gpd.read_file(glob(f"{PATH}/{shapefile}/*.shp")[0]).to_crs("EPSG:4326")

# set your limits
# box(xmin, ymin, xmax, ymax)
clipped = gdf.clip( box(100, -45, 160, -30))

regions = regionmask.from_geopandas(clipped, name= "clipped_shapefile", overlap=True) 

# need some data
filename = "/g/data/ia39/ncra/extratropical_storms/5km/GWLs/lows_AGCD-05i_ACCESS-CM2_ssp370_r4i1p1f1_BOM_BARPA-R_v1-r1_GWL12.nc"
ds = xr.open_dataset(filename, use_cftime = True,)

mask = regions.mask_3D(ds)

# then calculate the stats for this clipped region
dims = ("lat", "lon",)
var="low_freq"
df_summary = acs_regional_stats(ds=ds,var=var, mask=mask, dims = dims, how = ["min", "median", "max"])
df_summary
```
### Stats FAQs...?
tba

## Development principles
This code has been developed to make consistent plotting and statistical analysis quick and easy across ACS hazard teams. These teams regularly get information requests with quick turn around times (~days), so having easy to use and flexible code to produce report ready plots is critical for devilering high quality figures and data summaries.

We want to enable scientists to focus on producing data and interpreting that information.

These plotting and stats tools should reduce the duplication of work for scientists who spend time adjusting their plots. We want these tools to be easy to use and flexible for all purposes across the hazard teams so all the plots are presented consistently and all teams are supported.

Using these functions should be the easiest way to produce a nice plot. We want this to be the path of least resistance. The easier the functions are to use, the more people will use them. This means supporting with good documentation, good example notebooks, and adding new functionalities according to user feedback.

These plotting and stats tools are optimised for Australian hazard maps, but flexible so that they can be adjusted for any data, anywhere in the world, with any shape files, any color schemes, or any projection etc.

To test updates in [acs_area_stats.py](https://github.com/AusClimateService/plotting_maps/blob/main/acs_area_stats.py)  or [acs_plotting_maps.py](https://github.com/AusClimateService/plotting_maps/blob/main/acs_plotting_maps.py), we rerun the example notebooks to ensure the functions still perform as expected. (This could be automated to run when changes are pushed to git.)

This code is in active use by a range of teams. Take care to maintain backward compatibility while adding features. If that is not practical, then, comunicate the changes with the users. Ideally, I would like to “release” this version of the code. Eg see https://github.com/AusClimateService/plotting_maps/releases/new

## TODO
**Figures to make:**
-	SLR observations with station data. For the Climate Hazards report, recreate the sea level observations of gridded ocean data and station data in the same plot.
-	Lightning plot. For the Climate Hazards report, recreate the lightning observations plot using the plot_acs_hazards function so that it is in the consistent format.

**Documentation:**
-	Clear tutorials for how to use different functionalities
-	“release” this version of the code as v1.0.0. It is in wide use and care is being taken to keep backward compatibility while adding features. Eg see https://github.com/AusClimateService/plotting_maps/releases/new

**Improve plotting function and axillaries:**
-	Improve the aesthetics and proportions of plotting, especially with dataset/date_range/baseline annotations. Design aesthetics where focused on vertical orientations for 4 panel plots without these annotations for a particular report.
-	Improve the aesthetics of plotting select_area. Eg remove boundaries of neighbouring regions (if desired)
-	Try plotting Antarctica and see what settings work. Especially xlim, ylim, projection. I think this code should work for any region of the world, with the right settings.
-	Forest mask for forested areas. For example, FFDI is not useful in places where there is not connected vegetation/fuel. This is probably particularly for arid desert areas of central Australia. Changes in climate and land use may cause changes over time.
-	Improve colormap for fire climate classes. This colour scheme is not completely colourblind friendly. Perhaps modify the colours to increase the contrast. 
-	Enable rotating all labels and tick labels so that they are horizontal (easier to read). May need to reduce the labels to every second increment. Eg for temperature.
-	Create dictionaries for each hazard to enable automation of figures. Eg, use one keyword to select titles, colormaps and ticks.
-	Possibly automate the scaling of the colourbar to the data limits of the plot. (I am personally against this idea. Let's come up with standard colormaps and colourscales so that all figures of that one variable or hazard have a standard and comparable scale.)\
-	Possibly automate the arrows of the colourbar. (I don’t think the arrows on the colorbar should be determined by the data in the plot, I think they should be only limited by possible physical values of that metric so that all colourbars of that metric are comparable. Determine if you want the arrows to be determined by the plotting data or the metric’s possible physical values.)
-	If hazard data had consistent file naming practices (DRS) and consistent attribute labels, then the plotting functions could be further automated. At the moment, the data files are named in different patterns, the files might have different names for coordinates (eg “time”, “lat”, “lon”)
-	Use a keyword to make plots appropriate for different uses eg journal, report, powerpoint, poster etc similar to https://seaborn.pydata.org/generated/seaborn.set_context.html

**New plotting function:**
-	Fully flexible custom n x m grid of plots. At the moment, minor modifications within multiplot are needed to make a custom plot for new layouts. It may be possible to make a function that can take in dimensions and a list of dataarrays to make a figure of many plots. This should use a similar format to the existing multi-panel plots and allow plotting gridded data, station data, stippling, ocean data, etc.

**Stats functions:**
-	Optimise workflow to enable area-averaged time series (stats or just area mean). This function can be very memory intensive. Need to apply a strategy or strategies to reduce memory use. A possible option may be to calculate and save area averages for every year. Saving outputs in annual files is a common practice for climate models.
-	Calculate statistics along streamlines. Similar to area averages, but for a custom transect. Eg for rivers instead of catchments. Eg issue https://github.com/AusClimateService/plotting_maps/issues/23


