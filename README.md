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

 - Basic usage: Single plot of Australia eg temperature
 - Plot ocean data: Single plot of ocean data eg marine heat waves
 - Plot stations data: Single plot of station data eg coastal flooding
 - Plot multiple data types. Gridded data and station data can be plotted on the same plot: eg ocean data and station data (station and gridded data on the same plot)
 - Plot categorical data: Single plot of categorical data eg aridity
 - Plot categorical data with stippling: Single plot of hazard data with stippling eg Fire climate classes
 - Plot a selected region: Single plot of single state/region
 - Plot a selection of regions: Single plot of multiple selected regions
 - Plot multi-paneled plots for a range of future global warming levels (GWLs). For example, `plot_acs_hazard_4pp` and `plot_acs_hazard_1plus3` and both four panel plots for gwl1.2, gwl1.5, gwl2.0, and gwl3.0.  `plot_acs_hazard_1plus3` plots the first (gwl1.2) panel as the baseline and the subsequent 3 gwls as anomalies from this baseline.
 - All the above functionality is available in multi-panelled plots. Functions exist for 1, 2, 3, and 4-panelled plots in vertical or horizontal orientations. (also 2-by-2 “square” for 4pp) The hazard plotting function eg plot_acs_hazard_4pp for four-panelled-plots is constructed using helper functions

### Colours and design
Using suggested colormaps and scales will improve the consistency across teams producing similar variables. This will support comparison across different plots.

Most colours have been tested for common red-green colorblindness eg Deuteranopia. [Coblis](https://www.color-blindness.com/coblis-color-blindness-simulator/) is a handy tool to understand what your plots look like with a range of colorblind types.

Colorscales follow [IPCC design principles](https://www.ipcc.ch/site/assets/uploads/2022/09/IPCC_AR6_WGI_VisualStyleGuide_2022.pdf) and [ACS design guide (internal BOM document)](https://bom365-my.sharepoint.com/:w:/g/personal/amy_walsh_bom_gov_au/EU0i7YY8nlNHrFo3shk35nwBbl-0A4gFqG9QyxKajo2l1A). Subject matter experts gave guidance on common colourscales used in their field.
ACS has specific guidelines on figure layout and text label sizes etc.

We have provided dictionaries with suggested region shapefiles, cmap colormaps, and tick intervals. Using the recommended items may help make plotting between users more consistent, but if they are not fit for your purpose, you may specify whatever you like.

## Suggested regions, colormaps and scales for plotting

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

```




### acs_area_stats.py
This module enables calculating a range of statistics for areas defined by shapefiles, including area averages. It is best used for reducing 2D maps into a table of summary statistics for each region or state. The function can be used for more dimensions (eg, lat, lon, time, model) but may be slow and memory intensive depending on the data.
 - The function works for continuous and numerical variables
 - The function also works for calculating stats for categorical data, including calculating mode, median (if ordinal), and each category's proportions.
 - The function can calculate the states for many models individually or across the multi member ensemble.
 - The function can be used for time series extraction for regions, but it can be very memory intensive (TODO set up a workflow to cope with large data input)

### masks
Shapefiles and masks that define regions can be at /g/data/ia39/shapefiles/... . 

These shapefiles and masks can be used to outline some selected regions, calculate area statistics, or any other use you like. 

[More information on the shapefiles](https://github.com/aus-ref-clim-data-nci/shapefiles) is in the readme and example notebooks.

You may apply your own shapefiles or masks. You may need to rename some columns so that functions work as intended.

### other

See the github “issues” https://github.com/AusClimateService/plotting_maps/issues?q=is%3Aissue for some history of added functionality etc.

## Getting started:

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

Include the projects you need when you start an ARE session. Eg, storage: "gdata/ia39+gdata/hh5+gdata/mn51"

![image](https://github.com/user-attachments/assets/97b5b23d-4d21-45ab-bbc0-feeff5d74388)


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
from acs_plotting_maps import plot_acs_hazard, regions_dict, cmap_dict, tick_dict, plot_acs_hazard_3pp
import xarray as xr
```

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
                cmap=cmap_dict["pr"],
                ticks=tick_dict['pr_annual'],
                cbar_label="annual rainfall [mm]",
                cbar_extend="max",
                title = "Rainfall",
                watermark="EXPERIMENTAL IMAGE ONLY", # When you are making your final figures, remove watermark using: watermark=""
                dataset_name = ds.source_id,
                date_range=f"1 January 1960 to 31 December 1960",
                outfile = "~/figures/out.png");
```

![ann_pr_plot](https://github.com/user-attachments/assets/0791c5e8-c756-427a-9122-eb1d670e4410)

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
                    watermark="EXPERIMENTAL IMAGE ONLY", # When you are making your final figures, remove watermark using: watermark=""
                    watermark_color="k",
                    cmap = cmap_dict["tasmax"],
                    ticks = np.arange(18,53,2),)
```

![three-panel-plot](https://github.com/user-attachments/assets/9338a639-da39-4e48-ab7d-f38bc01d6cfa)

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


## Time Series extraction
For time series extraction of point locations see https://github.com/AusClimateService/TimeSeriesExtraction

