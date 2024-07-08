#!/usr/bin/env python
# coding: utf-8

# # Plotting maps example
# 
# This notebook gives a minimal example of using acs_plotting_maps.py
# 
# Code is available publically here [https://github.com/AusClimateService/plotting_maps/blob/main/acs_plotting_maps.py]
# 
# The README [https://github.com/AusClimateService/plotting_maps/tree/main] gives cloning instructions and virtual environment requirements to ensure it runs smoothly.
# 
# For example, in your working directory (eg navigate to your home, scratch, or user directory in a project using cd), clone this repository to access this code 
# ```
# $ git clone https://github.com/AusClimateService/plotting_maps.git plotting_maps
# ```
# 
# This code is designed to work with hh5 analysis3-24.04 virtual environment. Eg:
# ```
# $ module use /g/data/hh5/public/modules
# $ module load conda/analysis3-24.04
# ```
# 

# # Step 1 - access plotting package
# Navigate to the directory that you have cloned the plotting_maps repo to. eg
# ```
# cd ~/plotting_maps
# ```

# In[1]:


cd /g/data/mn51/users/gt3409/plotting_maps/


# Then import the plotting function ```plot_acs_hazard``` and helpful dictionaries ```regions_dict, cmap_dict, tick_dict```

# In[2]:


# import ACS plotting maps and Xarray.
from acs_plotting_maps import plot_acs_hazard, regions_dict, cmap_dict, tick_dict
import xarray as xr


# In[3]:


# If you want to define your own color maps or ticks,
# you may want to use these packages
import matplotlib.colors as colors
import numpy as np


# # Step 2 - Load and prepare hazard data
# Use xarray to load hazard data.\
# If this data is not a two-dimensional array, perform your desired selection or calculation (eg mean, min, max, percentile) to reduce data to 2D.

# In[4]:


# load some dataset
ds = xr.open_dataset("/g/data/ia39/ncra/bushfire/fire_climate_classes_AGCD-05i_MM_ssp370_v1-r1-ACS-NRNBC_GWL12.nc")
# ds = ds.rename({"latitude":"lat", "longitude": "lon"})
ds


# # Step 3 - Plot
# Use ```plot_acs_hazard``` to visualise the hazard on a map of Australia.
# 
# There are quite a few options to modify this plotting. At a minimum, you will need:
#  - **data**, a 2D xarray.DataArray of your hazard
#  - **regions**, use the regions_dict to access region or state boundary shape data
#  - **title**, title of plot naming the index or hazard you are plotting
#  - **date_range**, date range of the data you have plotted, appears as a subtitle under the title
#  - **cmap**, use the cmap_dict to access a range of recommended colormaps
#  - **ticks**, use the tick_dict to access a range of useful ranges or input your own list or array
#  - **cbar_label**, is the label for the colorbar. Give name and unit
#  - **cbar_extend**, controls the arrows of the colorbar. Indicates that values beyond the colorbar are possible. Use 'neither' for finite ranges eg deciles. Use 'both' for anomalies or temperatures. Use 'max' for total rainfall, where negative values are not possible, but very large positive values are possible. "min" is also an option.
#  - **dataset_name**, name of the data source eg "AGCD v2", "BARPA-R ACCESS-CM2"
#  - (**baseline**, If plotting anomalies, give the base period as a string, eg "1961-1990")

# In[5]:


# example of categorical data
plot_acs_hazard(data = ds.fire_climate_class,
                regions = regions_dict['ncra_regions'],
                title = "Fire Climate Classes",
                date_range = "GWL1.2",
                cmap = colors.ListedColormap(["#84a19b", "#e0d7c6", "#486136", "#737932", "#a18a6e",]),
                ticks = [100, 101, 102, 103, 104,],
                tick_labels = ["Tropical Savanna", "Arid grass \nand woodland", "Wet Forest", "Dry Forest", "Grassland",],
                cbar_label = "classes",
                dataset_name = "MRNBC bias corrected, multi-model",
                outfile = "figures/out.png",
                );


# In[6]:


# example of two level categorical data
ds = xr.open_dataset("/g/data/ia39/ncra/bushfire/fire_climate_classes_shift_AGCD-05i_MM_ssp370_v1-r1-ACS-NRNBC_GWL12_to_GWL15.nc")
da = ds.fire_climate_class
plot_acs_hazard(data = da,
                regions = regions_dict['ncra_regions'],
                title = "Fire Climate Classes",
                date_range = "GWL1.2 to GWL1.5",
                cmap = colors.ListedColormap(["palegreen", "grey",]),
                ticks = [0, 1,],
                tick_labels = ["no shift", "shifted",],
                cbar_label = "classes",
                dataset_name = "MRNBC bias corrected, multi-model",
                outfile = "figures/out.png",
                contour=False,
                );


# In[7]:


# Categorical example with non-linear index values and NaN values over land

ds_ai = xr.open_dataset("/g/data/ia39/ncra/drought_aridity/ai/AI-atmospheric-based_NHP1-AUS-5_rcp85_bias-adjusted_2D_GWL20_percentiles10-50-90.nc")

da = ds_ai.sel(quantile = 0.5).AI

plot_acs_hazard(data =  da.where(da<0.65),
                regions = regions_dict['ncra_regions'],
                title = "Aridity Index - 50th percentile / median",
                date_range = "GWL2.0",
                cmap = cmap_dict["aridity"],
                ticks = tick_dict['aridity_index_ticks'],
                tick_labels = tick_dict['aridity_index_labels'],
                cbar_label = "Aridity Index",
                dataset_name = "NHP1-AUS-5_rcp85_bias-adjusted",
                outfile = "figures/out.png",
                );


# In[ ]:





# In[8]:


filename = "/g/data/ia39/ncra/extratropical_storms/5km/GWLs/lows_AGCD-05i_ACCESS-CM2_ssp370_r4i1p1f1_BOM_BARPA-R_v1-r1_GWL12.nc"
ds_pr = xr.open_dataset(filename)
ds_pr


# In[9]:


plot_acs_hazard(data = ds_pr.low_freq,
                regions = regions_dict['ncra_regions'],
                title = "Lows",
                date_range = "GWL1.2",
                cmap = cmap_dict["xts_freq"],
                ticks = tick_dict['xts_freq'],
                cbar_label = "low freq",
                dataset_name = "",
                outfile = "figures/out.png",
                );


# ## Step 4: Calculate NCRA region statistics

# In[10]:


# import needed packages
from acs_area_statistics import acs_regional_stats, regions


# In[11]:


acs_regional_stats(ds=ds_pr, 
                   infile = filename.split("/")[-1], 
                   mask = "fractional", 
                   how = ["mean", "min", "max"],)


# In[12]:


path = "/g/data/ia39/ncra/drought_aridity/ai/"
filename = "AI-atmospheric-based_NHP1-AUS-5_rcp85_bias-adjusted_2D_GWL20_percentiles10-50-90.nc"
ds_ai = xr.open_dataset(path+filename)

da_summary = acs_regional_stats(ds=ds_ai.sel(quantile = 0.5), infile = filename, mask = "fractional", dims=("lat", "lon"), how = ["mean", "min", "max"],)
da_summary


# In[13]:


ds_ai = xr.open_dataset("/g/data/ia39/ncra/drought_aridity/ai/AI-atmospheric-based_NHP1-AUS-5_rcp85_bias-adjusted_2D_GWL20_percentiles10-50-90.nc")

da_summary = acs_regional_stats(ds=ds_ai.sel(quantile = 0.5), mask = "fractional", dims=("lat", "lon"), how = ["mean", "min", "max"], outfile="out.csv")
da_summary


# In[14]:


get_ipython().run_cell_magic('time', '', 'mask = regions.mask_3D_frac_approx(ds)\n')


# In[15]:


da_summary = acs_regional_stats(ds=ds, var = "fire_climate_class",  mask = "fractional", dims=("lat", "lon"), how = ["mean", "min", "max", "mode"], outfile="out.csv")
da_summary


# # Access the docstring for more info

# In[16]:


get_ipython().run_line_magic('pinfo', 'plot_acs_hazard')


# In[17]:


get_ipython().run_line_magic('pinfo', 'acs_regional_stats')


# In[ ]:




