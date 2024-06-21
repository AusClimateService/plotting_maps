# # Cookie cutter

# Use this code to produce summary statistics for Hazards over NCRA regions, with the flexibility to apply the method to any shapefile region.
# Typical statistics include median, mean, min, max, 10th, 90th percentiles

# This method has used guidance from [https://github.com/aus-ref-clim-data-nci/shapefiles/blob/master/python_tutorial.ipynb]

import xarray as xr
import geopandas as gpd
import regionmask
import numpy as np
import pandas as pd

# read in the shapefile with regions you will use Australia and NCRA regions
# from acs_plotting_maps import regions_dict
ncra_gdf = gpd.read_file(f'/g/data/ia39/aus-ref-clim-data-nci/shapefiles/data/ncra_regions/ncra_regions.shp')
ncra_gdf["abbrevs"]=['VIC', 'NT','TAS', 'SA', 'NSW', 'WAN', 'WAS', 'SQ', 'NQ']

aus_gdf = gpd.read_file(f'/g/data/ia39/aus-ref-clim-data-nci/shapefiles/data/australia/australia.shp')
aus_gdf["abbrevs"]=['AUS']
aus_gdf["NAME"]=['Australia']

gdf = pd.concat([ncra_gdf[["NAME", "abbrevs", "geometry"]], aus_gdf[["NAME", "abbrevs", "geometry"]]],)
gdf.index = np.arange(0,10)
regions = regionmask.from_geopandas(gdf, names="NAME", abbrevs="abbrevs")
regions

def acs_regional_stats(ds,
                        var = None, 
                        mask = None, 
                        start = None, 
                        end = None, 
                        dims = ("time", "lat", "lon"), 
                        how = "mean",
                        quantile = None, 
                        select_abbr = None,
                        select_name = None,
                        overlap_threshold = None):
    """
    This function takes an Xarray dataset (ds) with variable (var) and multiple dimensions (eg time, lat, and lon), 
    then selects the time range between two years (start and end), 
    and applies regions.mask_3D_frac_approx fractional mask (frac) to compute a regional statistic (how, eg "mean") over two or three dimensions.

    Parameters
    ----------
    ds: xr.Dataset or xr.DataArray
        expects an xr.Dataset with variable var and dimensions time, lat, and lon.
        
    var: str
        name of variable in ds, eg "pr" or "tas".

    mask:  xarray.DataArray 'mask' or ["fractional", "centred", "min_overlap"] 
        expects a precalculated mask, or will calculate one of the three options. If "min_overlap" selected, you must provide an overlap_threshold value.
        fractional mask from regionmask.from_geopandas(ncra_gdf, names="NAME", abbrevs="abbrevs").mask_3D_frac_approx(ds). 
        If None is provided, then will calculate a fractional mask based on NCRA regions, but this will take about one minute to calculate. 
        This is best avoided by calculating frac outside the function.
        
    start: string or int
        start year to slice data array. eg start year of global warming level GWL. If either start or end is None, then time selection is not performed. 
        
    end: string or int
        end year to slice data array. eg end year of global warming level GWL. If either start or end is None, then time selection is not performed. 
        
    dims: tuple of dimensions
        Dimensions to reduce data. Suggestion to reduce to one value per region: dims= ("time", "lat", "lon"). Or to reduce to 1D time series per region ("lat", "lon").
        
    how: str
        statistic to reduce data. One of ['mean', 'median', 'min', 'max', 'sum', 'std', 'var']. Default "mean".
    
    quantile: float [0.0,1.], optional
        Percentile to calculate. Will overwrite any "how" method

    select_abbr: list of str
        list of regions by abbreviation to perform statistic on. eg ["VIC", "NSW"]

    select_name: list of str
        list of regions by name to perform statistic on.  eg ["Victoria", "New South Wales & ACT"]

    overlap_threshold: float between 0.0 and 1.
        If mask =  "min_overlap", you must provide an overlap_threshold value.
        If mask != "min_overlap", this value will be ignored.
    
    Returns
    -------
    An xarray.DataArray (1D or 2D) with a dimension called "region" and values of the regional statistics.
    If no valid statistic is called through how or quantile, then the DataArrayWeighted is returned


    Example usage
    --------------
    # open data 
    ds = xr.open_dataset(filename, use_cftime = True,)
    # define mask
    mask_frac = regions.mask_3D_frac_approx(ds)

    # apply function
    da_mean = acs_regional_stats(ds=ds, var="pr", mask =mask_frac, start="1991", end="2010", dims = ("time", "lat", "lon"), how = "mean")
    """
    if isinstance(mask, xr.DataArray):
        # Prefered method: use precalculated xarray.DataArray 'mask'
        mask = mask
    elif mask == "fractional":
        # !warning very slow!
        print("!warning very slow! Calculating fractional mask every time is very slow. \
        \nPlease calculate ```mask = regions.mask_3D_frac_approx(ds)``` before this function.")
        mask = regions.mask_3D_frac_approx(ds)
    elif mask == "centred":
        # !warning slow!
        print("!warning slow! Calculating mask every time is slow. \
        \nPlease calculate ```mask = regions.mask_3D(ds)``` before this function.")
        mask = regions.mask_3D(ds)
    elif mask == "min_overlap":
        # !warning very slow!
        print("!warning very slow! Calculating fractional mask with minimum overlap every time is very slow. \
        \nPlease calculate ```mask = regions.mask_3D_frac_approx(ds) >= overlap_threshold``` before this function.")
        assert (0. <= overlap_threshold <= 1.), "You have selected min_overlap mask. Please specify overlap_threshold between [0.,1.]"
        mask = regions.mask_3D_frac_approx(ds) >= overlap_threshold
    else:
        print("ERROR: missing mask. mask must be xarray.DataArray 'mask' or ['fractional', 'centred', 'min_overlap']. Aborting ...")
        return
        
    if select_abbr is not None:
        mask = mask[np.isin(mask.abbrevs,list(select_abbr))]
    elif select_name is not None:
        mask = mask[np.isin(mask.names,list(select_name))]
    else:
        # use all regions
        mask = mask
    
    if start is not None and end is not None:
        # select GWL time slice
        ds = ds.sel(time = slice(str(start), str(end)))
    
    # calculate weights due to latitude
    lat_weights = np.cos(np.deg2rad(ds['lat']))
    # create your weighted 3D xr.Dataset
    if var is None:
        try:
            var = list(ds.keys())[0]
            ds_weighted = ds[var].weighted(mask * lat_weights)
        except:
            print(f"Please enter var. One of {list(ds.keys())}")
            return 
    else:
        ds_weighted = ds[var].weighted(mask * lat_weights)
        
    # perform a statistic calculation
    if quantile is not None:
        return ds_weighted.quantile(quantile, dim=dims).drop_vars("quantile").rename(f"{var}_p{quantile*100:.0f}")
    if how == "mean":
        return ds_weighted.mean(dim=dims).rename(f"{var}_mean")
    if how == "sum":
        return ds_weighted.sum(dim=dims).rename(f"{var}_sum")
    if how == "std":
        return ds_weighted.std(dim=dims).rename(f"{var}_std")
    if how == "var":
        return ds_weighted.var(dim=dims).rename(f"{var}_var")
    if how == "min":
        return ds_weighted.quantile(0.0, dim=dims).drop_vars("quantile").rename(f"{var}_min")
    if how == "median":
        return ds_weighted.quantile(0.5, dim=dims).drop_vars("quantile").rename(f"{var}_median")
    if how == "max":
        return ds_weighted.quantile(1., dim=dims).drop_vars("quantile").rename(f"{var}_max")
        
    print("No statistic calculated. Please provide valid how, one of: ['mean', 'median', 'min', 'max', 'sum', 'std', 'var']")
    return ds_weighted

