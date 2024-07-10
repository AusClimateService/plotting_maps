""" Cookie cutter - ACS area statistics
Use this module to produce summary statistics for Hazards over NCRA regions,
with the flexibility to apply the method to any shapefile region.
Typical statistics include median, mean, min, max, 10th, 90th percentiles.

This method has used guidance from
[https://github.com/aus-ref-clim-data-nci/shapefiles/blob/master/python_tutorial.ipynb]"""

import xarray as xr
import geopandas as gpd
import regionmask
import numpy as np
import pandas as pd

# read in the shapefile with regions you will use Australia and NCRA regions
# from acs_plotting_maps import regions_dict
ncra_gdf = gpd.read_file(
    "/g/data/ia39/aus-ref-clim-data-nci/shapefiles/data/ncra_regions/ncra_regions.shp"
)
ncra_gdf["abbrevs"] = ncra_gdf["short_name"]
ncra_gdf["NAME"] = ncra_gdf["regionname"]

aus_gdf = gpd.read_file(
    "/g/data/ia39/aus-ref-clim-data-nci/shapefiles/data/australia/australia.shp"
)
aus_gdf["abbrevs"] = ["AUS"]
aus_gdf["NAME"] = ["Australia"]

# make sure the projections are the same
ncra_gdf = ncra_gdf.to_crs(aus_gdf.crs)

gdf = pd.concat(
    [
        ncra_gdf[["NAME", "abbrevs", "geometry"]],
        aus_gdf[["NAME", "abbrevs", "geometry"]],
    ],
)
gdf.index = np.arange(0, 10)
regions = regionmask.from_geopandas(gdf, names="NAME", abbrevs="abbrevs", overlap=True)


def acs_regional_stats(
    ds=None,
    infile=None,
    var=None,
    mask=None,
    start=None,
    end=None,
    dims=None,
    how=None,
    outfile=None,
    select_abbr=None,
    select_name=None,
    overlap_threshold=None,
    bins=None,
    bin_labels=None,
):
    """
    This function takes an Xarray dataset (ds) with variable (var)
    and multiple dimensions (eg time, lat, and lon),
    then selects the time range between two years (start and end),
    and applies regions.mask_3D_frac_approx fractional mask (frac)
    to compute a regional statistic (how, eg "mean") over two or three dimensions.
    Best used with numerical data without nans.

    Parameters
    ----------
    ds: xr.Dataset or xr.DataArray
        expects an xr.Dataset with variable var and dimensions time, lat, and lon.

    infile: str
        NetCDF file to read in as xr.Dataset

    var: str
        name of variable in ds, eg "pr" or "tas".
        If None, then tries to infer the var name from the data

    mask:  xarray.DataArray 'mask' or ["fractional", "centred", "min_overlap"]
        expects a precalculated mask, or will calculate one of the three options. 
        If "min_overlap" selected, you must provide an overlap_threshold value.
        fractional mask from regionmask.from_geopandas(ncra_gdf, 
                                                       names="NAME",
                                                       abbrevs="abbrevs"
                                                       ).mask_3D_frac_approx(ds).
        If None is provided, then will calculate a fractional mask based on NCRA regions,
        but this will take about one minute to calculate.
        This is best avoided by calculating frac outside the function.

    start: string or int
        start year to slice data array. eg start year of global warming level GWL. 
        If either start or end is None, then time selection is not performed.

    end: string or int
        end year to slice data array. eg end year of global warming level GWL.
        If either start or end is None, then time selection is not performed.

    dims: tuple of dimensions
        Dimensions to reduce data. If None, will reduce to one value per region. 
        Suggestion to reduce to one value per region: dims= ("time", "lat", "lon"). 
        Or to reduce to 1D time series per region ("lat", "lon").

    how: list of str
        List of statistics to reduce data. 
        List of ['mean', 'median', 'min', 'max', 'sum', 'std', 'var', 'p10', 'p90']. 
        (any pxx where xx is between 0 and 100)

    quantile: float [0.0,1.], optional
        REPLACED with how = "pxx" where xx is number between 0 and 100. Percentile to calculate.

    outfile: str
        csv filename to save dataframe.

    select_abbr: list of str
        list of regions by abbreviation to perform statistic on. eg ["VIC", "NSW"]

    select_name: list of str
        list of regions by name to perform statistic on.  eg ["Victoria", "New South Wales & ACT"]

    overlap_threshold: float between 0.0 and 1.
        If mask =  "min_overlap", you must provide an overlap_threshold value.
        If mask != "min_overlap", this value will be ignored.

    Returns
    -------
    An xarray.DataArray (1D or 2D) with dimension called "region" and values of regional statistics.
    If no valid statistic is called through how, then the DataArrayWeighted is returned


    Example usage
    --------------
    # open data
    ds = xr.open_dataset(filename, use_cftime = True,)
    # define mask
    mask_frac = regions.mask_3D_frac_approx(ds)

    # apply function
    da_mean = acs_regional_stats(ds=ds, 
                                 var="pr",
                                 mask=mask_frac,
                                 start="1991",
                                 end="2010", 
                                 dims=("time", "lat", "lon"),
                                 how="mean")
    """
    # Open file if ds not provided and infile is provided
    if ds is None and infile is not None:
        ds = xr.open_dataset(
            infile,
            use_cftime=True,
        )

    # --- start calculating mask ---
    if isinstance(mask, xr.DataArray):
        # Prefered method: use precalculated xarray.DataArray 'mask'
        pass
    elif mask == "fractional":
        # !warning very slow!
        print(
            "!warning very slow! Calculating fractional mask every time is very slow. \
        \nPlease consider calculating `mask = regions.mask_3D_frac_approx(ds)` before function."
        )
        # mask = regions.mask_3D_frac_approx(ds)
        try:
            mask = regions.mask_3D_frac_approx(ds)
        except:
            # This loop deals with rounding errors in lat lon that make the fractional mask fail
            mask = None
            i = 6
            while mask is None and i > 0:
                try:
                    mask = regions.mask_3D_frac_approx(ds)
                    print(f"rounded lat and lon to {i} decimal places")
                except:
                    # Adjust lat and lon to correct for float problems!
                    ds = ds.assign_coords(
                        lat=ds.lat.astype("double").round(i),
                        lon=ds.lon.astype("double").round(i),
                    )
                    i -= 1
    elif mask == "centred":
        # !warning slow!
        print(
            "!warning slow! Calculating mask every time is slow. \
        \nPlease consider calculating ```mask = regions.mask_3D(ds)``` before this function."
        )
        mask = regions.mask_3D(ds)
    elif mask == "min_overlap":
        # !warning very slow!
        print(
            "!warning very slow! Calculating fractional mask every time is very slow.\
        \nPlease consider calculating `mask = regions.mask_3D_frac_approx(ds) >= overlap_threshold`\
         before this function."
        )
        assert (
            0.0 <= overlap_threshold <= 1.0
        ), "You have selected min_overlap mask. Please specify overlap_threshold between [0.,1.]"
        mask = regions.mask_3D_frac_approx(ds) >= overlap_threshold
    else:
        print(
            "ERROR: missing mask.\
            \nMask must be xarray.DataArray 'mask' or ['fractional', 'centred', 'min_overlap'].\
            \nAborting ..."
        )
        return None
    # --- end calculating mask ---

    if select_abbr is not None:
        mask = mask[np.isin(mask.abbrevs, list(select_abbr))]
    elif select_name is not None:
        mask = mask[np.isin(mask.names, list(select_name))]
    else:
        # use all regions
        pass

    if start is not None and end is not None:
        # select GWL time slice
        ds = ds.sel(time=slice(str(start), str(end)))

    # calculate weights due to latitude
    lat_weights = np.cos(np.deg2rad(ds["lat"]))
    # drop redundant coords
    redundant_coords = set(lat_weights.coords) - set(lat_weights.dims)
    lat_weights = lat_weights.drop_vars(redundant_coords)
    mask = mask.drop_vars(redundant_coords)

    ds_masked = ds.where(mask)[var]

    # create your weighted 3D xr.Dataset
    if var is None:
        try:
            # find the variable that has the full dimensions of dataset (excluding time_bnds etc)
            while var is None:
                for v in list(ds.variables):
                    if list(ds[v].coords) == list(ds.coords):
                        var = v
        except:
            print(f"Please enter var. One of {list(ds.keys())}")
            return None
    ds_weighted = ds[var].weighted(mask * lat_weights)

    # if calculating stats over all dims:
    if dims is None:
        dims = list(ds.coords)

    # for every stat in the how list, add to summary list
    summary_list = []
    for stat in how:
        # perform a statistic calculation
        if stat.replace("p", "").isnumeric() and (int(stat.replace("p", "")) <= 100):
            q = int(stat.replace("p", "")) / 100
            summary_list.append(
                ds_weighted.quantile(q, dim=dims, )
                .drop_vars("quantile")
                .rename(f"{var}_{stat}")
            )
        elif stat == "mean":
            summary_list.append(ds_weighted.mean(dim=dims,).rename(f"{var}_{stat}"))
        elif stat == "sum":
            summary_list.append(ds_weighted.sum(dim=dims,).rename(f"{var}_{stat}"))
        elif stat == "std":
            summary_list.append(ds_weighted.std(dim=dims,).rename(f"{var}_{stat}"))
        elif stat == "var":
            summary_list.append(ds_weighted.var(dim=dims,).rename(f"{var}_{stat}"))
        elif stat == "min":
            summary_list.append(
                ds_masked.groupby("region").min(dim=dims,).rename(f"{var}_{stat}")
            )
            if bins is not None:
                df_masked = ds_masked.to_dataframe()
                df_masked["category"] = pd.cut(df_masked[var], bins, labels=bin_labels, ordered=True)
                summary_list.append(
                    df_masked.groupby("region")["category"]
                    .min()
                    .to_xarray()
                    .rename(f"{var}_cat_{stat}")
                )
        elif stat == "median":
            summary_list.append(
                ds_weighted.quantile(0.5, dim=dims,)
                .drop_vars("quantile")
                .rename(f"{var}_{stat}")
            )
        elif stat == "max":
            summary_list.append(
                ds_masked.groupby("region").max(dim=dims,).rename(f"{var}_{stat}")
            )
            if bins is not None:
                df_masked = ds_masked.to_dataframe()
                df_masked["category"] = pd.cut(df_masked[var], bins, labels=bin_labels, ordered=True)
                summary_list.append(
                    df_masked.groupby("region")["category"]
                    .max()
                    .to_xarray()
                    .rename(f"{var}_cat_{stat}")
                )

        elif stat == "mode":
            if bins is not None:
                df_masked = ds_masked.to_dataframe()
                df_masked["category"] = pd.cut(df_masked[var], bins, labels=bin_labels, ordered=True)
                summary_list.append(
                    df_masked.groupby("region")["category"]
                    .agg(pd.Series.mode)
                    .to_xarray()
                    .rename(f"{var}_cat_{stat}")
                )
            else:
                # mode cannot use the fractional masking in the same way as other statistics
                summary_list.append(
                    ds_masked.to_dataframe().groupby(["region"])[var]
                    .agg(pd.Series.mode)
                    .to_xarray()
                    .rename(f"{var}_{stat}")
                )
        else:
            print(
                f"{stat} statistic not calculated. \
Please provide valid how as a list including, one of: \
['mean', 'median', 'min', 'max','mode', 'sum', 'std', 'var', 'p10', 'p90']"
            )
    ds_summary = xr.merge(summary_list)
    df = ds_summary.to_dataframe()
    # drop columns with constant value
    df_summary = df[[col for col in df.columns if df[col].nunique() != 1]]

    if (
        outfile is None
        and infile is not None
        and select_abbr is None
        and select_name is None
    ):
        outfile = infile.replace(".nc", f"_summary-{'-'.join(how)}_ncra-regions.csv")
    if outfile is not None:
        try:
            df_summary.to_csv(outfile)
        except PermissionError:
            print(f"Could not save to {outfile}")
    return ds_summary
