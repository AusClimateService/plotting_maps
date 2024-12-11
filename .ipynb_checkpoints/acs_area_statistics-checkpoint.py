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
from glob import glob
import shapely

shapefile_list = ["aus_local_gov",
                  "aus_states_territories",
                  "australia", 
                  "broadacre_regions", 
                  "NCRA_Marine_region",
                  "ncra_regions", 
                  "NCRA_regions_coastal_waters_GDA94", 
                  "nrm_regions",
                  "plantations"]


name_dict = {"aus_local_gov":"LGA_NAME22", 
             "aus_states_territories":"STE_NAME21",
             "australia": "AUS_NAME21",
             "broadacre_regions": "name",
             "NCRA_Marine_region":"Label",
             "ncra_regions": "regionname", 
             "NCRA_regions_coastal_waters_GDA94": "regionname",
             "nrm_regions":"SubClusNm",
             "plantations":"CodeName",
            }


abbr_dict = {"aus_local_gov":"LGA_CODE22", 
             "aus_states_territories":"ABBREV",
             "australia": "AUS_CODE21",
             "broadacre_regions": "aagis",
             "NCRA_Marine_region":"RegionID",
             "ncra_regions": "short_name", 
             "NCRA_regions_coastal_waters_GDA94": "NCRA",
             "nrm_regions": "SubClusAb",
             "plantations":"REGCODE",
            }

def get_regions(shapefiles):
    """
    This function takes a list of names of shapefiles from ia39 and
    returns a combined regionmask. Renames columns so that there are
    columns named "NAME" and "abbrevs".  
    
    Parameters
    -----------
    name: str
        one of "aus_local_gov", "aus_states_territories", "australia", 
        "nrm_regions", "ncra_regions","broadacre_regions",
        "NCRA_regions_coastal_waters_GDA94"

    Returns
    -------
    geopandas dataframe
    
    """
    gdfs = {}
    PATH = "/g/data/ia39/aus-ref-clim-data-nci/shapefiles/data"
    
    for i, shapefile in enumerate(shapefiles):
        regions = gpd.read_file(glob(f"{PATH}/{shapefile}/*.shp")[0]).rename(columns = {name_dict[shapefile]:"NAME", abbr_dict[shapefile]:"abbrevs"}).to_crs(crs = "GDA2020")
        # simplify regions to ~110m resolution
        regions[["geometry"]] =shapely.simplify(regions[["geometry"]], 0.001)
        gdfs[i] = regions
    gdf = pd.concat(gdfs)
    gdf.index = np.arange(0, len(gdf))
    return regionmask.from_geopandas(gdf, names="NAME", abbrevs="abbrevs", name= "-".join(shapefiles), overlap=True) 


# read in the shapefile with regions you will use Australia and NCRA regions
# from acs_plotting_maps import regions_dict
# regions = get_regions(["ncra_regions", "australia"])

def calc_mask(mask=None, regions=None, ds=None, overlap_threshold=None):
    """Helper function for calculating masks.
    If mask is already an xr.dataarray, then just return this mask.
    Else if, mask is a string ["fractional", "centred", "min_overlap"], then calculate the relevant mask.
    Else, raise an exception.

    Parameters
    ----------
    mask:  xarray.DataArray 'mask' or ["fractional", "centred", "min_overlap"]
        expects a precalculated mask, or will calculate one of the three options. 
        There are some mask files in /g/data/ia39/aus-ref-clim-data-nci/shapefiles/masks/AGCD-05i/
        which provide ~5km gridded masks for shapefiles in the above list. eg use
        /g/data/ia39/aus-ref-clim-data-nci/shapefiles/masks/AGCD-05i/mask-3D-frac-approx_ncra-regions.nc 
        for fractional mask for NCRA regions.
        If "min_overlap" selected, you must provide an overlap_threshold value.
        fractional mask from regionmask.from_geopandas(ncra_gdf, 
                                                       names="NAME",
                                                       abbrevs="abbrevs"
                                                       ).mask_3D_frac_approx(ds).
        If None is provided, then will calculate a fractional mask based on NCRA regions,
        but this will take about one minute to calculate.
        This is best avoided by calculating the fractional mask outside the function.

    regions:
        regions to use to calculate mask if no mask is provided.
        use get_regions function for regions in shapefile list
        eg get_regions(["ncra_regions", "australia"])  

    ds: xr.Dataset or xr.DataArray
        expects an xr.Dataset with variable var and dimensions time, lat, and lon.

    overlap_threshold: float between 0.0 and 1.
        If mask =  "min_overlap", you must provide an overlap_threshold value.
        If mask != "min_overlap", this value will be ignored.
    
    Returns
    -------
    mask as xr.dataarray
    """
    if isinstance(mask, xr.DataArray):
        # Prefered method: use precalculated xarray.DataArray 'mask'
        pass
    elif mask == "fractional":
        # !warning very slow!
        print(
            "!warning very slow! Calculating fractional mask every time is very slow. \
        \nPlease consider calculating `mask = regions.mask_3D_frac_approx(ds)` before function."
        )
        if regions is None:
            regions = get_regions(["ncra_regions", "australia"])

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
        if regions is None:
            regions = get_regions(["ncra_regions", "australia"])
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
        if regions is None:
            regions = get_regions(["ncra_regions", "australia"])
        mask = regions.mask_3D_frac_approx(ds) >= overlap_threshold
    else:
        raise Exception(
            "Missing mask. Mask must be xarray.DataArray 'mask' or ['fractional', 'centred', 'min_overlap'].\
            \nAborting ..."
        )
    return mask

def calc_stats(ds_weighted=None,
               ds_masked=None,
               var=None,
               dims=None,
               how=None,
               bins=None,
               bin_labels=None):
    """
    Helper function for acs_regional_stats. 
    Takes prepared xr.datasets and calculates regional statistics into a summarised xr.dataset.
    
    Parameters
    ----------
    ds_weighted: xr.dataset
    
    ds_masked: xr.dataset

    var: str
        name of variable in ds, eg "pr" or "tas".
        If None, then tries to infer the var name from the data

    dims: tuple of dimensions
        Dimensions to reduce data. If None, will reduce to one value per region. 
        Suggestion to reduce to one value per region: dims= ("time", "lat", "lon"). 
        Or to reduce to 1D time series per region ("lat", "lon").

    how: list of str
        List of statistics to reduce data. 
        List of ['mean', 'median', 'min', 'max', 'mode', 'sum', 'std', 'var', 'proportions', 'p10', 'p90', ]. 
        (any pxx where xx is between 0 and 100)

    bins: int, sequence of scalars, or IntervalIndex
        For categorical data. The criteria to bin by using pandas.cut 'bins'

    bin_labels: array or False, default None
        For labeled categorical data.
        passed to pandas.cut 'labels'
        Specifies the labels for the returned bins. Must be the same length as the resulting bins. 
    
    Returns
    -------
    xr.dataset of calculated summary statistics of ds_weighted and/or ds_masked
    """
    
    # for every stat in the how list, add to summary list
    summary_list = []
    for stat in how:
        # perform a statistic calculation
        if stat.replace("p", "").isnumeric() and (int(stat.replace("p", "")) <= 100):
            q = int(stat.replace("p", "")) / 100
            summary_list.append(
                ds_weighted.quantile(q, dim=dims, keep_attrs=True)
                .drop_vars("quantile")
                .rename(f"{var}_{stat}")
            )
        elif stat == "mean":
            summary_list.append(ds_weighted.mean(dim=dims, keep_attrs=True).rename(f"{var}_{stat}"))
        elif stat == "sum":
            summary_list.append(ds_weighted.sum(dim=dims, keep_attrs=True).rename(f"{var}_{stat}"))
        elif stat == "std":
            summary_list.append(ds_weighted.std(dim=dims, keep_attrs=True).rename(f"{var}_{stat}"))
        elif stat == "var":
            summary_list.append(ds_weighted.var(dim=dims, keep_attrs=True).rename(f"{var}_{stat}"))
        elif stat == "min":
            summary_list.append(
                ds_masked.groupby("region").min(dim=dims, keep_attrs=True).rename(f"{var}_{stat}")
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
                ds_weighted.quantile(0.5, dim=dims, keep_attrs=True)
                .drop_vars("quantile")
                .rename(f"{var}_{stat}")
            )
        elif stat == "max":
            summary_list.append(
                ds_masked.groupby("region").max(dim=dims, keep_attrs=True).rename(f"{var}_{stat}")
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
        elif stat == "proportions":
            if bins is not None:
                df_masked = ds_masked.to_dataframe()
                df_masked["category"] = pd.cut(df_masked[var], bins, labels=bin_labels, ordered=True)
                proportion_dict = {}
                props = df_masked.groupby("region").value_counts(["category"], normalize=True).round(4)
                for i in range(len(ds_masked.region)):
                    proportion_dict[i] = props[i].to_dict() 
                    # limit to first 20 categories
                    proportion_dict[i]  = dict(list(proportion_dict[i].items())[0: 20]) 
                    
                d = {"coords": {
                        "region": {"dims": "region", "data": list(proportion_dict.keys()),}
                    },
                     "dims": "region",
                    "data":list(proportion_dict.values()),
                    "name": stat}
                summary_list.append(
                    xr.DataArray.from_dict(d).rename(f"{var}_cat_{stat}")
                )
            else:
                df_masked = ds_masked.to_dataframe()
                proportion_dict = {}
                props = df_masked.groupby("region").value_counts([var], normalize=True).round(4)
                for i in range(len(ds_masked.region)):
                    proportion_dict[i] = props[i].to_dict() 
                    # limit to first 20 categories
                    proportion_dict[i]  = dict(list(proportion_dict[i].items())[0: 20]) 
                    
                d = {"coords": {
                        "region": {"dims": "region", "data": list(proportion_dict.keys()),}
                    },
                     "dims": "region",
                    "data":list(proportion_dict.values()),
                    "name": stat}
                summary_list.append(
                    xr.DataArray.from_dict(d).rename(f"{var}_{stat}")
                )

        else:
            raise Exception(f"{stat} statistic not calculated. \
Please provide valid how as a list including, one of: \
['mean', 'median', 'min', 'max', 'mode', 'sum', 'std', 'var', 'p10', 'p90', 'proportions']"
            )
    ds_summary = xr.merge(summary_list)
    return ds_summary

def acs_regional_stats(
    ds=None,
    infile=None,
    var=None,
    mask=None,
    regions=None,
    dims=None,
    how=None,
    outfile=None,
    select_abbr=None,
    select_name=None,
    overlap_threshold=None,
    bins=None,
    bin_labels=None,
    chunks=None,
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
        There are some mask files in /g/data/ia39/aus-ref-clim-data-nci/shapefiles/masks/AGCD-05i/
        which provide ~5km gridded masks for shapefiles in the above list. eg use
        /g/data/ia39/aus-ref-clim-data-nci/shapefiles/masks/AGCD-05i/mask-3D-frac-approx_ncra-regions.nc 
        for fractional mask for NCRA regions.
        If "min_overlap" selected, you must provide an overlap_threshold value.
        fractional mask from regionmask.from_geopandas(ncra_gdf, 
                                                       names="NAME",
                                                       abbrevs="abbrevs"
                                                       ).mask_3D_frac_approx(ds).
        If None is provided, then will calculate a fractional mask based on NCRA regions,
        but this will take about one minute to calculate.
        This is best avoided by calculating the fractional mask outside the function.

    regions:
        regions to use to calculate mask if no mask is provided.
        use get_regions function for regions in shapefile list
        eg get_regions(["ncra_regions", "australia"])

    dims: tuple of dimensions
        Dimensions to reduce data. If None, will reduce to one value per region. 
        Suggestion to reduce to one value per region: dims= ("time", "lat", "lon"). 
        Or to reduce to 1D time series per region ("lat", "lon").

    how: list of str
        List of statistics to reduce data. 
        List of ['mean', 'median', 'min', 'max', 'mode', 'sum', 'std', 'var', 'proportions', 'p10', 'p90', ]. 
        (any pxx where xx is between 0 and 100)

    outfile: str
        csv filename to save dataframe.

    select_abbr: list of str
        list of regions by abbreviation to perform statistic on. eg ["VIC", "NSW"]
        Prefered over select_name.

    select_name: list of str
        list of regions by name to perform statistic on.  eg ["Victoria", "New South Wales & ACT"]

    overlap_threshold: float between 0.0 and 1.
        If mask =  "min_overlap", you must provide an overlap_threshold value.
        If mask != "min_overlap", this value will be ignored.

    bins: int, sequence of scalars, or IntervalIndex
        For categorical data. The criteria to bin by using pandas.cut 'bins'

    bin_labels: array or False, default None
        For labeled categorical data.
        passed to pandas.cut 'labels'
        Specifies the labels for the returned bins. Must be the same length as the resulting bins. 

    chunks: int, default None
        Chunk size to control memory resources.
        If not None, the loop over chunks of "time" dim of size "chunks". 
        Eg, for daily data, try chunks=365.
        

    Returns
    -------
    An xarray.DataArray (1D or 2D) with dimension called "region" and values of regional statistics.
    If no valid statistic is called through how, then the DataArrayWeighted is returned


    Example usage
    --------------
    # open data
    ds = xr.open_dataset(filename, use_cftime = True,)
    # define mask
    regions = get_regions(["ncra_regions", "australia"])
    mask_frac = regions.mask_3D_frac_approx(ds)

    # apply function
    da_mean = acs_regional_stats(ds=ds, 
                                 var="pr",
                                 mask=mask_frac,
                                 dims=("lat", "lon"),
                                 how=["mean"])
    """
    # Open file if ds not provided and infile is provided
    if ds is None and infile is not None:
        ds = xr.open_dataset(
            infile,
            use_cftime=True,
        )

    # --- start calculating mask ---
    mask = calc_mask(mask=mask, regions=regions, ds=ds, overlap_threshold=overlap_threshold)
    # --- end calculating mask ---

    # --- start selecting regions ---
    if select_abbr is not None:
        mask = mask[np.isin(mask.abbrevs, list(select_abbr))]
    elif select_name is not None:
        mask = mask[np.isin(mask.names, list(select_name))]
    else:
        # use all regions
        pass
    # --- end selecting regions ---

    # --- start assign var ---
    if var is None:
        try:
            # find the variable that has the full dimensions of dataset (excluding time_bnds etc)
            while var is None:
                for v in list(ds.variables):
                    if list(ds[v].coords) == list(ds.coords):
                        var = v
        except:
            raise Exception(f"Please enter var. One of {list(ds.keys())}")
    # --- end assign var ---

    # mask invalid nan and inf etc
    ds[var].values = np.ma.masked_invalid(ds[var].values)

    # --- start assign dims ---
    # if calculating stats over all dims:
    if dims is None:
        dims = list(ds.coords)
    # --- end assign dims

    # --- start apply weight and mask to ds ---
    # calculate weights due to latitude
    lat_weights = np.cos(np.deg2rad(ds["lat"]))
    # drop redundant coords
    redundant_coords = set(lat_weights.coords) - set(lat_weights.dims)
    lat_weights = lat_weights.drop_vars(redundant_coords)
    try:
        mask = mask.drop_vars(redundant_coords)
    except:
        pass
    mask_x_lat_weights = mask * lat_weights
    # create your weighted 3D xr.Dataset witihin chunks condition
    # ds_weighted = ds[var].weighted(mask_x_lat_weights)
    if any([stat in ["proportions", "mode", "max", "min"] for stat in how]):
        #only calculate this if needed
        ds_masked = ds[var].where(mask)
    else:
        ds_masked=None
    # --- end apply weight and mask to ds ---

    # --- start calc stats ---
    if chunks is None:
        ds_weighted = ds[var].weighted(mask_x_lat_weights)
        ds_summary = calc_stats(ds_weighted=ds_weighted,
                                ds_masked=ds_masked,
                                var=var,
                                dims=dims,
                                how=how,
                                bins=bins,
                                bin_labels=bin_labels,
                               )
    else:
        assert "time" in set(ds.coords), "'time' must be in ds.coords to use 'chunks' keyword"
        summaries = []
        for i in np.arange(len(ds["time"])//chunks+1):
            ds_weighted_chunk = ds[var].isel(time=slice(chunks*i, chunks*(i+1))).weighted(mask_x_lat_weights)
            if ds_masked is not None:
                ds_masked_chunk = ds_masked.isel(time=slice(chunks*i, chunks*(i+1)))
            else:
                ds_masked_chunk=None
            summaries.append(calc_stats(ds_weighted=ds_weighted_chunk,
                                        ds_masked=ds_masked_chunk,
                                        var=var,
                                        dims=dims,
                                        how=how,
                                        bins=bins,
                                        bin_labels=bin_labels,
                                       )
                            )
        ds_summary = xr.concat(summaries, dim = "time",)
    # --- end calc stats ---

    # --- start save to csv ---
    if (
        outfile is None
        and infile is not None
        and select_abbr is None
        and select_name is None
    ):
        outfile = infile.replace(".nc", f"_summary-{'-'.join(how)}_ncra-regions.csv")
    if outfile is not None:
        try:
            ds_summary.to_dataframe().to_csv(outfile)
        except PermissionError:
            print(f"Could not save to {outfile}")
    # --- end save to csv ---
    return ds_summary
