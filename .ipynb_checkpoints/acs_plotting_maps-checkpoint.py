"""Standardising Australia Hazard Maps
This module plots maps to consistently present climate hazards for Australia.
It is code is designed to work with hh5 analysis3-24.04 venv"""

import datetime
import os

# import packages used in this workflow
import geopandas as gpd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import image, cm
import xarray as xr
import cartopy.crs as ccrs
from glob import glob

# import colormap packages
import cmaps
from matplotlib.colors import ListedColormap, BoundaryNorm, LinearSegmentedColormap

from shapely.geometry import box

import matplotlib as mpl
mpl.rcParams['hatch.linewidth'] = 0.3 
plt.rcParams['savefig.facecolor']='white'

# define some standard imput for the maps
projection = ccrs.LambertConformal(
    central_latitude=-24.75,
    central_longitude=134.0,
    cutoff=30,
    standard_parallels=(-10, -40),
)

from pathlib import Path
logo = image.imread(f"{Path(__file__).parent}/ACS_Logo_Blue_on_white_Stacked.png")

# # Suggested colormaps and scales
# Using suggested colormaps and scales will improve the consistency across teams
# producing similar variables. This will support comparison across different plots.
# - see many colormaps here:
# https://www.ncl.ucar.edu/Document/Graphics/color_table_gallery.shtml
# This suggested colormaps are matched with possible variables to plot.
# This includes color maps for the total amount and anomalies

#ipcc colormaps from github.com/IPCC-WG1/colormaps/

cmap_mustard = LinearSegmentedColormap.from_list(
    "mustard",
    # [(195 / 235, 152 / 235, 21 / 235), (229 / 235, 208 / 235, 147 / 235)],
    # ["#f4ea94", "#e2ce17", "#61580a"]
    # [ "#3c320a", "#5a4c10", "#977f1b","#d3b125", "#e2c85e",  "#eddd9b",]
    ["#5a4c10", "#977f1b","#d3b125", "#e2c85e",  "#eddd9b",]
)
cmap_mustard.set_bad(color="lightgrey")

cmap_BuGnPi = LinearSegmentedColormap.from_list('BuGnPi', np.vstack((cm.BrBG(np.linspace(0.55, 1, 30)[::-1]),[0.9,0.9,0.9,0.9], cm.pink_r(np.linspace(0.25, 1, 30)) )))
cmap_BuGnPi_r = LinearSegmentedColormap.from_list('BuGnPi_r',cmap_BuGnPi(np.linspace(0,1,256))[::-1])

cmap_dir = f"{Path(__file__).parent}/continuous_colormaps_rgb_0-1"

cmap_dict = {
    "sst": cmaps.cmocean_tempo,
    "sst_anom": cmaps.cmocean_balance_r,
    "mhw_days": cm.YlOrRd,
    "mhw_intensity": cm.hot_r,
    "hot_r": cm.hot_r,
    "surface_pH": cm.YlGnBu,
    "surface_aragonite_sat": cmaps.cmocean_delta,
    "tas": cm.Spectral_r,
    "tas_anom": cm.RdBu_r,
    "tas_anom_1": cm.seismic,
    "tas_deciles_bwr": cm.bwr,
    "EHF_days": cm.YlOrRd,
    "EHF_days_1": cm.YlOrBr,
    "EHF_duration": cm.hot_r,
    "AFDRS_category": ListedColormap(["white", "green", "orange", "red", "darkred"]),
    "ffdi_category": ListedColormap(
        ["green", "blue", "yellow", "orange", "red", "darkred"]
    ),
    "fire_climate": ListedColormap(
        [ "#84a19b", "#e0d7c6", "#486136", "#737932", "#a18a6e", ]
    ),
    "fire_climate_alternative": ListedColormap(
        [ "#355834", "#F1F5F2", "#14281D", "#6E633D", "#C2A878", ]
    ),
    'tasmax_bom': ListedColormap(
        [ '#E3F4FB','#C8DEE8','#91C4EA','#56B6DC','#00A2AC','#30996C',
         '#7FC69A','#B9DA88','#DCE799', '#FCE850','#EACD44','#FED98E',
         '#F89E64','#E67754','#D24241', '#AD283B','#832D57','#A2667A','#AB9487']
    ), #not colourblind safe
    'tasmax' : ListedColormap(
        ['#014636','#016c59','#02818a','#3690c0','#67a9cf','#a6bddb',
         '#d0d1e6','#ece2f0','#fff7fb','#ffffcc','#ffeda0','#fed976',
         '#feb24c','#fd8d3c','#fc4e2a','#e31a1c','#bd0026','#800026',
         '#510019','#2E000E']
    ),
    "pr": cm.YlGnBu,
    "pr_1": cmaps.cmocean_deep,
    "pr_days": cm.Blues,
    "pr_GMT_drywet": cmaps.GMT_drywet,
    "pr_anom": cm.BrBG,
    "pr_anom_1": cmaps.cmocean_curl,
    "pr_anom_12lev": cmaps.precip_diff_12lev,
    "pr_chance_extremes": cmaps.cmocean_solar_r,
    "tc_days": cm.RdPu,
    "tc_intensity": cm.PuRd,
    "tc_days_anom": cm.PuOr,
    "tc_intensity_anom": cm.PiYG_r,
    "xts_freq": cmaps.cmocean_dense,
    "xts_intensity": cmaps.cmocean_matter,
    "xts_freq_anom": cmaps.cmocean_balance_r,
    "xts_intensity_anom": cmaps.cmocean_curl_r,
    "drought_severity": cm.RdYlGn, # not colorblind friendly
    "drought_severity_r": cm.RdYlGn_r, # not colorblind friendly
    "drought_duration": cmaps.hotres,
    "drought_duration_r": cmaps.hotres_r,
    "aridity": cmap_mustard,
    # "aridity_anom": cmaps.NEO_div_vegetation_a, # not colorblind friendly
    # "aridity_anom_r": cmaps.NEO_div_vegetation_a_r, # not colorblind friendly
    "aridity_anom": cmap_BuGnPi_r,
    "aridity_anom_r": cmap_BuGnPi,
    "BrBu": LinearSegmentedColormap.from_list("BrBu", ["#3c320a", "#d3b125", "lightgrey", "royalblue", "navy"]),
    "BuBr": LinearSegmentedColormap.from_list("BuBr", ["navy", "royalblue", "lightgrey", "#d3b125", "#3c320a"]),
    "anom_BlueYellowRed": cmaps.BlueYellowRed,
    "anom_BlueYellowRed_r": cmaps.BlueYellowRed_r,
    "anom": cmaps.BlueWhiteOrangeRed,
    "anom_r": cmaps.BlueWhiteOrangeRed_r,
    "anom_b2r": cmaps.cmp_b2r,
    "anom_b2r_r": cmaps.cmp_b2r_r,
    "anom_coolwarm": cmaps.MPL_coolwarm,
    "anom_coolwarm_r": cmaps.MPL_coolwarm_r,
    "anom_deciles": cm.bwr,
    "anom_deciles_r": cm.bwr_r,
    "anom_veg_1": cmaps.NEO_div_vegetation_a, # not colorblind friendly
    "anom_veg_1_r": cmaps.NEO_div_vegetation_a_r, # not colorblind friendly
    "BuGnRd": cmaps.temp1,
    "rh_19lev": cmaps.rh_19lev,
    "sunshine_9lev": cmaps.sunshine_9lev,
    "sunshine_diff_12lev": cmaps.sunshine_diff_12lev,
    "inferno": cm.inferno,
    "Oranges": cm.Oranges,
    "Oranges_r": cm.Oranges_r,
    "OrRd": cm.OrRd,
    "Greens": cm.Greens,
    "topo": cmaps.OceanLakeLandSnow,
    "gmt_relief": cmaps.GMT_relief,
    "ipcc_chem_div": LinearSegmentedColormap.from_list('colormap', np.loadtxt(f"{cmap_dir}/chem_div.txt")),
    "ipcc_chem_seq": LinearSegmentedColormap.from_list('colormap', np.loadtxt(f"{cmap_dir}/chem_seq.txt")),
    "ipcc_cryo_div": LinearSegmentedColormap.from_list('colormap', np.loadtxt(f"{cmap_dir}/cryo_div.txt")),
    "ipcc_cryo_seq": LinearSegmentedColormap.from_list('colormap', np.loadtxt(f"{cmap_dir}/cryo_seq.txt")),
    "ipcc_misc_div": LinearSegmentedColormap.from_list('colormap', np.loadtxt(f"{cmap_dir}/misc_div.txt")),
    "ipcc_misc_seq_1": LinearSegmentedColormap.from_list('colormap', np.loadtxt(f"{cmap_dir}/misc_seq_1.txt")),
    "ipcc_misc_seq_2": LinearSegmentedColormap.from_list('colormap', np.loadtxt(f"{cmap_dir}/misc_seq_2.txt")),
    "ipcc_misc_seq_3": LinearSegmentedColormap.from_list('colormap', np.loadtxt(f"{cmap_dir}/misc_seq_3.txt")),
    "ipcc_prec_div": LinearSegmentedColormap.from_list('colormap', np.loadtxt(f"{cmap_dir}/prec_div.txt")),
    "ipcc_prec_seq": LinearSegmentedColormap.from_list('colormap', np.loadtxt(f"{cmap_dir}/prec_seq.txt")),
    "ipcc_slev_div": LinearSegmentedColormap.from_list('colormap', np.loadtxt(f"{cmap_dir}/slev_div.txt")),
    "ipcc_slev_seq": LinearSegmentedColormap.from_list('colormap', np.loadtxt(f"{cmap_dir}/slev_seq.txt")),
    "ipcc_temp_div": LinearSegmentedColormap.from_list('colormap', np.loadtxt(f"{cmap_dir}/temp_div.txt")),
    "ipcc_temp_seq": LinearSegmentedColormap.from_list('colormap', np.loadtxt(f"{cmap_dir}/temp_seq.txt")),
    "ipcc_wind_div": LinearSegmentedColormap.from_list('colormap', np.loadtxt(f"{cmap_dir}/wind_div.txt")),
    "ipcc_wind_seq": LinearSegmentedColormap.from_list('colormap', np.loadtxt(f"{cmap_dir}/wind_seq.txt"))
}


# Here are some suggestions for the ticks/ scale for some variables.
# Some scales are taken from climate maps on bom.gov.au/climate
tick_dict = {
    "pr_annual": [0, 50, 100, 200, 300, 400, 600, 1000, 1500, 2000, 3000, 6000],
    "pr_6mon": [0, 50, 100, 200, 300, 400, 600, 900, 1200, 1800, 2400, 6000],
    "pr_3mon": [0, 10, 25, 50, 100, 200, 300, 400, 600, 800, 1200, 2500],
    "pr_mon": [0, 1, 5, 10, 25, 50, 100, 200, 300, 400, 600, 1200],
    "pr_hour": [0, 1, 2, 5, 10, 15, 20, 30, 50, 75, 100, 200,],
    "pr_days": [0, 2, 3, 5, 10, 20, 30, 40, 50, 75, 100, 125, 150, 175],
    "pr_anom_mon": [ -1000, -400, -200, -100, -50, -25, -10, 0, 10, 25, 50, 100, 200, 400, 1000,],
    "pr_anom_3mon": [ -2000, -600, -400, -200, -100, -50, -25, 0, 25, 50, 100, 200, 400, 600, 2000,],
    "pr_anom_6mon": [ -3000, -1200, -800, -400, -200, -100, -50, 0, 50, 100, 200, 400, 800, 1200, 3000,],
    "pr_anom_ann": [ -4000, -1800, -1200, -800, -400, -200, -100, 0, 100, 200, 400, 800, 1200, 1800, 4000,],
    "pr_diff_mon": [ -1000, -400, -200, -100, -50, -25, -10, 10, 25, 50, 100, 200, 400, 1000,],
    "pr_diff_ann": [ -3000, -1800, -1200, -800, -400, -200, -100, 100, 200, 400, 800, 1200, 1800, 3000,],
    "frost_days": [0, 10, 20, 30, 40, 50, 75, 100, 150, 300],
    "frost_days_mon": [0, 2, 5, 10, 15, 20, 25, 31],
    "tas": np.arange(-9, 52, 3),
    "tas_anom_day": np.arange(-14, 14.1, 2),
    "tas_anom_mon": np.arange(-7, 7.1, 1),
    "tas_anom_ann": np.arange(-3.5, 3.6, 0.5),
    "apparent_tas": np.arange(-6, 42, 3),
    "percent": np.arange(0, 101, 10),
    "xts_freq": [0.00, 0.005, 0.01, 0.02, 0.03, 0.05, 0.07, 0.10, 0.12, 0.15],
    "fire_climate_ticks": [ 100, 101, 102, 103, 104, ],
    "fire_climate_labels": [
        "Tropical Savanna",
        "Arid grass \nand woodland",
        "Wet Forest",
        "Dry Forest",
        "Grassland",
    ],
    "aridity_index_ticks": [0.0, 0.05, 0.2, 0.5, 0.65],
    "aridity_index_labels": ["Hyper-arid", "Arid", "Semi-arid", "Dry sub-humid"],
}

# # Load the State and Region shape files
class RegionShapefiles:
    """Load and return a shapefile based on its name."""

    def __init__(self, path, shapefiles):
        """Create an instance of the RegionShapefiles class.
        Parameters
        ----------
        path : str
            The path to the shapefiles directory.
        shapefiles : list
            A list of shapefile names to load.
        """
        self.path = path
        self.shapefiles = shapefiles
        self._regions_dict = {}

    def __call__(self, name):
        """Retrieve the shapefile for the given region name.
        Parameters
        ----------
        name : str
            The name of the region to retrieve.
        Returns
        -------
        GeoDataFrame or GeoSeries
            The shapefile data for the specified region.
        """
        if name not in self._regions_dict:
            if name in self.shapefiles:
                self._regions_dict[name] = gpd.read_file(glob(f"{self.path}/{name}/*.shp")[0])
            elif name == "not_australia":
                # Define a white mask for the area outside of Australian land
                # We will use this to hide data outside of the Australian land borders.
                # note that this is not a data mask,
                # the data under the masked area is still loaded and computed, but not visualised
                base_name = name[4:]  # Remove 'not_' prefix
                base_region = self(base_name).copy()
                base_region.crs = crs
                # This mask is a rectangular box around the maximum land extent of Australia
                # with a buffer of 10 degrees on every side,
                # with the Australian land area cut out so only the ocean is hidden.
                not_region = gpd.GeoSeries(
                    data=[
                        box(*box(*base_region.total_bounds).buffer(20).bounds
                        ).difference(base_region["geometry"].values[0])
                    ],
                    crs=ccrs.PlateCarree(),
                )
                self._regions_dict[name] = not_region
            else:
                raise ValueError(f"Shapefile for region '{name}' not found.")
        return self._regions_dict[name]

    def __getitem__(self, name):
        if name not in self._regions_dict:
            self(name)
        return self._regions_dict[name]

    def __setitem__(self, name, value):
        self._regions_dict[name] = value

    def keys(self):
        return self._regions_dict.keys()

    def __len__(self):
        return len(self._regions_dict)

    def __repr__(self):
        return repr(self._regions_dict)

    def update(self, *args, **kwargs):
        self._regions_dict.update(*args, **kwargs)


# Define the path and shapefiles
# These will be used for state boundaries, LGAs, NRM, etc
PATH = "/g/data/ia39/aus-ref-clim-data-nci/shapefiles/data"
shapefile_list = ["aus_local_gov",
                  "aus_states_territories",
                  "australia", 
                  "broadacre_regions", 
                  "NCRA_Marine_region",
                  "ncra_regions", 
                  "NCRA_regions_coastal_waters_GDA94", 
                  "nrm_regions",
                  "plantations"]

# Create an instance of the RegionShapefiles class
regions_dict = RegionShapefiles(PATH, shapefile_list)

# define a white mask for the area outside of Australian land
# We will use this to hide data outside of the Australian land borders.
# note that this is not a data mask,
# the data under the masked area is still loaded and computed, but not visualised
australia = regions_dict["australia"].copy()

# Define the CRS of the shapefile manually
australia = australia.to_crs(crs = "GDA2020")

# This mask is a rectangular box around the maximum land extent of Australia
# with a buffer of 10 degrees on every side,
# with the Australian land area cut out so only the ocean is hidden.
not_australia = gpd.GeoSeries(
    data=[
        box(*box(*australia.total_bounds).buffer(20).bounds).difference(
            australia["geometry"].values[0]
        )
    ],
    crs=ccrs.PlateCarree(),
)


def crop_cmap_center(cmap, ticks, mid, extend=None):
    """This function is used to align the centre of a colormap (cmap)
    to a specified midpoint (mid). Allows the cmap to be normalised on 
    the specified ticks with cbar_extend arrows taken into account.
    Intended for divergent colormaps that show anomalies that are mostly
    all positive (or negative), for example, temperature anomalies in 
    future climate projections.
    The shorter side of the colormap is cropped and not stretched - unlike 
    matplotlib.colors.TwoSlopeNorm.

    Parameters
    -------
    cmap: matplotlib colormap
        original colormap to crop

    ticks: list or np.array
        list or array of numerical ticks

    mid: float
        Set the midpoint value of the colormap. 
        For example, 0.0

    extend: {'neither', 'both', 'min', 'max'}
        Make pointed end(s) for out-of-range values (unless 'neither').
        These are set for a given colormap using the colormap set_under
        and set_over methods.
    
    Returns
    -------
    matplotlib.colors.LinearSegmentedColormap
    """
    ticks=np.array(ticks)
    # number of color segments:
    below = (ticks<mid).sum()
    above = (ticks>mid).sum()
    if extend =="both" or extend == "max":
        above=above+1
    if extend =="both" or extend =="min":
        below=below+1
    
    #total segments
    N = below+above
    
    #porportion below mid point 
    prop_below_mid = below/(max(below,above))
    # propotion above mid point
    prop_above_mid = above/(max(below,above))
    
    new_cmap_list = cmap(np.linspace(0.5*(1-prop_below_mid),
                                     0.5*(1+prop_above_mid),
                                     N)
                        )
    new_cmap = LinearSegmentedColormap.from_list("new_cmap",new_cmap_list)
    return new_cmap


# Define subfunctions for different parts of the plotting 
# so that they can be reused for single panel and multi panel plots
def plot_data(regions=None,
              data=None, 
              station_df = None,
              markersize=None,
              stippling=None,
              xlim=(114, 162),
              ylim=(-43.5, -7.5),
              cmap=cm.Greens,
              cbar_extend="both",
              ticks=None,
              tick_labels=None,
              contourf=False,
              contour=False,
              ax=None,
              subtitle = "",
              subtitle_xy = None,
              facecolor="none",
              edgecolor="k",
              area_linewidth=0.3,
              mask_not_australia = True,
              mask_australia=False,
              agcd_mask=False,
              select_area = None,
              vcentre=None,
             ):
    """This function takes one axis and plots the hazard data to one map of Australia. 
    This function takes gridded "data" and/or a "station_df" dataframe and "regions" shapefiles 
    to visualise hazard data from a 2D Xarray data array and plots the data on a map
    of Australia with the regions outlines.

    Parameters
    ----------
    regions: geopandas.GeoDataFrame
        region geometries for regions/states/catchments etc

    data: xr.DataArray
        a 2D xarray DataArray which has already computed the 
        average, sum, anomaly, metric or index you wish to visualise.
        This function is resolution agnostic.

    station_df: pd.DataFrame, optional
        a pandas.DataFrame with columns ["lon", "lat", variable]. 
        If station_df is given, then variable values are represented as dots on 
        the mapaccordingg to the lat lon coordinates and coloured according to
        cmap colors and ticks.

    markersize: int, optional
        default None. If None then the markersize will adjust to the size of the
        figure and the number of markers in the plot such that when there are
        many markers and the figure is small, the markersize is smaller.

    stippling: xr.DataArray
        a True/False mask to define regions of stippling hatching. 
        Intended to show information such as "model agreement for direction of change".

    xlim: tuple, optional
        longitude min and max of the plot area.
        default is cropped to Australian continent xlim=(114, 162)

    ylim: tuple, optional
        latitude min and max of the plot area.
        default is cropped to Australian continent ylim=(-43.5, -7.5),

    cmap:
        defines the colormap used for the data.
        See cmap_dict for suggested colormaps.
        Default cmap set to cm.Greens.
        Please choose appropriate colormaps for your data.

    cbar_extend: one of {'neither', 'both', 'min', 'max'}.
        eg "both" changes the ends of the colorbar to arrows to indicate that
        values are possible outside the scale shown.
        If contour or contourf is True, then cbar_extend will be overridden to "none".

    ticks: list or arraylike
        Define the ticks on the colorbar. Define any number of intervals. 
        This will make the color for each interval one discrete color, 
        instead of a smooth color gradient.
        If None, linear ticks will be auto-generated to fit the provided data.

    tick_labels: list
        Labels for categorical data. 
        If tick_labels is used, then pcolormesh is used to plot data 
        and does not allow contour or contourf to be used.
        Tick labels will correspond to the ticks.

    contourf: bool
        if True then the gridded data is visualised as smoothed filled contours. 
        Default is False.
        Use with caution when plotting data with negative and positive values;
        Check output for NaNs and misaligned values.  
        High resolution data is slow to compute.

    contour: bool
        if True then the gridded data is visualised as smoothed unfilled grey contours.
        Default is False.
        High resolution data is slow to compute.
        Using both contourf and contour results in smooth filled contours
        with grey outlines between the color levels.

    ax: matplotlib.axes.Axes
        Axes object of existing figure to put the plotting.

    subtitle: str
        default ""
        Intended to label global warming levels for subplots eg. "GWL 1.2"
        
    subtitle_xy: tuple, optional
        (x, y) location of subtitle relative to transAxes.
        defines the top left location for the subtitle. 

    facecolor: color
        color of land when plotting the regions without climate data and select_area is None. 
        facecolor recommendations include "white", "lightgrey", "none".
        default is "none"

    edgecolor: color
        defines the color of the state/region borders. 
        edgecolor recommendations include "black" and "white".

    mask_not_australia: boolean
        decides whether or not the area outside of Australian land is hidden 
        under white shape.
        Default is True.

    mask_australia: boolean
        decides whether or not Australian land is hidden under white shape.
        Eg, use when plotting ocean only.
        Default is False.

    agcd_mask: boolean
        If True, applies a ~5km mask for data-sparse inland areas of Australia.
        Default is False.

    area_linewidth: float, optional
        the width of state/region borders only. All other linewidths are hardcoded.

    select_area: list
        If None, then don't add region borders geometries.
        
    vcentre: float, eg 0
        default is None.
        Align centre of colormap to this value. 
        Intended for using a divergent colormap with uneven number of ticks 
        around the centre, eg for future temperature anomalies with a larger
        positive range compared to the negative range.

    Returns
    -------
    ax, norm, cont, middle_ticks
    ax, the matplotlib.axes.Axes with the plot
    norm, the normalisation for the colormap and plotted contours according to ticks
    cont, the QuadContourSet or QuadMesh of the plotted gridded data
    middle_ticks, the location to label categorical tick labels
    """

    if vcentre is not None:
        cmap = crop_cmap_center(cmap, ticks, vcentre, extend=cbar_extend)

    ax.set_extent([xlim[0], xlim[1], ylim[0], ylim[1]])
    
    middle_ticks=[]
    
    if data is None:
        norm=None
        cont=None
    else:
        data = data.squeeze()

        if agcd_mask:
            # mask data where observations are sparse
            directory = "/g/data/ia39/aus-ref-clim-data-nci/shapefiles/masks/AGCD-05i"
            agcd = xr.open_dataset(f"{directory}/mask-fraction_agcd_v1-0-2_precip_weight_1960_2022.nc").fraction
            data = data.where(agcd>=0.8)
        
        if ticks is None:
            norm = None
        else:
            # if ticks are labelled or if there is one more tick than tick labels,
            # do the usual normalisation
            if tick_labels is None or (len(tick_labels) == len(ticks) - 1):
                norm = BoundaryNorm(ticks, cmap.N, extend = cbar_extend)
                if tick_labels is not None:
                    middle_ticks = [
                        (ticks[i + 1] + ticks[i]) / 2 for i in range(len(ticks) - 1)
                    ]
                else:
                    middle_ticks = []
            else:
                middle_ticks = [
                    (ticks[i + 1] + ticks[i]) / 2 for i in range(len(ticks) - 1)
                ]
                outside_bound_first = [ticks[0] - (ticks[1] - ticks[0]) / 2]
                outside_bound_last = [ticks[-1] + (ticks[-1] - ticks[-2]) / 2]
                bounds = outside_bound_first + middle_ticks + outside_bound_last
                norm = BoundaryNorm(bounds, cmap.N, extend = "neither")
    
        # plot the hazard data
        if contourf and tick_labels is None:
            if data.max()>=0 and data.min()<=0: 
                print("Using contourf to plot data. Use with caution and check output for data crossing zero")
            cont = ax.contourf(
                data.lon,
                data.lat,
                data,
                cmap=cmap,
                norm=norm,
                levels=ticks,
                extend=cbar_extend,
                zorder=2,
                transform=ccrs.PlateCarree(),
            )
        else:
            cont = ax.pcolormesh(
                data.lon,
                data.lat,
                data,
                cmap=cmap,
                norm=norm,
                zorder=2,
                transform=ccrs.PlateCarree(),
            )
       
        if contour and tick_labels is None:
            ax.contour(
                data.lon,
                data.lat,
                data,
                colors="grey",
                norm=norm,
                levels=ticks,
                extend=cbar_extend,
                linewidths=0.2,
                zorder=3,
                transform=ccrs.PlateCarree(),
            )

    # for station data
    if station_df is not None:
        # assuming columns are named "lon", "lat", variable,
        gdf = gpd.GeoDataFrame(
            station_df, geometry=gpd.points_from_xy(station_df.lon, station_df.lat), crs=ccrs.PlateCarree()
            )
        var = gdf.columns[[2]][0]
        norm = BoundaryNorm(ticks, cmap.N, extend=cbar_extend)
        cont = ax.scatter(x=station_df.lon,
                          y=station_df.lat,
                          s=markersize, 
                          c=station_df[var],
                          # edgecolors="k", 
                          alpha = 0.8,
                          zorder=7,
                          transform=ccrs.PlateCarree(), 
                          cmap= cmap,
                          norm = norm)
    
    if stippling is not None:
        ax.contourf(stippling.lon,
                    stippling.lat,
                    stippling,
                    alpha=0,
                    hatches = ["","xxxxxx"],
                    zorder=10,
                    transform=ccrs.PlateCarree(),
                   )

    # cover area outside australia land area eg mask ocean
    if mask_not_australia:
        # inside the shape, fill white
        ax.add_geometries(
            not_australia,
            crs=ccrs.PlateCarree(),
            facecolor="white",
            linewidth=0,
            zorder=5,
        )

    # cover australia land area eg for ocean data
    if mask_australia:
        # inside the shape, fill white
        ax.add_geometries(
            australia["geometry"],
            crs=ccrs.PlateCarree(),
            facecolor="lightgrey",
            linewidth=0.3,
            edgecolor="k",
            zorder=4,
        )

    if select_area is None:
        # add region borders unless you have selected area
        ax.add_geometries(
            regions["geometry"],
            crs=ccrs.PlateCarree(),
            facecolor=facecolor,
            edgecolor=edgecolor,
            linewidth=area_linewidth,
            zorder=6,
        )

    # subtitle
    if subtitle_xy is None:
        subtitle_xy = (0.03, 0.14)
    ax.text(
        x=subtitle_xy[0],
        y=subtitle_xy[1],
        s=subtitle,
        fontsize=12,
        horizontalalignment="left",
        verticalalignment="top",
        transform=ax.transAxes,
        zorder=10,
        wrap=True,
    )

    return ax, norm, cont, middle_ticks

def plot_cbar(cont=None,
              norm=None,
              ax=None,
              cbar_extend=None, 
              cbar_label=None,
              ticks=None, 
              tick_labels=None,
              middle_ticks=[], 
              cax_bounds =None,
              contour=False,
              location=None,
              rotation=None,
             ):
    """This function defines and plots the colorbar. 
    It takes crucial information from the plot_data function.
    
    Parameters
    ----------
    cont:
        output from matplotlib plotting
        
    norm:
        normalisation
        
    ax: matplotlib.axes.Axes
        Axes to put the colorbar
        
    cbar_extend: one of {'neither', 'both', 'min', 'max'}.
        eg "both" changes the ends of the colorbar to arrows to indicate that
        values are possible outside the scale shown.
         
    cbar_label: str
        Title for colorbar. Include name of metric and [units]
        
    ticks: list or array
        numerical location of ticks
        
    tick_labels: list
        If categorical data, these labels go inbetween the numerical bounds set by ticks
        
    middle_ticks: list
        If categorical data, this specifies the location of the tick labels.
        Typically in the middle of the bounds set by ticks
        
    cax_bounds: [left, bottom, width, height]
        Colorbar axes bounds relative to ax
        
    contour: bool
        draw lines on colorbar if True
        Default is False
        
    location: {"top", "bottom", "left", "right"}
        location of the colorbar. Defaults to right.
        
    rotation: [-360,360]
        rotation of tick labels in degrees. Set to 0 for horizontal.

    Returns
    -------
    matplotlib.colorbar    
    
    """

    if cax_bounds is not None:
        cax = ax.inset_axes(cax_bounds)
    else:
        cax=None
    
    cbar = None
    
    if norm is None:
        return cbar
    
    if tick_labels is None:
        cbar = plt.colorbar(
            cont,
            ax=ax,
            extend=cbar_extend,
            cax=cax,
            ticks=ticks,
            norm=norm,
            location=location,
            fraction=0.046, 
            pad=0.04
        )
    else:
        # for categorical data
        cbar = plt.colorbar(
            cont,
            ax=ax,
            extend='neither',
            cax=cax,
            ticks=ticks,
            norm=norm,
            drawedges=True,
            location=location,
            fraction=0.046, 
            pad=0.04
        )
        if location=="bottom":
            if len(ticks) == len(tick_labels):
                cbar.ax.set_xticks(ticks, tick_labels, wrap=True, verticalalignment="top")
            elif len(middle_ticks) == len(tick_labels):
                cbar.ax.set_xticks(middle_ticks, tick_labels, wrap=True, verticalalignment="top")

        else:
            if len(ticks) == len(tick_labels):
                cbar.ax.set_yticks(ticks, tick_labels, wrap=True)
            elif len(middle_ticks) == len(tick_labels):
                cbar.ax.set_yticks(middle_ticks, tick_labels, wrap=True)

    cbar.ax.tick_params(labelsize=8)
    if contour and tick_labels is None:
        cbar.add_lines(cont)
    
    # Label colorbar
    if cbar is not None:
        if len(cbar_label)>35 and "\n" not in cbar_label:
            fontstretch = "condensed"
        else:
            fontstretch = "normal"
        cbar.ax.set_title(cbar_label, 
                          zorder=10,
                          loc="center",
                          fontsize=10,
                          fontstretch=fontstretch,
                          verticalalignment="baseline",
                         wrap=True)

    cbar.ax.tick_params(rotation=rotation) 
    return cbar

def plot_select_area(select_area=None,
                     ax=None, 
                     xlim=None,
                     ylim=None,
                     regions=None,
                     land_shadow=False,
                     area_linewidth=0.3,
                    ):
    """This function takes a list of named areas to plot and adjusts 
    the limits of the axis.   
    
    Parameters
    ----------
    select_area: list
        list of selected areas to plot. Must be name of area in regions.

    ax: matplotlib.axes.Axes
        axis 

    xlim:
        longitude limits to use if select_area is None

    ylim:
        latitude limits to use if select_area is None  

    regions: geopandas.GeoDataFrame
        region border data, must contain a column name with "NAME" in it
        to select those areas.

    land_shadow: bool
        whether or not to shade in the non-selected areas. Can help 
        visualise land-sea borders.

    area_linewidth: float
        default 0.3
        linewidth of area edges. Larger values have thicker borders.

    Returns
    -------
    matplotlib.axes.Axes
    """

    if select_area is None:
        ax.set_extent([xlim[0], xlim[1], ylim[0], ylim[1]], crs=ccrs.PlateCarree())
        pass
    else:
        assert isinstance(select_area, list), "select_area must be a list"
        # select state
        name_column = [name for name in regions.columns if "NAME" in name.upper()][0]
        area = regions.loc[regions[name_column].isin(select_area)]
        area= area.to_crs(crs = "GDA2020")
        map_total_bounds = area.total_bounds
        minx, miny, maxx, maxy = map_total_bounds
        mid_x = (minx + maxx) / 2
        mid_y = (miny + maxy) / 2
        max_range = np.max([(maxy - miny), (maxx - minx)])
        buffer = 0.1 * max_range        
    
        not_area = gpd.GeoSeries(
            data=[
                box(*box(*map_total_bounds).buffer(10 * buffer).bounds).difference(
                    area.dissolve()["geometry"].values[0]
                )
            ],
            crs=ccrs.PlateCarree(),
        )
    
        # mask outside selected area
        if land_shadow:
            # show land as light grey
            ax.add_geometries(not_area,
                              crs=ccrs.PlateCarree(),
                              facecolor="lightgrey",
                              linewidth=area_linewidth,
                              edgecolor="k", 
                              zorder=4)
        else:
            # mask white
            ax.add_geometries(not_area,
                              crs=ccrs.PlateCarree(),
                              facecolor="white",
                              linewidth=area_linewidth, 
                              edgecolor="k", 
                              zorder=4)
    
        ax.set_extent([mid_x - 0.6 * max_range,
                       mid_x + 0.8 * max_range,
                       mid_y - 0.7 * max_range,
                       mid_y + 0.7 * max_range],
                      crs=ccrs.PlateCarree())
    return ax

def plot_titles(title="title",
                date_range = "DD Mon YYYY to DD Mon YYYY", 
                baseline = None, 
                dataset_name= None,
                issued_date=None,
                watermark= None, 
                watermark_color="r",
                ax=None,
                text_xy = None,
                title_ha = "left",
               orientation = "none",):
    """
    Set the plot title and axis labels
    
    Parameters
    ----------
    title: str
        Main text. Size 14, bold. Location set by text_xy["title"]. 
        Horizontal alignment set by title_ha.
        
    date_range: str
        Text under title. Size 12. Horizontal alignment set by title_ha.
        Intended for data date range. Also can be used for any subtitle.
        Default "DD Mon YYYY to DD Mon YYYY" to indicate prefered date format.
        
    baseline: str
        Text in bottom left corner. Size 8. 
        Intended to describe the baseline period of the data.
        If None, then no text is printed.
        Default is None. 
        
    dataset_name: str
        Text inside bottom right. Size 8.
        Intended to describe data source.
        If None, then no text is printed.
        Default is None.
        
    issued_date: str
        Text on bottom right under the border.
        The date that the plot is current/valid/produced.
        If None (default), then today's date is used.
        To suppress any text, use issued_date="".
        Default is None.
        
    watermark: str
        Large text over plot. Use to indicate draft figures etc.
        Default is None. 
        
    watermark_color: color
        Option to change watermark text colour from red,
        eg if figure colours are red and then you can't read the watermark.
        default is "r" for red text.
        
    ax: matplotlib.axes.Axes
        axes
        
    text_xy: dictionary 
        Expects a dictionary with "title", "date_range" and "watermark" keys
        dictionary items with tuples describing text locations.
        Can omit "watermark" if watermark is None.

    title_ha: {"left", "center", "right"}
        Title horizontal alignment.
        Default "left"
    
    orientation: {"vertical", "horizontal"}
        Orientation of figures, used to control layout of Copywrite text.
        I.e. vertically stacked, or horizontal side-by-side
        
    Returns
    -------
    matplotlib.axes.Axes with text for titles etc.
    
    """
    
    ax.text(
        x=text_xy["title"][0],
        y=text_xy["title"][1],
        s=f"{title}",
        fontsize=14,
        weight="bold",
        horizontalalignment=title_ha,
        verticalalignment="baseline",
        transform=ax.transAxes,
        zorder=10,
        wrap=True,
    )

    ax.text(
        x=text_xy["date_range"][0],
        y=text_xy["date_range"][1],
        s=f"{date_range}",
        fontsize=12,
        horizontalalignment=title_ha,
        verticalalignment="top",
        transform=ax.transAxes,
        zorder=10,
        wrap=True,
    )
    
    if baseline is not None:
        # print base period inside bottom left corner
        ax.text(
            x=0.01,
            y=0.01,
            s=f"Base period: {baseline}",
            fontsize=8,
            verticalalignment="baseline",
            transform=ax.transAxes,
            zorder=10,
        )
    # print copyright outside bottom left corner
    if orientation == "vertical":
        string = f"\u00A9 Commonwealth of Australia {datetime.datetime.now().year}, \nAustralian Climate Service"
    else:
        string = f"\u00A9 Commonwealth of Australia {datetime.datetime.now().year}, Australian Climate Service"
    ax.text(
        x=0.01,
        y=0.00,
        s=string,
        fontsize=8,
        transform=ax.transAxes,
        verticalalignment="bottom",
        zorder=10,
        wrap=True,
    )
    if dataset_name is not None:
        # print data source inside bottom right
        ax.text(
            x=0.99,
            y=0.01,
            s=f"Dataset: {dataset_name}",
            fontsize=8,
            transform=ax.transAxes,
            horizontalalignment="right",
            verticalalignment="baseline",
            zorder=10,
        )
    # print issued date on bottom right under the border.
    # Set to today's date if None supplied
    # Suppress this by issued_date=""
    if issued_date is None:
        issued_date = datetime.datetime.today().date().strftime("%d %B %Y")
    if len(issued_date)>=1:
        ax.text(
            x=0.99,
            y=-0.01,
            s=f"Issued: {issued_date}",
            fontsize=8,
            transform=ax.transAxes,
            horizontalalignment="right",
            verticalalignment="top",
            zorder=10,
        )
    
    if watermark is not None:
        ax.text(
            x=text_xy["watermark"][0],
            y=text_xy["watermark"][1],
            s=watermark.upper(),
            fontsize=36,
            transform=ax.transAxes,
            horizontalalignment="center",
            verticalalignment="center",
            zorder=10,
            wrap=True,
            alpha=0.5,
            color=watermark_color,
        )
    ax.axis('off')
    return ax



# # Define a function for plotting maps
# This is the function you call to plot all the graphs
def plot_acs_hazard(
    name='ncra_regions',
    regions=regions_dict['ncra_regions'],
    data=None,
    station_df=None,
    markersize=None,
    stippling=None,
    mask_not_australia=True,
    mask_australia=False,
    agcd_mask=False,
    facecolor="none",
    edgecolor="black",
    figsize=(4.5,3),
    title=None,
    date_range="",
    projection=None,
    area_linewidth=0.3,
    xlim=(114, 162),
    ylim=(-43, -8),
    cmap=cm.Greens,
    cmap_bad="lightgrey",
    cbar_extend="both",
    ticks=None,
    tick_labels=None,
    cbar_label="",
    baseline=None,
    dataset_name=None,
    issued_date=None,
    label_states=False,
    contourf=False,
    contour=False,
    select_area=None,
    land_shadow=False,
    watermark="EXPERIMENTAL\nIMAGE ONLY",
    watermark_color = "r",
    show_logo = False,
    infile=None,
    outfile=None,
    savefig=True,
    tick_rotation=None,
    vcentre=None,
):
    """This function takes gridded data or station_df dataframe and shapefiles 
    to visualise hazard data from a 2D Xarray data array
    and plots the data on a map of Australia with the shape outlines.

    Parameters
    ----------
    name: str
        name of a shapefile collection in 
        /g/data/ia39/aus-ref-clim-data-nci/shapefiles/data/
        to get regions from.

    regions: geopandas.GeoDataFrame
        if None, then will try to read from regions_dict[{name}].

    data: xr.DataArray
        a 2D xarray DataArray that has already computed the 
        average, sum, anomaly, metric or index you wish to visualise.
        This function is resolution agnostic.

    station_df: pd.DataFrame, optional
        a pandas.DataFrame with columns ["lon", "lat", variable]. 
        If station_df is given, then variable values are represented as dots on 
        the map according to the lat lon coordinates and coloured according to
        cmap colours and ticks.

    markersize: int, optional
        default None. If None then the markersize will adjust to the size of the
        figure and the number of markers in the plot such that when there are
        many markers and the figure is small, the markersize is smaller.

    stippling: xr.DataArray
        a True/False mask to define regions of stippling hatching. 
        Intended to show model agreement, eg for the direction of change.

    mask_not_australia: boolean
        decides whether or not the area outside of Australian land is hidden 
        under white shape.
        Default is True.

    mask_australia: boolean
        decides whether or not Australian land is hidden under white shape.
        Eg, use when plotting ocean only.
        Default is False.

    agcd_mask: boolean
        If True, applies a ~5km mask for data-sparse inland areas of Australia.
        Default is False.

    facecolor: color
        color of land when plotting the regions without climate data. 
        facecolor recommendations include "white", "lightgrey", "none".

    edgecolor: color
        defines the color of the state/region borders. 
        edgecolor recommendations include "black" and "white".

    figsize: tuple
        defines the width and height of the figure in inches.
        Reccommend (8,6) for Australia-wide plots and (6,6) for individual states

    title: str
        defines the text inside the plot.
        If none is given, then will print the name of the shape file.

    date_range: str
        decribes the start and end date of the data analysed. 
        This is printed under the title. 
        format: dd Month yyyy to dd Month yyy.

    projection: optional
        formally spuriously named "crs". Defines the projection of the plots.
        defaults to a LambertConformal projection unless using "select_area".

    area_linewidth: float, optional
        the width of state/region borders only. All other linewidths are hardcoded.

    xlim: tuple, optional
        longitude min and max of the plot area.

    ylim: tuple, optional
        latitude min and max of the plot area.

    cmap:
        defines the colormap used for the data.
        See cmap_dict for suggested colormaps.
        If none, cmap set to cm.Greens.
        Please choose appropriate colormaps for your data.

    cmap_bad: color
        define the color to set for "bad" or missing values
        default "lightgrey"

    cbar_extend: one of {'neither', 'both', 'min', 'max'}.
        eg "both" changes the ends of the colorbar to arrows to indicate that
        values are possible outside the scale show.
        If contour or contourf is True, then cbar_extend will be overridden to "none".

    ticks: list or arraylike
        Define the ticks on the colorbar. Define any number of intervals. 
        This will make the color for each interval one discrete color, 
        instead of a smooth color gradient.
        If None, linear ticks will be auto-generated to fit the provided data.

    tick_labels: list
        Labels for categorical data. 
        If tick_labels is used, then pcolormesh is used to plot data 
        and does not allow contour or contourf to be used.
        Tick labels will correspond to the ticks.

    cbar_label: string
        defines the title for the color bar. 
        This should indicate the variable name and the units eg 
        "daily rainfall [mm]",
        "annual rainfall [mm]", 
        "monthly rainfall anomaly [mm]",
        "tas [\N{DEGREE SIGN}C]".

    baseline: string
        the baseline period for anomalies, eg "1961 - 1990".

    dataset_name: string
        describes the source of the data eg "AGCD v2" or "BARPA-R ACCESS-CM2"

    issued_date: string
        The date of issue. If None is supplied, then today's date is printed.

    label_states: bool
        if set to True and Australian states are the shapefile selected,
        then each state is labelled with its three-letter abbreviation. 
        Default False.

    contourf: bool
        if True then the gridded data is visualised as smoothed filled contours. 
        Default is False.
        Use with caution when plotting data with negative and positive values;
        Check output for NaNs and misaligned values.  

    contour: bool
        if True then the gridded data is visualised as smoothed unfilled grey contours.
        Default is True.
        Using both contourf and contour results in smooth filled contours
        with grey outlines between the color levels.

    select_area: list
        A list of areas (eg states) that are in the geopandas.GeoDataFrame.
        Inspect the regions gdf for area names. eg ["Victoria", "New South Wales"]

    land_shadow: bool
        Used when select_area is not None. 
        This option controls whether to show Australian land area that is outside 
        the select area in grey for visual context.
        Default False.

    watermark: str
        text over the plot for images not in their final form. 
        If the plot is in final form, set to None. 
        Suggestions include "PRELIMINARY DATA", "DRAFT ONLY", 
        "SAMPLE ONLY (NOT A FORECAST)", "EXPERIMENTAL IMAGE ONLY"

    watermark_color: default "r"
        for the watermark, this changes the colour of the text.
        The default is red. Only change color if red is not visible. 

    show_logo: bool, default False
        show ACS logo in top left corner.

    infile: str
        Not yet tested. 
        The idea is to read in 2D netCDF data and use this as the mappable data.

    outfile: str
        The location to save the figure. 
        If None, then figure is saved here f"figures/{title.replace(' ', '-')}.png"

    savefig: bool
        default is True
        If set to False, then fig is not saved.

    tick_rotation: int [-360,360]
        Angle to rotate colorbar tick labels.
        Default is None. Tick labels will be horizontal if colorbar is vertical,
        or vertical if colorbar is horizontal.
        Easiest to read id tick_rotation = 0
        
    vcentre: float, eg 0
        default is None.
        Align centre of colormap to this value. 
        Intended for using a divergent colormap with uneven number of ticks 
        around the centre, eg for future temperature anomalies with a larger
        positive range compared to the negative range.

    Returns
    -------
    The map is saved as a png in a "figures" file in your working directory.
    This function returns fig and ax.
    """
    
    if regions is None:
        try:
            regions = regions_dict[name]
        except:
            print(f"Could not read regions_dict[{name}]")

    regions = regions.to_crs(crs = "GDA2020")

    # Set default crs for Australia maps and selection maps
    if projection is None:
        if select_area is None:
            # Default for Australian map
            projection = ccrs.LambertConformal(
                central_latitude=-24.75,
                central_longitude=134.0,
                cutoff=30,
                standard_parallels=(-10, -40),
            )
        else:
            projection = ccrs.PlateCarree()

    # Set up the plot
    fig = plt.figure(
        figsize=figsize,
        zorder=1,
        layout="constrained",
        frameon=False,
    )
    ax = plt.axes(
        projection=projection,
        frameon=False,
    )
    ax.set_global()

    if infile is not None:
        data = xr.open_dataset(infile)


    if contour or contourf:
        cbar_extend = "neither"

    # plot hazard data ------------------------
    cmap.set_bad(cmap_bad)
    if station_df is not None and markersize is None:
        markersize=(100 - 80*len(station_df)/5000)*(figsize[0]*figsize[1])/48
    ax, norm, cont, middle_ticks =plot_data(regions=regions,
                                            data=data, 
                                            station_df = station_df,
                                            markersize=markersize,
                                            xlim=xlim,
                                            ylim=ylim,
                                            cmap=cmap,
                                            cbar_extend=cbar_extend,
                                            ticks=ticks,
                                            tick_labels=tick_labels,
                                            contourf=contourf,
                                            contour=contour,
                                            ax=ax,
                                            subtitle="",
                                            facecolor=facecolor,
                                            mask_not_australia = mask_not_australia,
                                            mask_australia=mask_australia,
                                            agcd_mask=agcd_mask,
                                            area_linewidth=area_linewidth,
                                            stippling=stippling,
                                            vcentre=vcentre,
                                            )
                    
    # ---------------------------------

    # if select a specific area -----------
    ax = plot_select_area(select_area=select_area, 
                          ax=ax,
                          xlim=xlim,
                          ylim=ylim,
                          regions=regions,
                          land_shadow=land_shadow,
                          area_linewidth=area_linewidth,
                         )
    # ---------------------------------------------

    # colorbar------------------------
    if data is not None:
        cbar = plot_cbar(cont=cont,
                      norm=norm,
                      ax=ax,
                      cbar_extend=cbar_extend, 
                      cbar_label=cbar_label,
                      location = "right",
                      ticks=ticks, 
                      tick_labels=tick_labels,
                      middle_ticks=middle_ticks,
                      cax_bounds = [0.82, 0.15, 0.03, 0.7],
                         rotation = tick_rotation,
                      )
    # ---------------------------------------

    # Annotations and titles ---------------------

    #plot border and annotations
    ax111 = fig.add_axes([0.,0.,1,1], xticks=[], yticks=[]) #(left, bottom, width, height)

    # text annotation xy locations for 1-panel plot
    text_xy_1pp = {"title": (0.03, 0.12),  
                   "date_range": (0.03, 0.11),
                   "watermark": (0.4, 0.5),}
    
    ax111 = plot_titles(title=title,
                        date_range = date_range, 
                        baseline = baseline, 
                        dataset_name= dataset_name,
                        issued_date=issued_date,
                        watermark=watermark, 
                        watermark_color=watermark_color,
                        ax=ax111,
                        text_xy = text_xy_1pp,
                        title_ha = "left",
                        )
    # draw border
    # ax111.axis(True)
    ax111.axis(False)

    # -----------------------------------------------


    if outfile is None:
        PATH = os.path.abspath(os.getcwd())
        outfile = f"{PATH}/figures/{title.replace(' ', '-')}.png"
        os.makedirs(os.path.dirname(outfile), exist_ok=True)

    if savefig:
        plt.savefig(outfile, dpi=300,)
    return fig, ax


# # Define a function for plotting maps
# This is the function you call to plot all the graphs
def plot_acs_hazard_3pp(
    name='ncra_regions',
    regions=None,
    ds_gwl15=None,
    ds_gwl20=None,
    ds_gwl30=None,
    station_df_gwl15=None,
    station_df_gwl20=None,
    station_df_gwl30=None,
    stippling_gwl15=None,
    stippling_gwl20=None,
    stippling_gwl30=None,
    mask_not_australia=True,
    mask_australia=False,
    agcd_mask=False,
    facecolor="none",
    edgecolor="black",
    figsize=None,
    markersize=None,
    title=None,
    date_range="",
    subplot_titles=None,
    projection=None,
    area_linewidth=0.3,
    xlim=(114,154),
    ylim=(-43.5, -7.5),
    cmap=cm.Greens,
    cmap_bad="lightgrey",
    cbar_extend="both",
    ticks=None,
    tick_labels=None,
    cbar_label="",
    baseline=None,
    dataset_name=None,
    issued_date=None,
    label_states=False,
    contourf=False,
    contour=False,
    select_area=None,
    land_shadow=False,
    watermark="EXPERIMENTAL IMAGE ONLY",
    watermark_color = "r",
    show_logo = False,
    infile=None,
    outfile=None,
    savefig=True,
    orientation="horizontal",
    tick_rotation=None,
    vcentre=None,
):
    """
    Three panel plot with shared cmap and titles etc. 
    As with plot_acs_hazard, but takes three xarray data arrays:
    ds_gwl15, ds_gwl20, ds_gwl30. (left, middle and right) or 
    (top, middle and bottom) 
    This function is intended for plotting multiple Global Warming 
    Levels, but it will plot any valid data (xr.DataArrays or pd.DataFrames)

    Parameters
    ----------
    name: str
        name of a shapefile collection in 
        /g/data/ia39/aus-ref-clim-data-nci/shapefiles/data/
        to get regions from.
        
    regions: geopandas.GeoDataFrame
        if None, then will try to read from regions_dict[{name}].
        
    ds_gwl15: xr.DataArray
        The first DataArray to plot.
        Expects a 2D xarray DataArray that has already computed the 
        average, sum, anomaly, metric or index you wish to visualise.
        This function is resolution agnostic.
        
    ds_gwl20: xr.DataArray
        The second DataArray to plot.
        
    ds_gwl30: xr.DataArray
        The third DataArray to plot.
        
    station_df_gwl15: pd.DataFrame, optional
        The first pandas.DataFrame with columns ["lon", "lat", variable]. 
        If station_df_gwl15 is given, then variable values are represented as dots on 
        the map according to the lat lon coordinates and coloured according to
        cmap colours and ticks. Use markersize to change dot size.

    station_df_gwl20: pd.DataFrame, optional
        The second pandas.DataFrame with columns ["lon", "lat", variable].
        
    station_df_gwl30: pd.DataFrame, optional
        The third pandas.DataFrame with columns ["lon", "lat", variable]. 
        
    stippling_gwl15: xr.DataArray
        a True/False mask to define regions of stippling hatching 
        for the first subplot. 
        Intended to show model agreement, eg for the direction of change.
        
    stippling_gwl20: xr.DataArray
        a True/False mask to define regions of stippling hatching 
        for the second subplot. 
        
    stippling_gwl30: xr.DataArray
        a True/False mask to define regions of stippling hatching 
        for the third subplot. 

    mask_not_australia: boolean
        decides whether or not the area outside of Australian land is hidden 
        under white shape.
        Default is True.

    mask_australia: boolean
        decides whether or not Australian land is hidden under white shape.
        Eg, use when plotting ocean only.
        Default is False.

    agcd_mask: boolean
        If True, applies a ~5km mask for data-sparse inland areas of Australia.
        Default is False.

    facecolor: color
        color of land when plotting the regions without climate data. 
        facecolor recommendations include "white", "lightgrey", "none".
        Default is "none"

    edgecolor: color
        defines the color of the state/region borders. 
        edgecolor recommendations include "black" and "white".
        Default is "black"

    figsize: tuple
        defines the width and height of the figure in inches.
        ACS recommends a maximum width of 6.7" (17cm) and 
        maximum height of ~7.5" (19cm)

    markersize: optional 
        Markersize for station_df dots.
        default None. If None then the markersize will adjust to the size of the
        figure and the number of markers in the plot such that when there are
        many markers and the figure is small, the markersize is smaller.
        
    title: str
        A title should describe what is shown in the map. 
        The title should be written in plain English and 
        centred at the top of the visualization.
        If title is None, then defaults to the name of the shapefile.
        
    date_range: str
        date_range (or subtitle)
        Expected to decribe the start and end date of the data analysed. 
        This is printed under the title. 
        format: dd Month yyyy to dd Month yyy.
        Default=""
        
    subplot_titles: list of strings
        subplot_titles for labeling each subplot title
        if None, then subtitles are ["GWL1.5", "GWL2.0", "GWL3.0"]
        otherwise specify a list of three strings.
        
    projection:
        Specify projection of the maps. The default suits Australia.
        Formally "crs".
        If None, defaults to
        ccrs.LambertConformal(central_latitude=-24.75,
                              central_longitude=134.0,
                              cutoff=30,
                              standard_parallels=(-10, -40),
        unless select_area is not None, then defaults to
        ccrs.PlateCarree()
                            
    area_linewidth: float
        linewidth of state/region borders.
        Default =0.3
        
    xlim: tuple of floats
        longitude limits
        Default = (114,154)
        
    ylim: tuple of floats
        latitude limits
        Default = (-43.5, -7.5)
        
    cmap: matplotlib colormap
        color map for gridded and/or station data
        See cmap_dict for suggested colormaps.
        If none, cmap set to cm.Greens.
        Please choose appropriate colormaps for your data.

    cmap_bad: color
        define the color to set for "bad" or missing values
        default "lightgrey"
        
    cbar_extend: one of {'neither', 'both', 'min', 'max'}.
        eg "both" changes the ends of the colorbar to arrows to indicate that
        values are possible outside the scale show.
        If contour or contourf is True, then cbar_extend will be overridden to "none".
        
    ticks: list or arraylike
        Define the ticks on the colorbar. Define any number of intervals. 
        This will make the color for each interval one discrete color, 
        instead of a smooth color gradient.
        If None, linear ticks will be auto-generated to fit the provided data.

    tick_labels: list
        Labels for categorical data. 
        If tick_labels is used, then pcolormesh is used to plot data 
        and does not allow contour or contourf to be used.
        Tick labels will correspond to the ticks.
        
    cbar_label: string
        defines the title for the color bar. 
        This should indicate the variable name and the units eg 
        "daily rainfall [mm]",
        "annual rainfall [mm]", 
        "monthly rainfall anomaly [mm]",
        "tas [\N{DEGREE SIGN}C]".
        Default is ""
        
    baseline: string
        the baseline period for anomalies, eg "1961 - 1990".
        
    dataset_name: string
        describes the source of the data eg "AGCD v2" or "BARPA-R ACCESS-CM2"
        
    issued_date: string
        The date of issue. If None is supplied, then today's date is printed.
        To supress, set to ""
        
    label_states: bool
        if set to True and Australian states are the shapefile selected,
        then each state is labelled with its three-letter abbreviation. 
        Default False.
        
    contourf: bool
        if True then the gridded data is visualised as smoothed filled contours. 
        Default is False.
        Use with caution when plotting data with negative and positive values;
        Check output for NaNs and misaligned values.  

    contour: bool
        if True then the gridded data is visualised as smoothed unfilled grey contours.
        Default is True.
        Using both contourf and contour results in smooth filled contours
        with grey outlines between the color levels.

    select_area: list
        A list of areas (eg states) that are in the geopandas.GeoDataFrame.
        Inspect the regions gdf for area names. eg ["Victoria", "New South Wales"]

    land_shadow: bool
        Used when select_area is not None. 
        This option controls whether to show Australian land area that is outside 
        the select area in grey for visual context.
        Default False.

    watermark: str
        text over the plot for images not in their final form. 
        If the plot is in final form, set to None. 
        Suggestions include "PRELIMINARY DATA", "DRAFT ONLY", 
        "SAMPLE ONLY (NOT A FORECAST)", "EXPERIMENTAL IMAGE ONLY"

    watermark_color: default "r"
        for the watermark, this changes the colour of the text.
        The default is red. Only change color if red is not visible. 

    show_logo: bool, default False
        show ACS logo in top left corner.

    infile: str
        Not yet tested. 
        The idea is to read in 2D netCDF data and use this as the mappable data.

    outfile: str
        The location to save the figure. 
        If None, then figure is saved here f"figures/{title.replace(' ', '-')}.png"

    savefig: bool
        default is True
        If set to False, then fig is not saved.

    orientation: {"horizontal", "vertical"}
        whether the three plots are orientatied in a vertical stack or
        horizontally.
        Default "horizontal"
 
    tick_rotation: int [-360,360]
        Angle to rotate colorbar tick labels.
        Default is None. Tick labels will be horizontal if colorbar is vertical,
        or vertical if colorbar is horizontal.
        Easiest to read if tick_rotation = 0
        
    vcentre: float, eg 0
        default is None.
        Align centre of colormap to this value. 
        Intended for using a divergent colormap with uneven number of ticks 
        around the centre, eg for future temperature anomalies with a larger
        positive range compared to the negative range.

    Returns
    -------
    A three panel plot saved as a png in a "figures" file in your working directory.
    This function returns fig and ax.
    """

    if orientation=="horizontal":
        cax_bounds = [1.05,0,0.1,1]
        if tick_rotation is None:
            tick_rotation = 0
        nrows = 1
        ncols = 3
        cbar_location = "right"
        plots_rect = (0.01,0.05,0.98,0.85) #left bottom width height
        # text annotation xy locations
        text_xy = {"title": (0.5, 0.9),
                   "date_range": (0.5, 0.87),
                   "watermark": (0.45, 0.41),}
        subtitle_xy = None
        if figsize is None:
            figsize=(10, 3)
            
    elif orientation=="vertical":
        cax_bounds = [-0.4,-0.3,1.6,0.1]
        if tick_rotation is None:
            tick_rotation = -90
        nrows = 3
        ncols = 1
        cbar_location = "bottom"
        plots_rect = (0.01,0.07,0.98,0.85) #left bottom width height
        # text annotation xy locations
        text_xy = {"title": (0.5, 0.94),
               "date_range": (0.5, 0.93),
               "watermark": (0.45, 0.41),}
        subtitle_xy = (-0.55, 0.6)
        if figsize is None:
            figsize=(3, 6.5)
        
    else:
        print('orientation must be one of ["horizontal", "vertical",]')

    
    if regions is None:
        try:
            regions = regions_dict[name]
        except:
            print(f"Could not read regions_dict[{name}]")

    regions = regions.to_crs(crs = "GDA2020")

    # Set default projection for Australia maps and selection maps
    if projection is None:
        if select_area is None:
            # Default for Australian map
            projection = ccrs.LambertConformal(
                central_latitude=-24.75,
                central_longitude=134.0,
                cutoff=30,
                standard_parallels=(-10, -40),
            )
        else:
            projection = ccrs.PlateCarree()
        
    fig, axs = plt.subplots(nrows=nrows, ncols=ncols, 
                            sharey=True, sharex=True,
                            figsize=figsize, 
                            layout="constrained",
                            subplot_kw={'projection': projection, "frame_on":False},)

    cmap.set_bad(cmap_bad)
    
    station_dfs = [station_df_gwl15, station_df_gwl20, station_df_gwl30]
    if any(df is not None for df in station_dfs) and markersize is None:
        markersize=(100 - 80*len(station_dfs[0])/5000)*(figsize[0]*figsize[1])/48/3

    if subplot_titles is None:
        subplot_titles = [f"GWL{[1.5,2.0,3.0][i]}" for i in range(3)]
        
    for i, ds in enumerate([ds_gwl15, ds_gwl20, ds_gwl30]):
        station_df = station_dfs[i]
        stippling = [stippling_gwl15, stippling_gwl20, stippling_gwl30][i]
        subtitle = subplot_titles[i]
        ax, norm, cont, middle_ticks = plot_data(regions=regions,
                                              data=ds, 
                                              station_df = station_df,
                                              markersize=markersize,
                                              xlim=xlim,
                                              ylim=ylim,
                                              cmap=cmap,
                                              cbar_extend=cbar_extend,
                                              ticks=ticks,
                                              tick_labels=tick_labels,
                                              contourf=contourf,
                                              contour=contour,
                                              ax=axs[i],
                                              figsize=figsize,
                                              subtitle=subtitle,
                                                 subtitle_xy=subtitle_xy,
                                              facecolor=facecolor,
                                              mask_not_australia = mask_not_australia,
                                              mask_australia=mask_australia,
                                              agcd_mask=agcd_mask,
                                              area_linewidth=area_linewidth,
                                              stippling=stippling,
                                                vcentre=vcentre,
                                                )
        
        # if select a specific area -----------
        ax = plot_select_area(select_area=select_area, 
                              ax=ax,
                              xlim=xlim,
                              ylim=ylim,
                              regions=regions,
                              land_shadow=land_shadow,
                              area_linewidth=area_linewidth,
                              )
        # ---------------------------------------------

                    
        ax.axis('off')

    # colorbar -----------------------------------------------------------
    # fig.subplots_adjust(left=0.05, bottom=0, right=0.85, top=0.95, wspace=0.05, hspace=0.05)
    # cbar_ax = fig.add_axes([0.87, 0.2, 0.03, 0.5]) #left bottom width height
    # cbar_ax.axis('off')
    fig.get_layout_engine().set(rect=plots_rect)
    cbar_ax = axs.flatten()[-1]

    cbar = plot_cbar(cont=cont,
                  norm=norm,
                  ax=cbar_ax,
                  cbar_extend=cbar_extend, 
                  cbar_label=cbar_label,
                  ticks=ticks, 
                  tick_labels=tick_labels,
                  middle_ticks=middle_ticks,
                  cax_bounds = cax_bounds,
                     rotation=tick_rotation,
                     location = cbar_location,
                  )
    #------------------------------------------


    # plot border and annotations -----------------
    ax111 = fig.add_axes([0.01,0.01,0.98,0.98], xticks=[], yticks=[]) #(left, bottom, width, height)
    
    ax111 = plot_titles(title=title,
                        date_range = date_range, 
                        baseline = baseline, 
                        dataset_name= dataset_name,
                        issued_date=issued_date,
                        watermark=watermark, 
                        watermark_color=watermark_color,
                        ax=ax111,
                        text_xy = text_xy,
                        title_ha = "center",
                        orientation = orientation,
                   )
    # draw border
    # ax111.axis(True)
    ax111.axis(False)
    # --------------------------------------------

    if outfile is None:
        PATH = os.path.abspath(os.getcwd())
        outfile = f"{PATH}/figures/{title.replace(' ', '-')}.png"
        os.makedirs(os.path.dirname(outfile), exist_ok=True)

    if savefig:
        plt.savefig(outfile, dpi=300,)
    return fig, ax

def plot_acs_hazard_4pp(
                name='ncra_regions',
                regions=None,
                ds_gwl12=None,
                ds_gwl15=None,
                ds_gwl20=None,
                ds_gwl30=None,
                station_df_gwl12=None,                    
                station_df_gwl15=None,
                station_df_gwl20=None,
                station_df_gwl30=None,
                stippling_gwl12=None,
                stippling_gwl15=None,
                stippling_gwl20=None,
                stippling_gwl30=None,
                mask_not_australia=True,
                mask_australia=False,
                agcd_mask=False,
                facecolor="none",
                edgecolor="black",
                figsize=None,
                markersize=None,
                title=None,
                date_range="",
                subplot_titles=None,
                projection=None,
                area_linewidth=0.3,
                xlim=(113, 154),
                ylim=(-43.5, -9.5),
                cmap=cm.Greens,
                cmap_bad="lightgrey",
                cbar_extend="both",
                ticks=None,
                tick_labels=None,
                cbar_label="",
                baseline=None,
                dataset_name=None,
                issued_date=None,
                label_states=False,
                contourf=False,
                contour=False,
                select_area=None,
                land_shadow=False,
                watermark="EXPERIMENTAL\nIMAGE ONLY",
                watermark_color = "r",
                show_logo = False,
                infile=None,
                outfile=None,
                savefig=True,
                orientation="vertical",
                tick_rotation=None,
                vcentre=None,
            ):
    """
    Four panel plot with shared cmap and titles etc. 
    As with plot_acs_hazard, but takes four xarray data arrays:
    ds_gwl12, ds_gwl15, ds_gwl20, ds_gwl30. (left to right), 
    (top-left, top-right, bottom-left, bottom-right), or 
    (top to bottom).
    This function is intended for plotting multiple Global Warming 
    Levels, but it will plot any valid data (xr.DataArrays or pd.DataFrames)

    Parameters
    ----------
    name: str
        name of a shapefile collection in 
        /g/data/ia39/aus-ref-clim-data-nci/shapefiles/data/
        to get regions from.
        
    regions: geopandas.GeoDataFrame
        if None, then will try to read from regions_dict[{name}].

    ds_gwl12: xr.DataArray
        The first DataArray to plot.
        Expects a 2D xarray DataArray that has already computed the 
        average, sum, anomaly, metric or index you wish to visualise.
        This function is resolution agnostic.
        
    ds_gwl15: xr.DataArray
        The second DataArray to plot.

    ds_gwl20: xr.DataArray
        The third DataArray to plot.
        
    ds_gwl30: xr.DataArray
        The fourth DataArray to plot.

    station_df_gwl12: pd.DataFrame, optional
        The first pandas.DataFrame with columns ["lon", "lat", variable]. 
        If station_df_gwl12 is given, then variable values are represented as dots on 
        the map according to the lat lon coordinates and coloured according to
        cmap colours and ticks. Use markersize to change dot size.
        
    station_df_gwl15: pd.DataFrame, optional
        The second pandas.DataFrame with columns ["lon", "lat", variable]. 

    station_df_gwl20: pd.DataFrame, optional
        The third pandas.DataFrame with columns ["lon", "lat", variable].
        
    station_df_gwl30: pd.DataFrame, optional
        The fourth pandas.DataFrame with columns ["lon", "lat", variable]. 

    stippling_gwl12: xr.DataArray
        a True/False mask to define regions of stippling hatching 
        for the first subplot. 
        Intended to show model agreement, eg for the direction of change.
        
    stippling_gwl15: xr.DataArray
        a True/False mask to define regions of stippling hatching 
        for the second subplot. 
        
    stippling_gwl20: xr.DataArray
        a True/False mask to define regions of stippling hatching 
        for the third subplot. 
        
    stippling_gwl30: xr.DataArray
        a True/False mask to define regions of stippling hatching 
        for the fourth subplot. 

    mask_not_australia: boolean
        decides whether or not the area outside of Australian land is hidden 
        under white shape.
        Default is True.

    mask_australia: boolean
        decides whether or not Australian land is hidden under white shape.
        Eg, use when plotting ocean only.
        Default is False.

    agcd_mask: boolean
        If True, applies a ~5km mask for data-sparse inland areas of Australia.
        Default is False.

    facecolor: color
        color of land when plotting the regions without climate data. 
        facecolor recommendations include "white", "lightgrey", "none".
        Default is "none"

    edgecolor: color
        defines the color of the state/region borders. 
        edgecolor recommendations include "black" and "white".
        Default is "black"

    figsize: tuple
        defines the width and height of the figure in inches.
        ACS recommends a maximum width of 6.7" (17cm) and 
        maximum height of ~7.5" (19cm)

    markersize: optional 
        Markersize for station_df dots.
        default None. If None then the markersize will adjust to the size of the
        figure and the number of markers in the plot such that when there are
        many markers and the figure is small, the markersize is smaller.
        
    title: str
        A title should describe what is shown in the map. 
        The title should be written in plain English and 
        centred at the top of the visualization.
        If title is None, then defaults to the name of the shapefile.
        
    date_range: str
        date_range (or subtitle)
        Expected to decribe the start and end date of the data analysed. 
        This is printed under the title. 
        format: dd Month yyyy to dd Month yyy.
        Default=""
        
    subplot_titles: list of strings
        subplot_titles for labeling each subplot title
        if None, then subtitles are ["GWL1.5", "GWL2.0", "GWL3.0"]
        otherwise specify a list of three strings.
        
    projection:
        Specify projection of the maps. The default suits Australia.
        Formally "crs".
        If None, defaults to
        ccrs.LambertConformal(central_latitude=-24.75,
                              central_longitude=134.0,
                              cutoff=30,
                              standard_parallels=(-10, -40),
        unless select_area is not None, then defaults to
        ccrs.PlateCarree()
                            
    area_linewidth: float
        linewidth of state/region borders.
        Default =0.3
        
    xlim: tuple of floats
        longitude limits
        Default = (113,154)
        
    ylim: tuple of floats
        latitude limits
        Default = (-43.5, -9.5)
        
    cmap: matplotlib colormap
        color map for gridded and/or station data
        See cmap_dict for suggested colormaps.
        Default cmap set to cm.Greens.
        Please choose appropriate colormaps for your data.

    cmap_bad: color
        define the color to set for "bad" or missing values
        default "lightgrey"
        
    cbar_extend: one of {'neither', 'both', 'min', 'max'}.
        eg "both" changes the ends of the colorbar to arrows to indicate that
        values are possible outside the scale show.
        If contour or contourf is True, then cbar_extend will be overridden to "none".
        Default is "both"
        
    ticks: list or arraylike
        Define the ticks on the colorbar. Define any number of intervals. 
        This will make the color for each interval one discrete color, 
        instead of a smooth color gradient.
        If None, linear ticks will be auto-generated to fit the provided data.

    tick_labels: list
        Labels for categorical data. 
        If tick_labels is used, then pcolormesh is used to plot data 
        and does not allow contour or contourf to be used.
        Tick labels will correspond to the ticks.
        
    cbar_label: string
        defines the title for the color bar. 
        This should indicate the variable name and the units eg 
        "daily rainfall [mm]",
        "annual rainfall [mm]", 
        "monthly rainfall anomaly [mm]",
        "tas [\N{DEGREE SIGN}C]".
        Default is ""
        
    baseline: string
        the baseline period for anomalies, eg "1961 - 1990".
        
    dataset_name: string
        describes the source of the data eg "AGCD v2" or "BARPA-R ACCESS-CM2"
        
    issued_date: string
        The date of issue. If None is supplied, then today's date is printed.
        To supress, set to ""
        
    label_states: bool
        if set to True and Australian states are the shapefile selected,
        then each state is labelled with its three-letter abbreviation. 
        Default False.
        
    contourf: bool
        if True then the gridded data is visualised as smoothed filled contours. 
        Default is False.
        Use with caution when plotting data with negative and positive values;
        Check output for NaNs and misaligned values.  

    contour: bool
        if True then the gridded data is visualised as smoothed unfilled grey contours.
        Default is True.
        Using both contourf and contour results in smooth filled contours
        with grey outlines between the color levels.

    select_area: list
        A list of areas (eg states) that are in the geopandas.GeoDataFrame.
        Inspect the regions gdf for area names. eg ["Victoria", "New South Wales"]

    land_shadow: bool
        Used when select_area is not None. 
        This option controls whether to show Australian land area that is outside 
        the select area in grey for visual context.
        Default False.

    watermark: str
        text over the plot for images not in their final form. 
        If the plot is in final form, set to None. 
        Suggestions include "PRELIMINARY DATA", "DRAFT ONLY", 
        "SAMPLE ONLY (NOT A FORECAST)", "EXPERIMENTAL IMAGE ONLY"
        default "EXPERIMENTAL\nIMAGE ONLY"

    watermark_color: default "r"
        for the watermark, this changes the colour of the text.
        The default is red. Only change color if red is not visible. 

    show_logo: bool, default False
        show ACS logo in top left corner.

    infile: str
        Not yet tested. 
        The idea is to read in 2D netCDF data and use this as the mappable data.

    outfile: str
        The location to save the figure. 
        If None, then figure is saved here f"figures/{title.replace(' ', '-')}.png"

    savefig: bool
        default is True
        If set to False, then fig is not saved.

    orientation: {"horizontal", "vertical", "square"}
        whether the four plots are orientatied in a vertical stack,
        horizontally, or in a 2-by-2 grid (left to right).
        Default "horizontal"
 
    tick_rotation: int [-360,360]
        Angle to rotate colorbar tick labels.
        Default is None. Tick labels will be horizontal if colorbar is vertical,
        or vertical if colorbar is horizontal.
        Easiest to read if tick_rotation = 0
        
    vcentre: float, eg 0
        default is None.
        Align centre of colormap to this value. 
        Intended for using a divergent colormap with uneven number of ticks 
        around the centre, eg for future temperature anomalies with a larger
        positive range compared to the negative range.

    Returns
    -------
    A four panel plot saved as a png in a "figures" file in your working directory.
    This function returns fig and ax.
    """

    if orientation=="horizontal":
        cax_bounds = [1.05,0,0.1,1]
        if tick_rotation is None:
            tick_rotation = 0
        nrows = 1
        ncols = 4
        cbar_location = "right"
        plots_rect = (0.01,0.05,0.98,0.85) #left bottom width height
        # text annotation xy locations
        text_xy = {"title": (0.5, 0.9),
                   "date_range": (0.5, 0.87),
                   "watermark": (0.45, 0.41),}
        subtitle_xy = None
        if figsize is None:
            figsize=(10, 3)
            
    elif orientation=="vertical":
        cax_bounds = [-0.4,-0.3,1.6,0.1]
        if tick_rotation is None:
            tick_rotation = -90
        nrows = 4
        ncols = 1
        cbar_location = "bottom"
        plots_rect = (0.01,0.05,0.98,0.85) #left bottom width height
        # text annotation xy locations
        text_xy = {"title": (0.5, 0.94),
               "date_range": (0.5, 0.93),
               "watermark": (0.45, 0.41),}
        subtitle_xy = (-0.54, 0.6)
        if figsize is None:
            figsize=(3, 8)       
        
    elif orientation=="square":   
        if tick_rotation is None:
            tick_rotation = -90
        ncols=2
        nrows=2
        cbar_location = "right"
        figsize=(6, 4.5)
        plots_rect = (0.01,0.05,0.84,0.85) #left bottom width height
        cax_bounds = [0.1,0,0.5,1]
        cbar_rect = [0.85, 0.2, 0.03, 0.5]

        # text annotation xy locations for 3-panel plot
        text_xy = {"title": (0.5, 0.93),
                   "date_range": (0.5, 0.92),
                   "watermark": (0.45, 0.41),}
        subtitle_xy = None

        if figsize is None:
            figsize=(6, 4.5)
    else:
        print('orientation must be one of ["horizontal", "vertical", "square"]')

    if regions is None:
        try:
            regions = regions_dict[name]
        except:
            print(f"Could not read regions_dict[{name}]")

    regions = regions.to_crs(crs = "GDA2020")

    # Set default projection for Australia maps and selection maps
    if projection is None:
        if select_area is None:
            # Default for Australian map
            projection = ccrs.LambertConformal(
                central_latitude=-24.75,
                central_longitude=134.0,
                cutoff=30,
                standard_parallels=(-10, -40),
            )
        else:
            projection = ccrs.PlateCarree()

    fig, axs = plt.subplots(nrows=nrows, ncols=ncols,
                            sharey=True, sharex=True, 
                            figsize=figsize,
                            layout="constrained",
                            subplot_kw={'projection': projection, "frame_on":False},)
    
    cmap.set_bad(cmap_bad)
    
    station_dfs = [station_df_gwl12, station_df_gwl15, station_df_gwl20, station_df_gwl30]
    if any(df is not None for df in station_dfs) and markersize is None:
        markersize=(100 - 80*len(station_dfs[0])/5000)*(figsize[0]*figsize[1])/48/4

    if subplot_titles is None:
        subplot_titles = [f"GWL{[1.2, 1.5, 2.0, 3.0][i]}" for i in range(4)]
        
    for i, ds in enumerate([ds_gwl12, ds_gwl15, ds_gwl20, ds_gwl30]):
        station_df = station_dfs[i]
        stippling = [stippling_gwl12, stippling_gwl15, stippling_gwl20, stippling_gwl30][i]
        subtitle = subplot_titles[i]
        ax, norm, cont, middle_ticks = plot_data(regions=regions,
                                              data=ds,
                                                 station_df=station_df,
                                                 markersize=markersize,
                                              xlim=xlim,
                                              ylim=ylim,
                                              cmap=cmap,
                                              cbar_extend=cbar_extend,
                                              ticks=ticks,
                                              tick_labels=tick_labels,
                                              contourf=contourf,
                                              contour=contour,
                                              ax=axs.flatten()[i],
                                              subtitle=subtitle,
                                                 subtitle_xy=subtitle_xy,
                                              facecolor=facecolor,
                                              mask_not_australia = mask_not_australia,
                                              mask_australia=mask_australia,
                                              agcd_mask=agcd_mask,
                                              area_linewidth=area_linewidth,
                                              stippling=stippling,
                                                vcentre=vcentre,)
        
        # if select a specific area -----------
        ax = plot_select_area(select_area=select_area, 
                              ax=ax,
                              xlim=xlim,
                              ylim=ylim,
                              regions=regions,
                              land_shadow=land_shadow,
                              area_linewidth=area_linewidth,
                              )
        # ---------------------------------------------
    
                    
        ax.axis('off')
    
    # colorbar -----------------------------------------------------------
    fig.get_layout_engine().set(rect=plots_rect)
    
    if orientation == "square":
        cbar_ax = fig.add_axes(cbar_rect) #left bottom width height
        cbar_ax.axis('off')
    else:
        cbar_ax = axs.flatten()[-1]
    
    cbar = plot_cbar(cont=cont,
                     norm=norm,
                     ax=cbar_ax,
                     cbar_extend=cbar_extend, 
                     cbar_label=cbar_label,
                     location=cbar_location,
                     ticks=ticks, 
                     tick_labels=tick_labels,
                     middle_ticks=middle_ticks,
                     cax_bounds=cax_bounds,
                     rotation = tick_rotation,
                     )    
    #------------------------------------------
    
    
    # plot border and annotations -----------------
    ax111 = fig.add_axes([0.01,0.01,0.98,0.98], xticks=[], yticks=[]) #(left, bottom, width, height)
    
    
    ax111 = plot_titles(title=title,
                        date_range = date_range, 
                        baseline = baseline, 
                        dataset_name= dataset_name,
                        issued_date=issued_date,
                        watermark=watermark, 
                        watermark_color=watermark_color,
                        ax=ax111,
                        text_xy = text_xy,
                        title_ha = "center",
                        orientation=orientation,
                   )
    # draw border
    # ax111.axis(True)
    ax111.axis(False)
    # --------------------------------------------
    
    if outfile is None:
        PATH = os.path.abspath(os.getcwd())
        outfile = f"{PATH}/figures/{title.replace(' ', '-')}.png"
        os.makedirs(os.path.dirname(outfile), exist_ok=True)
    
    if savefig:
        plt.savefig(outfile, dpi=300,)
    return fig, ax

def plot_acs_hazard_1plus3(
                name='ncra_regions',
                regions=None,
                ds_gwl12=None,
                station_df_gwl12=None, 
                stippling_gwl12=None,
                gwl12_cmap=cm.Greens,
                gwl12_cbar_extend="both",
                gwl12_cbar_label=None,
                gwl12_ticks=None,
                gwl12_tick_labels=None,
                gwl12_tick_rotation=None,
                gwl12_vcentre=None,
                ds_gwl15=None,
                ds_gwl20=None,
                ds_gwl30=None,                      
                station_df_gwl15=None,
                station_df_gwl20=None,
                station_df_gwl30=None,
                stippling_gwl15=None,
                stippling_gwl20=None,
                stippling_gwl30=None,
                mask_not_australia=True,
                mask_australia=False,
                agcd_mask=False,
                facecolor="none",
                edgecolor="black",
                figsize=None,
                markersize=None,
                title=None,
                date_range="",
                subplot_titles=None,
                projection=None,
                area_linewidth=0.3,
                xlim=(113, 154),
                ylim=(-43.5, -9.5),
                cmap=cm.Greens,
                cmap_bad="lightgrey",
                cbar_extend="both",
                ticks=None,
                tick_labels=None,
                cbar_label="",
                baseline=None,
                dataset_name=None,
                issued_date=None,
                label_states=False,
                contourf=False,
                contour=False,
                select_area=None,
                land_shadow=False,
                watermark="EXPERIMENTAL\nIMAGE ONLY",
                watermark_color = "r",
                show_logo = False,
                infile=None,
                outfile=None,
                savefig=True,
                orientation="vertical",
                tick_rotation=None,
                vcentre=None,
            ):
    """
    Four panel plot with 1 baseline plot and 3 future scenario 
    anomaly plots. The first plot has its own cmap and normalisation.
    The last three  plots share colormap and normalisation.
    
    As with plot_acs_hazard, but takes four xarray data arrays:
    ds_gwl12, ds_gwl15, ds_gwl20, ds_gwl30. (left to right), 
    (top-left, top-right, bottom-left, bottom-right), or 
    (top to bottom).
    This function is intended for plotting multiple Global Warming 
    Levels, but it will plot any valid data (xr.DataArrays or pd.DataFrames)

    Parameters
    ----------
    name: str
        name of a shapefile collection in 
        /g/data/ia39/aus-ref-clim-data-nci/shapefiles/data/
        to get regions from.
        default "ncra_regions"
        
    regions: geopandas.GeoDataFrame
        if None, then will try to read from regions_dict[{name}].

    ds_gwl12: xr.DataArray
        The first DataArray to plot.
        Expects a 2D xarray DataArray that has already computed the 
        average, sum, anomaly, metric or index you wish to visualise.
        This function is resolution agnostic.

    station_df_gwl12: pd.DataFrame, optional
        The first pandas.DataFrame with columns ["lon", "lat", variable]. 
        If station_df_gwl12 is given, then variable values are represented as dots on 
        the map according to the lat lon coordinates and coloured according to
        cmap colours and ticks. Use markersize to change dot size.
    
    stippling_gwl12: xr.DataArray
        a True/False mask to define regions of stippling hatching 
        for the first subplot. 
        Intended to show model agreement, eg for the direction of change.

    gwl12_cmap: matplotlib.colors.Colormap
        colormap for baseline plot
        default cm.Greens
        
    gwl12_cbar_extend: {"both", "neither", "min", "max"}
        arrows for colorbar for first (gwl12) plot.
        default "both".
        
    gwl12_cbar_label: str
        title for first (gwl12) colorbar.
        
    gwl12_ticks: list or np.array
        ticks for normalising colorbar for gwl12 data. 
        Used for both gridded (ds_gwl12) and station data (gwl12_station_df)
        
    gwl12_tick_labels: list
        tick labels for gwl12 data if categorical data
        
    gwl12_tick_rotation: [-360,360]
        degrees to rotate tick labels for gwl12 colorbar
        0 is horizontal
         
    gwl12_vcentre: float, optional
        Align centre of colormap to this value. 
        Intended for using a divergent colormap with uneven number of ticks 
        around the centre, eg for future temperature anomalies with a larger
        positive range compared to the negative range.
        
    ds_gwl15: xr.DataArray
        The second DataArray to plot.

    ds_gwl20: xr.DataArray
        The third DataArray to plot.
        
    ds_gwl30: xr.DataArray
        The fourth DataArray to plot.
        
    station_df_gwl15: pd.DataFrame, optional
        The second pandas.DataFrame with columns ["lon", "lat", variable]. 

    station_df_gwl20: pd.DataFrame, optional
        The third pandas.DataFrame with columns ["lon", "lat", variable].
        
    station_df_gwl30: pd.DataFrame, optional
        The fourth pandas.DataFrame with columns ["lon", "lat", variable]. 
        
    stippling_gwl15: xr.DataArray
        a True/False mask to define regions of stippling hatching 
        for the second subplot. 
        
    stippling_gwl20: xr.DataArray
        a True/False mask to define regions of stippling hatching 
        for the third subplot. 
        
    stippling_gwl30: xr.DataArray
        a True/False mask to define regions of stippling hatching 
        for the fourth subplot. 

    mask_not_australia: boolean
        decides whether or not the area outside of Australian land is hidden 
        under white shape.
        Default is True.

    mask_australia: boolean
        decides whether or not Australian land is hidden under white shape.
        Eg, use when plotting ocean only.
        Default is False.

    agcd_mask: boolean
        If True, applies a ~5km mask for data-sparse inland areas of Australia.
        Default is False.

    facecolor: color
        color of land when plotting the regions without climate data. 
        facecolor recommendations include "white", "lightgrey", "none".
        Default is "none"

    edgecolor: color
        defines the color of the state/region borders. 
        edgecolor recommendations include "black" and "white".
        Default is "black"

    figsize: tuple
        defines the width and height of the figure in inches.
        ACS recommends a maximum width of 6.7" (17cm) and 
        maximum height of ~7.5" (19cm)
        Defaults depend on "orientation"

    markersize: optional 
        Markersize for station_df dots.
        default None. If None then the markersize will adjust to the size of the
        figure and the number of markers in the plot such that when there are
        many markers and the figure is small, the markersize is smaller.
        
    title: str
        A title should describe what is shown in the map. 
        The title should be written in plain English and 
        centred at the top of the visualization.
        If title is None, then defaults to the name of the shapefile.
        
    date_range: str
        date_range (or subtitle)
        Expected to decribe the start and end date of the data analysed. 
        This is printed under the title. 
        format: dd Month yyyy to dd Month yyy.
        Default=""
        
    subplot_titles: list of strings
        subplot_titles for labeling each subplot title
        if None, then subtitles are ["GWL1.5", "GWL2.0", "GWL3.0"]
        otherwise specify a list of three strings.
        
    projection:
        Specify projection of the maps. The default suits Australia.
        Formally "crs".
        If None, defaults to
        ccrs.LambertConformal(central_latitude=-24.75,
                              central_longitude=134.0,
                              cutoff=30,
                              standard_parallels=(-10, -40),
        unless select_area is not None, then defaults to
        ccrs.PlateCarree()
                            
    area_linewidth: float
        linewidth of state/region borders.
        Default =0.3
        
    xlim: tuple of floats
        longitude limits
        Default = (113, 154)
        
    ylim: tuple of floats
        latitude limits
        Default = (-43.5, -9.5)
        
    cmap: matplotlib colormap
        color map for gridded and/or station data
        See cmap_dict for suggested colormaps.
        Default cmap set to cm.Greens.
        Please choose appropriate colormaps for your data.

    cmap_bad: color
        define the color to set for "bad" or missing values
        default "lightgrey"
        
    cbar_extend: one of {'neither', 'both', 'min', 'max'}.
        eg "both" changes the ends of the colorbar to arrows to indicate that
        values are possible outside the scale show.
        If contour or contourf is True, then cbar_extend will be overridden to "none".
        Default is "both"
        
    ticks: list or arraylike
        Define the ticks on the colorbar. Define any number of intervals. 
        This will make the color for each interval one discrete color, 
        instead of a smooth color gradient.
        If None, linear ticks will be auto-generated to fit the provided data.

    tick_labels: list
        Labels for categorical data. 
        If tick_labels is used, then pcolormesh is used to plot data 
        and does not allow contour or contourf to be used.
        Tick labels will correspond to the ticks.
        
    cbar_label: string
        defines the title for the color bar. 
        This should indicate the variable name and the units eg 
        "daily rainfall [mm]",
        "annual rainfall [mm]", 
        "monthly rainfall anomaly [mm]",
        "tas [\N{DEGREE SIGN}C]".
        Default is ""
        
    baseline: string
        the baseline period for anomalies, eg "1961 - 1990".
        
    dataset_name: string
        describes the source of the data eg "AGCD v2" or "BARPA-R ACCESS-CM2"
        
    issued_date: string
        The date of issue. If None is supplied, then today's date is printed.
        To supress, set to ""
        
    label_states: bool
        if set to True and Australian states are the shapefile selected,
        then each state is labelled with its three-letter abbreviation. 
        Default False.
        
    contourf: bool
        if True then the gridded data is visualised as smoothed filled contours. 
        Default is False.
        Use with caution when plotting data with negative and positive values;
        Check output for NaNs and misaligned values.  

    contour: bool
        if True then the gridded data is visualised as smoothed unfilled grey contours.
        Default is True.
        Using both contourf and contour results in smooth filled contours
        with grey outlines between the color levels.

    select_area: list
        A list of areas (eg states) that are in the geopandas.GeoDataFrame.
        Inspect the regions gdf for area names. eg ["Victoria", "New South Wales"]

    land_shadow: bool
        Used when select_area is not None. 
        This option controls whether to show Australian land area that is outside 
        the select area in grey for visual context.
        Default False.

    watermark: str
        text over the plot for images not in their final form. 
        If the plot is in final form, set to None. 
        Suggestions include "PRELIMINARY DATA", "DRAFT ONLY", 
        "SAMPLE ONLY (NOT A FORECAST)", "EXPERIMENTAL IMAGE ONLY"
        default "EXPERIMENTAL\nIMAGE ONLY"

    watermark_color: default "r"
        for the watermark, this changes the colour of the text.
        The default is red. Only change color if red is not visible. 

    show_logo: bool, default False
        show ACS logo in top left corner.

    infile: str
        Not yet tested. 
        The idea is to read in 2D netCDF data and use this as the mappable data.

    outfile: str
        The location to save the figure. 
        If None, then figure is saved here f"figures/{title.replace(' ', '-')}.png"

    savefig: bool
        default is True
        If set to False, then fig is not saved.

    orientation: {"horizontal", "vertical", "square"}
        whether the four plots are orientatied in a vertical stack,
        horizontally, or in a 2-by-2 grid (left to right).
        Default "horizontal"
 
    tick_rotation: int [-360,360]
        Angle to rotate colorbar tick labels.
        Default is None. Tick labels will be horizontal if colorbar is vertical,
        or vertical if colorbar is horizontal.
        Easiest to read if tick_rotation = 0
        
    vcentre: float, eg 0
        default is None.
        Align centre of colormap to this value. 
        Intended for using a divergent colormap with uneven number of ticks 
        around the centre, eg for future temperature anomalies with a larger
        positive range compared to the negative range.

    Returns
    -------
    A four panel plot (one baseline and three anomalies) saved as a png 
    in a "figures" file in your working directory.
    This function returns fig and ax.
    """

    
    if orientation=="horizontal":
        cax_bounds = [1.05,0,0.1,1]
        if tick_rotation is None:
            tick_rotation = 0
        nrows = 1
        ncols = 4
        cbar_location = "right"
        plots_rect = (0.01,0.05,0.98,0.85) #left bottom width height
        # text annotation xy locations
        text_xy = {"title": (0.5, 0.9),
                   "date_range": (0.5, 0.87),
                   "watermark": (0.45, 0.41),}
        subtitle_xy = None
        if figsize is None:
            figsize=(10, 3)
        
    elif orientation=="vertical":
        cax_bounds = [-0.4,-0.3,1.6,0.1]
        if tick_rotation is None:
            tick_rotation = -90
        nrows = 4
        ncols = 1
        cbar_location = "bottom"
        plots_rect = (0,0.05,1,0.85) #left bottom width height
        # text annotation xy locations
        text_xy = {"title": (0.5, 0.94),
               "date_range": (0.5, 0.93),
               "watermark": (0.45, 0.41),}
        subtitle_xy = (-0.55, 0.6)
        if figsize is None:
            figsize=(3, 8)
        
    elif orientation=="square":
        cax_bounds = [1.05,0,0.1,1]
        if tick_rotation is None:
            tick_rotation = 0
        nrows = 2
        ncols = 2
        cbar_location = "right"
        plots_rect = (0.01,0.05,0.98,0.85) #left bottom width height
        # text annotation xy locations
        text_xy = {"title": (0.5, 0.93),
                   "date_range": (0.5, 0.92),
                   "watermark": (0.45, 0.41),}
        subtitle_xy = None
        if figsize is None:
            figsize=(6,4.5)
    else:
        print('orientation must be one of ["horizontal", "vertical", "square"]')

    if gwl12_tick_rotation is None:
        gwl12_tick_rotation = tick_rotation
    
    if regions is None:
        try:
            regions = regions_dict[name]
        except:
            print(f"Could not read regions_dict[{name}]")

    regions = regions.to_crs(crs = "GDA2020")

    # Set default projection for Australia maps and selection maps
    if projection is None:
        if select_area is None:
            # Default for Australian map
            projection = ccrs.LambertConformal(
                central_latitude=-24.75,
                central_longitude=134.0,
                cutoff=30,
                standard_parallels=(-10, -40),
            )
        else:
            projection = ccrs.PlateCarree()

    fig, axs = plt.subplots(nrows=nrows,
                            ncols=ncols,  
                            sharey=True,
                            sharex=True,
                            figsize=figsize, 
                            layout="constrained",
                            subplot_kw={'projection': projection, "frame_on":False},)

    gwl12_cmap.set_bad(cmap_bad)
    cmap.set_bad(cmap_bad)

    station_dfs = [station_df_gwl12, station_df_gwl15, station_df_gwl20, station_df_gwl30]
    if any(df is not None for df in station_dfs) and markersize is None:
        markersize=(100 - 80*len(station_dfs[0])/5000)*(figsize[0]*figsize[1])/48/4

    if subplot_titles is None:
        subplot_titles = [f"GWL{[1.2, 1.5, 2.0, 3.0][i]}" for i in range(4)]

    # -------- plot baseline plot and its colorbar ---------------------
    ax, norm, cont, middle_ticks = plot_data(regions=regions,
                                             data=ds_gwl12, 
                                             station_df = station_df_gwl12,
                                             markersize=markersize,
                                             xlim=xlim,
                                             ylim=ylim,
                                             cmap=gwl12_cmap,
                                             cbar_extend=gwl12_cbar_extend,
                                             ticks=gwl12_ticks,
                                             tick_labels=gwl12_tick_labels,
                                             contourf=contourf,
                                             contour=contour,
                                             ax=axs.flatten()[0],
                                             subtitle=subplot_titles[0],
                                             subtitle_xy=subtitle_xy,
                                             facecolor=facecolor,
                                             mask_not_australia = mask_not_australia,
                                             mask_australia=mask_australia,
                                             agcd_mask=agcd_mask,
                                             area_linewidth=area_linewidth,
                                             stippling=stippling_gwl12,
                                             vcentre=gwl12_vcentre,
                                            )
    cbar = plot_cbar(cont=cont,
                  norm=norm,
                  ax=axs.flatten()[0],
                  cbar_extend=gwl12_cbar_extend, 
                  cbar_label=gwl12_cbar_label,
                     location=cbar_location,
                  ticks=gwl12_ticks, 
                  tick_labels=gwl12_tick_labels,
                  middle_ticks=middle_ticks,
                  cax_bounds=cax_bounds,
                  rotation=gwl12_tick_rotation,
                  )
    # ------- end plot baseline plot and its colorbar ---------------------

    # ------- plot three scenarios as anomalies from baseline--------------
    for i, ds in enumerate([ds_gwl15, ds_gwl20, ds_gwl30]):
        station_df = station_dfs[i+1]
        stippling = [stippling_gwl15, stippling_gwl20, stippling_gwl30][i]
        subtitle = subplot_titles[i+1]
        ax, norm, cont, middle_ticks = plot_data(regions=regions,
                                                 data=ds, 
                                                 station_df = station_df,
                                                 markersize=markersize,
                                                 xlim=xlim,
                                                 ylim=ylim,
                                                 cmap=cmap,
                                                 cbar_extend=cbar_extend,
                                                 ticks=ticks,
                                                 tick_labels=tick_labels,
                                                 contourf=contourf,
                                                 contour=contour,
                                                 ax=axs.flatten()[i+1],
                                                 subtitle=subtitle,
                                                 subtitle_xy=subtitle_xy,
                                                 facecolor=facecolor,
                                                 mask_not_australia = mask_not_australia,
                                                 mask_australia=mask_australia,
                                                 agcd_mask=agcd_mask,
                                                 area_linewidth=area_linewidth,
                                                 stippling=stippling,
                                                vcentre=vcentre,)
        
        # if select a specific area -----------
        ax = plot_select_area(select_area=select_area, 
                              ax=ax,
                              xlim=xlim,
                              ylim=ylim,
                              regions=regions,
                              land_shadow=land_shadow,
                              area_linewidth=area_linewidth,
                              )
        # ---------------------------------------------                    
        ax.axis('off')
    
    # colorbar -----------------------------------------------------------
    cbar = plot_cbar(cont=cont,
                  norm=norm,
                  ax=axs.flatten()[-1],
                  cbar_extend=cbar_extend, 
                  cbar_label=cbar_label,
                  ticks=ticks, 
                  tick_labels=tick_labels,
                  middle_ticks=middle_ticks,
                  cax_bounds =cax_bounds,
                  location=cbar_location,
                    rotation=tick_rotation,
                    )
    
    #------------------------------------------
    
    
    # plot border and annotations -----------------
    fig.get_layout_engine().set(rect=plots_rect)
    
    ax111 = fig.add_axes([0.01,0.01,0.98,0.98], xticks=[], yticks=[]) #(left, bottom, width, height)
        
    ax111 = plot_titles(title=title,
                        date_range = date_range, 
                        baseline = baseline, 
                        dataset_name= dataset_name,
                        issued_date=issued_date,
                        watermark=watermark, 
                        watermark_color=watermark_color,
                        ax=ax111,
                        text_xy = text_xy,
                        title_ha = "center",
                        orientation = orientation,
                   )
    # draw border
    # ax111.axis(True)
    ax111.axis(False)
    # --------------------------------------------
    
    if outfile is None:
        PATH = os.path.abspath(os.getcwd())
        outfile = f"{PATH}/figures/{title.replace(' ', '-')}.png"
        os.makedirs(os.path.dirname(outfile), exist_ok=True)
    
    if savefig:
        plt.savefig(outfile, dpi=300,)
    return fig, ax


