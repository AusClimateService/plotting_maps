# # Standardising Australia Hazard Maps
# The goal of this script is to map climate hazards for Australia so that data can be consistently presented.
# this code is designed to work with hh5 analysis3-24.04 venv

import datetime

# import packages used in this workflow
import geopandas as gpd
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as image
import xarray as xr
import cartopy.crs as ccrs

# import colormap packages
import matplotlib.cm as cm
import cmaps
from matplotlib.colors import ListedColormap, BoundaryNorm

from shapely.geometry import box

# define some standard imput for the maps
crs = ccrs.LambertConformal(central_latitude = -24.75, central_longitude = 134.0, cutoff=30, standard_parallels=(-10,-40))
xlim = (110,170)
ylim = (-45, -5)
logo = image.imread("ACS_Logo_Blue_on_white_Stacked.png")


# # Suggested colormaps and scales
# Using suggested colormaps and scales will improve the consistency across teams producing similar variables. This will support comparison across different plots.
# - see many colormaps here: https://www.ncl.ucar.edu/Document/Graphics/color_table_gallery.shtml
# This are a dictionary of suggested colormaps matched with possible variables to plot.  This includes color maps for total amount and for anomalies

cmap_dict = {"sst":cmaps.cmocean_tempo,
             "sst_anom": cmaps.cmocean_balance_r,
             "mhw_days": cm.YlOrRd,
             "mhw_intensity": cm.hot_r,
             "surface_pH":cm.YlGnBu,
             "surface_pH_1":cmaps.cmocean_phase,
             "surface_pH_2":cm.cool,
             "surface_aragonite_sat": cmaps.cmocean_delta,
             
             "tas":cm.Spectral_r,
             "tas_anom":cm.RdBu_r,
             "tas_anom_1":cm.seismic,
             "tas_deciles_bwr":cm.bwr,
             "EHF_days": cm.YlOrRd,
             "EHF_days_1": cm.YlOrBr,
             "EHF_duration":cm.hot_r,
             
             "AFDRS_category":ListedColormap(["white","green","orange", "red", "darkred"]),
             "ffdi_category":ListedColormap(["green", "blue", "yellow","orange", "red", "darkred"]),

             "pr":cm.YlGnBu,
             "pr_1":cmaps.cmocean_deep,
             "pr_days":cm.Blues,
             "pr_GMT_drywet": cmaps.GMT_drywet,
             "pr_rainbow_1":cmaps.prcp_1,
             "pr_rainbow_2":cmaps.prcp_2,
             "pr_rainbow_3":cmaps.prcp_3,         
             "pr_anom":cm.BrBG,
             "pr_anom_1":cmaps.cmocean_curl,
             "pr_anom_12lev":cmaps.precip_diff_12lev,
             "pr_chance_extremes": cmaps.cmocean_solar_r,
             
             "tc_days": cm.RdPu,
             "tc_intensity":cm.PuRd,
             "tc_days_anom":cm.PuOr,
             "tc_intensity_anom":cm.PiYG_r,

             "xts_freq":cmaps.cmocean_dense,
             "xts_intensity":cmaps.cmocean_matter,
             "xts_freq_anom":cmaps.cmocean_balance_r,
             "xts_intensity_anom": cmaps.cmocean_curl_r,
             
             "drought_severity": cm.RdYlGn_r,   
             "drought_duration":cmaps.hotres,
             
             "anom_BlueYellowRed":cmaps.BlueYellowRed,
             "anom_BlueYellowRed_r":cmaps.BlueYellowRed_r,
             "anom":cmaps.BlueWhiteOrangeRed,
             "anom_r":cmaps.BlueWhiteOrangeRed_r,
             "anom_b2r":cmaps.cmp_b2r,
             "anom_b2r_r":cmaps.cmp_b2r_r,
             "anom_coolwarm":cmaps.MPL_coolwarm,
             "anom_coolwarm_r":cmaps.MPL_coolwarm_r,
             "anom_deciles":cm.bwr,
             "anom_deciles_r":cm.bwr_r,
             "anom_veg_1":cmaps.NEO_div_vegetation_a, 
             "anom_veg_1_r":cmaps.NEO_div_vegetation_a_r,          

             "BuGnRd":cmaps.temp1,

             "rh_19lev": cmaps.rh_19lev,
             "sunshine_9lev":cmaps.sunshine_9lev,
             "sunshine_diff_12lev":cmaps.sunshine_diff_12lev,
             "inferno":cm.inferno,
             "Oranges":cm.Oranges,
             "OrRd":cm.OrRd,
             "Greens": cm.Greens,

             "topo":cmaps.OceanLakeLandSnow,
             "gmt_relief":cmaps.GMT_relief,
            }


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
     
# # Load the State and Region shape files
# write a dictionary of the shapefile geopandas dataframes. These will be used for state boundaries, LGAs, NRM, etc
regions_dict = {}

shape_files = ["aus_local_gov", "aus_states_territories", "australia", "nrm_regions", "river_regions", "ncra_regions"]
for name in shape_files:
    regions_dict.update({name: gpd.read_file(f'/g/data/ia39/aus-ref-clim-data-nci/shapefiles/data/{name}/{name}.shp')})
regions_dict.update({"broadacre_regions" : gpd.read_file('/g/data/ia39/aus-ref-clim-data-nci/shapefiles/data/broadacre_regions/aagis_asgs16v1_g5a.shp')})

# define a white mask for the area outside of Australian land
# We will use this to hide data outside of the Australian land borders. 
# note that this is not a data mask, the data under the masked area is still loaded and computed, but not visualised 

australia = regions_dict["australia"].copy()

# Define the CRS of the shapefile manually
australia.crs = crs

# This mask is a rectangular box around the maximum land extent of Australia with a buffer of 10 degrees on every side, with the Australian land area cut out so only the ocean is hidden.
not_australia = gpd.GeoSeries(data = [box(*box(*australia.total_bounds).buffer(20).bounds).difference(australia["geometry"].values[0])],
                              crs=ccrs.PlateCarree())

# # Define a function for plotting maps
# This is the function you call to plot all the graphs
def plot_acs_hazard(name = "aus_states_territories",
                    regions = None,
                    data = None,
                    mask_not_australia = True,
                    facecolor = None, 
                    edgecolor = "black",
                    figsize = (8,6),
                    title = None,
                    date_range = "",
                    crs = None,
                    area_linewidth = 0.3,
                    xlim = (110,170),
                    ylim = (-45, -5),
                    cmap = cm.Greens,
                    cbar_extend = "neither",
                    ticks = None, 
                    tick_labels = None,
                    cbar_label = "",
                    baseline = None,
                    dataset_name = None,
                    issued_date = None,                      
                    label_states = False,
                    contourf = True,
                    contour = True,
                    select_area = None,
                    land_shadow = False,
                    watermark = "EXPERIMENTAL\nIMAGE ONLY",
                    infile = None,
                    outfile = None,
                    ):
    """This function takes a name of an Australian shapefile collection for data in /g/data/ia39/aus-ref-clim-data-nci/shapefiles/data/ and hazard data from a 2D Xarray data array
    and plots the data on a map of Australia with the shape outlines. 

    Parameters
    ----------
    name: str
        name of a shapefile collection in /g/data/ia39/aus-ref-clim-data-nci/shapefiles/data/ to get regions from.
        
    regions: geopandas.GeoDataFrame
        if None, then will try to read from regions_dict[{name}].
        
    data: xr.DataArray
        a 2D xarray DataArray which has already computed the average, sum, anomaly, metric or index you wish to visualise. This function is resolution agnostic. 
        
    mask_not_australia: boolean
        decides whether or not the area outside of Australian land is hidden. Default is True.
        
    facecolor: color
        color of land when you plot the regions without climate data. facecolorrecommendationss include "white", "lightgrey", "none".
        
    edgecolor: color
        defines the color of the state/region borders. edgecolor recommendations include "black" and "white".
        
    figsize: tuple
        defines the width and height of the figure in inches.I  Reccommend (8,6) for Australia-wide plots and (6,6) for individual states etc..
        
    title: str
        defines the text inside the plot. If none is given, then will print the name of the shape file.
        
    date_range: str
        decribes the start and end date of the data analysed. This is printed under the title. format: dd Month yyyy to dd Month yyy.
        
    crs: optional
        defines the coordinate reference system. Similar to transform or projection.
        
    area_linewidth: float, optional
        the width of state/region borders only. All other linewidths are hardcoded.
        
    xlim: tuple, optional
        longitude min and max of the plot area.
        
    ylim: tuple, optional
        latitude min and max of the plot area.
        
    cmap:
        defines the colormap used for the data. See cmap_dict for suggested colormaps. If none, cmap set to cm.Greens. Please choose appropriate colormap for your data.
        
    cbar_extend: one of {'neither', 'both', 'min', 'max'}.
        eg "both" changes the ends of the colorbar to arrows to indicate that values are possible outside the scale show. If contour or contourfise True, then cbar_extend will be overridden to "none".
        
    ticks: list or arraylike
        Define the ticks on the colorbar. Define any number of intervals. This will make the color for each interval one discrete coloh, instead of a smooth color gradient. If None, linear ticks will be auto-generated to fit the provided data. 

    tick_labels: list
        Labels for categorical data. If tick_labels is used, then pcolormesh is used to plot data and does not allow contour or contourf to be used. Tick labels will correspond to the ticks. Ticks are assumed to be linear eg [101, 102, 103, 104] for four categories.
        
    cbar_label: string
        defines the title for the color bar. This should indicate the variable name and the units eg "daily rainfall [mm]", "annual rainfall [mm]", "monthly rainfall anomaly [mm]", "tas [\N{DEGREE SIGN}C]".
        
    baseline: string
        the baseline period for anomalies, eg "1961 - 1990".
        
    dataset_name: string
        describes the source of the data eg "AGCD v2" or "BARPA-R ACCESS-CM2"
        
    issued_date: string
        The date of issue. If None is supplied, then today's date is printed.
        
    label_states: bool
        if set to True and Australian states are the shapefile selected, then each state is labelled with its three-letter abbreviation. Default False.
        
    contourf: bool
        if True then the gridded data is visualised as smoothed filled contours. Default is True.
        
    contour: bool
        if True then the gridded data is visualised as smoothed unfilled grey contours. Default is True. Using both contourf and contour results in smooth filled contours with grey outlines between the color levels.
        
    select_area: list
        A list of areas (eg states) that are in the geopandas.GeoDataFrame. Inspect the regions gdf for area names. eg ["Victoria", "New South Wales"]
        
    land_shadow: bool
        Used when select_area is not None. This option controls whether to show Australian land area that is outside the select area in grey for visual context. Default False.
        
    watermark: str
        red text over the plot for images not in their final form. If the plot is in final form, set to None. Suggestions include "PRELIMINARY DATA", "DRAFT ONLY", "SAMPLE ONLY (NOT A FORECAST)", "EXPERIMENTAL IMAGE ONLY"
        
    infile: str
        Not yet tested. The idea is to read in 2D netCDF data and use this as the mappable data.
        
    outfile: str
        The location to save the figure. If None, then figure is saved here f"figures/{title.replace(' ', '_')}.png"

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

    # Set default crs for Australia maps and selction maps
    if crs is None:
        if select_area is None:
            # Default for Australian map
            crs = ccrs.LambertConformal(central_latitude = -24.75, central_longitude = 134.0, cutoff=30, standard_parallels=(-10,-40))
        else:
            crs = ccrs.PlateCarree()
        
    # Define the CRS of the shapefile manually
    regions = regions.to_crs(crs.proj4_init)

    # Set up the plot
    fig = plt.figure(figsize=figsize, zorder=1, layout="constrained", )
    ax = plt.axes(projection = crs, )
    ax.set_global()

    if infile is not None:
        data = xr.open_dataset(infile)
    
    if data is not None:
        data = data.squeeze()
        
        facecolor = "none"
        
        if ticks is None:
            norm=None
        else:   
            if tick_labels is None:
                norm = BoundaryNorm(ticks, cmap.N)
            else:
                step = ticks[1]-ticks[0] 
                bounds = np.arange(ticks[0]-0.5*step, ticks[-1]+0.5*step+1,step)
                norm = BoundaryNorm(bounds, cmap.N)
        
        # plot the hazard data
        # rename coords if necessary:
        if ("latitude" and "longitude") in list(data.coords):
            data = data.rename({"latitude":"lat", "longitude": "lon"})
        
        if contourf and tick_labels is None:
            cont = ax.contourf(data.lon, data.lat, data,
                           cmap = cmap,
                           norm=norm,
                           levels=ticks, 
                           extend = cbar_extend,
                           zorder=2, 
                           transform= ccrs.PlateCarree())
        else:
            cont = ax.pcolormesh(data.lon, data.lat, data, 
                                 cmap = cmap,
                                 norm=norm,
                                 zorder=2,
                                 transform= ccrs.PlateCarree())

        if tick_labels is None:
            cbar = plt.colorbar(cont,
                                ax=ax, 
                                extend = cbar_extend,
                                cax=ax.inset_axes([0.8,0.2, 0.03, 0.6]), 
                                ticks=ticks,
                                norm=norm,)
        else:
            # for categorical data
            cbar = plt.colorbar(cont,
                                ax=ax, 
                                extend = None,
                                cax=ax.inset_axes([0.8,0.2, 0.03, 0.6]), 
                                norm=norm,
                                drawedges=True,
                                ticks=ticks,
                               )
            cbar.ax.set_yticklabels(tick_labels)
            
        cbar.ax.set_title(cbar_label, zorder=8, y=1.1, loc="center")
        if contour and tick_labels is None:
            cont = plt.contour(data.lon, data.lat, data,
                               colors="grey", norm=norm, levels=ticks,
                               extend = cbar_extend, linewidths=0.2, zorder=3,
                               transform= ccrs.PlateCarree())
            cbar.add_lines(cont)

        if mask_not_australia:
            # outside the shape, fill white
            ax.add_geometries(not_australia, crs= not_australia.crs, facecolor="white", linewidth=0, zorder=5,)

    if label_states and name=="aus_states_territories":
        # label the states with their name in the centre of the state
        states = ["WA", "NT", "SA", "QLD", "NSW", "VIC", "TAS"]
        lats = [-26.06, -19.7, -29.0, -24.2, -32.3, -37.1, -42.4]
        lons = [119.9, 133,  133, 143.6, 143.7, 142.0, 145.3]
        for lon, lat, state in zip(lons, lats, states):
            ax.text(lon, lat, state, size=12, zorder=10, transform = ccrs.PlateCarree())

    ax.add_geometries(regions["geometry"],
                  crs=crs,
                  facecolor=facecolor,
                  edgecolor=edgecolor, 
                  linewidth=area_linewidth,
                  zorder=4,)

    # if select a specific area
    if select_area is None:
        ax.set_extent([xlim[0], xlim[1], ylim[0], ylim[1]])
    else:
        assert isinstance(select_area, list), "select_area must be a list"
        # select state
        name_column = [name for name in regions.columns if 'NAME' in name.upper()][0]
        area =regions.loc[regions[name_column].isin(select_area)] 
        map_total_bounds = area.total_bounds
        minx, miny, maxx, maxy = map_total_bounds
        mid_x = (minx+maxx)/2
        mid_y = (miny+maxy)/2
        max_range = np.max([(maxy-miny), (maxx-minx)])
        buffer=0.1*max_range
        
        not_area =  gpd.GeoSeries(data = [box(*box(*map_total_bounds).buffer(10*buffer).bounds).difference(area.dissolve()["geometry"].values[0])],
                                      crs=crs, )

        # mask outside selected area
        if land_shadow:
            # show land as light grey
            not_area.plot(ax=ax, facecolor="lightgrey", linewidth=0, zorder = 4)
        else:
            # mask white
            not_area.plot(ax=ax, facecolor="white", linewidth=0, zorder = 4)
               
        plt.xlim(mid_x - 0.6*max_range, mid_x + 0.8*max_range)
        plt.ylim(mid_y - 0.7*max_range, mid_y + 0.7*max_range)      
        
    
    # Set the plot title and axis labels
    if title is None:
        if data is None:
            title=name
        else:
            title = f"{data.name} {name}"

    ax.text(0.1,0.06,f"{title}\n{date_range}",
            fontsize = "large",
            horizontalalignment='left',
            transform=ax.transAxes,
            zorder=10)

    if baseline is not None:
        # print base period inside bottom left corner
        ax.text(x=0.01, y=0.01, s= f"Base period: {baseline}", fontsize=6, transform=ax.transAxes, zorder = 10)
    # print copyright outside bottom left corner
    ax.text(x=0.01, y=-0.03, s= f"\u00A9 Commonwealth of Australia {datetime.datetime.now().year}, Australian Climate Service",
            fontsize=6, transform=ax.transAxes, zorder = 10)
    # print data source inside bottom right
    ax.text(x=0.99, y=0.01, s= f"Dataset: {dataset_name}", fontsize=6, transform=ax.transAxes, horizontalalignment="right", zorder = 10)
    # print issued date on bottom right under the border. Set to today's date if None supplied
    if issued_date is None:
        issued_date = datetime.datetime.today().date().strftime('%d %B %Y')
    ax.text(x=0.99, y=-0.03, s= f"Issued: {issued_date}", fontsize=6, transform=ax.transAxes, horizontalalignment="right", zorder = 10)

    if watermark is not None:
         ax.text(x=0.45, y=0.5, s=watermark.upper(), 
                 fontsize=36,
                 transform=ax.transAxes,
                 horizontalalignment="center",
                 verticalalignment = "center",
                 zorder = 10, 
                 wrap=True,
                 alpha = 0.5,
                 color='r')

    # remove any lat lon labels and ticks
    plt.xlabel('')
    plt.ylabel('')
    plt.tick_params(axis="both", which = "both", bottom =False, left=False, labelbottom = False, labelleft=False)

    ax.patch.set_edgecolor('black')  
    ax.patch.set_linewidth(2)
    
    fig.set_figheight(figsize[1])
    fig.set_figwidth(figsize[0])

    # Place logo in top left
    ins = ax.inset_axes([0.,0.78,0.3,0.3],)
    ins.set_xticks([])
    ins.set_yticks([])
    ins.imshow(logo,  ) 

    if outfile is None:
        outfile = f"figures/{title.replace(' ', '_')}.png"

    plt.savefig(outfile, dpi = 300)
    return fig, ax

