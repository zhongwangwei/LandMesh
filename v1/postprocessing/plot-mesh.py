from netCDF4 import Dataset
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import matplotlib.patches as mpatches
from matplotlib import cm, colors
import cartopy as cart
import cartopy.feature as cfeature




def isout(xs, nlo):
    for t in range(0, nlo - 1):
        for s in range(t + 1, nlo):
            if abs(xs[t] - xs[s]) > 100:
                return 1
    return 0


ncfile = "gridfile_NXP036_lbx.nc4"
nc_obj = Dataset(ncfile)
lons = nc_obj['GLONM'][:]
lats = nc_obj['GLATM'][:]
wlons = nc_obj['GLONW'][:]
ngrwm = nc_obj['itab_w%im'][:, :]
ref_pl = nc_obj['ref_pl'][:, :]

ref = np.zeros(wlons.shape[0])
ref_s = np.zeros(16)

for i in range(0, wlons.shape[0]):
    if ref_pl[0, i] == 1:
        ref[i] += 1
        ref_s[0] += 1
    if ref_pl[2, i] == 1:
        ref[i] += 2
        ref_s[3] += 1
    if ref_pl[4, i] == 1:
        ref[i] += 4
        ref_s[5] += 1
    if ref_pl[7, i] == 1:
        ref[i] += 8
        ref_s[7] += 1


norm = colors.Normalize(vmin=0, vmax=ref.max())
k = cm.Set2(norm(ref))

# Modify the figure creation and aspect ratio
lon_range = 107 - 66  # = 41
lat_range = 43 - 22   # = 21
aspect_ratio = lon_range / lat_range

fig = plt.figure(figsize=(24, 24/aspect_ratio))  # This will maintain the geographical ratio
ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=0.0, globe=None))
plt.axis([66, 107, 22, 43])
plt.rcParams['font.family'] = 'MicroSoft YaHei'

# Add DEM data reading
dem_file = "glofas_dem.nc"
dem_data = Dataset(dem_file)
dem_lons = dem_data['lon'][:]  # Add longitude coordinates
dem_lats = dem_data['lat'][:]   # Add latitude coordinates
dem_elv = dem_data['elv'][:]         # Read elevation data

# Find indices for our region of interest
lon_mask = (dem_lons >= 66) & (dem_lons <= 107)
lat_mask = (dem_lats >= 22) & (dem_lats <= 43)

# Slice the DEM data to our region
dem_elv_subset = dem_elv[lat_mask][:, lon_mask]

# Plot the subset of DEM data
plt.imshow(dem_elv_subset, 
           extent=[66, 107, 22, 43],
           cmap='terrain',
           transform=ccrs.PlateCarree(),
           alpha=0.7)

for i in range(0, wlons.shape[0]):
    x = np.zeros(7)
    y = np.zeros(7)

    nlon = 7

    for j in range(0, 7):
        x[j] = lons[int(ngrwm[i, j]) - 1]
        y[j] = lats[int(ngrwm[i, j]) - 1]
        if ngrwm[i, j] == 1 or ngrwm[i, j] == 0:
            nlon -= 1

    if isout(x, nlon) == 0:
        if ref[i] == 0:
            if nlon == 7:
                ax.add_patch(mpatches.Polygon(
                    np.array(
                        [[x[0], y[0]], [x[1], y[1]], [x[2], y[2]], [x[3], y[3]], [x[4], y[4]], [x[5], y[5]],
                         [x[6], y[6]]]),
                    ec='black',
                    fc='none', linewidth=1.0, alpha=0.3))
            elif nlon == 6:
                ax.add_patch(mpatches.Polygon(
                    np.array([[x[0], y[0]], [x[1], y[1]], [x[2], y[2]], [x[3], y[3]], [x[4], y[4]], [x[5], y[5]]]),
                    ec='black',
                    fc='none', linewidth=1.0, alpha=0.3))
            elif nlon == 5:
                ax.add_patch(mpatches.Polygon(
                    np.array([[x[0], y[0]], [x[1], y[1]], [x[2], y[2]], [x[3], y[3]], [x[4], y[4]]]),
                    ec='black',
                    fc='none', linewidth=1.0, alpha=0.3))
        else:
            if nlon == 7:
                ax.add_patch(mpatches.Polygon(
                    np.array(
                        [[x[0], y[0]], [x[1], y[1]], [x[2], y[2]], [x[3], y[3]], [x[4], y[4]], [x[5], y[5]],
                         [x[6], y[6]]]),
                    ec='black',
                    fc=k[i], linewidth=1.0, alpha=0.3))
            elif nlon == 6:
                ax.add_patch(mpatches.Polygon(
                    np.array([[x[0], y[0]], [x[1], y[1]], [x[2], y[2]], [x[3], y[3]], [x[4], y[4]], [x[5], y[5]]]),
                    ec='black',
                    fc=k[i], linewidth=1.0, alpha=0.3))
            elif nlon == 5:
                ax.add_patch(mpatches.Polygon(
                    np.array([[x[0], y[0]], [x[1], y[1]], [x[2], y[2]], [x[3], y[3]], [x[4], y[4]]]),
                    ec='black',
                    fc=k[i], linewidth=1.0, alpha=0.3))

    elif isout(x, nlon) == 1:
        for j in range(0, 7):
            if x[j] < 0:
                x[j] += 360
        if nlon == 7:
            ax.add_patch(mpatches.Polygon(
                np.array(
                    [[x[0], y[0]], [x[1], y[1]], [x[2], y[2]], [x[3], y[3]], [x[4], y[4]], [x[5], y[5]], [x[6], y[6]]]),
                ec='black',
                fc='none', linewidth=0.3, alpha=1))
        elif nlon == 6:
            ax.add_patch(mpatches.Polygon(
                np.array([[x[0], y[0]], [x[1], y[1]], [x[2], y[2]], [x[3], y[3]], [x[4], y[4]], [x[5], y[5]]]),
                ec='black',
                fc='none', linewidth=0.3, alpha=1))
        elif nlon == 5:
            ax.add_patch(mpatches.Polygon(
                np.array([[x[0], y[0]], [x[1], y[1]], [x[2], y[2]], [x[3], y[3]], [x[4], y[4]]]),
                ec='black',
                fc='none', linewidth=0.3, alpha=1))

        for j in range(0, 7):
            x[j] -= 360

        if nlon == 7:
            ax.add_patch(mpatches.Polygon(
                np.array(
                    [[x[0], y[0]], [x[1], y[1]], [x[2], y[2]], [x[3], y[3]], [x[4], y[4]], [x[5], y[5]], [x[6], y[6]]]),
                ec='black',
                fc='none', linewidth=0.3, alpha=1))
        elif nlon == 6:
            ax.add_patch(mpatches.Polygon(
                np.array([[x[0], y[0]], [x[1], y[1]], [x[2], y[2]], [x[3], y[3]], [x[4], y[4]], [x[5], y[5]]]),
                ec='black',
                fc='none', linewidth=0.3, alpha=1))
        elif nlon == 5:
            ax.add_patch(mpatches.Polygon(
                np.array([[x[0], y[0]], [x[1], y[1]], [x[2], y[2]], [x[3], y[3]], [x[4], y[4]]]),
                ec='black',
                fc='none', linewidth=0.3, alpha=1))

ax.add_feature(cfeature.COASTLINE, linewidth=1, color="gray", zorder=21)
ax.add_feature(cfeature.OCEAN, color="white", zorder=20)

ax.add_patch(
    mpatches.Polygon(np.array([[-180, 80], [-180, -80], [180, -80], [180, 80]]), ec='gray', fc='none',
                     linewidth=1,
                     alpha=0.9, zorder=220))

# Before saving/showing, modify these lines:
plt.tight_layout()
plt.margins(0, 0)
ax.set_aspect(1)  # This forces the plot to maintain the correct geographical proportions

plt.savefig("threshold.png", dpi=300, bbox_inches='tight', pad_inches=0)
plt.show()
