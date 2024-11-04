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


ncfile = "../data/gridfile_NXP036_lbx.nc4"
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

fig = plt.figure(figsize=(24, 16))

ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=0.0, globe=None))
plt.axis([-180, 180, -80, 80])
plt.rcParams['font.family'] = 'MicroSoft YaHei'

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
                    fc='none', linewidth=0.3, alpha=0.8))
            elif nlon == 6:
                ax.add_patch(mpatches.Polygon(
                    np.array([[x[0], y[0]], [x[1], y[1]], [x[2], y[2]], [x[3], y[3]], [x[4], y[4]], [x[5], y[5]]]),
                    ec='black',
                    fc='none', linewidth=0.3, alpha=0.8))
            elif nlon == 5:
                ax.add_patch(mpatches.Polygon(
                    np.array([[x[0], y[0]], [x[1], y[1]], [x[2], y[2]], [x[3], y[3]], [x[4], y[4]]]),
                    ec='black',
                    fc='none', linewidth=0.3, alpha=0.8))
        else:
            if nlon == 7:
                ax.add_patch(mpatches.Polygon(
                    np.array(
                        [[x[0], y[0]], [x[1], y[1]], [x[2], y[2]], [x[3], y[3]], [x[4], y[4]], [x[5], y[5]],
                         [x[6], y[6]]]),
                    ec='black',
                    fc=k[i], linewidth=0.3, alpha=0.8))
            elif nlon == 6:
                ax.add_patch(mpatches.Polygon(
                    np.array([[x[0], y[0]], [x[1], y[1]], [x[2], y[2]], [x[3], y[3]], [x[4], y[4]], [x[5], y[5]]]),
                    ec='black',
                    fc=k[i], linewidth=0.3, alpha=0.8))
            elif nlon == 5:
                ax.add_patch(mpatches.Polygon(
                    np.array([[x[0], y[0]], [x[1], y[1]], [x[2], y[2]], [x[3], y[3]], [x[4], y[4]]]),
                    ec='black',
                    fc=k[i], linewidth=0.3, alpha=1))

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

plt.savefig("threshold.png", dpi=300)
plt.show()
