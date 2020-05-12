import numpy as np
from pathlib import Path
from scipy.interpolate import RectBivariateSpline
from tkinter.messagebox import showerror, showwarning

# folder name where hgt files are located:
hgtfolder = 'hgt'
hgtpath = Path.joinpath(Path.cwd(), hgtfolder)
x = np.linspace(0, 1, 3601, dtype=np.dtype('float32'))

# rewrited with showerror instead of print


class MissingHgtError(Exception):
    pass


def get_path_profile(phi_t, psi_t, phi_r, psi_r, coord_samples=500, fill_missing=False):
    [dist, atrdeg, di, delta, phi_values, psi_values] = get_path_geometry(phi_t, psi_t, phi_r, psi_r, coord_samples)
    subpaths, ufnames, flagsmissing = get_path_sections(phi_values, psi_values, fill_missing)
    hi = get_elevations(subpaths, ufnames, flagsmissing)
    if not fill_missing:
        if flagsmissing.any():
            missing = np.array2string(ufnames[np.nonzero(flagsmissing)])
            showwarning('Missing HGT files', f"Missing hgt files and filled with 0 values:\n{missing}")
    # print(f'Distance: {dist:.5f}km\nBearing angle: {atrdeg:.1f}°\nMax elevation:{np.max(hi):.0f}m')
    return dist, atrdeg, di, hi, delta


def get_elevations(subpaths, ufnames, flagsmissing):
    interpvalues = np.empty(0, dtype=int)
    for i in range(len(subpaths)):
        # open hgt filename according to subpath
        path = Path.joinpath(hgtpath, ufnames[i])
        if not flagsmissing[i]:
            with open(path, 'rb') as data:
                elevationdata = np.fromfile(data, np.dtype('>i2'), 3601 ** 2).reshape(3601, 3601)

            # interpolate elevation values with scipy.interpolate.griddata
            f1 = RectBivariateSpline(x, x, elevationdata, kx=1, ky=1)
            latcorr = int(ufnames[i][1:3]) + 1
            loncorr = int(ufnames[i][4:7])
            interpvalues = np.hstack((interpvalues, f1.ev(latcorr - subpaths[i][0], subpaths[i][1] - loncorr)))
        else:
            interpvalues = np.hstack((interpvalues, np.zeros_like(subpaths[i][0])))

    return interpvalues


def get_path_sections(phi_values, psi_values, fill_missing):
    fnames = get_hgt_names(phi_values, psi_values)

    # check if all hgt files exist in hgt directory
    _, indices = np.unique(fnames, return_index=True)
    indices = np.sort(indices)
    ufnames = fnames[indices]

    flagsmissing = np.full_like(indices, 0)
    for i in range(len(indices)):
        file = fnames[indices[i]]
        path = Path.joinpath(hgtpath, file)
        flagsmissing[i] = not path.exists()
    if not fill_missing:
        if flagsmissing.any():
            missing = np.array2string(ufnames[np.nonzero(flagsmissing)])
            showerror('Missing HGT files', f"Missing hgt files from hgt directory: {missing}")
            raise MissingHgtError(f"Missing hgt files from hgt directory:\n{missing}")

    # create sub-paths for each difference hgt file:
    subpaths = [None] * len(indices)
    for i in range(len(indices) - 1):
        startidx = indices[i]
        endidx = indices[i+1]
        subpaths[i] = (phi_values[startidx:endidx], psi_values[startidx:endidx])
        # subpaths[i] = list(zip(phi_values[startidx:endidx], psi_values[startidx:endidx]))
    subpaths[-1] = (phi_values[indices[-1]:], psi_values[indices[-1]:])
    # subpaths[-1] = list(zip(phi_values[indices[-1]:], psi_values[indices[-1]:]))
    return subpaths, ufnames, flagsmissing


def get_hgt_names(phi_values, psi_values):
    n_s = np.where(phi_values >= 0, 'N', 'S')
    e_w = np.where(psi_values >= 0, 'E', 'W')
    lat = abs(phi_values).astype(int).astype(str)
    lat = np.char.zfill(lat, 2)
    lon = abs(psi_values).astype(int).astype(str)
    lon = np.char.zfill(lon, 3)
    fnames = np.char.add(np.char.add(n_s, lat), np.char.add(e_w, lon))
    fnames = np.char.add(fnames, '.hgt')
    return fnames


def get_path_geometry(phi_t, psi_t, phi_r, psi_r, coord_samples):
    phi_t = np.deg2rad(phi_t)
    psi_t = np.deg2rad(psi_t)
    phi_r = np.deg2rad(phi_r)
    psi_r = np.deg2rad(psi_r)

    # The angle subtended by the path at the centre of the Earth, δ, from the stations’ geographic
    # coordinates using eq (65):
    delta = np.arccos(np.sin(phi_t) * np.sin(phi_r) +
                      np.cos(phi_t) * np.cos(phi_r) *
                      np.cos(psi_t - psi_r))

    # great circle distance:
    R = 6371.009
    dist = R * delta  # km

    # bearing (azimuthal direction clockwise from true North) from station t to station r
    # using eq 67:
    atr = np.arccos((np.sin(phi_r) - np.sin(phi_t) * np.cos(delta)) /
                    (np.sin(delta) * np.cos(phi_t)))

    if psi_t > psi_r:
        atr = 2 * np.pi - atr
    atrdeg = np.rad2deg(atr)

    if coord_samples == 0:
        return dist, atrdeg, None, None, None, None

    di = np.linspace(0, dist, coord_samples)

    phi_values = np.arcsin(np.sin(phi_t) * np.cos(di / R) + np.cos(phi_t) * np.sin(di / R) * np.cos(atr))
    psi_values = psi_t + np.arctan2(np.sin(atr) * np.sin(di / R) * np.cos(phi_t),
                                    np.cos(di / R) - np.sin(phi_t) * np.sin(phi_values))
    psi_values = np.remainder(psi_values + 3 * np.pi, (2 * np.pi) - np.pi)

    phi_values = np.rad2deg(phi_values)
    psi_values = np.rad2deg(psi_values)

    return [dist, atrdeg, di, delta, phi_values, psi_values]
