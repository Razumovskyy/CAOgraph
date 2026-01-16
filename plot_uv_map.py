#!/usr/bin/env python3
"""Plot UV index grid from E_<case>.UV using PyGMT."""

import argparse
from pathlib import Path

import numpy as np
import pygmt
import xarray as xr


def read_grid(grid_path: Path):
    """Read latitude/longitude vectors from Input/grid."""
    with grid_path.open() as f:
        nlat, dlat, lat_min, lat_max = map(float, f.readline().split())
        nlon, dlon, lon_min, lon_max = map(float, f.readline().split())
        nlat, nlon = int(nlat), int(nlon)
        lats = [float(f.readline().strip()) for _ in range(nlat)]
        lons = [float(f.readline().strip()) for _ in range(nlon)]
    return np.array(lats), np.array(lons), (lat_min, lat_max, dlat), (lon_min, lon_max, dlon)


def read_uv_field(data_path: Path, nlat: int, nlon: int, lats_expected):
    """Parse E_<case>.UV values into a (nlat, nlon) array."""
    values = np.zeros((nlat, nlon), dtype=float)
    with data_path.open() as f:
        lat_idx = 0
        while lat_idx < nlat:
            lat_line = f.readline()
            if not lat_line:
                raise ValueError(f"Reached EOF after {lat_idx} latitudes; expected {nlat}.")
            lat_val = float(lat_line.strip())
            if abs(lat_val - lats_expected[lat_idx]) > 1e-3:
                raise ValueError(
                    f"Latitude mismatch at row {lat_idx}: file has {lat_val}, grid expects {lats_expected[lat_idx]}"
                )

            row = []
            while len(row) < nlon:
                data_line = f.readline()
                if not data_line:
                    raise ValueError(f"Latitude {lat_val}: expected {nlon} values, got {len(row)} before EOF.")
                row.extend(float(v) for v in data_line.split())
            values[lat_idx, :] = row[:nlon]
            lat_idx += 1
    return values


def plot_uv(grid: xr.DataArray, output_path: Path):
    """Plot the UV index grid to a PNG using PyGMT."""
    grid = grid.sel(lat=slice(-67, None))
    
    # w = 9  # 7–13 обычно достаточно для таких полос
    # grid = (
    #     grid.pad(lon=(w//2, w//2), mode="wrap")
    #         .rolling(lon=w, center=True).mean()
    #         .isel(lon=slice(w//2, -w//2))
    # )
    # # опционально: легкое сглаживание по широте
    # grid = grid.rolling(lat=3, center=True).mean()

    w = 3  # half-window; итоговое окно = 2*w+1 (тут 7)
    gpad = xr.concat([grid.isel(lon=slice(-w, None)), grid, grid.isel(lon=slice(0, w))], dim="lon")
    grid = gpad.rolling(lon=2*w+1, center=True).mean().isel(lon=slice(w, w + grid.lon.size))

    
    vmax = float(np.nanmax(grid))
    cmax = max(1.0, np.ceil(vmax))

    pygmt.makecpt(cmap="turbo", series=[0, cmax, max(0.5, cmax / 20)], continuous=True)

    fig = pygmt.Figure()
    region = [
        float(grid.lon.min()), 
        float(grid.lon.max()), 
        -67.0, 
        float(grid.lat.max())
    ]
    fig.shift_origin(xshift="1c", yshift="1c")
    fig.grdimage(grid, region=region, projection="W0/16c", cmap=True, frame=["xaf", "yaf"])
    # Draw coastlines on top without opaque land/water fills so the data stays visible.
    fig.coast(shorelines="1/0.25p,black", borders="1/0.2p,gray", frame="WSen")
    fig.colorbar(frame="af")
    fig.savefig(str(output_path), crop=True, resize="+m1c")
    return output_path


def main():
    parser = argparse.ArgumentParser(description="Plot UV index map from E_<case>.UV output.")
    parser.add_argument("--input", default="Output/E_Calc_12Dec.UV", type=Path, help="Path to E_<case>.UV file")
    parser.add_argument("--grid", default="Input/grid", type=Path, help="Path to grid definition (Input/grid)")
    parser.add_argument("--output", default="uvi_map.png", type=Path, help="Output PNG filename")
    args = parser.parse_args()

    lats, lons, _, _ = read_grid(args.grid)
    data = read_uv_field(args.input, nlat=len(lats), nlon=len(lons), lats_expected=lats)

    uv_grid = xr.DataArray(data, coords={"lat": lats, "lon": lons}, dims=("lat", "lon"), name="uv_index")
    plot_uv(uv_grid, args.output)
    print(f"Saved map to {args.output}")


if __name__ == "__main__":
    main()
