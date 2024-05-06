import pytest
import sys
if sys.version_info < (3, 9):
    from importlib_resources import files
else:
    from importlib.resources import files


def test_macrostrat_from_bounds():
    from macrostratpy import macrostrat_from_bounds
    import geopandas as gpd
    import platformdirs

    fp = files("macrostratpy.data") / "AlaskaTungstenSkarn.gpkg"

    data = gpd.read_file(fp)
    data_4326 = data.to_crs(epsg=4326, inplace=False)

    bounds = {
        'n': data_4326.bounds['maxy'].max(),
        's': data_4326.bounds['miny'].min(),
        'e': data_4326.bounds['maxx'].max(),
        'w': data_4326.bounds['minx'].min()
    }

    output_filename = platformdirs.user_data_path() / "AlaskaTungstenSkarn" / "macrostrat_alaska.json"

    print(f"Saving data to {output_filename}")
    # output_path = "/home/efvega/data/macrostrat/upper_midwest_packaged/macrostrat_upper_midwest.json"
    zoom_level = 0


    macrostrat_from_bounds(bounds = bounds, output_path=output_filename, zoom_level=zoom_level)