import pytest
import sys
if sys.version_info < (3, 9):
    from importlib_resources import files
else:
    from importlib.resources import files
from pathlib import Path


def test_macrostrat_from_bounds():
    from macrostratpy import macrostrat_from_bounds
    import geopandas as gpd
    import platformdirs

    # fp = files("macrostratpy.data") / "AlaskaTungstenSkarn.gpkg"
    fp = Path("/home/efvega/Downloads/h3_mvt_reg_aoi.gpkg")

    data = gpd.read_file(fp)
    data_4326 = data.to_crs(epsg=4326, inplace=False)
    data_4326 = data.geometry.to_crs(epsg=4326)

    bounds = {
        'n': data_4326.bounds['maxy'].max(),
        's': data_4326.bounds['miny'].min(),
        'e': data_4326.bounds['maxx'].max(),
        'w': data_4326.bounds['minx'].min()
    }

    # output_filename = platformdirs.user_data_path() / "hackathon_month9" / "macrostrat_lines.json"
    # # print(f"Saving data to {output_filename}")
    # # # output_path = "/home/efvega/data/macrostrat/upper_midwest_packaged/macrostrat_upper_midwest.json"
    # zoom_level = 6
    # #
    # #
    #
    # # layername = "units"
    # layername = "lines"
    #
    # macrostrat_from_bounds(bounds = bounds, output_path=output_filename, layername=layername, zoom_level=zoom_level)

    output_filename = platformdirs.user_data_path() / "hackathon_month9" / "ta1_map_units.json"
    layername="units"
    zoom_level=8

    macrostrat_from_bounds(
        bounds = bounds,
        output_path=output_filename,
        layername=layername,
        macrostrat_server='https://tileserver.development.svc.macrostrat.org/',
        service = 'map',#'maps',
        params = 'source_id=1724',#'lithology=siliciclastic&lithology=gravel',
        zoom_level=zoom_level,
    )

    #macrostrat_from_bounds(bounds=bounds, output_path=output_filename, layername=layername, zoom_level=zoom_level)
if __name__ == "__main__":
    test_macrostrat_from_bounds()