import geopandas as gpd
import numpy as np

# Filter the lines
# df = gpd.read_file("/home/efvega/.local/share/hackathon_month9/macrostrat_lines.json")
# filtered_df = df[df['type'].str.contains('fault')]
# filtered_df.to_file("/home/efvega/.local/share/hackathon_month9/macrostrat_faults.geojson", driver='GeoJSON')

df = gpd.read_file("/home/efvega/.local/share/hackathon_month9/macrostrat.json")


lithologies = ['limestone', 'dolostone', 'dolomite', 'carbonate']
filtered_df = df[df['lith'].str.contains('|'.join(lithologies))]

age_subset = filtered_df[np.logical_and(filtered_df['best_age_bottom'] <= 538.8, filtered_df['best_age_top'] >= 458.4)]
# Then age filter

age_subset.to_file("/home/efvega/.local/share/hackathon_month9/macrostrat_units_subset.geojson", driver='GeoJSON')