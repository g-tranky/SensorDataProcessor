# Script modified from Phil Stepanian script to linkNARR to bird tracks
# Script optimized by Adam Gernes


# Import packages
from datetime import date, datetime, time

from netCDF4 import Dataset
import numpy as np
import pandas as pd

# Define the csv file containing locations and times
# filename = './grace_lightData_edit.csv'
filename = "./df_all.csv"


# Open and read the flight track data into a pandas data frame
df = pd.read_csv(filename)
df["utc_time"] = pd.to_datetime(df["utc_time"])

# Convert datetime to the NARR time units of "hours since 1800 01 01 00:00:00"
print("Converting utc time to NARR time...")
base_narr_time = datetime.combine(date(1800, 1, 1), time(0, 0, 0))
seconds_to_hours = (df["utc_time"] - base_narr_time).dt.seconds / 3600
days_to_hours = (df["utc_time"] - base_narr_time).dt.days * 24
df["narr_time"] = days_to_hours + seconds_to_hours
print("Done converting time.")


def find_lat_long(row, narr_lat, narr_long):
    # this could be more efficient, left as is
    track_lat = row["lat"]
    track_long = row["long"]
    NS_distance = (narr_lat - track_lat) * 111
    WE_distance = (narr_long - track_long) * (
        0.001
        * (
            np.pi
            * 6378137
            * np.cos(np.pi / 180 * narr_lat)
            / 180
            / np.sqrt(1 - 0.006694380004261 * np.sin(np.pi / 180 * narr_lat) ** 2)
        )
    )
    column_yNARRIdx_now, column_xNARRIdx_now = np.unravel_index(
        np.argmin(np.sqrt(NS_distance**2 + WE_distance**2)), narr_lat.shape
    )
    return column_yNARRIdx_now, column_xNARRIdx_now


def get_cloud_data_for_year(df: pd.DataFrame, year: int):
    # only get rows from this year
    df_year = df[df["utc_time"].dt.year == year]

    # read nc file for this year
    fname = f"tcdc.{year}.nc"
    print(f"Reading {fname}...")
    ncid = Dataset(fname, "r")

    # Get the data we want from the netCDF4 dataset
    narr_timeMono = ncid.variables["time"][
        :
    ]  # (time) [hours since 1800 01 01 00:00:00]
    narr_lat = ncid.variables["lat"][:]  # (y, x)
    narr_long = ncid.variables["lon"][:]  # (y, x)
    narr_cloud = ncid.variables["tcdc"][:]  # (time, y, x) [%]
    ncid.close()

    # find the closest time to the cloud data
    print("Finding closest times for cloud data...")
    time_mono = df_year["narr_time"]
    min_time_diff_idx = np.argmin(
        np.abs(time_mono.values[:, np.newaxis] - narr_timeMono), axis=1
    )

    print("Finding lat and long indexes for cloud data...")
    narr_idx = df_year.apply(find_lat_long, args=(narr_lat, narr_long), axis=1)

    cloud_data = narr_cloud[min_time_diff_idx, narr_idx.str[0], narr_idx.str[1]]
    return pd.DataFrame({"cloud": cloud_data}, index=df_year.index)


years = [2022, 2023]
df["cloud"] = pd.NA
for year in years:
    print(f"Processing Year {year}...")
    df.update(get_cloud_data_for_year(df, year))

save_file = f"df_all_cloud-{datetime.now().isoformat()}.csv"
print(f"Saving results as {save_file}...")
df.to_csv(save_file)
print("Done.")
