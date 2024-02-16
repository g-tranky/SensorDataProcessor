from sys import argv
import sqlite3 as sql
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import List, Union

from pandas import DataFrame, concat, read_csv


@dataclass
class DatMetaData:
    lat: float = 0
    long: float = 0
    elev_m: float = 0
    local_timezone: str = ""
    serial_num: str = ""
    hardware_id: str = ""


def pair_dat_csv_files(data_root: Path) -> List[Union[Path, Path]]:
    csv_paths = list(data_root.glob("**/*.csv"))
    dat_paths = list(data_root.glob("**/*.dat"))
    dat_dict = {path.name: path for path in dat_paths}
    pairs = []
    for path in csv_paths:
        csv_dat_key = path.name[:-4] + ".dat"
        if csv_dat_key in dat_dict:
            pairs.append((path, dat_dict[csv_dat_key]))
    return pairs


def parse_dat_comments(dat_file: Path) -> DatMetaData:
    print("Processing... ", dat_file)
    keys = {
        "# Local timezone": 0,
        "# SQM serial number": 1,
        "# SQM hardware identity": 2,
        "# Position (lat, lon, elev(m))": 3,
    }
    meta = DatMetaData()
    with open(dat_file) as f:
        lines = f.readlines()
    for line in [line for line in lines if line[0] == "#"]:
        i = line.find(":")
        key = line[:i]
        value = line[i + 1 :].strip()
        if key in keys:
            if keys[key] == 0:
                meta.local_timezone = value
            elif keys[key] == 1:
                meta.serial_num = value
            elif keys[key] == 2:
                meta.hardware_id = value
            elif keys[key] == 3:
                values = [float(s.strip()) for s in value.split(",")]
                meta.lat = values[0]
                meta.long = values[1]
                meta.elev_m = values[2]
    return meta


def clean_csv(csv_file: Path) -> None:
    with open(csv_file, "r+") as f:
        content = f.read()
        f.seek(0)
        content = content.replace(";", ",").replace(",\n", "\n")
        f.write(content)
        f.truncate()


def csv_to_df(csv_file: Path) -> DataFrame:
    print("Processing... ", csv_file)
    clean_csv(csv_file)
    new_data = {"location": csv_file.name[: csv_file.name.find("_")]}
    df = read_csv(csv_file)
    return concat([df, DataFrame(new_data, index=df.index)], axis=1)


def build_df(meta: DatMetaData, df: DataFrame) -> DataFrame:
    new_df = concat([df, DataFrame(asdict(meta), index=df.index)], axis=1)
    new_df.columns = new_df.columns.str.strip()
    new_df.rename(
        columns={
            "UTC Date & Time": "utc_time",
            "Local Date & Time": "local_time",
            "Temperature": "temperature",
            "Voltage": "voltage",
            "MSAS": "msas",
            "Record type": "record_type",
            "MoonPhaseDeg": "moon_phase_deg",
            "MoonElevDeg": "moon_elev_deg",
            "MoonIllum%": "moon_illum_perc",
            "SunElevDeg": "sun_elev_deg",
        },
        inplace=True,
    )
    return new_df


def insert_df_to_db(df: DataFrame, con: sql.Connection) -> None:
    df.to_sql("sensor_data", con, if_exists="replace", index=False)


def build_db() -> sql.Connection:
    con = sql.connect("Data.db")
    cur = con.cursor()
    cur.execute(
        """
            CREATE TABLE IF NOT EXISTS sensor_data
            (ID INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
            utc_time TIMESTAMP NOT NULL,
            local_time TIMESTAMP NOT NULL,
            location TEXT NOT NULL,
            temperature FLOAT NOT NULL,
            voltage FLOAT NOT NULL,
            msas FLOAT NOT NULL,
            record_type INT NOT NULL,
            moon_phase_deg FLOAT NOT NULL,
            moon_elev_deg FLOAT NOT NULL,
            moon_illum_perc FLOAT NOT NULL,
            sun_elev_deg FLOAT NOT NULL,
            lat FLOAT NOT NULL,
            long FLOAT NOT NULL,
            elev_m FLOAT NOT NULL,
            local_timezone TIMESTAMP NOT NULL,
            serial_num VARCHAR(10) NOT NULL,
            hardware_id VARCHAR(10) NOT NULL)
            """
    )
    return con


def insert_data_to_db(root_dir: Path, con: sql.Connection) -> None:
    paths = pair_dat_csv_files(root_dir)
    data = (
        build_df(parse_dat_comments(dat_path), csv_to_df(csv_path))
        for csv_path, dat_path in paths
    )
    insert_df_to_db(concat(data), con)


if __name__ == "__main__":
    if len(argv) <= 1:
        print("Usage: python main.py <data directory>")
        exit()
    insert_data_to_db(Path(argv[1]), build_db())
