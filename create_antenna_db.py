import pandas as pd
from pathlib import Path
import numpy as np

# define folder containg csv files.
# Each csv file contains horizontal and vertical radiation pattern in columns (0, 1) respectively, and follow the
# filename format "ANTENNANAME_ELTILT_MAXGAIN.csv"
# example file patterns are located in '\pattern' folder
folder = Path.joinpath(Path.cwd(), "patterns")
files = folder.glob("*.csv")

# create dataframe:
df = pd.DataFrame(columns=['Antenna', 'Tilt', 'Gain', 'Pattern'])
for file in files:
    fnamesplit = file.stem.split('_')
    antenna = fnamesplit[0]
    tiltval = float(fnamesplit[1])
    gainval = float(fnamesplit[2])
    pattern = np.genfromtxt(file, 'float32', delimiter=';')
    temp_df = pd.DataFrame({'Antenna': [antenna], 'Tilt': [tiltval], 'Gain': [gainval], 'Pattern': [pattern]})
    df = df.append(temp_df, ignore_index=True)

# set multindex to select antenna using name and electrical tilt:
df.set_index(['Antenna', 'Tilt'], inplace=True)
df.sort_index(inplace=True)

# store dataframe:
df.to_hdf('database.hdf', key='db')

# to extract antenna pattern and specs, use:
# row = antennas_db.loc[(antenna, eltilt)]
# maxgain = row['Gain'].item()
# horizontal = row['Pattern'].item()[:, 0]
# vertical = row['Pattern'].item()[:, 1]

