# Radar-Interference-Prediction
This a Python program to compute the basic transmission loss according to ITU-R P.452-16 recommendation.


# Usage
Run main.py to initialize the GUI. 


"The main module used is p452.py, and all the other  subfunctions of it are located in '/p452_modules' folder. 

The input parameters for p452 calculations are passed through the GUI, but regarding the path elevation profile and radio-climate zones they may be loaded either from a xls(x) file, like the test examples provided by ITU, or calculated from SRTM 30m data (heights along the path) as long as there are the necessary hgt files inside the /hgt folder.



The .xlsx path profile file should consist of 3 columns, describing the di, hi, and radio-climate zone arrays respectively. When computing automatically given the SRTM data files, radio-climate zone for each point can be set to a fixed 'A2' value (see p452 documentation for details) or you can set a rule to calculate the zone depending on elevation of each point (hi) (File -> Settings).\n\nRx station can be " \
          "selected from a csv file 'radars.csv' that should be located in the root directory of the program, " \
          "with columns: [name, lat, lon, height above ground, gain], or alternativley can be defined manually. " \
          "A plot with the path elevation profile vs distance is showed, with earth curvature and refraction " \
          "impact on line of sight, if this is selected in settings. SRTM hgt files can be downloaded from:\n" \
          "https://search.earthdata.nasa.gov/search\n (an earthdata account should be created first)\nIf any hgt " \
          "is missing from directory, there is the choice to fill with 0s, or a warning appears with the list of all " \
          "missing hgt files.\n" \
          "\nFinally, transmitter gain can be calculated using the selected antenna radiation pattern, from a pandas " \
          "dataframe stored as 'database.hdf' in the root directory. The dataframe has columns " \
          "[antenna name, el.tilt, gain, pattern], whre pattern is a (2, 361) shaped numpy array containing " \
          "the horizontal and vertical radiation pattern of antennas. An example database is given, but you can use " \
          "the module 'create_antenna_db.py' to use other radiation patterns (.csv files) to create your " \
          "own database (no gui for now)."
