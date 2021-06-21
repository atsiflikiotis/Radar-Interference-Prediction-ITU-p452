# Radar-Interference-Prediction
This a Python GUI program to compute the basic transmission loss according to ITU-R P.452-16 recommendation, compatible with NASA's SRTM data to generate the elevation profile path.

# Initilization
Run main.py to initialize the GUI. 


# Usage
"The main module used is p452.py (implementation of ITU P.452-16 recommendation) and all the other subfunctions of ITU are located in '/p452_modules' folder. 

The input parameters for p452 calculations are passed through the GUI, but regarding the path elevation profile and radio-climate zones they may be loaded either from a xls(x) file, like the test examples provided by ITU, or calculated from SRTM 30m data (heights along the path) as long as there are the necessary hgt files inside the /hgt folder.



The .xlsx path profile file should consist of 3 columns, describing the di, hi, and radio-climate zone arrays respectively. When computing automatically given the SRTM data files, radio-climate zone for each point can be set to a fixed 'A2' value (see p452 documentation for details) or you can set a rule to calculate the zone depending on elevation of each point (hi) (File ->  Settings).



Rx station can be selected from a csv file 'radars.csv' that should be located in the root directory of the program, with columns: [name, lat, lon, height above ground, gain], or alternativley can be defined manually. A plot with the path elevation profile vs distance is showed, with earth curvature and refraction impact on line of sight, if this is selected in settings. SRTM hgt files can be downloaded from:


https://search.earthdata.nasa.gov/search

(an earthdata account should be created first)

If any hgt file is missing from directory, there is the choice to fill with 0s, or a warning appears with the list of all missing hgt files.


Finally, transmitter gain can be calculated using the selected antenna radiation pattern, from a pandas dataframe stored as 'database.hdf' in the root directory. The dataframe has columns: [antenna name, el.tilt, gain, pattern], where pattern is a (2, 361) shaped numpy array containing the horizontal and vertical radiation pattern of antennas (relative gain loss (dB) at each angle). An example database and antennas patterns are given, but you can use the module 'create_antenna_db.py' to use other radiation patterns (.csv files) to create your own database (no gui toll for database as of now)."

Screenshot of the main GUI:

https://raw.githubusercontent.com/atsiflikiotis/Radar-Interference-Prediction/master/Screenshot1.png?token=ALAVZ555NSU2JNCNGOBULL3AZ3VRM

# More info:
ITU-R P.452: Prediction procedure for the evaluation of interference between stations on the surface of the Earth at frequencies above about 0.1 GHz

visit for more info: 

https://www.itu.int/rec/R-REC-P.452-16-201507-I/en


There is a matlab implementation already published by ITU-R Study Group 3, for the main calculation of transmission loss according to p452 (without SRTM path profile generation). You can download it from below url (validation examples also provided in the site):

https://www.itu.int/en/ITU-R/study-groups/rsg3/Pages/iono-tropo-spheric.aspx
