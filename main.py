import re
import tkinter as tk
import tkinter.filedialog
import tkinter.ttk as ttk
from pathlib import Path
from tkinter import font
from tkinter.messagebox import showerror, showinfo

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

import srtm_path as srtm
from p452 import p452_loss

version = 'version 1.0 2020'


def dmstodec(coords):
    if re.search('[SsWw]', coords):
        sign = -1
    else:
        sign = 1

    coords = coords.strip('SsWw ')

    splitlist = re.split(r'[°\s\'\"]', coords)
    splitlist = [*filter(len, splitlist)]

    if len(splitlist) == 1:
        dec = sign * float(splitlist[0])
    elif len(splitlist) == 3:
        dec = sign * float(splitlist[0]) + float(splitlist[1]) / 60 + float(splitlist[2]) / 3600
    else:
        raise ValueError

    return dec


def readme():
    msg = "This is a program to compute the basic transmission loss according to ITU-R P.452-16 " \
          "recommendation (07/2015).\n" \
          "\nThe main module used is p452.py, and all the other  subfunctions of it are located " \
          "in '/p452_modules' folder.\n\nThe input parameters for p452 calculations are passed through the GUI, but " \
          "regarding the path elevation profile and radio-climate zones they may be loaded either from " \
          "a xls(x) file, like the test examples provided by ITU, or calculated from SRTM 30m data (heights along the " \
          "path) as long as there are the necessary hgt files inside the /hgt folder.\n\n" \
          "The .xlsx path profile file should consist of 3 columns, describing the di, hi, and radio-climate zone " \
          "arrays respectively. When computing automatically given the SRTM data files, radio-climate zone for each " \
          "point can be set to a fixed 'A2' value (see p452 documentation for details) or you can set a rule to " \
          "calculate the zone depending on elevation of each point (hi) (File -> Settings).\n\nRx station can be " \
          "selected from a csv file 'radars.csv' that should be located in the root directory of the program, " \
          "with columns: [name, lat, lon, height above ground, gain], or alternativley can be defined manually. " \
          "A plot with the path elevation profile vs distance is showed, with earth curvature and refraction " \
          "impact on line of sight, if this is selected in settings. SRTM hgt files can be downloaded from:\n" \
          "https://search.earthdata.nasa.gov/search\n (an earthdata account should be created first)\nIf any hgt " \
          "is missing from directory, there is the choice to fill with 0s, or a warning appears with the list of all " \
          "missing hgt files.\n" \
          "\nFinally, transmitter gain can be calculated using the selected antenna radiation pattern, from a pandas " \
          "dataframe stored as 'database.hdf' in the root directory. The dataframe has columns " \
          "[antenna name, el.tilt, gain, pattern], where pattern is a (2, 361) shaped numpy array containing " \
          "the horizontal and vertical radiation pattern of antennas (relative gain loss (dB) at each angle). An " \
          "example database and antennas pattern are given, but you can use the module 'create_antenna_db.py' to use " \
          "other radiation patterns (.csv files) to create your own database (no gui tool for antenna database as of now)."
    showinfo('README', msg)

    # '- The main implementation of the recommendation is MATLAB function'
    # '  tl_p452.m placed in this folder that can be used independently'
    # '  of this Graphical User Interface but needs the functions '
    # '  defined in the folder ./src.'


def aboutcommand():
    try:
        with open('LICENSE.md', 'r') as file:
            text = file.read()
    except FileNotFoundError:
        text = "License file cannot be loaded."

    text = version + '\n\n' + text
    showinfo('About', text)


def openimage(fignumber, title):
    topwindow = tk.Toplevel()
    topwindow.title(title)
    fname = f'Fig{fignumber}.gif'
    canvas = tk.Canvas(topwindow, width=1200, height=700)
    canvas.grid(row=0, column=0)
    img = tk.PhotoImage(file=fname)
    canvas.create_image(10, 10, image=img, anchor=tk.NW)
    # keep reference
    canvas.img = img


class ValidEntry(ttk.Entry):

    def __init__(self, master, validtype, *args, **kwargs):
        super().__init__(master, *args, justify='center', **kwargs)
        vcmd = (self.register(self.onvalidate), '%d', '%P', '%V')
        if validtype == 'coords':
            self.configure(validate='focusout')
        else:
            self.configure(validate='all')
        self.configure(validatecommand=vcmd)
        self.validtype = validtype

    def onvalidate(self, action, value, event):
        if self.validtype == 'coords':
            if not value:
                return True
            try:
                dmstodec(value)
                return True
            except ValueError:
                showerror('Error', 'Invalid coordinates format, '
                                   'must be decimal or DMS in format: D Min Secs(opt: N/S or E/W)')
                self.delete(0, tk.END)
                self.insert(0, 0)
                return False
        else:
            if event == 'key':
                if action != '1':
                    return True
                try:
                    float(value)
                    return True
                except ValueError:
                    self.bell()
                    return False
            elif event == 'focusout':
                try:
                    float(value)
                    return True
                except ValueError:
                    self.delete(0, tk.END)
                    self.insert(0, 0)
                    showerror('Error', 'You have to enter a number to coordinates entry, temporary filled with 0 value')
                    return False
            else:
                return True


class MainGUI(ttk.Frame):

    def __init__(self, master, **kwargs):
        super().__init__(master, **kwargs)

        master.title('Radar Interference Prediction (ITU P.472)')
        master.geometry('1400x700')
        # master.iconbitmap('icon.ico')

        try:
            self.antennas = pd.read_hdf('database.hdf', key='db')
            antennaslist = self.antennas.index.unique(level='Antenna').tolist()
            self.antennasloaded = True
        except FileNotFoundError:
            self.antennasloaded = False
            antennaslist = ['No dB file found']

        try:
            self.radars = pd.read_csv('radars.csv', delimiter=';')
            radarslist = self.radars['Name'].tolist()
            radarsloaded = True
        except FileNotFoundError:
            radarsloaded = False
            radarslist = ['No .csv list found']

        menubar = tk.Menu(master)
        filemenu = tk.Menu(menubar, tearoff=0)
        radiomapsmenu = tk.Menu(menubar, tearoff=0)
        aboutmenu = tk.Menu(menubar, tearoff=0)

        filemenu.add_command(label='Settings', command=self.opensettings)

        filemenu.add_command(label='Exit', command=quit_me)
        radiomapsmenu.add_command(label='Median Annual values of ΔΝ',
                                  command=lambda: openimage(1, 'Median Annual values of ΔΝ'))
        radiomapsmenu.add_command(label='Average Annual values of ΔΝ',
                                  command=lambda: openimage(2, 'Average Annual values of ΔΝ'))
        radiomapsmenu.add_command(label='Median sea-level surface refractivity',
                                  command=lambda: openimage(3, 'Median sea-level surface refractivity'))
        radiomapsmenu.add_command(label='Sea-level surface refractivity',
                                  command=lambda: openimage(4, 'sea-level surface refractivity'))

        aboutmenu.add_command(label='About', command=aboutcommand)
        aboutmenu.add_command(label='Help', command=readme)

        menubar.add_cascade(label='File', menu=filemenu)
        menubar.add_cascade(label='Radio-meteo maps', menu=radiomapsmenu)
        menubar.add_cascade(label='Help', menu=aboutmenu)
        master.config(menu=menubar)

        master.rowconfigure(0, weight=0)
        master.rowconfigure(1, weight=1)
        master.rowconfigure(2, weight=0)

        master.columnconfigure(0, weight=0)
        master.columnconfigure(1, weight=0)
        master.columnconfigure(2, weight=1)
        master.columnconfigure(3, weight=0)

        # default settings and settings variables, stores as class instances
        self.zoneselvar = tk.StringVar()
        self.earthcurvvar = tk.StringVar()
        self.fillvoidsrtmvar = tk.IntVar()

        self.zonesel = 'auto'  # as opposed to 'fixed' (=A1, inland)
        self.earthcurv = 'curved'  # as opposed to 'flat' earth in path profile plot
        self.fillvoidsrtm = 0  # don't fill with 0 missing srtm data, raise an error, set 1 to fill with 0s
        self.seathreshvalue = 0  # when 'auto' fill radio-climate zones, set sea zone when h <= value
        self.coastlandthreshvalue = 50  # when 'auto' fill climate zones, set coastland zone when h <= value
        self.coordsamples = 700  # maximum number of samples, when generating path profile from SRTM files
        self.validpathloaded = False  # turned to valid if manual path is selected and valid xls file loaded

        # Build frames
        paramframe = ttk.LabelFrame(master, text='Parameters')
        paramframe.grid(row=1, column=1, sticky='nsew', pady=(5, 20), padx=15)

        frame2 = tk.Frame(master)
        frame2.grid(row=1, column=2, sticky='nsew', pady=5, padx=15)
        frame2.rowconfigure(0, weight=0)  # transmitter frame
        frame2.rowconfigure(1, weight=0)  # receiver frame
        frame2.rowconfigure(2, weight=0)  # button frame
        frame2.rowconfigure(3, weight=1)  # results frame
        frame2.rowconfigure(4, weight=0)
        frame2.columnconfigure(0, weight=1)

        # parametersframe
        freqlabel = tk.Label(paramframe, text='Frequency (GHz):')
        freqlabel.grid(row=0, column=0, padx=5, pady=(7, 4), sticky='nsw')
        self.freqentry = ValidEntry(paramframe, 'float', width=6)
        self.freqentry.insert(0, 2.6)
        self.freqentry.grid(row=0, column=1, padx=7, pady=(7, 4), sticky='nse')
        timelabel = tk.Label(paramframe, text='Time percentage (%):')
        timelabel.grid(row=1, column=0, padx=5, pady=4, sticky='nsw')
        self.timeentry = ValidEntry(paramframe, 'float', width=6)
        self.timeentry.insert(0, 0.001)
        self.timeentry.grid(row=1, column=1, padx=7, pady=4, sticky='nse')
        self.poldict = {'Vertical polarization': 2, 'Horizontal polarization': 1, 'Slant polarization': 3}
        self.polarbox = ttk.Combobox(paramframe, values=[*self.poldict.keys()], width=21, state='readonly')
        self.polarbox.grid(row=2, column=0, columnspan=2, padx=5, pady=6, sticky='ns')
        self.polarbox.current(0)

        pathframe = tk.LabelFrame(paramframe, text='Path profile')
        pathframe.grid(row=3, padx=5, pady=8, sticky='nsew', columnspan=2)
        self.pathselectvar = tk.IntVar(0)
        pathrbtn1 = ttk.Radiobutton(pathframe, text='Generate from SRTM data', variable=self.pathselectvar, value=0,
                                    command=self.pathrbtn)
        pathrbtn2 = ttk.Radiobutton(pathframe, text='Load from file:', variable=self.pathselectvar, value=1,
                                    command=self.pathrbtn)
        pathrbtn1.grid(row=0, column=0, sticky='nsw', padx=5, pady=(5, 3))
        pathrbtn2.grid(row=1, column=0, sticky='nsw', padx=(5, 25), pady=(3, 5))
        self.pathloadbtn = ttk.Button(pathframe, text='Load', command=self.loadpathfromfile, width=8, state='disabled')
        self.pathloadbtn.grid(row=1, column=2, sticky='nsw', padx=(8, 5), pady=(3, 5))

        meteoframe = tk.LabelFrame(paramframe, text='Meteorological data')
        meteoframe.grid(row=4, column=0, columnspan=2, padx=5, pady=8, sticky='nsew')
        meteoframe.columnconfigure(0, weight=0)
        meteoframe.columnconfigure(1, weight=1)
        pressurelabel = tk.Label(meteoframe, text='Pressure (hPa):')
        pressurelabel.grid(row=0, column=0, sticky='nsw', padx=5, pady=4)
        self.pressureentry = ValidEntry(meteoframe, 'float', width=7)
        self.pressureentry.insert(0, 1013)
        self.pressureentry.grid(row=0, column=1, sticky='nse', padx=7, pady=4)
        templabel = tk.Label(meteoframe, text='Temperature (°C):')
        templabel.grid(row=1, column=0, sticky='nsw', padx=5, pady=4)
        self.tempentry = ValidEntry(meteoframe, 'float', width=7)
        self.tempentry.insert(0, 25)
        self.tempentry.grid(row=1, column=1, sticky='nse', padx=7, pady=4)
        dnlabel = tk.Label(meteoframe, text='ΔΝ (N-units/km):')
        dnlabel.grid(row=2, column=0, sticky='nsw', padx=5, pady=4)
        self.dnentry = ValidEntry(meteoframe, 'float', width=7)
        self.dnentry.insert(0, 50)
        self.dnentry.grid(row=2, column=1, sticky='nse', padx=7, pady=4)
        n0label = tk.Label(meteoframe, text='N0 (N-units/km):')
        n0label.grid(row=3, column=0, sticky='nsw', padx=5, pady=4)
        self.n0entry = ValidEntry(meteoframe, 'float', width=7)
        self.n0entry.insert(0, 350)
        self.n0entry.grid(row=3, column=1, sticky='nse', padx=7, pady=4)

        othersframe = tk.LabelFrame(paramframe, text='Other parameters')
        othersframe.grid(row=5, column=0, columnspan=2, sticky='nsew', padx=5, pady=7)
        othersframe.columnconfigure(0, weight=0)
        othersframe.columnconfigure(1, weight=1)
        dctlabel = tk.Label(othersframe, text='dct (km):')
        dctlabel.grid(row=0, column=0, sticky='nsw', padx=5, pady=4)
        self.dctentry = ValidEntry(othersframe, 'float', width=7)
        self.dctentry.insert(0, 500)
        self.dctentry.grid(row=0, column=2, sticky='nse', padx=5, pady=4)
        dcrlabel = tk.Label(othersframe, text='dcr (km):')
        dcrlabel.grid(row=1, column=0, sticky='nsw', padx=5, pady=4)
        self.dcrentry = ValidEntry(othersframe, 'float', width=7)
        self.dcrentry.insert(0, 500)
        self.dcrentry.grid(row=1, column=2, sticky='nse', padx=5, pady=4)

        self.clutterdict = {'None': (0, 0), 'High crop fields': (4, 0.1), 'Park land': (4, 0.1),
                            'Irregularly spaced sparse trees': (4, 0.1), 'Orchard (regularly spaced)': (4, 0.1),
                            'Sparse houses': (4, 0.1), 'Village center': (5, 0.07), 'Deciduous trees': (15, 0.05),
                            'Mixed tree forest': (15, 0.05), 'Coniferous trees': (20, 0.05),
                            'Tropical rain forest': (20, 0.03), 'Suburban': (9, 0.025), 'Dense suburban': (12, 0.02),
                            'Urban': (20, 0.02), 'Dense Urban': (25, 0.02), 'High-rise urban': (35, 0.02),
                            'Industrial zone': (20, 0.05)}

        txclutterlbl = tk.Label(othersframe, text='Tx clutter type:')
        rxclutterlbl = tk.Label(othersframe, text='Rx clutter type:')
        self.txclutterbox = ttk.Combobox(othersframe, values=[*self.clutterdict.keys()], width=27,
                                         state='readonly', justify='center')
        self.rxclutterbox = ttk.Combobox(othersframe, values=[*self.clutterdict.keys()], width=27,
                                         state='readonly', justify='center')

        txclutterlbl.grid(row=2, column=0, padx=5, pady=(6, 1), sticky='nsw')
        self.txclutterbox.grid(row=3, column=0, columnspan=3, padx=5, pady=1, sticky='nsw')
        rxclutterlbl.grid(row=4, column=0, padx=5, pady=(6, 1), sticky='nsw')
        self.rxclutterbox.grid(row=5, column=0, columnspan=3, padx=5, pady=1, sticky='nsw')
        self.txclutterbox.current(0)
        self.rxclutterbox.current(0)

        # frame2 subframes:
        txframe = ttk.LabelFrame(frame2, text='Transmitter')
        txframe.grid(row=0, column=0, padx=5, pady=(0, 6), sticky='nsew')
        rxframe = ttk.LabelFrame(frame2, text='Receiver station')
        rxframe.grid(row=1, column=0, padx=5, pady=6, sticky='nsew')
        calcbtnframe = tk.Frame(frame2)
        calcbtnframe.grid(row=2, column=0, pady=2, padx=5, sticky='nsew')
        calcbtnframe.columnconfigure(0, weight=1)
        resultsframe = ttk.LabelFrame(frame2, text='Results')
        resultsframe.grid(row=3, column=0, pady=(2, 20), padx=5, sticky='nsew')

        # rx frame
        txposframe = tk.Frame(txframe)
        txposframe.grid(row=0, column=0, columnspan=9, padx=5, pady=5, sticky='nsew')
        txlatlbl = tk.Label(txposframe, text='Lat (φ):')
        txlatlbl.grid(row=0, column=0, padx=(5, 2), pady=5, sticky='nsw')
        self.txlatentry = ValidEntry(txposframe, 'coords', width=9)
        self.txlatentry.grid(row=0, column=1, padx=2, pady=5, sticky='nsw')
        txlonlbl = tk.Label(txposframe, text='Lon (λ):')
        txlonlbl.grid(row=0, column=2, padx=(30, 2), pady=5, sticky='nsw')
        self.txlonentry = ValidEntry(txposframe, 'coords', width=9)
        self.txlonentry.grid(row=0, column=3, padx=2, pady=5, sticky='nsw')
        txheightlbl = tk.Label(txposframe, text='Height (agl, m):')
        txheightlbl.grid(row=0, column=4, padx=(35, 2), pady=5, sticky='nsw')
        self.txheightentry = ValidEntry(txposframe, 'float', width=5)
        self.txheightentry.grid(row=0, column=5, padx=2, pady=5, sticky='nsw')
        txpowerlbl = tk.Label(txposframe, text='Tx Power (dBm):')
        txpowerlbl.grid(row=0, column=6, padx=(35, 2), pady=5, sticky='nsw')
        self.txpowerentry = ValidEntry(txposframe, 'float', width=6)
        self.txpowerentry.grid(row=0, column=7, padx=2, pady=5, sticky='nsw')

        self.txgainselectvar = tk.IntVar()
        if self.antennasloaded:
            self.txgainselectvar.set(0)
        else:
            self.txgainselectvar.set(1)

        self.txgainrbtn1 = ttk.Radiobutton(txframe, text='Calculate Gain:', variable=self.txgainselectvar, value=0,
                                           command=self.txgainrbtn)
        self.txgainrbtn1.grid(row=2, column=0, padx=5, pady=3, rowspan=1, sticky='nsw')
        self.txgainrbtn2 = ttk.Radiobutton(txframe, text='Fixed Gain (dBi):', variable=self.txgainselectvar, value=1,
                                           command=self.txgainrbtn)
        self.txgainrbtn2.grid(row=3, column=0, padx=5, pady=(10, 3), sticky='nsw')
        antennalbl = tk.Label(txframe, text='Antenna type')
        antennalbl.grid(row=1, column=1, sticky='ns', padx=(10, 2), pady=1)
        eltiltlbl = tk.Label(txframe, text='El. Tilt (deg)')
        eltiltlbl.grid(row=1, column=2, sticky='ns', padx=8, pady=1)
        azimlbl = tk.Label(txframe, text='Azimuth (deg)')
        azimlbl.grid(row=1, column=3, sticky='ns', padx=8, pady=1)
        mechtiltlbl = tk.Label(txframe, text='Mech. Tilt (deg)')
        mechtiltlbl.grid(row=1, column=4, sticky='ns', padx=8, pady=1)
        self.txgainentry = ValidEntry(txframe, 'float', width=7)
        self.txgainentry.grid(row=3, column=1, padx=15, pady=(10, 3), sticky='sw')

        self.antennabox = ttk.Combobox(txframe, values=antennaslist, state='readonly', width=30, justify='center')
        self.antennabox.grid(row=2, column=1, padx=(10, 2))
        self.antennabox.current(0)
        self.antennabox.bind("<<ComboboxSelected>>", self.updatetilts)

        self.tiltsbox = ttk.Combobox(txframe, values=[], state='readonly', width=4, justify='center')
        self.tiltsbox.grid(row=2, column=2, padx=8)
        self.updatetilts("<<ComboboxSelected>>")

        self.azimentry = ValidEntry(txframe, 'float', width=6)
        self.azimentry.grid(row=2, column=3, sticky='ns', padx=8, pady=2)

        self.mechtiltentry = ValidEntry(txframe, 'float', width=5)
        self.mechtiltentry.grid(row=2, column=4, sticky='ns', padx=8, pady=2)

        if self.antennasloaded:
            self.txgainselectvar.set(0)
        else:
            self.txgainselectvar.set(1)

        # rxframe:
        self.rxposvar = tk.IntVar()
        rxposbtn1 = ttk.Radiobutton(rxframe, text='Rx from list:', variable=self.rxposvar, value=0,
                                    command=self.rxselection)
        rxposbtn1.grid(row=0, column=0, padx=5, pady=5, sticky='nsw')
        rxposbtn2 = ttk.Radiobutton(rxframe, text='Manual Rx position:', variable=self.rxposvar, value=1,
                                    command=self.rxselection)
        rxposbtn2.grid(row=1, column=0, padx=5, pady=(5, 5), sticky='nsw')
        self.rxcombobox = ttk.Combobox(rxframe, values=radarslist, state='readonly', width=25, justify='center')
        self.rxcombobox.grid(row=0, column=1, padx=15, pady=5, sticky='nsw')
        self.rxcombobox.current(0)
        if not radarsloaded:
            self.rxposvar.set(1)

        rxposframe = tk.Frame(rxframe)
        rxposframe.grid(row=1, column=1, pady=(5, 5))
        rxlatlbl = tk.Label(rxposframe, text='Lat (φ):')
        rxlatlbl.grid(row=0, column=0, padx=(12, 2), pady=5, sticky='nsw')
        self.rxlatentry = ValidEntry(rxposframe, 'coords', width=9, state='disabled')
        self.rxlatentry.grid(row=0, column=1, padx=2, pady=5, sticky='nsw')
        rxlonlbl = tk.Label(rxposframe, text='Lon (λ):')
        rxlonlbl.grid(row=0, column=2, padx=(17, 2), pady=5, sticky='nsw')
        self.rxlonentry = ValidEntry(rxposframe, 'coords', width=9, state='disabled')
        self.rxlonentry.grid(row=0, column=3, padx=2, pady=5, sticky='nsw')
        rxheightlbl = tk.Label(rxposframe, text='Height (agl, m):')
        rxheightlbl.grid(row=0, column=4, padx=(20, 2), pady=5, sticky='nsw')
        self.rxheightentry = ValidEntry(rxposframe, 'float', width=5, state='disabled')
        self.rxheightentry.grid(row=0, column=5, padx=2, pady=5, sticky='nsw')
        rxpowerlbl = tk.Label(rxposframe, text='Rx Gain (dBm):')
        rxpowerlbl.grid(row=0, column=6, padx=(20, 2), pady=5, sticky='nsw')
        self.rxgainentry = ValidEntry(rxposframe, 'float', width=6, state='disabled')
        self.rxgainentry.grid(row=0, column=7, padx=(2, 5), pady=5, sticky='nsw')

        # initialize radiobutton selections, just limit selections if preset lists for antenna or rx stations are loaded
        self.rxselection()
        self.txgainrbtn()

        # button frame
        self.calcbutton = ttk.Button(calcbtnframe, text='Run', command=self.calculateloss, width=10)
        self.calcbutton.grid(row=0, column=0)

        # results frame
        resultsframe.rowconfigure(0, weight=0)
        resultsframe.rowconfigure(1, weight=0)
        resultsframe.rowconfigure(2, weight=1)
        resultsframe.rowconfigure(3, weight=0)

        resultsframe.columnconfigure(0, weight=1)
        resultsframe.columnconfigure(1, weight=0)
        resultsframe.columnconfigure(2, weight=0)

        boxframe1 = ttk.LabelFrame(resultsframe, text='Transmission', width=200, height=140)
        boxframe1.grid(row=0, column=1, padx=5, pady=5, sticky='nw')
        boxframe1.grid_propagate(0)
        boxframe2 = ttk.LabelFrame(resultsframe, text='Path profile', width=200, height=150)
        boxframe2.grid(row=1, column=1, padx=5, pady=5, sticky='nw')
        boxframe2.grid_propagate(0)
        infoframe = tk.Frame(resultsframe)
        infoframe.grid(row=2, column=1, padx=5, pady=5, sticky='nsew')

        # plot frame
        plt.style.use('seaborn')
        fig, self.ax = plt.subplots()
        fig.patch.set_facecolor((0.97, 0.97, 0.97))
        self.ax.set_xlabel('d (km)')
        self.ax.set_ylabel('Elevation (m)')
        self.ax.set_title('Path elevation profile')
        fig.subplots_adjust(left=0.1, bottom=0.15, right=0.95, top=0.92)

        self.canvas = FigureCanvasTkAgg(fig, resultsframe)
        self.canvas.get_tk_widget().grid(row=0, column=0, rowspan=3, sticky='nsew')

        # transmission results
        text1 = tk.Label(boxframe1, text='Rx level (dBm):')
        text1.grid(row=0, column=0, sticky='nsw')
        text2 = tk.Label(boxframe1, text='Transmission loss (dB):')
        text2.grid(row=1, column=0, sticky='nsw')
        text3 = tk.Label(boxframe1, text='FSP & gaseous atten.(dB):')
        text3.grid(row=2, column=0, sticky='nsw')
        text4 = tk.Label(boxframe1, text='Diff. loss (Ldp) (dB):')
        text4.grid(row=3, column=0, sticky='nsw')
        text5 = tk.Label(boxframe1, text='Tx Gain (dBi):')
        text5.grid(row=4, column=0, sticky='nsw')

        self.tlossvar = tk.StringVar()
        self.rxlevelvar = tk.StringVar()
        self.fspgaslossvar = tk.StringVar()
        self.difflossvar = tk.StringVar()
        self.txgainvar = tk.StringVar()
        self.rxlevelvar.set(f'{np.NAN:9.2f}')
        self.tlossvar.set(f'{np.NAN:9.2f}')
        self.fspgaslossvar.set(f'{np.NAN:9.2f}')
        self.difflossvar.set(f'{np.NAN:9.2f}')
        self.txgainvar.set(f'{np.NAN:9.2f}')

        rxlevellbl = tk.Label(boxframe1, textvariable=self.rxlevelvar)
        rxlevellbl.grid(row=0, column=1, padx=2, sticky='e')
        tlosslbl = tk.Label(boxframe1, textvariable=self.tlossvar)
        tlosslbl.grid(row=1, column=1, padx=2, sticky='e')
        fspgaslbl = tk.Label(boxframe1, textvariable=self.fspgaslossvar)
        fspgaslbl.grid(row=2, column=1, padx=2, sticky='e')
        difflosslbl = tk.Label(boxframe1, textvariable=self.difflossvar)
        difflosslbl.grid(row=3, column=1, padx=2, sticky='e')
        txgainlbl = tk.Label(boxframe1, textvariable=self.txgainvar)
        txgainlbl.grid(row=4, column=1, padx=2, sticky='e')

        # path profile frame
        text4 = tk.Label(boxframe2, text='Distance (km): ')
        text4.grid(row=0, column=0, sticky='nsw')
        text5 = tk.Label(boxframe2, text='Bearing (deg): ')
        text5.grid(row=1, column=0, sticky='nsw')
        text6 = tk.Label(boxframe2, text='Path type:')
        text6.grid(row=2, column=0, sticky='nsw')
        text6 = tk.Label(boxframe2, text='Tx height (m, amsl):')
        text6.grid(row=3, column=0, sticky='nsw')
        text7 = tk.Label(boxframe2, text='Rx height (m, amsl):')
        text7.grid(row=4, column=0, sticky='nsw')
        text8 = tk.Label(boxframe2, text='theta_t (deg):')
        text8.grid(row=5, column=0, sticky='nsw')

        self.distvar = tk.StringVar()
        self.bearingvar = tk.StringVar()
        self.pathtypevar = tk.StringVar()
        self.htsvar = tk.StringVar()
        self.hrsvar = tk.StringVar()
        self.thetatvar = tk.StringVar()
        self.distvar.set(f'{np.NAN:20.2f}')
        self.bearingvar.set(f'{np.NAN:20.2f}')
        self.pathtypevar.set('NLOS')
        self.htsvar.set(f'{np.NAN:20.2f}')
        self.hrsvar.set(f'{np.NAN:20.2f}')
        self.thetatvar.set(f'{np.NAN:20.2f}')

        distlbl = tk.Label(boxframe2, textvariable=self.distvar)
        distlbl.grid(row=0, column=1, padx=2, sticky='e')
        bearinglbl = tk.Label(boxframe2, textvariable=self.bearingvar)
        bearinglbl.grid(row=1, column=1, padx=2, sticky='e')
        pathtypelbl = tk.Label(boxframe2, textvariable=self.pathtypevar)
        pathtypelbl.grid(row=2, column=1, padx=2, sticky='e')
        htslbl = tk.Label(boxframe2, textvariable=self.htsvar)
        htslbl.grid(row=3, column=1, padx=2, sticky='e')
        hrslbl = tk.Label(boxframe2, textvariable=self.hrsvar)
        hrslbl.grid(row=4, column=1, padx=2, sticky='e')
        thetatlbl = tk.Label(boxframe2, textvariable=self.thetatvar)
        thetatlbl.grid(row=5, column=1, padx=2, sticky='e')

        # clickable labels for plot info and detailed p452 output
        plotinfo = tk.Label(infoframe, text='Plot info', cursor='hand2', fg='blue')
        plotinfo.grid(row=0, column=1, padx=5, pady=2)
        plotinfo.bind("<Button-1>", self.showplotinfo)
        detailedoutput = tk.Label(infoframe, text='p452 output', cursor='hand2', fg='blue')
        detailedoutput.grid(row=0, column=2, padx=5, pady=2)
        detailedoutput.bind("<Button-1>", self.showdetailedoutput)
        f = font.Font(plotinfo, plotinfo.cget("font"))
        f.configure(underline=True)
        plotinfo.configure(font=f)
        detailedoutput.configure(font=f)

        self.line1 = None
        self.line2 = None
        self.line3 = None
        self.fill1 = None
        self.fill2 = None
        self.results = None

    def pathrbtn(self):
        if self.pathselectvar.get():
            # manual path
            self.pathloadbtn.configure(state='normal')
            self.txgainselectvar.set(1)
            self.txgainrbtn1.configure(state='disabled')
            self.txgainrbtn2.configure(state='normal')
            self.txgainrbtn()
            # self.txlonentry.configure(state='disabled')
            # self.rxlonentry.configure(state='disabled')
        else:
            # auto path selected:
            self.pathloadbtn.configure(state='disabled')
            self.txgainselectvar.set(0)
            self.txgainrbtn1.configure(state='normal')
            self.txgainrbtn()
            # self.txlonentry.configure(state='normal')
            if self.rxposvar.get():
                self.rxlonentry.configure(state='normal')

    def txgainrbtn(self):
        if self.txgainselectvar.get():
            self.txgainentry.configure(state='normal')
            self.antennabox.configure(state='disabled')
            self.tiltsbox.configure(state='disabled')
            self.azimentry.configure(state='disabled')
            self.mechtiltentry.configure(state='disabled')
        else:
            self.txgainentry.configure(state='disabled')
            self.antennabox.configure(state='readonly')
            self.tiltsbox.configure(state='readonly')
            self.azimentry.configure(state='normal')
            self.mechtiltentry.configure(state='normal')

    def rxselection(self):
        if self.rxposvar.get():
            self.rxcombobox.configure(state='disabled')
            self.rxlatentry.configure(state='normal')
            self.rxlonentry.configure(state='normal')
            self.rxheightentry.configure(state='normal')
            self.rxgainentry.configure(state='normal')
        else:
            self.rxcombobox.configure(state='readonly')
            self.rxlatentry.configure(state='disabled')
            self.rxlonentry.configure(state='disabled')
            self.rxheightentry.configure(state='disabled')
            self.rxgainentry.configure(state='disabled')

    def loadpathfromfile(self):

        filepath = tk.filedialog.askopenfilename(initialdir=Path.cwd(), title="Select file",
                                                 filetypes=[("Excel file", ".xlsx .xls")])
        pathdf = pd.read_excel(filepath, names=['d', 'h', 'zone'], header=None)

        # fast check if path profile file contains header or not:
        try:
            float(pathdf.iloc[0, 0])
        except ValueError:
            pathdf.drop(0, inplace=True)

        try:
            data = pathdf.to_numpy()
            d = data[:, 0].astype('float64')
            h = data[:, 1].astype('float64')
            zone = data[:, 2].astype('str')
        except ValueError:
            showerror('Invalid path profile', 'Invalid path profile. File must include 3 columns: '
                                              '[d: number, h: number, zone: [A1/A2/B]')
            return
        if data.shape[0] < 4 or data.shape[1] != 3:
            showerror('Invalid path profile', 'Invalid path profile. It should contain more than 3 rows-data '
                                              'points and exactly 3 columns: [d, h, zonetype]')
            return
        if not set(zone).issubset({'A1', 'A2', 'B'}):
            showerror('Error in climate zones', 'Error in climate zones. Zone column must contain '
                                                'only "A1" (Coastal land), "A2 (inland)" or "B" (sea) values.')
            return

        self.di = d
        self.hi = h
        self.zonei = zone
        self.validpathloaded = True

    def getpathsrtm(self, phi_t, psi_t, phi_r, psi_r):
        dist, atrdeg, di, hi, delta = srtm.get_path_profile(phi_t, psi_t, phi_r, psi_r,
                                                            coord_samples=self.coordsamples,
                                                            fill_missing=self.fillvoidsrtm)

        zonei = np.asarray(['A2'] * self.coordsamples)
        zonei[hi <= self.coastlandthreshvalue] = 'A1'
        zonei[hi <= self.seathreshvalue] = 'B'

        return dist, atrdeg, di, hi, zonei, delta

    def calculateloss(self):
        *response, = self.validateparam()
        if response[0]:
            (f, p, pressure, temp, DN, N0, dct, dcr, phi_t, psi_t, txheight, txpower, azim, mechtilt,
             fixedtxgain, phi_r, psi_r, rxheight, rxgain) = response[2]

            pol = self.poldict[self.polarbox.get()]
            ha_t = self.clutterdict[self.txclutterbox.get()][0]
            dk_t = self.clutterdict[self.txclutterbox.get()][1]
            ha_r = self.clutterdict[self.rxclutterbox.get()][0]
            dk_r = self.clutterdict[self.rxclutterbox.get()][1]

            if self.pathselectvar.get():
                if not self.validpathloaded:
                    self.loadpathfromfile()  # sets di, hi, zonei arrays

                dist = self.di[-1] - self.di[0]
                atrdeg = srtm.get_path_geometry(phi_t, psi_t, phi_r, psi_r, 0)[1]
            else:
                dist, atrdeg, di, hi, zonei, delta = self.getpathsrtm(phi_t, psi_t, phi_r, psi_r)
                self.di = di
                self.hi = hi
                self.zonei = zonei

            if self.txgainselectvar.get():
                self.results = p452_loss(f, p, self.di, self.hi, self.zonei, txheight, rxheight, phi_t, phi_r, rxgain,
                                         pol, dct, dcr, DN, N0, pressure, temp, Gt=fixedtxgain,
                                         ha_t=ha_t, dk_t=dk_t, ha_r=ha_r, dk_r=dk_r)

                (Lb, Lbfsg, Lb0p, Lb0b, Ld50, Ldp, Lbs, Lba, theta_t, Gt, pathtype) = self.results
            else:
                # calculate gain by antenna pattern inside p452 function:
                eltilt = float(self.tiltsbox.get())
                antenna = self.antennabox.get()
                self.results = p452_loss(f, p, self.di, self.hi, self.zonei, txheight, rxheight, phi_t, phi_r, rxgain,
                                         pol, dct, dcr, DN, N0, pressure, temp, psi_t=psi_t, psi_r=psi_r,
                                         antennasdb=self.antennas, antennaname=antenna, azim=azim, eltilt=eltilt,
                                         mechtilt=mechtilt, ha_t=ha_t, dk_t=dk_t, ha_r=ha_r, dk_r=dk_r)
                (Lb, Lbfsg, Lb0p, Lb0b, Ld50, Ldp, Lbs, Lba, theta_t, Gt, pathtype) = self.results

            rxlevel = Gt + rxgain + txpower - Lb

            # convert theta_t mrad to degrees
            theta_t = np.rad2deg(theta_t / 1000)

            # transmitter and receiver heights amsl:
            hts = txheight + self.hi[0]
            hrs = rxheight + self.hi[-1]

            # plot path profile and line of sight
            if self.line1:
                self.line1.remove()
                self.line2.remove()
                self.fill1.remove()
                self.fill2.remove()
            if self.line3:
                self.line3.remove()

            ec = np.full_like(self.di, 0)
            rc = ec
            # median effective Earth radius factor k50 for the path
            k50 = 157 / (157 - DN)
            Ro = 6371  # km
            ae = k50 * Ro

            # if self.pathselectvar.get() == 0 and self.earthcurv == 'curved':
            if self.earthcurv == 'curved':
                # earth curvature
                # setting reference point x0 as the center of path (x0 = dist/2)
                ec = -1000 * (self.di - dist / 2) ** 2 / (2 * Ro)
                ec -= min(ec)

                # line of sight curvature due to refraction
                rc = ec * (1 - 1 / k50)
                rc -= min(rc)

            # stragiht line los without accounting refraction
            los = self.di * (hrs - hts) / dist + hts

            # modify path heights and line of sight curvature:
            self.h = self.hi + ec
            los += rc

            # plot elevation along path:
            self.line1, = self.ax.plot(self.di, self.h, color='grey')

            if pathtype == 'los':
                pathtype = 'Line of Sight'
                # if los, blue color to line connecting tx and rx
                self.line2, = self.ax.plot(self.di, los, color='blue')
                self.line3 = None
            else:
                # no los:
                self.line2, = self.ax.plot(self.di, los, color='red')
                # find obstruction point (point with maximum theta_tx):
                dii = self.di[1:-1]
                hii = self.h[1:-1] - ec[1:-1]
                thetatx = np.arctan((hii - hts) / (1000 * dii) - dii / (2 * ae))
                idx = np.argmax(thetatx)
                dmax = self.di[idx]
                self.line3, = self.ax.plot((0, dmax), (hts, self.h[idx + 1]), color='cyan', linewidth=1)

            self.fill1 = self.ax.fill_between(self.di, self.h, ec, color='darkgoldenrod', alpha=0.6)
            self.fill2 = self.ax.fill_between(self.di, ec, 0, color='silver')

            self.ax.set_xlim(0, dist)
            self.ax.set_ylim(0, max(2.5 * max(self.h), 1.7 * hts, 1.7 * hrs))
            self.canvas.draw()

            # update transmission results
            self.rxlevelvar.set(f'{rxlevel:9.2f}')
            self.tlossvar.set(f'{Lb:9.2f}')
            self.fspgaslossvar.set(f'{Lb0p:9.2f}')
            self.difflossvar.set(f'{Ldp:9.2f}')
            self.txgainvar.set(f'{Gt:9.2f}')

            # update path parameters results
            self.distvar.set(f'{dist:20.2f}')
            self.bearingvar.set(f'{atrdeg:20.2f}')
            self.htsvar.set(f'{hts:20.2f}')
            self.hrsvar.set(f'{hrs:20.2f}')
            self.thetatvar.set(f'{theta_t:20.2f}')
            self.pathtypevar.set(f'{pathtype}')

        else:
            showerror('Error', response[1])

    def validateparam(self):
        freq = float(self.freqentry.get())
        if freq < 0.1:
            return False, 'Frequency should be greater than 0.1GHz'

        timeper = float(self.timeentry.get())
        if timeper > 50:
            return False, 'Time percentage should be lower than 50%'

        pressure = float(self.pressureentry.get())
        if pressure < 0:
            return False, 'Pressure can\'t be negative'

        temp = float(self.tempentry.get())

        DN = float(self.dnentry.get())
        N0 = float(self.n0entry.get())
        if N0 <= DN:
            return False, 'N0 should be greater than ΔΝ value'

        dct = float(self.dctentry.get())
        dcr = float(self.dcrentry.get())

        if len(self.txlatentry.get()) == 0:
            return False, 'Please enter Tx latitude (decimal or DMS)'
        else:
            phi_t = dmstodec(self.txlatentry.get())

        if len(self.txlonentry.get()) == 0:
            return False, 'You need to enter Tx longitude value (decimal or DMS)'
        else:
            psi_t = dmstodec(self.txlonentry.get())

        if len(self.txheightentry.get()) == 0 or len(self.txpowerentry.get()) == 0:
            return False, 'You need to enter Tx antenna height (above ground) and power.'
        else:
            txheight = float(self.txheightentry.get())
            if txheight <= 0:
                return False, "Transmitter antenna height must be can't be 0"

            txpower = float(self.txpowerentry.get())

        if not self.txgainselectvar.get():
            if len(self.azimentry.get()) == 0 or len(self.mechtiltentry.get()) == 0:
                return False, 'You need to enter Tx antenna azimuth and mech. tilt.'
            else:
                azim = float(self.azimentry.get())
                mechtilt = float(self.mechtiltentry.get())
                fixedtxgain = None
        else:
            if len(self.txgainentry.get()) == 0:
                return False, 'You need to enter a number for Tx antenna Gain (dBi).'
            fixedtxgain = float(self.txgainentry.get())
            azim = None
            mechtilt = None

        if self.rxposvar.get():
            if len(self.rxlatentry.get()) == 0:
                return False, 'You need to enter a valid Rx station latitude (decimal or DMS).'
            phi_r = dmstodec(self.rxlatentry.get())

            if len(self.rxlonentry.get()) == 0:
                return False, 'You need to enter a valid Rx station longitude (decimal or DMS).'
            psi_r = dmstodec(self.rxlonentry.get())

            if len(self.rxheightentry.get()) == 0 or len(self.rxgainentry.get()) == 0:
                return False, 'You need to enter Rx station height (above ground) and gain (dBi).'
            rxheight = float(self.rxheightentry.get())
            if rxheight <= 0:
                return False, "Rx antenna height above ground must be greater than 0."
            rxgain = float(self.rxgainentry.get())
        else:
            # get receiver station from preset list
            try:
                rxidx = self.rxcombobox.current()
                phi_r = self.radars.loc[rxidx, 'lat']
                psi_r = self.radars.loc[rxidx, 'lon']
                rxheight = self.radars.loc[rxidx, 'h_agl']
                rxgain = self.radars.loc[rxidx, 'gain']
            except KeyError:
                return False, 'Error when trying to load Rx station attributes'

        # check if transmitter and receiver have the same coords:
        if (phi_t, psi_t) == (phi_r, psi_r):
            return False, 'Tx position is identical to Rx position.'

        return True, None, (freq, timeper, pressure, temp, DN, N0, dct, dcr, phi_t, psi_t,
                            txheight, txpower, azim, mechtilt, fixedtxgain, phi_r, psi_r, rxheight, rxgain)

    def updatetilts(self, evnt):
        if not self.antennasloaded:
            return
        antennasel = self.antennabox.get()
        tiltslist = self.antennas.loc[antennasel].index.tolist()
        self.tiltsbox.configure(values=tiltslist)
        self.tiltsbox.current(0)

    def opensettings(self):
        self.settingswindow = tk.Toplevel()
        self.settingswindow.title('Settings')
        x = self.master.winfo_x()
        y = self.master.winfo_y()
        self.settingswindow.geometry("+%d+%d" % (x + 100, y + 70))

        # self.settingswindow.iconbitmap('icon.ico')

        # read saved variables
        self.zoneselvar.set(self.zonesel)
        self.earthcurvvar.set(self.earthcurv)
        self.fillvoidsrtmvar.set(self.fillvoidsrtm)

        self.settingswindow.columnconfigure(0, weight=0)
        self.settingswindow.columnconfigure(1, weight=1)

        zoneframe = ttk.LabelFrame(self.settingswindow, text='Radio-climate zones')
        text = tk.Label(self.settingswindow, text='Radio-climate zone calculation options, '
                                                  'when the path is not loaded from file:')
        text.grid(row=0, column=0, columnspan=2, padx=(10, 70), pady=(10, 2), sticky='nsew')
        zoneframe.grid(row=1, column=0, columnspan=2, padx=10, pady=(2, 10), sticky='nsw')
        zonerb1 = ttk.Radiobutton(zoneframe, text='Calculate zone based on heights:',
                                  variable=self.zoneselvar, value='auto')
        zonerb2 = ttk.Radiobutton(zoneframe, text='Fixed A2 (inland) zone', variable=self.zoneselvar, value='fixed')
        seathreshlbl = tk.Label(zoneframe, text="Sea ('B') zone when h<=:")
        coastlandthreshlbl = tk.Label(zoneframe, text="Coastland ('A2') zone when h<=:")
        self.seathreshentry = ValidEntry(zoneframe, 'float', width=5)
        self.coastlandthreshentry = ValidEntry(zoneframe, 'float', width=5)

        self.seathreshentry.delete(0, tk.END)
        self.coastlandthreshentry.delete(0, tk.END)
        self.seathreshentry.insert(0, self.seathreshvalue)
        self.coastlandthreshentry.insert(0, self.coastlandthreshvalue)

        zonerb1.grid(row=0, column=0, padx=3, pady=1, columnspan=4, sticky='nsw')
        seathreshlbl.grid(row=1, column=2, padx=(10, 5), pady=1, sticky='nsw')
        self.seathreshentry.grid(row=1, column=3, padx=(2, 10), pady=1, sticky='nsw')
        coastlandthreshlbl.grid(row=1, column=0, padx=(20, 2), pady=1, sticky='nsw')
        self.coastlandthreshentry.grid(row=1, column=1, padx=2, pady=1, sticky='nsw')
        zonerb2.grid(row=2, column=0, padx=3, pady=5, columnspan=4, sticky='nsw')

        earthcurvbtn = ttk.Checkbutton(self.settingswindow, text='Consider earth curvature and refraction when '
                                                                 'plotting path',
                                       variable=self.earthcurvvar, onvalue='curved', offvalue='flat')
        earthcurvbtn.grid(row=2, column=0, columnspan=2, pady=6, padx=15, sticky='nsw')

        fillvoidsrtm = ttk.Checkbutton(self.settingswindow, text='Fill missing SRTM data with 0 heights',
                                       variable=self.fillvoidsrtmvar)
        fillvoidsrtm.grid(row=3, column=0, columnspan=2, pady=6, padx=15, sticky='nsw')

        coordsampleslbl = tk.Label(self.settingswindow, text='Max number of SRTM samples in path profile generator:')
        coordsampleslbl.grid(row=4, column=0, padx=(15, 2), pady=6, sticky='nsw')
        self.coordsamplesentry = ValidEntry(self.settingswindow, 'float', width=6)
        self.coordsamplesentry.grid(row=4, column=1, padx=2, pady=6, sticky='nsw')
        self.coordsamplesentry.delete(0, tk.END)
        self.coordsamplesentry.insert(0, self.coordsamples)

        submitsettingsbtn = ttk.Button(self.settingswindow, text='Save', command=self.submitsettings)
        submitsettingsbtn.grid(row=5, column=0, columnspan=2, pady=20)

    def submitsettings(self):
        sea = float(self.seathreshentry.get())
        coast = float(self.coastlandthreshentry.get())
        samples = int(self.coordsamplesentry.get())

        if samples == 0:
            self.coordsamples = 700
        elif samples > 1500:
            self.coordsamples = 1500
        else:
            self.coordsamples = samples

        if self.zoneselvar.get() == 'auto':
            if sea >= coast:
                showerror('Error', 'Sea threshold value must be smaller than coast-land threshold value.')
            else:
                self.zonesel = 'auto'
                self.seathreshvalue = sea
                self.coastlandthreshvalue = coast
                self.earthcurv = self.earthcurvvar.get()
                self.fillvoidsrtm = self.fillvoidsrtmvar.get()
                self.settingswindow.destroy()
        else:
            self.zonesel = 'fixed'
            self.earthcurv = self.earthcurvvar.get()
            self.fillvoidsrtm = self.fillvoidsrtmvar.get()
            self.settingswindow.destroy()

    def showdetailedoutput(self, e):
        if self.results:
            resultswd = tk.Toplevel()
            resultswd.transient(self.master)
            resultswd.grab_set()
            resultswd.title('p452 output parameters')
            x = self.master.winfo_x()
            y = self.master.winfo_y()
            resultswd.geometry(f"+{x + 500:d}+{y + 200:d}")
            paramslbl = ['Lb', 'Lbfsg', 'Lb0p', 'Lb0b', 'Ld50', 'Ldp', 'Lbs', 'Lba']
            paramsvl = self.results[:8]
            text = ''
            for i, lbl in enumerate(paramslbl):
                text += f'{paramslbl[i]:5s} = {paramsvl[i]:8.4f} dB\n'
            textwidget = tk.Text(resultswd, padx=10, pady=10, width=35, height=8)
            textwidget.grid(row=0, column=0, sticky='nsew')
            textwidget.insert(tk.END, text)
            textwidget.config(state='disabled')
            closebtn = ttk.Button(resultswd, text='Exit', command=lambda: resultswd.destroy())
            closebtn.grid(row=1, column=0, pady=6, padx=5)
        else:
            showerror('No results', 'First run a prediction to get results')

    def showplotinfo(self, e):
        showinfo('Plot explained', 'Path elavation profile is plotted, traveling along a great circle arc, from '
                                   'transmitter to receiver station. WGS84 reference ellipsoid is used.\n'
                                   'When enabled in settings, earth curvature (Ro=6371km) and refraction is taken into '
                                   'consideration, so what you actually see as elevation profile is the distance above '
                                   'the chord running tx to rx station.\nRefraction is also considered, thus the '
                                   'line-of-sight isn\'t actually a straight line, but a curved line with a '
                                   'curvature depending on ΔΝ value (average radio-refractive index lapse-rate through '
                                   'the lowest 1 km of the atmosphere). Refraction and curvature of Earth have opposite'
                                   ' impact of tx-rx line of sight determination, since the curvature of Earth makes '
                                   'distant objects look lower, and refraction makes them look higher.')


def quit_me():
    root.quit()
    root.destroy()

if __name__ == '__main__':
    root = tk.Tk()
    root.protocol("WM_DELETE_WINDOW", quit_me)
    s = ttk.Style()
    s.theme_use('vista')
    mainapp = MainGUI(root)
    root.mainloop()
