#!/usr/bin/env python
import tkinter
import tkinter.ttk
from sys import version_info
if version_info[0] == 2:
    import Tkinter as Tk
else:
    import tkinter as Tk

from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import photoncount
import pandas as pd
import os
import csv

######################################################################
### CONSTANTS
# units conversion
NM_TO_M = 1e-9
M_TO_NM = 1e9
NM_TO_PM = 1e3
KM_TO_M = 1e3
# plot window sizes
WINDOW_WIDTH = 9
WINDOW_HEIGHT = 5

# Init value parameters
D_INIT = 4.2 # Telescope diameter in m
R_INIT = 8e4 # Resolving power

T_T = 0.3 # Telescope throughput
T_I = 0.4 # Instrument throughput
QE = 0.9 # Detector QE
CLD = 1.0 # CLD Focal station BS transmission or reflection

SN_INIT = 1e3 #SNR
V_INIT = 7.0 #Velocity

# entry field width
FWIDTH = 8

class photongui():
    #
    # Computes and displays the integration to reach a desired level
    # of S/N per diffraction limited spatial pixel and given spectral
    # resolution.
    #
    # written by Jorrit Leenaarts
    def __init__(self, master):
        self.ph = photoncount.Photocount()

        self.master = master
        self.master.title("Generic solar telescope integration time estimator              J. Leenaarts - Institute for Solar Physics")

        self.title_frame = None
        self.wavelength_frame = None

        ### Widgets & Widgets parameters###
        self.init_parameters_widgets()
        ## Widgets
        self.init_widgets()

        self.init_parameters()

        #init rowcounter
        self.rowcounter = 0

        self.init_gui()

        #init dictionary:
        self.properties_dict = dict
        self.update_properties_dict()

        self.init_plot()
        self.redraw()

    def title_widget(self):
        self.title_frame = tkinter.Frame(self.master, bd=1)
        input_label = Tk.Label(self.title_frame, text='Set model parameters here', font=("Helvetica", 20, "bold"))
        input_label.grid(row=0, column=0, columnspan=4, sticky='n')

    def throughput_widget(self):
        self.throughput_frame = tkinter.Frame(self.master, bd=1, padx=4, pady=4, relief="sunken")
        Tk.Label(self.throughput_frame, text='Total Throughput [0,1]:', font=("Helvetica", 16, "bold")).grid(
            row=0, column=0, sticky='w')
        T_widget = Tk.Label(self.throughput_frame, width=5, textvariable=self.T, font=("Helvetica", 15, "bold"))
        #T_widget.bind('<Return>', self.redraw_from_event)
        T_widget.grid(row=0, column=1, sticky='w')

        self.T_T_widget = Tk.Spinbox(self.throughput_frame, from_=0.2, to=1.0, increment=0.02, width=5, textvariable=self.T_T,
                                command=self.redraw, font=("Helvetica", 15))
        self.T_I_widget = Tk.Spinbox(self.throughput_frame, from_=0.2, to=1.0, increment=0.02,width=5, textvariable=self.T_I,
                                command=self.redraw, font=("Helvetica", 15))
        self.QE_widget = Tk.Spinbox(self.throughput_frame, from_=0.02, to=1.0, increment=0.02,
                               command=self.redraw, width=5, textvariable=self.QE, font=("Helvetica", 15))

        self.CLD_widget = Tk.Spinbox(self.throughput_frame, from_=0.02, to=1.0, increment=0.02,
                               command=self.redraw, width=5, textvariable=self.CLD, font=("Helvetica", 15))

        Tk.Label(self.throughput_frame, text='Telescope:', font=("Helvetica", 15)).grid(
            row=1, column=0, sticky='w')
        self.T_T_widget.grid(row=1, column=1)
        self.T_T_widget.bind('<Return>', self.redraw_from_event)

        Tk.Label(self.throughput_frame, text='Instrument:', font=("Helvetica", 15)).grid(
            row=1, column=2, sticky='w')
        self.T_I_widget.grid(row=1, column=3, sticky='w')
        self.T_I_widget.bind('<Return>', self.redraw_from_event)

        Tk.Label(self.throughput_frame, text='QE:', font=("Helvetica", 15)).grid(
            row=1, column=4, sticky='w')
        self.QE_widget.grid(row=1, column=5, sticky='w')
        self.QE_widget.bind('<Return>', self.redraw_from_event)

        Tk.Label(self.throughput_frame, text='CLD Focal station:', font=("Helvetica", 15)).grid(
            row=1, column=6, sticky='w')
        self.CLD_widget.grid(row=1, column=7, sticky='w')
        self.CLD_widget.bind('<Return>', self.redraw_from_event)

    def performance_widget(self):
        self.performance_frame = tkinter.Frame(self.master, bd=1, padx=4, pady=4, relief="sunken")
        Tk.Label(self.performance_frame, text='Required Performances', font=("Helvetica", 16, "bold")).grid(
            row=0, column=0, sticky='w')
        # S/N
        Tk.Label(self.performance_frame, text='Desired S/N:', font=("Helvetica", 16)).grid(
            row=1, column=0, sticky='w')
        self.SN_widget = Tk.Spinbox(self.performance_frame, from_=1, to=1000000, increment=1,
                             width=FWIDTH, textvariable=self.SN, command=self.redraw, font=("Helvetica", 16))
        self.SN_widget.bind('<Return>', self.redraw_from_event)
        self.SN_widget.grid(row=1, column=1, sticky='w')

        #Resolving Power
        Tk.Label(self.performance_frame, text='Resolving Power (R):', font=("Helvetica", 16)).grid(
            row=2, column=0, sticky='w')
        self.R_entry = Tk.Spinbox(self.performance_frame, from_=1, to=1000000, increment=1, width=FWIDTH,
                             textvariable=self.R, command=self.redraw, font=("Helvetica", 16))
        self.R_entry.grid(row=2, column=1, sticky='w')
        self.R_entry.bind('<Return>', self.redraw_from_event)

        #Speed evolution
        Tk.Label(self.performance_frame, text='Evolution speed (km/s):', font=("Helvetica", 16)).grid(
            row=3, column=0, sticky='w')
        self.v_entry = Tk.Spinbox(self.performance_frame,  from_=0, to=1000, increment=0.1,
                             width=FWIDTH, textvariable=self.v, command=self.redraw, font=("Helvetica", 16))
        self.v_entry.bind('<Return>', self.redraw_from_event)
        self.v_entry.grid(row=3, column=1, sticky='w')

        # polarimetry radio button
        Tk.Label(self.performance_frame, text='Polarimetry:', font=("Helvetica", 16)).grid(
            row=4, column=0, sticky='w')
        for i, option in zip((0, 1), ('no', 'yes')):
            r = Tk.Radiobutton(self.performance_frame, text=option, variable=self.polarimetry,
                            value=i, command=self.redraw, font=("Helvetica", 16))
            r.grid(row=4, column=i+1, sticky='w')

    def wavelength_widget(self):
        self.wavelength_frame = tkinter.Frame(self.master, bd=1, padx=4, pady=4, relief="sunken")
        values = [('Ca II K', 393.33), ('H I beta', 486.10), ('Mg I b', 517.30), ('Fe I @525', 525.0),
            ('Fe I @543', 543.0), ('Na I D @589', 589.60), ('Na I D @590', 590.00), ('Fe I @630', 630.2),
            ('H I alpha', 656.3), ('K I @769', 769.9), ('Fe I@846', 846.80), ('CaII @849', 849.80),
            ('Fe I @851', 851.40), ('Ca II @854', 854.2)]
        if (self.ph.file == 1):
            values.append(('He I', 1083.0))
            values.append(('Fe I IR', 1564.8))
        self.cwl = Tk.ttk.Combobox(self.wavelength_frame, values = values)
        self.bp = Tk.Spinbox(self.wavelength_frame, from_=0, to=10, increment=0.02, width=5, textvariable=self.bandpass,
                             command=self.redraw, font=("Helvetica", 15))
        self.cwls = Tk.Spinbox(self.wavelength_frame, from_=-1000, to=1000, increment=0.01, width=5,
                               textvariable=self.cwl_shift, command=self.redraw, font=("Helvetica", 15))
        Tk.Label(self.wavelength_frame, text='Target Wavelength Parameters:', font=("Helvetica", 16, "bold")).grid(
            row=0, column=0, sticky='w')
        Tk.Label(self.wavelength_frame, text='Target Wavelength:', font=("Helvetica", 15)).grid(
            row=1, column=0, sticky='w')
        self.cwl.bind('<<ComboboxSelected>>', self.redraw_from_event)
        self.cwl.grid(row=1, column=1, sticky='w')

        Tk.Label(self.wavelength_frame, text='Bandpass(nm):', font=("Helvetica", 15)).grid(
            row=2, column=0, sticky='w')
        self.bp.grid(row=2, column=1, sticky='w')
        self.bp.bind('<Return>', self.redraw_from_event)

        Tk.Label(self.wavelength_frame, text='CWL shift(nm):', font=("Helvetica", 15)).grid(
            row=2, column=2, sticky='w')
        self.cwls.grid(row=2, column=3, sticky='w')
        self.cwls.bind('<Return>', self.redraw_from_event)

    def pixel_widget(self):
        self.pixel_frame = tkinter.Frame(self.master, bd=1, padx=4, pady=4, relief="sunken")

        Tk.Label(self.pixel_frame, text='Pixel Parameters', font=("Helvetica", 16, "bold")).grid(
            row=0, column=0, sticky='w')
        """
        spatial binning radio buttons
        """
        Tk.Label(self.pixel_frame, text='Spatial binning', font=("Helvetica", 15)).grid(
            row=1, column=0, sticky='w')
        j = 1
        for i, option in zip((1., 2., 3., 4.), ('1x1', '2x2', '3x3', ' 4x4')):
            r = Tk.Radiobutton(self.pixel_frame, text=option, variable=self.binning,
                            value=i, command=self.redraw, font=("Helvetica", 15))
            r.grid(row=1, column=j, sticky='w')
            j += 1
        # resulting spatial resolution
        Tk.Label(self.pixel_frame, text='Spatial pixel sampling (arcsec)', font=("Helvetica", 15)).grid(
            row=2, column=0, sticky='w')
        self.resmin_entry = Tk.Entry(self.pixel_frame, width=FWIDTH, textvariable=self.resmin, font=("Helvetica", 15))
        self.resmin_entry.grid(row=2, column=1, sticky='w')
        Tk.Label(self.pixel_frame, text='-', font=("Helvetica", 15)).grid(row=2, column=2, sticky='w')
        self.resmax_entry = Tk.Entry(self.pixel_frame, width=FWIDTH, textvariable=self.resmax, font=("Helvetica", 15))
        self.resmax_entry.grid(row=2, column=3, sticky='w')
        self.pixel_frame.grid_columnconfigure(2, weight=1)
        self.pixel_frame.grid_rowconfigure(2, weight=1)
        self.pixel_frame.grid_columnconfigure(3, weight=1)
        #self.pixel_frame.grid_rowconfigure(2, weight=1)

        # resulting spectral resolution
        Tk.Label(self.pixel_frame, text='Spectral pixel sampling (pm)', font=("Helvetica", 15)).grid(
            row=3, column=0, sticky='w')
        self.spresmin_entry = Tk.Entry(self.pixel_frame, width=FWIDTH, textvariable=self.spresmin, font=("Helvetica", 15))
        self.spresmin_entry.grid(row=3, column=1, sticky='w')
        Tk.Label(self.pixel_frame, text='-', font=("Helvetica", 15)).grid(row=3, column=2, sticky='w')
        self.spresmax_entry = Tk.Entry(self.pixel_frame, width=FWIDTH, textvariable=self.spresmax, font=("Helvetica", 15))
        self.spresmax_entry.grid(row=3, column=3, sticky='w')

    def telescope_widget(self):
        self.telescope_frame = tkinter.Frame(self.master, bd=1, padx=4, pady=4, relief="sunken")
        self.telescope_frame.grid_columnconfigure(0, weight=1)
        Tk.Label(self.telescope_frame, text='Telescope Parameters:', font=("Helvetica", 16, 'bold')).grid(
            row=0, column=0, sticky='w')
        Tk.Label(self.telescope_frame, text='Aperture diameter (m):', font=("Helvetica", 15)).grid(
            row=1, column=0, sticky='w')
        self.D_widget = Tk.Entry(self.telescope_frame, width=FWIDTH, textvariable=self.D, font=("Helvetica", 15))
        self.D_widget.bind('<Return>', self.redraw_from_event)
        self.D_widget.grid(row=1, column=1, sticky='w')

        Tk.Label(self.telescope_frame, text='Strehl:', font=("Helvetica", 15)).grid(
            row=2, column=0, sticky='w')
        self.strehl_widget = Tk.Entry(self.telescope_frame, width=FWIDTH, textvariable=self.strehl, font=("Helvetica", 15))
        self.strehl_widget.grid(row=2, column=1, sticky='w')

    def init_widgets(self):
        self.title_widget()
        self.wavelength_widget()
        self.throughput_widget()
        self.performance_widget()
        self.pixel_widget()
        self.telescope_widget()

    def init_gui(self):
        #window label
        self.title_frame.grid(row=self.rowcounter, column=0, sticky="n")
        self.rowcounter += 1

        # telescope diameter
        self.telescope_frame.grid(row=self.rowcounter, column=0, sticky='w')
        self.rowcounter += 1

        # target wavelengths combobox
        self.wavelength_frame.grid(row=self.rowcounter, column=0, sticky='w')
        self.rowcounter += 1

        # Throughput
        self.throughput_frame.grid(row=self.rowcounter, column=0, sticky='w')
        self.rowcounter += 1
        
        # Performances
        self.performance_frame.grid(row=self.rowcounter, column=0, sticky='w')
        self.rowcounter += 1

        self.pixel_frame.grid(row=self.rowcounter, column=0, sticky='w')
        self.rowcounter += 1

        # export button
        export_button = Tk.Button(self.master, text="Export graphs", command=self.exportGraphs, font=("Helvetica", 18))
        export_button.grid(row=self.rowcounter, column=0, sticky='ew')
        self.rowcounter += 1

        # quit button
        quit_button = Tk.Button(self.master, text="Quit", command=self.quit, font=("Helvetica", 18))
        quit_button.grid(row=self.rowcounter, column=0, sticky='ew')

    def init_plot(self):
        # create Figure for the plot 
        self.fig = Figure(figsize=(WINDOW_WIDTH, WINDOW_HEIGHT), dpi=100)
        self.axes = []

        nrow = 2
        ncol = 2
        for i in range(nrow*ncol):
            self.axes.append(self.fig.add_subplot(nrow, ncol, i+1))

        # link figure to gui
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.master)
        self.canvas.get_tk_widget().grid(row=0, column=4, rowspan=self.rowcounter)
        self.canvas.draw()

    def init_parameters_widgets(self):
        # Widgets Parameters
        self.cwl_shift = Tk.DoubleVar(value=0.0)
        self.bandpass = Tk.DoubleVar(value=2.0)
        self.lmin = Tk.DoubleVar()
        self.lmax = Tk.DoubleVar()
        self.D = Tk.DoubleVar(value=D_INIT)  # T elescope diameter in [m]
        self.strehl = Tk.DoubleVar(value=1) # Telescope Strehl ratio
        self.R = Tk.IntVar(value=R_INIT)  # Desired resolving power
        self.polarimetry = Tk.IntVar(value=0)  # Polarimetric mode
        self.T = Tk.DoubleVar()  # Overall transmission
        self.T_T = Tk.DoubleVar(value=T_T) # Telescope transmission
        self.T_I = Tk.DoubleVar(value=T_I) # Instrument transmission (without QE)
        self.QE = Tk.DoubleVar(value=QE) # qe of the detector
        self.CLD = Tk.DoubleVar(value=CLD) # transmission of the Focal Station BS
        self.SN = Tk.IntVar(value=SN_INIT)  # Desired SNR
        self.v = Tk.DoubleVar(value=V_INIT)  # Velocity of the structure
        self.binning = Tk.DoubleVar(value=1.0)
        self.resmin = Tk.DoubleVar()
        self.resmax = Tk.DoubleVar()
        self.spresmin = Tk.DoubleVar()
        self.spresmax = Tk.DoubleVar()

    def init_parameters(self):
        self.cwl.current(13)

        cwl = self.cwl.get().split('}', 1)
        cwl = float(cwl[1])

        self.lmin.set(cwl-(self.bandpass.get())/2)
        self.lmax.set(cwl+(self.bandpass.get())/2)

        # throughput
        self.T.set(self.T_T.get()*self.T_I.get()*self.QE.get())

        self.set_spatres(self.lmin.get(), self.lmax.get())
        self.set_specres(self.lmin.get(), self.lmax.get())

        #lmin lmax:
        self.lmin, self.lmax = cwl-(self.bandpass.get())/2.0, cwl+(self.bandpass.get())/2.0

    def update_properties_dict(self):
        self.properties_dict = {
                'D' : self.D.get(), # telescope diameter
                'lmin' : self.lmin,
                'lmax' : self.lmax,
                'polarimetry': self.polarimetry.get(),
                'R' : self.R.get(), # spectral resolution
                'T' : self.T.get(), # telescope transmission
                'SN': self.SN.get(), # signal to noise ratio
                'v' : self.v.get() * KM_TO_M, # m/s
                'binning' : self.binning.get(), # spatial binning
                'strehl': self.strehl.get()
        }

    def quit(self):
        self.master.destroy()

    def set_spatres(self, lmin, lmax):
        self.resmin.set("{:.4f}".format(self.ph.spatres(lmin*NM_TO_M, self.D.get())))
        self.resmax.set("{:.4f}".format(self.ph.spatres(lmax*NM_TO_M, self.D.get())))

    def set_specres(self, lmin, lmax):
        self.spresmin.set("{:.4f}".format(self.ph.specres(lmin, self.R.get())))
        self.spresmax.set("{:.4f}".format(self.ph.specres(lmax, self.R.get())))

    def redraw_from_event(self, event):
        self.redraw()

    def redraw(self):
        #gets and converts the target wavelength value
        cwl = self.cwl.get().split('}', 1)
        cwl = float(cwl[1])

        cwl = cwl + float(self.cwls.get())
        self.T.set(self.T_T.get()*self.T_I.get()*self.QE.get()*self.CLD.get())

        #aux variables
        self.lmin, self.lmax = cwl-(self.bandpass.get())/2.0, cwl+(self.bandpass.get())/2.0
        # update of property dict
        self.update_properties_dict()

        self.ph.set_properties(self.properties_dict)
        self.ph.compute()
        self.set_spatres(self.lmin, self.lmax)
        self.set_specres(self.lmin, self.lmax)

        # plot
        self.plot()

    def exportGraphs(self):
        xax = self.ph.ll * M_TO_NM
        df = pd.DataFrame()
        #FIXME: I think you can simplify this with the tolist() method
        # No, as I modify the values with round() and int(), I need to select values one by one
        wavelengths = []
        nflux = []
        dx = []
        intTime = []
        optimalIntTime = []
        # Creation of df with data
        for i in range(len(xax)):
            wavelengths.append(round(xax[i], 5))
            nflux.append(int(self.ph.Ilambda[i]))
            dx.append(round(self.ph.dx[i], 2))
            intTime.append(round(self.ph.t[i], 2))
            optimalIntTime.append(round(self.ph.dt[i], 2))
        df['wavelength'] = wavelengths
        df['nflux'] = nflux
        df['dx'] = dx
        df['int_time'] = intTime
        df['optimal_int_time'] = optimalIntTime
        # Creation of csv file
        nameOfFile = self.openFile()
        if nameOfFile is None or nameOfFile == '':
            return
        with open(str(nameOfFile), 'w') as nameOfFile:
            writer = csv.writer(nameOfFile, delimiter=';')
        # Create the header of the csv file
            writer.writerow(['Telescope Diameter: ' + self.D_widget.get()])
            writer.writerow(['Telescope Sthrel: ' + self.strehl_widget.get()])
            writer.writerow(['TARGET WAVELENGTH: ' + self.cwl.get()])
            writer.writerow(['CWL shift: ' + self.cwls.get()])
            writer.writerow(['BP(nm): ' + self.bp.get()])
            writer.writerow(['TELESCOPE: ' + self.T_T_widget.get()])
            writer.writerow(['INSTRUMENT: ' + self.T_I_widget.get()])
            writer.writerow(['QE: ' + self.QE_widget.get()])
            writer.writerow(['CLD FOCAL STATION: ' + self.CLD_widget.get()])
            writer.writerow(['Desired S/N: ' + self.SN_widget.get()])
            writer.writerow(['Resolving Power (R): ' + self.R_entry.get()])
            writer.writerow(['Evolution Speed(km/h): ' + self.v_entry.get()])
            if (self.polarimetry.get() == 1):
                writer.writerow(['POLARIMETRY: ' + 'Yes'])
            else:
                writer.writerow(['POLARIMETRY: ' + 'No'])
            writer.writerow(['Spatial Binning: ' + str(self.binning.get())])
            writer.writerow(['Minimal spatial pixel sampling (arcsec): ' + self.resmin_entry.get()])
            writer.writerow(['Maximal spatial pixel sampling (arcsec): ' + self.resmax_entry.get()])
            writer.writerow(['Minimal spectral pixel sampling (pm): ' + self.spresmin_entry.get()])
            writer.writerow(['Maximal spectral pixel sampling (pm): ' + self.spresmax_entry.get()])
            writer.writerow(['wavelength(nm)', 'nflux(W/m^2/m/sr)', 'dx(arcsec)', 'int_time(s)', 'optimal_int_time(s)'])
        # Copy data from df and output to csv
            for i in range(len(df)):
                writer.writerow([
                    df['wavelength'][i],
                    df['nflux'][i],
                    df['dx'][i],
                    df['int_time'][i],
                    df['optimal_int_time'][i]])

    def openFile(self):
        route = './csv_data/'
        if not os.path.exists(route):
            os.makedirs(route)
        nameOfFile = route + str(self.cwl.get().split(' ')[len(self.cwl.get().split(' ')) - 1]) + '_graphs.csv'
        if (os.path.exists(nameOfFile)):
            answer = tkinter.messagebox.askyesno(
                message="The file will be deleted! Are you sure you want to overwrite it?",
                title="overwrite?")
            if answer is False:
                print('File not overwritten')
                return
            print('File already exists, overwriting...')
        else:
            print('Creating file ' + nameOfFile + '...')
        return nameOfFile


    def plot(self):
        xax = self.ph.ll * M_TO_NM

        for ax in self.axes:
            ax.clear()

        ax = self.axes[0]
        ax.plot(xax, self.ph.Ilambda, label='atlas')
        ax.plot(xax, self.ph.Ilambdac, label='smeared to R')

        ax.legend(loc='best')
        ax.get_xaxis().set_ticklabels([])
        ax.set_xlim(self.ph.ran )
        ax.set_ylabel(r'$I_\lambda$ [W m$^{-2}$ m$^{-1}$ sr$^{-1}$]')
        ax.set_title('spectrum')
        ax.grid(True)

        ax = self.axes[1]
        ax.plot(xax, self.ph.dx)
        ax.set_xlim(self.ph.ran)
        ax.get_xaxis().set_ticklabels([])
        ax.set_ylabel(r'$\Delta x$ [arcsec]')
        ax.set_title('optimal pixel size for given signal speed')
        ax.grid(True)

        ax = self.axes[2]
        ax.plot(xax, self.ph.t, label='with given transmission')
        ax.plot(xax, self.ph.tideal, label='perfect telescope')
        ax.legend(loc='best')
        ax.set_xlabel(r'$\lambda$ (nm)')
        ax.set_xlim(self.ph.ran)
        ax.set_ylabel(r' $\Delta t$ [s]')
        ax.set_title('integration time for given spatial and spectral pixel size')
        ax.grid(True)

        ax = self.axes[3]
        ax.plot(xax, self.ph.dt)
        ax.set_xlim(self.ph.ran)
        ax.set_ylabel(r'$\Delta t$ [s]')
        ax.set_title('optimal integration time for given signal speed')
        ax.set_xlabel(r'$\lambda$ (nm)')
        ax.grid(True)

        self.fig.tight_layout()

        self.canvas.draw()

def create_main_window():
    root = Tk.Tk()
    #root.grid_rowconfigure(0, weight=1)
    #root.grid_columnconfigure(0, weight=1)
    root.resizable(True, True)
    photongui(root)
    root.mainloop()

if __name__ == "__main__":
    create_main_window()