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
        self.init_plot()
        self.redraw()

    def title_widget(self):
        self.title_frame = tkinter.Frame(self.master, bd=1)
        input_label = Tk.Label(self.title_frame, text='Set model parameters here', font=("Helvetica", 20, "bold"))
        input_label.grid(row=0, column=0, columnspan=4, sticky='n')

    def throughput_widget(self):
        self.throughput_frame = tkinter.Frame(self.master, bd=1, padx=4, pady=4, relief="sunken")
        Tk.Label(self.throughput_frame, text='Total Throughput [0,1]', font=("Helvetica", 16, "bold")).grid(
            row=0, column=0, sticky='w')
        T_widget = Tk.Entry(self.throughput_frame, width=5, textvariable=self.T, font=("Helvetica", 15))
        #T_widget.bind('<Return>', self.redraw_from_event)
        T_widget.grid(row=0, column=1, sticky='w')

        T_T_widget = Tk.Entry(self.throughput_frame, width=5, textvariable=self.T_T, font=("Helvetica", 15))
        T_I_widget = Tk.Entry(self.throughput_frame, width=5, textvariable=self.T_I, font=("Helvetica", 15))
        QE_widget = Tk.Entry(self.throughput_frame, width=5, textvariable=self.QE, font=("Helvetica", 15))

        Tk.Label(self.throughput_frame, text='Telescope:', font=("Helvetica", 15)).grid(
            row=1, column=0, sticky='w')
        T_T_widget.grid(row=1, column=1, sticky='w')
        T_T_widget.bind('<Return>', self.redraw_from_event)

        Tk.Label(self.throughput_frame, text='Instrument:', font=("Helvetica", 15)).grid(
            row=1, column=2, sticky='w')
        T_I_widget.grid(row=1, column=3, sticky='w')
        T_I_widget.bind('<Return>', self.redraw_from_event)

        Tk.Label(self.throughput_frame, text='QE:', font=("Helvetica", 15)).grid(
            row=1, column=4, sticky='w')
        QE_widget.grid(row=1, column=5, sticky='w')
        QE_widget.bind('<Return>', self.redraw_from_event)

    def performance_widget(self):
        self.performance_frame = tkinter.Frame(self.master, bd=1, padx=4, pady=4, relief="sunken")
        Tk.Label(self.performance_frame, text='Required Performances', font=("Helvetica", 16, "bold")).grid(
            row=0, column=0, sticky='w')
        # S/N
        Tk.Label(self.performance_frame, text='Desired S/N:', font=("Helvetica", 16)).grid(
            row=1, column=0, sticky='w')
        SN_widget = Tk.Spinbox(self.performance_frame, from_=1, to=1000000, increment=1,
                             width=FWIDTH, textvariable=self.SN, command=self.redraw, font=("Helvetica", 16))
        SN_widget.bind('<Return>', self.redraw_from_event)
        SN_widget.grid(row=1, column=1, sticky='w')

        #Resolving Power
        Tk.Label(self.performance_frame, text='Resolving Power (R):', font=("Helvetica", 16)).grid(
            row=2, column=0, sticky='w')
        R_entry = Tk.Spinbox(self.performance_frame, from_=1, to=1000000, increment=1, width=FWIDTH,
                             textvariable=self.R, command=self.redraw, font=("Helvetica", 16))
        R_entry.grid(row=2, column=1, sticky='w')
        R_entry.bind('<Return>', self.redraw_from_event)

        #Speed evolution
        Tk.Label(self.performance_frame, text='Evolution speed (km/s):', font=("Helvetica", 16)).grid(
            row=3, column=0, sticky='w')
        v_entry = Tk.Spinbox(self.performance_frame,  from_=0, to=1000, increment=0.1,
                             width=FWIDTH, textvariable=self.v, command=self.redraw, font=("Helvetica", 16))
        v_entry.bind('<Return>', self.redraw_from_event)
        v_entry.grid(row=3, column=1, sticky='w')

        # polarimetry radio button
        Tk.Label(self.performance_frame, text='Polarimetry:', font=("Helvetica", 16)).grid(
            row=4, column=0, sticky='w')
        for i, option in zip((0, 1), ('no', 'yes')):
            r = Tk.Radiobutton(self.performance_frame, text=option, variable=self.polarimetry,
                            value=i, command=self.redraw, font=("Helvetica", 16))
            r.grid(row=4, column=i+1, sticky='w')

    def wavelength_widget(self):
        self.wavelength_frame = tkinter.Frame(self.master, bd=1, padx=4, pady=4, relief="sunken")
        self.cwl = Tk.ttk.Combobox(self.wavelength_frame, values=[('Ca II K', 393.33), ('H I beta', 486.10),
            ('Mg I b', 517.30), ('Fe I @525', 525.0), ('Fe I @543', 543.0), ('Na I D @589', 589.60), ('Na I D @590', 590.00),
            ('Fe I @630', 630.2), ('H I alpha', 656.3), ('K I @769', 769.9), ('Fe I@846', 846.80), ('CaII @849', 849.80),
            ('Fe I @851', 851.40), ('Ca II @854', 854.2), ('He I', 1083.0), ('Fe I IR', 1564.8)])
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
        resmin_entry = Tk.Entry(self.pixel_frame, width=FWIDTH, textvariable=self.resmin, font=("Helvetica", 15))
        resmin_entry.grid(row=2, column=1, sticky='w')
        Tk.Label(self.pixel_frame, text='-', font=("Helvetica", 15)).grid(row=2, column=2, sticky='w')
        resmax_entry = Tk.Entry(self.pixel_frame, width=FWIDTH, textvariable=self.resmax, font=("Helvetica", 15))
        resmax_entry.grid(row=2, column=3, sticky='w')
        self.pixel_frame.grid_columnconfigure(2, weight=1)
        self.pixel_frame.grid_rowconfigure(2, weight=1)
        self.pixel_frame.grid_columnconfigure(3, weight=1)
        #self.pixel_frame.grid_rowconfigure(2, weight=1)

        # resulting spectral resolution
        Tk.Label(self.pixel_frame, text='Spectral pixel sampling (pm)', font=("Helvetica", 15)).grid(
            row=3, column=0, sticky='w')
        spresmin_entry = Tk.Entry(self.pixel_frame, width=FWIDTH, textvariable=self.spresmin, font=("Helvetica", 15))
        spresmin_entry.grid(row=3, column=1, sticky='w')
        Tk.Label(self.pixel_frame, text='-', font=("Helvetica", 15)).grid(row=3, column=2, sticky='w')
        spresmax_entry = Tk.Entry(self.pixel_frame, width=FWIDTH, textvariable=self.spresmax, font=("Helvetica", 15))
        spresmax_entry.grid(row=3, column=3, sticky='w')

    def telescope_widget(self):
        self.telescope_frame = tkinter.Frame(self.master, bd=1, padx=4, pady=4, relief="sunken")
        self.telescope_frame.grid_columnconfigure(0, weight=1)
        Tk.Label(self.telescope_frame, text='Telescope Parameters:', font=("Helvetica", 16, 'bold')).grid(
            row=0, column=0, sticky='w')
        Tk.Label(self.telescope_frame, text='Aperture diameter (m):', font=("Helvetica", 15)).grid(
            row=1, column=0, sticky='w')
        D_widget = Tk.Entry(self.telescope_frame, width=FWIDTH, textvariable=self.D, font=("Helvetica", 15))
        D_widget.bind('<Return>', self.redraw_from_event)
        D_widget.grid(row=1, column=1, sticky='w')

        Tk.Label(self.telescope_frame, text='Strehl:', font=("Helvetica", 15)).grid(
            row=2, column=0, sticky='w')
        strehl_widget = Tk.Entry(self.telescope_frame, width=FWIDTH, textvariable=self.strehl, font=("Helvetica", 15))
        strehl_widget.grid(row=2, column=1, sticky='w')

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
        self.T.set(self.T_T.get()*self.T_I.get()*self.QE.get())

        #aux variables
        lmin, lmax = cwl-(self.bandpass.get())/2.0, cwl+(self.bandpass.get())/2.0
        # build dictionary
        properties_dict = {}
        properties_dict['D'] = self.D.get() # telescope diameter
        properties_dict['lmin'] = lmin
        properties_dict['lmax'] = lmax
        properties_dict['polarimetry'] = self.polarimetry.get()
        properties_dict['R'] = self.R.get() # spectral resolution
        properties_dict['T'] = self.T.get() # telescope transmission
        properties_dict['SN'] = self.SN.get() # signal to noise ration
        properties_dict['v'] = self.v.get() * KM_TO_M # m/s
        properties_dict['binning'] = self.binning.get() # spatial binning
        properties_dict['strehl'] = self.strehl.get()

        self.ph.set_properties(properties_dict)
        self.ph.compute()
        self.set_spatres(lmin, lmax)
        self.set_specres(lmin, lmax)

        # plot
        self.plot()

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
        ax.set_xlim(self.ph.ran )
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
    root.resizable(True, True)
    photongui(root)
    root.mainloop()

if __name__ == "__main__":
    create_main_window()