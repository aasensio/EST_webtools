#!/usr/bin/env python
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
WINDOW_WIDTH = 12
WINDOW_HEIGHT = 9

#

D_INIT = 4.0

# entry field width
FWIDTH = 8

class photongui:
    #
    # Computes and displays the integration to reach a desired level
    # of S/N per diffraction limited spatial pixel and given spectral
    # resolution.
    #
    # written by Jorrit Leenaarts
    def __init__(self, master):
        self.ph = photoncount.Photocount()

        # wavelength parameters
        self.cwl_shift = 0.0
        self.bandpass = 0.5

        self.cwl = None

        self.master = master
        self.master.title("Generic solar telescope integration time estimator              J. Leenaarts - Institute for Solar Physics")

        self.init_parameters()
        rowcounter = self.init_buttons()
        self.init_plot(rowcounter)

        self.redraw()

    ######################################################################

    def init_buttons(self):
 
        #start adding buttons and text
        rowcounter = 0

        #window label
        input_label = Tk.Label(self.master, text='Set model parameters here.', font=("Helvetica", 20))
        input_label.grid(row=rowcounter, column=0, columnspan=4, sticky='n')
        rowcounter+=1

        # telescope diameter
        Tk.Label(self.master, text='Aperture diameter (m)',font=("Helvetica", 17)).grid(row=rowcounter, column=0, sticky='w')
        D_entry = Tk.Entry(self.master, width=FWIDTH, textvariable=self.D, font=("Helvetica", 17))
        D_entry.bind('<Return>', self.redraw_from_event)
        D_entry.grid(row=rowcounter, column=1, sticky='w')
        rowcounter+=1

        # target wavelengths combobox
        Tk.Label(self.master, text='Target Wavelength', font=("Helvetica", 18)).grid(
            row=rowcounter, column=0, sticky='w')
        self.cwl = Tk.ttk.Combobox(self.master, values=[('Ca II K', 393.33), ('H I beta', 486.1), ('Mg I b', 517.30), ('Fe I', 525.0 ), ('Fe I V', 630.2 ),
                                                   ('H I alpha', 656.3 ), ('K I', 769.9), ('Ca II', 854.2), ('He I',1083.0),
                                                   ('Fe I IR', 1564.8)])
        self.cwl.current(7)
        self.cwl.bind('<<ComboboxSelected>>', self.redraw_from_event)
        self.cwl.grid(row=rowcounter,column=1)
        bp = Tk.Entry(self.master, width=FWIDTH, textvariable=self.bandpass, font=("Helvetica", 17))
        bp.grid(row=rowcounter, column=2)
        rowcounter+=1

        # lmin,
        Tk.Label(self.master, text='lambda_min/max (nm) [329-1250]',font=("Helvetica", 18)).grid(row=rowcounter, column=0,sticky='w')
        lmin_entry = Tk.Entry(self.master, width=FWIDTH, textvariable=self.lmin, font=("Helvetica", 17))
        lmin_entry.bind('<Return>',self.redraw_from_event)
        lmin_entry.grid(row=rowcounter, column=1)

        #lmax
        Tk.Label(self.master, text='-', font=("Helvetica", 17)).grid(row=rowcounter, column=2)
        lmax_entry = Tk.Entry(self.master, width=FWIDTH, textvariable=self.lmax, font=("Helvetica", 17))
        lmax_entry.bind('<Return>', self.redraw_from_event)
        lmax_entry.grid(row=rowcounter, column=3)
        rowcounter+=1

        # polarimetry radio button
        Tk.Label(self.master, text='polarimetry',font=("Helvetica", 17)).grid(row=rowcounter, column=0, sticky='w')
        for i,option in zip((0,1), ('no', 'yes')):
            r = Tk.Radiobutton(self.master, text=option, variable=self.polarimetry,
                            value=i, command=self.redraw, font=("Helvetica", 17))
            r.grid(row=rowcounter, column=i+1)
        rowcounter+=1

        """
        spatial binning radio buttons
        """
        Tk.Label(self.master, text='spatial binning',font=("Helvetica", 17)).grid(row=rowcounter, column=0,sticky='w')
        j=1
        for i,option in zip((1.,3.,4.), ('1x1','3x3','4x4')):
            r = Tk.Radiobutton(self.master,text=option, variable=self.binning, 
                            value=i, command=self.redraw,font=("Helvetica", 17))
            r.grid(row=rowcounter,column=j)
            j+=1
        rowcounter+=1


        # spectral resolution
        Tk.Label(self.master, text='spectral resolution R',font=("Helvetica", 17)).grid(row=rowcounter, column=0,sticky='w')
        R_entry = Tk.Entry(self.master, width=FWIDTH, textvariable=self.R, font=("Helvetica", 17))
        R_entry.bind('<Return>',self.redraw_from_event)
        R_entry.grid(row=rowcounter, column=1)
        rowcounter+=1

        # telescope transmission
        Tk.Label(self.master, text='total transmission [0,1]',font=("Helvetica", 17)).grid(row=rowcounter, column=0,sticky='w')
        T_entry = Tk.Entry(self.master, width=FWIDTH, textvariable=self.T, font=("Helvetica", 17))
        T_entry.bind('<Return>', self.redraw_from_event)
        T_entry.grid(row=rowcounter, column=1)
        rowcounter+=1
        
        # S/N
        Tk.Label(self.master, text='desired S/N',font=("Helvetica", 17)).grid(row=rowcounter, column=0,sticky='w')
        T_entry = Tk.Entry(self.master, width=FWIDTH, textvariable = self.SN, font=("Helvetica", 17))
        T_entry.bind('<Return>', self.redraw_from_event)
        T_entry.grid(row=rowcounter, column=1)
        rowcounter+=1

        # evolution speed
        Tk.Label(self.master, text='evolution speed (km/s)', font=("Helvetica", 17)).grid(row=rowcounter, column=0, sticky='w')
        v_entry = Tk.Entry(self.master, width=FWIDTH, textvariable=self.v, font=("Helvetica", 17))
        v_entry.bind('<Return>', self.redraw_from_event)
        v_entry.grid(row=rowcounter, column=1)
        rowcounter+=1

        # resulting spatial resolution
        Tk.Label(self.master, text='resulting spatial pixel size (arcsec)',font=("Helvetica", 17)).grid(row=rowcounter, column=0,sticky='w')
        resmin_entry = Tk.Entry(self.master, width=FWIDTH, textvariable=self.resmin, font=("Helvetica", 17))
        resmin_entry.grid(row=rowcounter, column=1)
        Tk.Label(self.master, text='-',font=("Helvetica", 17)).grid(row=rowcounter, column=2)
        resmax_entry = Tk.Entry(self.master, width=FWIDTH, textvariable=self.resmax, font=("Helvetica", 17))
        resmax_entry.grid(row=rowcounter, column=3)
        rowcounter+=1

       # resulting spectral resolution
        Tk.Label(self.master, text='resulting spectral pixel size (pm)',font=("Helvetica", 17)).grid(row=rowcounter, column=0,sticky='w')
        spresmin_entry = Tk.Entry(self.master, width=FWIDTH, textvariable=self.spresmin, font=("Helvetica", 17))
        spresmin_entry.grid(row=rowcounter, column=1)
        Tk.Label(self.master, text='-',font=("Helvetica", 17)).grid(row=rowcounter, column=2)
        spresmax_entry = Tk.Entry(self.master, width=FWIDTH, textvariable=self.spresmax, font=("Helvetica", 17))
        spresmax_entry.grid(row=rowcounter, column=3)
        rowcounter+=1

        
        # quit button
        quit_button = Tk.Button(self.master, text="Quit", command=self.quit, font=("Helvetica", 18))
        quit_button.grid(row=rowcounter, column=0, columnspan=4,sticky='ew')

        return rowcounter

    ######################################################################

    def init_plot(self, rowspan):

        # create Figure for the plot 
        self.fig = Figure(figsize=(WINDOW_WIDTH, WINDOW_HEIGHT), dpi=100)
        self.axes=[]

        nrow = 2
        ncol = 2
        for i in range(nrow*ncol):
            self.axes.append(self.fig.add_subplot(nrow, ncol, i+1))

        # link figure to gui
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.master)
        self.canvas.get_tk_widget().grid(row=0,column=4,rowspan=rowspan)
        self.canvas.draw()

    ######################################################################
        
    def init_parameters(self):

        self.D = Tk.DoubleVar()
        self.D.set(D_INIT)

        lmin_init=854.0
        self.lmin = Tk.DoubleVar()
        self.lmin.set(lmin_init)

        lmax_init=855.0
        self.lmax = Tk.DoubleVar()
        self.lmax.set(lmax_init)
    
        self.polarimetry = Tk.IntVar()
        self.polarimetry.set(0)

        R_init = 8e4
        self.R = Tk.DoubleVar()
        self.R.set(R_init)

        T_init = 0.1
        self.T  = Tk.DoubleVar()
        self.T.set(T_init)

        SN_init = 1e3
        self.SN = Tk.DoubleVar()
        self.SN.set(SN_init)

        v_init = 7.0
        self.v = Tk.DoubleVar()
        self.v.set(v_init)

        binning_init = 1.0
        self.binning = Tk.DoubleVar()
        self.binning.set(binning_init)

        self.resmin = Tk.DoubleVar()
        self.resmax = Tk.DoubleVar()
        self.set_spatres(lmin_init, lmax_init)

        self.spresmin = Tk.DoubleVar()
        self.spresmax = Tk.DoubleVar()
        self.set_specres(lmin_init, lmax_init)
 
    ######################################################################

    def quit(self):
        self.master.destroy()

    ######################################################################

    def set_spatres(self, lmin, lmax):

        self.resmin.set("{:.4f}".format(self.ph.spatres(lmin*NM_TO_M, self.D.get())))
        self.resmax.set("{:.4f}".format(self.ph.spatres(lmax*NM_TO_M, self.D.get())))

    ######################################################################

    def set_specres(self, lmin, lmax):
        self.spresmin.set("{:.4f}".format(self.ph.specres(lmin, self.R.get())))
        self.spresmax.set("{:.4f}".format(self.ph.specres(lmax, self.R.get())))

    ######################################################################

    def redraw_from_event(self, event):

        self.redraw()

    ######################################################################

    def redraw(self):
        print(self.cwl.get())
        # read the input parameters from buttons and fields
        D = self.D.get() # telescope diameter
        lmin = self.lmin.get() # wavelength range
        lmax = self.lmax.get() 
        R  = self.R.get() # spectral resolution 
        SN = self.SN.get() # signal to noise ration     
        telescopetrans = self.T.get() # telescope transmission
        binning = self.binning.get() # spatial binning
        v = self.v.get() * KM_TO_M # m/s
        # get polarimetry. should account for 4 polarization
        # states per wavelength in a balanced design with efficiency
        # 1/sqrt(3) = 0.577
        if (self.polarimetry.get() == 0):
            pfac = 1.0
        else:
            #            pfac = 1.0 / 3.0**0.5
            pfac = 1.0 / 3.0
        #           pfac = 0.16

        # build dictionary
        properties_dict = {}
        properties_dict['D'] = D
        properties_dict['lmin'] = lmin #cwl.value + cwl_shift.value - (bp.value / 2.0)
        properties_dict['lmax'] = lmax #cwl.value + cwl_shift.value + (bp.value / 2.0)
        properties_dict['polarimetry'] = pfac
        properties_dict['R'] = R# 8e4
        properties_dict['T'] = telescopetrans #0.1
        properties_dict['SN'] = SN #1e3
        properties_dict['v'] = v#7.0
        properties_dict['binning'] = binning # 1.0
        properties_dict['strehl'] = 1.0

        self.ph.set_properties(properties_dict)
        self.ph.compute()

        # plot

        self.plot()


    ######################################################################

    def plot(self):

        xax = self.ph.ll * M_TO_NM

        for ax in self.axes:
            ax.clear()

    
        ax = self.axes[0]
        ax.plot(xax, self.ph.Ilambda ,label = 'atlas')
        ax.plot(xax, self.ph.Ilambdac, label = 'smeared to R')

        ax.legend(loc='best')
        ax.get_xaxis().set_ticklabels([])
        ax.set_xlim( self.ph.ran )
        ax.set_ylabel(r'$I_\lambda$ [W m$^{-2}$ m$^{-1}$ sr$^{-1}$]')
        ax.set_title('spectrum')
        ax.grid(True)

        ax = self.axes[1]
        ax.plot(xax, self.ph.dx)
        ax.set_xlim( self.ph.ran )
        ax.get_xaxis().set_ticklabels([])
        ax.set_ylabel(r'$\Delta x$ [arcsec]')
        ax.set_title('optimal pixel size for given signal speed')
        ax.grid(True)

        ax=self.axes[2]
        ax.plot(xax,self.ph.t, label='with given transmission')
        ax.plot(xax,self.ph.tideal, label='perfect telescope')
        ax.legend(loc='best')
        ax.set_xlabel(r'$\lambda$ (nm)')
        ax.set_xlim( self.ph.ran )
        ax.set_ylabel(r' $\Delta t$ [s]')
        ax.set_title('integration time for given spatial and spectral pixel size')
        ax.grid(True)

        ax = self.axes[3]
        ax.plot(xax, self.ph.dt )
        ax.set_xlim( self.ph.ran )
        ax.set_ylabel(r'$\Delta t$ [s]')
        ax.set_title('optimal integration time for given signal speed')
        ax.set_xlabel(r'$\lambda$ (nm)')
        ax.grid(True)

        self.fig.tight_layout()

        self.canvas.draw()

######################################################################
        
root = Tk.Tk()
my_gui = photongui(root)
root.mainloop()
