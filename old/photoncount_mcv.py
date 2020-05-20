#!/usr/bin/env python 
from sys import version_info
if version_info[0] == 2:
    import Tkinter as Tk
else:
    import tkinter as Tk
    
import numpy as np
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import fts
from math import pi
from scipy.signal import fftconvolve

######################################################################

class photongui:
    #
    # Computes and displays the integration to reach a desired level
    # of S/N per diffraction limited spatial pixel and given spectral
    # resolution.
    #
    # written by Jorrit Leenaarts
    #

    # units
    hplanck = 6.626e-34 # m**2 kg /s
    clight=2.99792458e8 #m /s
    eV = 1.602e-19 # J
    mAA_to_m = 1e-13
    nm_to_m = 1e-9
    m_to_nm = 1e9
    nm_to_pm = 1e3
    rad_to_arcsec = 206265.
    km_to_m = 1e3
    m_to_arcsec = 1./(7.25e5)

    # plot window sizes 
    window_width = 9 #CQN 22-April-2020
    window_height= 7 #CQN 22-April-2020

    # entry field width
    fwidth=7 #CQN 22-April-2020

    ######################################################################

    def __init__(self, master):

        # load data only once

        # atmospheric transmission at La Palma
        self.al, self. at = self.get_atmostrans()

        # solar intensity W/ (m2 ster m)
        self.atlas = fts.fts()
               
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
        input_label = Tk.Label(self.master, text='Set model parameters here.',font=("Helvetica", 14))
        input_label.grid(row=rowcounter, column=0, columnspan=4,sticky='n')
        rowcounter+=1

        # telescope diameter
        Tk.Label(self.master, text='Aperture diameter (m)',font=("Helvetica", 14)).grid(row=rowcounter, column=0,sticky='w')
        D_entry = Tk.Entry(self.master, width=self.fwidth, textvariable = self.D,font=("Helvetica", 14))
        D_entry.bind('<Return>',self.redraw_from_event)
        D_entry.grid(row=rowcounter, column=1, sticky='w')
        rowcounter+=1

        # lmin,
        Tk.Label(self.master, text='lambda_min/max (nm) [329-1250]',font=("Helvetica", 14)).grid(row=rowcounter, column=0,sticky='w')
        lmin_entry = Tk.Entry(self.master, width=self.fwidth, textvariable = self.lmin,font=("Helvetica", 14))
        lmin_entry.bind('<Return>',self.redraw_from_event)
        lmin_entry.grid(row=rowcounter, column=1)

        #lmax
        Tk.Label(self.master, text='-',font=("Helvetica", 14)).grid(row=rowcounter, column=2)
        lmax_entry = Tk.Entry(self.master, width=self.fwidth, textvariable = self.lmax,font=("Helvetica", 14))
        lmax_entry.bind('<Return>',self.redraw_from_event)
        lmax_entry.grid(row=rowcounter, column=3)
        rowcounter+=1

        # polarimetry radio button
        Tk.Label(self.master, text='polarimetry',font=("Helvetica", 14)).grid(row=rowcounter, column=0,sticky='w')
        for i,option in zip((0,1), ('no','yes')):
            r = Tk.Radiobutton(self.master,text=option, variable=self.polarimetry, 
                            value=i, command=self.redraw,font=("Helvetica", 14))
            r.grid(row=rowcounter,column=i+1)     
        rowcounter+=1

        # spatial binning radio buttons
        Tk.Label(self.master, text='spatial binning',font=("Helvetica", 14)).grid(row=rowcounter, column=0,sticky='w')
        j=1
        for i,option in zip((1.,2.,4.), ('1x1','2x2','4x4')):
            r = Tk.Radiobutton(self.master,text=option, variable=self.binning, 
                            value=i, command=self.redraw,font=("Helvetica", 14))
            r.grid(row=rowcounter,column=j)
            j+=1
        rowcounter+=1

        # spectral resolution
        Tk.Label(self.master, text='spectral resolution R',font=("Helvetica", 14)).grid(row=rowcounter, column=0,sticky='w')
        R_entry = Tk.Entry(self.master, width= self.fwidth, textvariable = self.R,font=("Helvetica", 14))
        R_entry.bind('<Return>',self.redraw_from_event)
        R_entry.grid(row=rowcounter, column=1)
        rowcounter+=1

        # telescope transmission
        Tk.Label(self.master, text='total transmission [0,1]',font=("Helvetica", 14)).grid(row=rowcounter, column=0,sticky='w')
        T_entry = Tk.Entry(self.master, width= self.fwidth, textvariable = self.T,font=("Helvetica", 14))
        T_entry.bind('<Return>',self.redraw_from_event)
        T_entry.grid(row=rowcounter, column=1)
        rowcounter+=1
        
        # S/N
        Tk.Label(self.master, text='desired S/N',font=("Helvetica", 14)).grid(row=rowcounter, column=0,sticky='w')
        T_entry = Tk.Entry(self.master, width= self.fwidth, textvariable = self.SN,font=("Helvetica", 14))
        T_entry.bind('<Return>',self.redraw_from_event)
        T_entry.grid(row=rowcounter, column=1)
        rowcounter+=1

        # evolution speed
        Tk.Label(self.master, text='evolution speed (km/s)',font=("Helvetica", 14)).grid(row=rowcounter, column=0,sticky='w')
        v_entry = Tk.Entry(self.master, width= self.fwidth, textvariable = self.v,font=("Helvetica", 14))
        v_entry.bind('<Return>',self.redraw_from_event)
        v_entry.grid(row=rowcounter, column=1)
        rowcounter+=1

        # resulting spatial resolution
        Tk.Label(self.master, text='resulting spatial pixel size (arcsec)',font=("Helvetica", 14)).grid(row=rowcounter, column=0,sticky='w')
        resmin_entry = Tk.Entry(self.master, width=self.fwidth, textvariable = self.resmin,font=("Helvetica", 14))
        resmin_entry.grid(row=rowcounter, column=1)
        Tk.Label(self.master, text='-',font=("Helvetica", 14)).grid(row=rowcounter, column=2)
        resmax_entry = Tk.Entry(self.master, width=self.fwidth, textvariable = self.resmax,font=("Helvetica", 14))
        resmax_entry.grid(row=rowcounter, column=3)
        rowcounter+=1

       # resulting spectral resolution
        Tk.Label(self.master, text='resulting spectral pixel size (pm)',font=("Helvetica", 14)).grid(row=rowcounter, column=0,sticky='w')
        spresmin_entry = Tk.Entry(self.master, width=self.fwidth, textvariable = self.spresmin,font=("Helvetica", 14))
        spresmin_entry.grid(row=rowcounter, column=1)
        Tk.Label(self.master, text='-',font=("Helvetica", 14)).grid(row=rowcounter, column=2)
        spresmax_entry = Tk.Entry(self.master, width=self.fwidth, textvariable = self.spresmax,font=("Helvetica", 14))
        spresmax_entry.grid(row=rowcounter, column=3)
        rowcounter+=1

        
        # quit button
        quit_button = Tk.Button(self.master, text="Quit", command=self.quit,font=("Helvetica", 14))
        quit_button.grid(row=rowcounter, column=0, columnspan=4,sticky='ew')

        return rowcounter

    ######################################################################

    def init_plot(self, rowspan):

        # create Figure for the plot 
        self.fig = Figure(figsize=(self.window_width, self.window_height), dpi=100)
        self.axes=[]

        nrow = 2
        ncol = 2
        for i in range(nrow*ncol):
            self.axes.append(self.fig.add_subplot(nrow, ncol, i+1))

        # link figure to gui
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.master)
        self.canvas.get_tk_widget().grid(row=0,column=4,rowspan=rowspan)
        self.canvas.draw() #CQN 22-April-2020

    ######################################################################
        
    def init_parameters(self):

        D_init = 4.0
        self.D = Tk.DoubleVar()
        self.D.set(D_init)

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

    def spatres(self, l):
        r = 1.22 * l / self.D.get() * self.rad_to_arcsec /2.
        return r

    ######################################################################

    def set_spatres(self, lmin, lmax):

        self.resmin.set("{:.4f}".format(self.spatres(lmin*self.nm_to_m)))
        self.resmax.set("{:.4f}".format(self.spatres(lmax*self.nm_to_m)))

    ######################################################################

    def set_specres(self, lmin, lmax):

        a = lmin/self.R.get()/2.0*self.nm_to_pm
        self.spresmin.set("{:.4f}".format(a))
        a = lmax/self.R.get()/2.0*self.nm_to_pm
        self.spresmax.set("{:.4f}".format(a))

    ######################################################################

    def redraw_from_event(self, event):

        self.redraw()

    ######################################################################

    def redraw(self):

        # read the input parameters from buttons and fields
        D = self.D.get() # telescope diameter
        lmin = self.lmin.get() # wavelength range
        lmax = self.lmax.get() 
        R  = self.R.get() # spectral resolution 
        SN = self.SN.get() # signal to noise ration     
        telescopetrans = self.T.get() # telescope transmission
        binning = self.binning.get() # spatial binning
        v = self.v.get() * self.km_to_m # m/s

        #get polarimetry. should account for 4 polarization
        #states per wavelength in a balanced design with efficiency
        # 1/sqrt(3) = 0.577
        if (self.polarimetry.get() == 0):
            pfac = 1.0
        else:
#            pfac = 1.0 / 3.0**0.5
            pfac = 1.0 / 3.0
#           pfac = 0.16

        # sanity check on wavelengths
        if lmin >= lmax:
            lmin = lmax - 1.0

        # set wavelength range. Extend slightly to avoid boundary effects in convolution
        ran = (lmin, lmax)
        lmin -= 0.1
        lmax += 0.1

                
        # cut out selected  wavelength range
        iw= np.where( (self.atlas.ll > lmin * self.nm_to_m) &
                      (self.atlas.ll < lmax * self.nm_to_m))[0]
        Ilambda = self.atlas.ii[iw]
        ll = self.atlas.ll[iw]

        # interpolate atmospheric transmission
        atrans = np.interp(ll, self.al, self.at) 

        # wavelength spacing
        dl = ll / R
        
        # spectrally convolved spectrum, assume wavelength spacing is
        # equidistant.  the spacing changes slowly over the entire
        # wavelength range so this is approx OK.
        iw = int( len(ll) / 2)
        dl2 = ll[iw+1] - ll[iw]
        ff = self.get_filter(ll, ll[iw] / R)
        Ilambdac = fftconvolve(Ilambda, ff, mode = 'same')

        # energy per photon [J]
        ephot = self.hplanck * self.clight / ll

        
        # nphot = Ilambda / ephot * A              * T * (delta omega)**2               * (delta lambda)
        #       = Ilambda / ephot * f0 * pi * R**2 * T * 1.22**2 * lambda**2 / (4 R**2) * (delta lambda)
        #       = Ilambda / ephot * f0 * pi        * T * 1.22**2 * lambda**2 / 4        * (delta lambda)

        # f0 is the fraction of the aperture that is unobscured. I
        # assume 1.0 here

        # number of photons per second per m per diffraction limited
        # spatial resolution element of the light entering the
        # telescope
        nflux = Ilambda / ephot * pi * atrans * 1.22**2 / 4.0 * ll**2

        # number of photons per second per spatial pixel
        nflux /=  4.0 
         # per spectral pixel
        nflux *= (dl / 2.0)
         # factor in spatial binning
        nflux *= binning**2
         # and factor in polarimetry
        nflux *= pfac
 
        nfluxideal = np.copy(nflux) # ideal telescope
        nflux *= telescopetrans # for given total transmission

        # time to reach desired S/N, computed from: S/N = sqrt(nflux * t)
        t = SN**2 / nflux
        tideal = SN**2 / nfluxideal

        # display spatial resolution
        self.set_spatres(lmin * binning , lmax * binning)
        self.set_specres(lmin, lmax)

        # compute Alex Feller ideal dt,dx
        phi = Ilambda / ephot * pi / 4.0 * D**2 * atrans * telescopetrans * pfac # photons/ (s ster m)
        phi = phi * self.nm_to_m / self.rad_to_arcsec**2 # photons/ (s arcsec**2 nm)
        dt = SN**2 / phi / (dl/self.nm_to_m/2.0) / (v*self.m_to_arcsec)**2
        dt = dt**(1./3.)
        dx = (v*self.m_to_arcsec) * dt 

        # plot
        self.plot(ll, t, tideal, Ilambda, Ilambdac, ran, dx, dt)


    ######################################################################

    def plot(self, ll, t, tideal, Ilambda, Ilambdac, ran, dx, dt):

        xax = ll * self.m_to_nm

        for ax in self.axes:
            ax.clear()

    
        ax = self.axes[0]
        ax.plot(xax, Ilambda ,label = 'atlas')
        ax.plot(xax, Ilambdac, label = 'smeared to R')
        ax.legend(loc='best')
        ax.get_xaxis().set_ticklabels([])
        ax.set_xlim( ran )
        ax.set_ylabel(r'$I_\lambda$ [W m$^{-2}$ m$^{-1}$ sr$^{-1}$]')
        ax.set_title('spectrum')
        ax.grid(True)

        ax = self.axes[1]
        ax.plot(xax, dx)
        ax.set_xlim( ran )
        ax.get_xaxis().set_ticklabels([])
        ax.set_ylabel(r'$\Delta x$ [arcsec]')
        ax.set_title('optimal pixel size for given signal speed')
        ax.grid(True)

        ax=self.axes[2]
        ax.plot(xax,t, label='with given transmission')
        ax.plot(xax,tideal, label='perfect telescope')
        ax.legend(loc='best')
        ax.set_xlabel(r'$\lambda$ (nm)')
        ax.set_xlim( ran )
        ax.set_ylabel(r' $\Delta t$ [s]')
        ax.set_title('integration time for given spatial and spectral pixel size')
        ax.grid(True)

        ax = self.axes[3]
        ax.plot(xax, dt )
        ax.set_xlim( ran )
        ax.set_ylabel(r'$\Delta t$ [s]')
        ax.set_title('optimal integration time for given signal speed')
        ax.set_xlabel(r'$\lambda$ (nm)')
        ax.grid(True)



        self.fig.tight_layout()


        self.canvas.draw() #CQN 22-April-2020
      
    ######################################################################  

    def get_atmostrans(self, M = 2.0):
        #
        # Reads atmospheric transmission at La Palma from "D. L. King,
        # Atmospheric Extinction at the Roque de los Muchachos
        # Observatory, La Palma, RGO/La Palma technical note no 31"
        # available at http://www.ing.iac.es/Astronomy/observing/manuals/man_tn.html
        #
        # M sets the elevation. M=2 means 30 degrees, zenith is M=1
        # 

        a=np.loadtxt('atmostrans.dat')
        ll = a[:,0] *1e-10 # from AA to m
        T = a[:,1]
        TT = 10.**(-0.4 * T * M )
        return ll, TT

    ######################################################################

    def get_filter(self, ll, fwhm):
        #
        # just a normalized gaussian
        #
        llc= ll[int(len(ll)/2)]
        sigma = fwhm / (2.0 * (2.0 * np.log(2.0))**0.5)
        trans = np.exp(-0.5 * ((ll - llc) / sigma)**2)
        trans /= trans.sum()
        return trans

######################################################################
        
root = Tk.Tk()
my_gui = photongui(root)
root.mainloop()
