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
    rad_to_arcsec = 206265.

    # telescope diameter
    D=4.0 #m

    def __init__(self, master):

        # load data only once

        # atmospheric transmission at La Palma
        self.al, self. at = self.get_atmostrans()

        # solar intensity W/ (m2 ster m)
        self.atlas = fts.fts()
       
        # pot window sizes 
        window_width = 8
        window_height= 8
        
        self.master = master
        master.title("EST - integration time estimator")

        # define and initialize buttons and fields
        lmin_init=853.0
        self.lmin = Tk.DoubleVar()
        self.lmin.set(lmin_init)

        lmax_init=855.0
        self.lmax = Tk.DoubleVar()
        self.lmax.set(lmax_init)
    
        self.polarimetry = Tk.IntVar()
        self.polarimetry.set(0)

        R_init = 1e5
        self.R = Tk.DoubleVar()
        self.R.set(R_init)

        T_init = 0.1
        self.T  = Tk.DoubleVar()
        self.T.set(T_init)

        SN_init = 1e3
        self.SN = Tk.DoubleVar()
        self.SN.set(SN_init)

        binning_init = 1.0
        self.binning = Tk.DoubleVar()
        self.binning.set(binning_init)

        self.resmin = Tk.DoubleVar()
        self.resmax = Tk.DoubleVar()
        self.display_spatres(lmin_init, lmax_init)
        
        # create Figure for the plot 
        self.fig = Figure(figsize=(window_width, window_height), dpi=100)
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.master)
        self.canvas.get_tk_widget().grid(row=0,column=4,rowspan=8)
        self.canvas.show()
 
        #start adding buttons and text
        rowcounter = 0

        #window label
        input_label = Tk.Label(self.master, text='Set model parameters here.')
        input_label.grid(row=rowcounter, column=0, columnspan=4,sticky='ew')
        rowcounter+=1

        # lmin,
        Tk.Label(self.master, text='lambda_min (nm) [>329]').grid(row=rowcounter, column=0)
        lmin_entry = Tk.Entry(self.master, width=6, textvariable = self.lmin)
        lmin_entry.bind('<Return>',self.redraw_from_event)
        lmin_entry.grid(row=rowcounter, column=1)

        #lmax
        Tk.Label(self.master, text='lambda_max (nm) [<1250]').grid(row=rowcounter, column=2)
        lmax_entry = Tk.Entry(self.master, width=6, textvariable = self.lmax)
        lmax_entry.bind('<Return>',self.redraw_from_event)
        lmax_entry.grid(row=rowcounter, column=3)
        rowcounter+=1

        # polarimetry radio button
        Tk.Label(self.master, text='polarimetry').grid(row=rowcounter, column=0)
        for i,option in zip((0,1), ('no','yes')):
            r = Tk.Radiobutton(self.master,text=option, variable=self.polarimetry, 
                            value=i, command=self.redraw)
            r.grid(row=rowcounter,column=i+1)     
        rowcounter+=1

        # spatial binning radio buttons
        Tk.Label(self.master, text='spatial binning').grid(row=rowcounter, column=0)
        j=1
        for i,option in zip((1.,2.,4.), ('1x1','2x2','4x4')):
            r = Tk.Radiobutton(self.master,text=option, variable=self.binning, 
                            value=i, command=self.redraw)
            r.grid(row=rowcounter,column=j)
            j+=1
        rowcounter+=1

        # resulting spatial resolution
        Tk.Label(self.master, text='resulting spatial resolution (arcsec)').grid(row=rowcounter, column=0)
        resmin_entry = Tk.Entry(self.master, width=6, textvariable = self.resmin)
        resmin_entry.grid(row=rowcounter, column=1)
        Tk.Label(self.master, text='--').grid(row=rowcounter, column=2)
        resmax_entry = Tk.Entry(self.master, width=6, textvariable = self.resmax)
        resmax_entry.grid(row=rowcounter, column=3)
        rowcounter+=1

        # spectral resolution
        Tk.Label(self.master, text='spectral resolution R').grid(row=rowcounter, column=0)
        R_entry = Tk.Entry(self.master, width=8, textvariable = self.R)
        R_entry.bind('<Return>',self.redraw_from_event)
        R_entry.grid(row=rowcounter, column=1)
        rowcounter+=1

        # telescope transmission
        Tk.Label(self.master, text='total transmission [0,1]').grid(row=rowcounter, column=0)
        T_entry = Tk.Entry(self.master, width=8, textvariable = self.T)
        T_entry.bind('<Return>',self.redraw_from_event)
        T_entry.grid(row=rowcounter, column=1)
        rowcounter+=1
        
        # S/N
        Tk.Label(self.master, text='desired S/N').grid(row=rowcounter, column=0)
        T_entry = Tk.Entry(self.master, width=8, textvariable = self.SN)
        T_entry.bind('<Return>',self.redraw_from_event)
        T_entry.grid(row=rowcounter, column=1)
        rowcounter+=1
        
        # quit button
        quit_button = Tk.Button(master, text="Quit", command=self.quit)
        quit_button.grid(row=rowcounter, column=0, columnspan=4,sticky='ew')
        rowcounter+=1

        self.redraw()

    # quit button
    def quit(self):
        self.master.destroy()

    def spatres(self, l):
        r = 1.22 * l / self.D * self.rad_to_arcsec
        return r

    def display_spatres(self, lmin, lmax):

        self.resmin.set("{:.4f}".format(self.spatres(lmin*self.nm_to_m)))
        self.resmax.set("{:.4f}".format(self.spatres(lmax*self.nm_to_m)))

    def redraw_from_event(self, event):

        self.redraw()

    def redraw(self):

        # read the input parameters from buttons and fields
        lmin = self.lmin.get() # wavelength range
        lmax = self.lmax.get() 
        R  = self.R.get() # spectral resolution 
        SN = self.SN.get() # signal to noise ration     
        telescopetrans = self.T.get() # telescope transmission
        binning = self.binning.get() # spatial binning

        #get polarimetry. should account for 4 polarization
        #states per wavelength in a balanced design with efficiency
        # 1/sqrt(3) = 0.577
        if (self.polarimetry.get() == 0):
            pfac = 1.0
        else:
            pfac = 1.0 / 3.0**0.5

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
        self.display_spatres(lmin * binning , lmax * binning)

        # plot
        self.plot(ll, t, tideal, Ilambda, Ilambdac, ran)


    def plot(self, ll, t, tideal, Ilambda, Ilambdac, ran):

        xax = ll * self.m_to_nm
        self.fig.clear()
        ax = self.fig.add_subplot(212)
        ax.plot(xax,t, label='with given transmission')
        ax.plot(xax,tideal, label='perfect telescope')
        ax.legend(loc='best')
        ax.set_xlabel(r'$\lambda$ (nm)')
        ax.set_xlim( ran )
        ax.set_ylabel(r' $t_\mathrm{integration}$ [s]')
        ax.grid(True)

        ax = self.fig.add_subplot(211)
        ax.plot(xax, Ilambda ,label = 'atlas')
        ax.plot(xax, Ilambdac, label = 'smeared to R')
        ax.legend(loc='best')
        ax.set_xlabel(r'$\lambda$ (nm)')
        ax.set_xlim( ran )
        ax.set_ylabel(r'$I_\lambda$ [W m$^{-2}$ m$^{-1}$ sr$^{-1}$]')
        ax.grid(True)

        self.canvas.show()
        

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

    def get_filter(self, ll, fwhm):
        #
        # just a normalized gaussian
        #
        llc= ll[int(len(ll)/2)]
        sigma = fwhm / (2.0 * (2.0 * np.log(2.0))**0.5)
        trans = np.exp(-0.5 * ((ll - llc) / sigma)**2)
        trans /= trans.sum()
        return trans
        
root = Tk.Tk()
my_gui = photongui(root)
root.mainloop()
