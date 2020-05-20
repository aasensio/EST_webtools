# Hello, Flask!
from flask import Flask, render_template, request
import photoncount
from bokeh.plotting import figure
from bokeh.embed import components, server_document
from bokeh.layouts import gridplot, row, column
from bokeh.models import Slider, ColumnDataSource, RadioGroup, TextInput, LabelSet
from bokeh.themes import Theme
from bokeh.server.server import Server
from tornado.ioloop import IOLoop
from functools import partial

ph = photoncount.Photocount()

app = Flask(__name__)

properties_dict = {}
properties_dict['D'] = 4.0
properties_dict['lmin'] = 854.0
properties_dict['lmax'] = 855.0
properties_dict['polarimetry'] = 0
properties_dict['R'] = 8e4
properties_dict['T'] = 0.1
properties_dict['SN'] = 1e3
properties_dict['v'] = 7.0
properties_dict['binning'] = 1.0

def modify_doc(doc):

    ph.set_properties(properties_dict)
    ph.compute()
    xax = ph.ll * ph.m_to_nm

    source = ColumnDataSource(data={'x': xax, 'Ilambda': ph.Ilambda, \
        'Ilambdac': ph.Ilambdac, 'dx': ph.dx, 't': ph.t, 'tideal': ph.tideal,\
        'dt': ph.dt})
    
    plot1 = figure(plot_height=400, sizing_mode='scale_width', title="Spectrum")     
    plot1.line(x='x', y='Ilambda', legend_label="Atlas", source=source)
    plot1.line(x='x', y='Ilambdac', color="#FF0000", legend_label="Degraded atlas", source=source)
    plot1.xaxis.axis_label = "\u03BB [nm]"
    plot1.yaxis.axis_label = "I(\u03BB) [W/m2/m/sr]"
    plot1.legend.location = "bottom_right"
    plot1.legend.click_policy="hide"

    plot2 = figure(plot_height=400, sizing_mode='scale_width', title="Optimal pixel size")
    plot2.line(x='x', y='dx', source=source)
    plot2.xaxis.axis_label = "\u03BB [nm]"
    plot2.yaxis.axis_label = "\u0394 x [arcsec]"

    plot3 = figure(plot_height=400, sizing_mode='scale_width', title="Integration time")
    plot3.line(x='x', y='t', source=source, legend_label="With transmission")
    plot3.line(x='x', y='tideal', source=source, color="#FF0000", legend_label="Perfect telescope")
    plot3.xaxis.axis_label = "\u03BB [nm]"
    plot3.yaxis.axis_label = "\u0394 t [s]"
    plot3.legend.location = "top_right"
    plot3.legend.click_policy="hide"

    plot4 = figure(plot_height=400, sizing_mode='scale_width', title="Optimal integration time")     
    plot4.line(x='x', y='dt', source=source)
    plot4.xaxis.axis_label = "\u03BB [nm]"
    plot4.yaxis.axis_label = "\u0394 t [s]"

    grid = gridplot([[plot1, plot2], [plot3, plot4]], plot_width=350, plot_height=350)

    def update():
        ph.set_properties(properties_dict)
        ph.compute()

        data = {'x': ph.ll * ph.m_to_nm, 'Ilambda': ph.Ilambda, \
        'Ilambdac': ph.Ilambdac, 'dx': ph.dx, 't': ph.t, 'tideal': ph.tideal,\
        'dt': ph.dt}
        
        source.data = data

        source2.data = {'resmin': [ph.resmin, ph.resmax]}

    def callback(attr, old, new, parameter):
        if (parameter in ['D', 'lmin', 'lmax', 'R', 'T', 'SN', 'v', 'polarimetry']):
            properties_dict[parameter] = float(new)        
        if (parameter == 'binning'):
            properties_dict[parameter] = float(2**new)
        update()

        
    # slider = Slider(start=852, end=860, value=852, step=0.1, title="Smoothing by N Days")
    # slider.on_change('value', callback)

    widget_diameter = TextInput(value="4.0", title="Diameter [m]")
    widget_diameter.on_change('value', partial(callback, parameter='D'))

    widget_lmin = TextInput(value="854.0", title="Minimum wavelength")
    widget_lmin.on_change('value', partial(callback, parameter='lmin'))

    widget_lmax = TextInput(value="855.0", title="Maximum wavelength")
    widget_lmax.on_change('value', partial(callback, parameter='lmax'))

    widget_lambda = row(widget_lmin, widget_lmax)    

    widget_pol = RadioGroup(labels=["No", "Yes"], active=0, orientation='horizontal')
    widget_pol.on_change('active', partial(callback, parameter='polarimetry'))

    widget_bin = RadioGroup(labels=["1x1", "2x2", "4x4"], active=0)
    widget_bin.on_change('active', partial(callback, parameter='binning'))

    widget_R = TextInput(value="80000.0", title="Spectral resolution")
    widget_R.on_change('value', partial(callback, parameter='R'))

    widget_T = TextInput(value="0.1", title="Total transmission [0,1]")
    widget_T.on_change('value', partial(callback, parameter='T'))

    widget_SN = TextInput(value="1000.0", title="Desired S/N")
    widget_SN.on_change('value', partial(callback, parameter='SN'))

    widget_v = TextInput(value="7.0", title="Evolution speed [km/s]")
    widget_v.on_change('value', partial(callback, parameter='v'))

    widgets = column(widget_diameter, widget_lambda, widget_pol, widget_bin, \
        widget_R, widget_T, widget_SN, widget_v)

    doc.add_root(row(widgets, grid))

    doc.theme = Theme(filename="theme.yaml")

# Index page, no args
@app.route('/', methods=['GET'])
def bkapp_page():

    script = server_document('http://localhost:5006/bkapp')
    return render_template("index.html", script=script, template='Flask')

def bk_worker():

    server = Server({'/bkapp': modify_doc}, io_loop=IOLoop(), allow_websocket_origin=["localhost:8000"])
    server.start()
    server.io_loop.start()

from threading import Thread
Thread(target=bk_worker).start()
	    
    # plot = update_figure(properties_dict)

    # script, div = components(plot)

    # return render_template("index.html", script=script, div=div)

# With debug=True, Flask server will auto-reload 
# when there are code changes
if __name__ == '__main__':
	app.run(port=8000, debug=False)