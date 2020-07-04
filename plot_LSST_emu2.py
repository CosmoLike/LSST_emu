import getdist.plots as gplot
from getdist import MCSamples
from getdist import loadMCSamples
import os
import matplotlib
import subprocess
import matplotlib.pyplot as plt
import numpy as np
### ------------------------------------------------
### ------------------------------------------------
### ------------------------------------------------
# GENERAL PLOT OPTIONS
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams['mathtext.rm'] = 'Bitstream Vera Sans'
matplotlib.rcParams['mathtext.it'] = 'Bitstream Vera Sans:italic'
matplotlib.rcParams['mathtext.bf'] = 'Bitstream Vera Sans:bold'
matplotlib.rcParams['xtick.bottom'] = True
matplotlib.rcParams['xtick.top'] = False
matplotlib.rcParams['ytick.right'] = False
matplotlib.rcParams['axes.edgecolor'] = 'black'
matplotlib.rcParams['axes.linewidth'] = '1.0'
matplotlib.rcParams['axes.labelsize'] = 'medium'
matplotlib.rcParams['axes.grid'] = True
matplotlib.rcParams['grid.linewidth'] = '0.0'
matplotlib.rcParams['grid.alpha'] = '0.18'
matplotlib.rcParams['grid.color'] = 'lightgray'
matplotlib.rcParams['legend.labelspacing'] = 0.77
matplotlib.rcParams['savefig.bbox'] = 'tight'
matplotlib.rcParams['savefig.format'] = 'pdf'
### ------------------------------------------------
### ------------------------------------------------
### ------------------------------------------------
### ------------------------------------------------
### ------------------------------------------------
### ------------------------------------------------
analysissettings0={'ignore_rows': u'0.1'}

analysissettings={
  'smooth_scale_2D':0.5, 
  'smooth_scale_1D': 0.5,
  'fine_bins_2D': u'1000',
  'fine_bins': u'1000',
  'num_bins_2D': u'1000',
  'num_bins': u'1000',
  'range_confidence' : u'0.025',
  'ignore_rows': u'0.0'
}
### ------------------------------------------------
### ------------------------------------------------
### ------------------------------------------------
### ------------------------------------------------
### ------------------------------------------------
### ------------------------------------------------
zpivot=0.4
samples=loadMCSamples('./LSST_3x2pt_sys_model_5',settings=analysissettings0)
p = samples.getParams()
samples.addDerived(1e9*p.sigma8,name='AA',label='10^9 A_\mathrm{s}',
  range=[samples.getLower('sigma8')*1e9,samples.getUpper('sigma8')*1e9])
samples.addDerived(100*p.h0,name='H0',label='H_0',
  range=[samples.getLower('h0')*1e2,samples.getUpper('h0')*1e2])
samples.addDerived(10*p.ns,name='ns10',label='10 n_s',
  range=[10*samples.getLower('ns'),10*samples.getUpper('ns')])
samples.addDerived(100*p.omegab,name='omegab100',label='100 \Omega_b',
  range=[100*samples.getLower('omegab'),100*samples.getUpper('omegab')])
samples.addDerived(10*p.omegam,name='omegam10',label='10 \Omega_m',
  range=[10*samples.getLower('omegam'),10*samples.getUpper('omegam')])
samples.addDerived(p.w0+(zpivot/(1.+zpivot))*p.wa,name='wp',label='w_p',
  range=[samples.getLower('w0')+(zpivot/(1.+zpivot))*samples.getLower('wa'),samples.getUpper('w0')+(zpivot/(1.+zpivot))*samples.getUpper('wa')])
samples.saveAsText('.VM_plot_all_probes_multi_v2_TMP_1')

samples=loadMCSamples('./LSST_3x2pt_sys_model_6',settings=analysissettings0)
p = samples.getParams()
samples.addDerived(1e9*p.sigma8,name='AA',label='10^9 A_\mathrm{s}',
  range=[samples.getLower('sigma8')*1e9,samples.getUpper('sigma8')*1e9])
samples.addDerived(100*p.h0,name='H0',label='H_0',
  range=[samples.getLower('h0')*1e2,samples.getUpper('h0')*1e2])
samples.addDerived(10*p.ns,name='ns10',label='10 n_s',
  range=[10*samples.getLower('ns'),10*samples.getUpper('ns')])
samples.addDerived(100*p.omegab,name='omegab100',label='100 \Omega_b',
  range=[100*samples.getLower('omegab'),100*samples.getUpper('omegab')])
samples.addDerived(10*p.omegam,name='omegam10',label='10 \Omega_m',
  range=[10*samples.getLower('omegam'),10*samples.getUpper('omegam')])
samples.addDerived(p.w0+(zpivot/(1.+zpivot))*p.wa,name='wp',label='w_p',
  range=[samples.getLower('w0')+(zpivot/(1.+zpivot))*samples.getLower('wa'),samples.getUpper('w0')+(zpivot/(1.+zpivot))*samples.getUpper('wa')])
samples.saveAsText('.VM_plot_all_probes_multi_v2_TMP_2')

samples=loadMCSamples('./LSST_3x2pt_sys_model_7',settings=analysissettings0)
p = samples.getParams()
samples.addDerived(1e9*p.sigma8,name='AA',label='10^9 A_\mathrm{s}',
  range=[samples.getLower('sigma8')*1e9,samples.getUpper('sigma8')*1e9])
samples.addDerived(100*p.h0,name='H0',label='H_0',
  range=[samples.getLower('h0')*1e2,samples.getUpper('h0')*1e2])
samples.addDerived(10*p.ns,name='ns10',label='10 n_s',
  range=[10*samples.getLower('ns'),10*samples.getUpper('ns')])
samples.addDerived(100*p.omegab,name='omegab100',label='100 \Omega_b',
  range=[100*samples.getLower('omegab'),100*samples.getUpper('omegab')])
samples.addDerived(10*p.omegam,name='omegam10',label='10 \Omega_m',
  range=[10*samples.getLower('omegam'),10*samples.getUpper('omegam')])
samples.addDerived(p.w0+(zpivot/(1.+zpivot))*p.wa,name='wp',label='w_p',
  range=[samples.getLower('w0')+(zpivot/(1.+zpivot))*samples.getLower('wa'),samples.getUpper('w0')+(zpivot/(1.+zpivot))*samples.getUpper('wa')])
samples.saveAsText('.VM_plot_all_probes_multi_v2_TMP_3')


### ------------------------------------------------
### ------------------------------------------------
### ------------------------------------------------
### ------------------------------------------------
g=gplot.getSinglePlotter(
  chain_dir=r'./',
  analysis_settings=analysissettings,
  width_inch=4.5
)

g.settings.lw_contour = 1.0
g.settings.legend_rect_border = False
g.settings.figure_legend_frame = False
g.settings.axes_fontsize = 15.5
g.settings.legend_fontsize = 10.0
g.settings.alpha_filled_add = 0.7
g.settings.lab_fontsize=20
g.legend_labels=False

roots = [
  '.VM_plot_all_probes_multi_v2_TMP_1',
  '.VM_plot_all_probes_multi_v2_TMP_2',
  '.VM_plot_all_probes_multi_v2_TMP_3',  
]

params = [['wp','wa']]

g.plots_2d(
  roots, 
  param_pairs=params, 
  filled=True, 
  imax_shaded=2,
  line_args=[
    {'lw': 1.4,'ls': 'solid', 'color':'lightskyblue'},    
    {'lw': 1.0,'ls': 'solid', 'color':'maroon'},
    {'lw': 1.0,'ls': 'solid', 'color':'black'}
  ],
  contour_colors=[
    'lightskyblue',
    'maroon',
    'black',
  ]
  ,lims=[-2.,-0.,-2.5,2.5]
)
### ------------------------------------------------
### ------------------------------------------------
### ------------------------------------------------
# FINISH PLOT SETUP
g.settings.tight_layout=True
g.finish_plot(
  legend_labels=[
    'model 5',
    'model 6',
    'model 7',        
  ],
  no_extra_legend_space=True,
  legend_loc=(0.60,0.76)
)
### ------------------------------------------------
### ------------------------------------------------
### ------------------------------------------------
g.export()
### ------------------------------------------------
### ------------------------------------------------
### ------------------------------------------------
#DELETE TMP FILES
subprocess.Popen(
  "rm .VM_plot_all_probes_multi_v2_TMP_[0-9].*", 
  shell=True, 
  cwd="."
)
### ------------------------------------------------
### ------------------------------------------------
### ------------------------------------------------