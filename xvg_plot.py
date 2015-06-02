#generic python modules
import argparse
import operator
from operator import itemgetter
import sys, os, shutil
import os.path

##########################################################################################
# RETRIEVE USER INPUTS
##########################################################################################

#=========================================================================================
# create parser
#=========================================================================================
version_nb = "0.0.2"
parser = argparse.ArgumentParser(prog = 'xvg_plot', usage='', add_help = False, formatter_class = argparse.RawDescriptionHelpFormatter, description =\
'''
**********************************************
v''' + version_nb + '''
author: Jean Helie (jean.helie@bioch.ox.ac.uk)
git: https://github.com/jhelie/xvg_plot
**********************************************

[ DESCRIPTION ]

This script simply plots the content of an xvg file. In case the --std flag is used
there must an even number of data column (the second half of columns is assumed to contain
std values for the first half).

[ USAGE ]

Option	      Default  	Description                    
-----------------------------------------------------
-f			: xvg file for energy evolution
-o			: output file
--comments	@,#	: lines starting with these characters will be considered as comment
--ymax			: upper boundary of y axis
--ymin			: lower boundary of y axis
--xmax			: lower boundary of x axis
--xmin			: lower boundary of x axis
--nbx			: nb ticks on x axis
--nby			: nb ticks on y axis
--hline			: plot an horizontal line at this y value (use commas to specify several values)
--vline			: plot a vertical line at this x value (use commas to specify several values)
--std			: if xvg contains std data plot it as area around average values

Other options
-----------------------------------------------------
--version		: show version number and exit
-h, --help		: show this menu and exit
 
''')

#options
parser.add_argument('-f', nargs=1, dest='xvgfilename', help=argparse.SUPPRESS, required=True)
parser.add_argument('-o', nargs=1, dest='output_file', default=["auto"], help=argparse.SUPPRESS)
parser.add_argument('--ymax', nargs=1, dest='ymax', default=[-10000000], type=float, help=argparse.SUPPRESS)
parser.add_argument('--ymin', nargs=1, dest='ymin', default=[-10000000], type=float, help=argparse.SUPPRESS)
parser.add_argument('--xmax', nargs=1, dest='xmax', default=[-10000000], type=float, help=argparse.SUPPRESS)
parser.add_argument('--xmin', nargs=1, dest='xmin', default=[-10000000], type=float, help=argparse.SUPPRESS)
parser.add_argument('--nbx', nargs=1, dest='nbx', default=[-1], type=int, help=argparse.SUPPRESS)
parser.add_argument('--nby', nargs=1, dest='nby', default=[-1], type=int, help=argparse.SUPPRESS)
parser.add_argument('--hline', nargs=1, dest='hline', default=[-10000000], help=argparse.SUPPRESS)
parser.add_argument('--vline', nargs=1, dest='vline', default=[-10000000], help=argparse.SUPPRESS)
parser.add_argument('--std', dest='std', action='store_true', help=argparse.SUPPRESS)
parser.add_argument('--comments', nargs=1, dest='comments', default=['@,#'], help=argparse.SUPPRESS)

#other options
parser.add_argument('--version', action='version', version='%(prog)s v' + version_nb, help=argparse.SUPPRESS)
parser.add_argument('-h','--help', action='help', help=argparse.SUPPRESS)

#=========================================================================================
# store inputs
#=========================================================================================

args = parser.parse_args()
args.xvgfilename = args.xvgfilename[0]
args.output_file = args.output_file[0]
args.ymax = args.ymax[0]
args.ymin = args.ymin[0]
args.xmax = args.xmax[0]
args.xmin = args.xmin[0]
args.nbx = args.nbx[0]
args.nby = args.nby[0]
args.hline = args.hline[0]
args.vline = args.vline[0]
args.comments = args.comments[0].split(',')
if args.output_file == "auto":
	args.output_file = args.xvgfilename[:-4] + '.svg'

#=========================================================================================
# import modules (doing it now otherwise might crash before we can display the help menu!)
#=========================================================================================

#generic science modules
try:
	import numpy as np
except:
	print "Error: you need to install the numpy module."
	sys.exit(1)
try:
	import matplotlib as mpl
	mpl.use('Agg')
	import matplotlib.colors as mcolors
	mcolorconv = mcolors.ColorConverter()
	import matplotlib.cm as cm				#colours library
	import matplotlib.ticker
	from matplotlib.ticker import MaxNLocator
	from matplotlib.font_manager import FontProperties
	fontP=FontProperties()
	from matplotlib.collections import LineCollection

except:
	print "Error: you need to install the matplotlib module."
	sys.exit(1)
try:
	import pylab as plt
except:
	print "Error: you need to install the pylab module."
	sys.exit(1)

#=======================================================================
# sanity check
#=======================================================================

if not os.path.isfile(args.xvgfilename):
	print "Error: file " + str(args.xvgfilename) + " not found."
	sys.exit(1)

if os.path.isfile(args.output_file):
	print "Error: file " + str(args.output_file) + " already exists, specify a different name via -o."
	sys.exit(1)

try:
	args.hline = [float(args.hline)]
except:
	tmp_values = args.hline.split(',')
	args.hline = [float(h) for h in tmp_values]

try:
	args.vline = [float(args.vline)]
except:
	tmp_values = args.vline.split(',')
	args.vline = [float(v) for v in tmp_values]

##########################################################################################
# MAIN
##########################################################################################

#get file content
with open(args.xvgfilename) as f:
	lines = f.readlines()

#determine legends and nb of lines to skip
x_label = "x axis"
y_label = "y axis"
fig_title = ""
c_index = 1
c_labels = {}
tmp_nb_rows_to_skip = 0
for l_index in range(0,len(lines)):		
	line = lines[l_index]
	if line[0] in args.comments:
		tmp_nb_rows_to_skip += 1
		#retrieve title
		if "@ title" in line:
			fig_title = line.split("@ title \"")[-1][:-2]
		#retrieve x axis label
		if "@ xaxis label" in line:
			x_label = line.split("@ xaxis label \"")[-1][:-2]
		#retrieve y axis label
		if "@ yaxis label" in line:
			y_label = line.split("@ yaxis label \"")[-1][:-2]
		#retrieve captions
		if "legend \"" in line:
			c_labels[c_index] = line.split("legend \"")[-1][:-2]
			c_index += 1

#get data
data = np.loadtxt(args.xvgfilename, skiprows = tmp_nb_rows_to_skip)

#create figure
fig, ax = plt.subplots()
fig.suptitle(fig_title)
	
#get number of data columns
nb_col_data = int(np.shape(data)[1] - 1)
if args.std and nb_col_data % 2 != 0:
	print "Error: --std specified but " + str(args.xvgfilename) + " contains an odd number of data column, see xvg_plot -h"
	sys.exit(1)

#plot data
if args.std:
	for c_index in range(1, int(nb_col_data/float(2) + 1)):
		base_line, = plt.plot(data[:,0], data[:,c_index], label = c_labels[c_index])
		plt.fill_between(data[:,0], data[:,c_index] - data[:,c_index + int(nb_col_data/float(2))], data[:,c_index] + data[:,c_index + int(nb_col_data/float(2))], color = base_line.get_color(), edgecolor = base_line.get_color(), linewidth = 0, alpha = 0.2)
else:
	for c_index in range(1, np.shape(data)[1]):
		plt.plot(data[:,0], data[:,c_index], label = c_labels[c_index])

#format axes and legend
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
if args.xmax != -10000000:
	ax.set_xlim(xmax=args.xmax)
if args.xmin != -10000000:
	ax.set_xlim(xmin=args.xmin)
if args.ymax != -10000000:
	ax.set_ylim(ymax=args.ymax)
if args.ymin != -10000000:
	ax.set_ylim(ymin=args.ymin)
if args.nbx != -1:
	ax.xaxis.set_major_locator(MaxNLocator(nbins=args.nbx))
if args.nby != -1:
	ax.yaxis.set_major_locator(MaxNLocator(nbins=args.nby))
if args.vline != [-10000000]:
	ymin, ymax = ax.get_ylim()
	for v in args.vline:
		plt.vlines(v, ymin, ymax, linestyles = 'dashed')
if args.hline != [-10000000]:
	xmin, xmax = ax.get_xlim()
	for h in args.hline:
		plt.hlines(h, xmin, xmax, linestyles = 'dashed')

plt.setp(ax.xaxis.get_majorticklabels(), fontsize="small" )
plt.setp(ax.yaxis.get_majorticklabels(), fontsize="small" )
fontP.set_size("small")
ax.legend(prop=fontP)
plt.xlabel(x_label, fontsize="small")
plt.ylabel(y_label, fontsize="small")
fig.savefig(args.output_file)
plt.close()


#=========================================================================================
# exit
#=========================================================================================
print "\nFinished successfully! Check result in file " + str(args.output_file)
print ""
sys.exit(0)
