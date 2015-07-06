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
version_nb = "0.0.3"
parser = argparse.ArgumentParser(prog = 'xvg_plot', usage='', add_help = False, formatter_class = argparse.RawDescriptionHelpFormatter, description =\
'''
**********************************************
v''' + version_nb + '''
author: Jean Helie (jean.helie@bioch.ox.ac.uk)
git: https://github.com/jhelie/xvg_plot
**********************************************

DESCRIPTION:

 This script simply plots the content of an xvg file. In case the --std flag is used
 there must an even number of data column (the second half of columns is assumed to contain
 std values for the first half).

 A second axis can be specified but this is not compatible with the --std option (for now).
 When using --col2 indices should be 0-based and ignore the 1st column (ie follow the xvg legend
 numbering system). Lines belonging to the second axis are in bold.

USAGE:

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

2nd axis parameters
-----------------------------------------------------
--col2		: indices of columns to plot on second axis (format: i,j,...)
--nbx2			: nb ticks on x axis
--nby2			: nb ticks on y axis
--ymax2			: upper boundary of y axis
--ymin2			: lower boundary of y axis
--xmax2			: lower boundary of x axis
--xmin2			: lower boundary of x axis

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
parser.add_argument('--col2', nargs=1, dest='col2', default=['none'], help=argparse.SUPPRESS)
parser.add_argument('--ymax2', nargs=1, dest='ymax2', default=[-10000000], type=float, help=argparse.SUPPRESS)
parser.add_argument('--ymin2', nargs=1, dest='ymin2', default=[-10000000], type=float, help=argparse.SUPPRESS)
parser.add_argument('--xmax2', nargs=1, dest='xmax2', default=[-10000000], type=float, help=argparse.SUPPRESS)
parser.add_argument('--xmin2', nargs=1, dest='xmin2', default=[-10000000], type=float, help=argparse.SUPPRESS)
parser.add_argument('--nbx2', nargs=1, dest='nbx2', default=[-1], type=int, help=argparse.SUPPRESS)
parser.add_argument('--nby2', nargs=1, dest='nby2', default=[-1], type=int, help=argparse.SUPPRESS)
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

args.col2 = args.col2[0]
args.ymax2 = args.ymax2[0]
args.ymin2 = args.ymin2[0]
args.xmax2 = args.xmax2[0]
args.xmin2 = args.xmin2[0]
args.nbx2 = args.nbx2[0]
args.nby2 = args.nby2[0]

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

if args.col2 != "none":
	try:
		global col2_indices
		col2_indices = []
		content = args.col2.split(",")
		for c in content:
			col2_indices.append(int(int(c) + 1))
	except:
		print "Error: wrong format for --col2, see --help."
		sys.exit(1)


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

#display info about second axis
if args.col2 != "none":
	print "\nThe following columns will be plotted on the second axis:"
	for c_index in col2_indices:
		try:
			print " " + str(c_labels[c_index])
		except:
			print "Error: column " + str(c_index - 1) + " not found (index should be 0 based, see --help)."
			sys.exit(1)

#get data
data = np.loadtxt(args.xvgfilename, skiprows = tmp_nb_rows_to_skip)

#create figure
fig, ax1 = plt.subplots()
fig.suptitle(fig_title)
	
#create secondary axis
ax2 = ax1.twinx()

#get number of data columns
nb_col_data = int(np.shape(data)[1] - 1)
if args.std and nb_col_data % 2 != 0:
	print "Error: --std specified but " + str(args.xvgfilename) + " contains an odd number of data column, see xvg_plot -h"
	sys.exit(1)

#plot data
if args.std:
	for c_index in range(1, int(nb_col_data/float(2) + 1)):
		base_line, = ax1.plot(data[:,0], data[:,c_index], label = c_labels[c_index])
		ax1.fill_between(data[:,0], data[:,c_index] - data[:,c_index + int(nb_col_data/float(2))], data[:,c_index] + data[:,c_index + int(nb_col_data/float(2))], color = base_line.get_color(), edgecolor = base_line.get_color(), linewidth = 0, alpha = 0.2)
else:
	for c_index in range(1, np.shape(data)[1]):
		if c_index in col2_indices:
			ax2.plot(data[:,0], data[:,c_index], label = c_labels[c_index], linewidth = 3)
		else:
			ax1.plot(data[:,0], data[:,c_index], label = c_labels[c_index])

#format axes and legend: ax1
ax1.spines['top'].set_visible(False)
if args.col2 == "none":
	ax1.spines['right'].set_visible(False)
ax1.xaxis.set_ticks_position('bottom')
ax1.yaxis.set_ticks_position('left')
if args.xmax != -10000000:
	ax1.set_xlim(xmax=args.xmax)
if args.xmin != -10000000:
	ax1.set_xlim(xmin=args.xmin)
if args.ymax != -10000000:
	ax1.set_ylim(ymax=args.ymax)
if args.ymin != -10000000:
	ax1.set_ylim(ymin=args.ymin)
if args.nbx != -1:
	ax1.xaxis.set_major_locator(MaxNLocator(nbins=args.nbx))
if args.nby != -1:
	ax1.yaxis.set_major_locator(MaxNLocator(nbins=args.nby))
if args.vline != [-10000000]:
	ymin, ymax = ax1.get_ylim()
	for v in args.vline:
		plt.vlines(v, ymin, ymax, linestyles = 'dashed')
if args.hline != [-10000000]:
	xmin, xmax = ax1.get_xlim()
	for h in args.hline:
		plt.hlines(h, xmin, xmax, linestyles = 'dashed')
plt.setp(ax1.xaxis.get_majorticklabels(), fontsize="small" )
plt.setp(ax1.yaxis.get_majorticklabels(), fontsize="small" )
fontP.set_size("small")
if args.col2 == "none":
	ax1.legend(prop=fontP)
else:
	ax1.legend(prop=fontP, loc = 2)
plt.xlabel(x_label, fontsize="small")
ax1.set_ylabel(y_label, fontsize="small")

#format axes and legend: ax2
if args.col2 != "none":
	ax2.spines['top'].set_visible(False)
	ax2.spines['left'].set_visible(False)
	ax2.yaxis.set_ticks_position('right')
	if args.xmax2 != -10000000:
		ax2.set_xlim(xmax=args.xmax2)
	if args.xmin2 != -10000000:
		ax2.set_xlim(xmin=args.xmin2)
	if args.ymax2 != -10000000:
		ax2.set_ylim(ymax=args.ymax2)
	if args.ymin2 != -10000000:
		ax2.set_ylim(ymin=args.ymin2)
	if args.nbx2 != -1:
		ax2.xaxis.set_major_locator(MaxNLocator(nbins=args.nbx2))
	if args.nby2 != -1:
		ax2.yaxis.set_major_locator(MaxNLocator(nbins=args.nby2))
	plt.setp(ax2.yaxis.get_majorticklabels(), fontsize="small", weight = "bold")
	ax2.legend(prop=fontP, loc = 1)

#save file
fig.savefig(args.output_file)
plt.close()


#=========================================================================================
# exit
#=========================================================================================
print "\nFinished successfully! Check result in file " + str(args.output_file)
print ""
sys.exit(0)
