#!/usr/bin/python

# Plot given data (2 column) with given axes, optional logscale
# and range.

import sys
import getopt

# Usage/help message display
def usage():
    print sys.argv[0], " [options] data-file-name legend-key fields [data-file-name legend-key fields ...]"
    print "The available options are:"
    print "-x         X-axis label within quotes"
    print "-y         Y-axis label within quotes"
    print "--lx       Use log scale on X-axis"
    print "--ly       Use log scale on Y-axis"
    print "-t         Title of the graph in quotes"
    print "-f         Terminal format"
    print "-o         Output filename"
    print "--xrange   x1:x2 within quotes"
    print "--yrange   y1:y2 within quotes"
    print "--bmargin  bottom margin (in what units?)"
    print "--l        Plot using lines (default is linespoints)"
    print "-k         Key (legend) parameters"
    print "-p         Point size"
    print "-a         Tics"
    print "-h, --help Display this help message"

def main():
    # Get arguments
    args = sys.argv[1:]

    # get arguments
    try:
        options, real_args = getopt.gnu_getopt(args, "hlx:y:t:f:o:k:p:a:", ["lx", "ly", "xrange=", "yrange=", "bmargin=", "help"])
    except getopt.GetoptError, err:
        print str(err)
        usage()
        sys.exit(2)

    # Is the data file argument provided?
    if len(real_args) == 0:
        usage()
        sys.exit(0)
        
    if len(real_args) % 3 != 0:
        usage()
        sys.exit(0)

    lines = False # plotting using lines turned off by default. Usually, plot by linespoints.
    bmargin_set = False # for checking bottom margin setting
    for opt, val in options:
        if opt in ("-h", "--help"):
            usage()
            sys.exit(0)
        elif opt == "-x":
            print "set xlabel " + "\"" + val + "\" font \"Helvetica,28\""
        elif opt == "-y":
            print "set ylabel " + "\"" + val + "\" font \"Helvetica,28\""
        elif opt == "-t":
            print "set title  " + "\"" + val + "\""
        elif opt == "--lx":
            print "set logscale x"
        elif opt == "--ly":
            print "set logscale y"
        elif opt == "--xrange":
            print "set xrange " + val
        elif opt == "--yrange":
            print "set yrange " + val
        elif opt == "-l":
            lines = True
        elif opt == "-f":
            print "set terminal " + val
        elif opt == "-o":
            print "set output " + "\"" + val + "\""
        elif opt == "-k":
            print "set key " + val
        elif opt == "--bmargin":
            bmargin_set = True
            print "set bmargin " + val
        elif opt == "-p":
            print "set pointsize " + val
        elif opt == "-a":
            print "set tics " + val
            
    print "set lmargin 10"
    if bmargin_set == False:
        print "set bmargin 4"
    print "set tmargin 2"
    print "set rmargin 2.5"
            
    outstr = "plot "
    
    nsets = len(real_args)/3
            
    for i in range(0, nsets):

        filename, legend_key, fields = real_args[0:3]

        outstr += ("\"" + filename + "\" using " + fields + " with ")
        if lines:
            outstr += "lines"
        else:
            outstr += "linespoints"
        outstr += (" title \"" + legend_key + "\"")
        
        if i != (nsets - 1):
            outstr += ", "

        real_args = real_args[3:]

    print outstr

# Call main() if this is the main thread
if __name__ == "__main__":
    main()

