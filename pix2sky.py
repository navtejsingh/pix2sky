#!/usr/bin/env python

"""
    ------------------------------------------------------------------------
    Routine to translate X and Y pixel coordinates to RA, DEC values in
    multiprocessing MPI mode

    Usage: mpiexec -np ncpus python pix2sky.py [options] image pixfile

    Input:
        ncpus: Number of physical cores or processors
        image: Input image with basic header keywords
        pixfile: Input file with (x, y) pixel value in column 1 and 2
      
      
    [Options]:
        --help: help
        --version: program version
        --verbose: show result messages
        --quiet: don't show result messages
        --filename: output file name (default is pix2sky.dat)
        --degree: (ra, dec) in degrees? (default is yes)        
        
    Output:
        pix2sky.dat: Output file with (x, y) and corresponding (ra, dec)
        
    Author:
        Navtej Singh

    Organization:
        Centre for Astronomy, National University of Ireland, Galway, Ireland
        
    Version:
        20 December 2011    1.0   Original version
    ------------------------------------------------------------------------    
"""

from optparse import OptionParser
import pix2sky_mpi as pm
    
    
# Entry point for the routine
if __name__ == '__main__':
    usage = "Usage: mpiexec -np ncpus python %prog [options] image pixfile"
    description = "Description. Utility to convert X/Y pixel image coordinates to RA/DEC sky coordinates in multiprocessing mode."
    parser = OptionParser(usage = usage, version = "%prog 1.0", description = description)
    parser.add_option("-v", "--verbose",
                      action="store_true", dest="verbose", default = False,
                      help = "print result messages to stdout"
                      )
    parser.add_option("-q", "--quiet",
                    action="store_false", dest="verbose", default = True,
                    help = "don't print result messages to stdout"
                    )
    parser.add_option("-d", "--degree", dest = "degree", metavar="DEGREE",
                    action="store", help = "ra/dec in degree? [default is yes]",
                    choices=['yes', 'no'], default = 'yes'
                    )    
    parser.add_option("-f", "--filename", dest = "filename",
                    action='store', metavar="FILE", help = "output file name [default is pix2sky.out]"
                    )
    (options, args) = parser.parse_args()
    
    # Check for number of input arguments
    if len(args) != 2:
        parser.error("Incorrect number of arguments")

    pm.pix2sky(args[0], args[1], options.degree, options.filename, options.verbose)
    
