#!/usr/bin/env python

"""
    ------------------------------------------------------------------------
    Routine to translate X and Y pixel coordinates to RA, DEC values in
    serial mode.

    Usage: python pix2sky_serial.py [options] image pixfile

    Input:
        image: Input image with basic header keywords
        pixfile: Input file with (x, y) pixel value in column 1 and 2
      
      
    [Options]:
        --help: help
        --version: program version
        --verbose: show result messages
        --quiet: don't show result messages
        --filename: output file name (default is pix2sky.dat)
        --degree: (ra, dec) in degrees? (default is yes)
        
        
    Author:
        Navtej Singh

    Organization:
        Centre for Astronomy, National University of Ireland, Galway, Ireland

    Version:
        15 December 2011     1.0     Original version  
    ------------------------------------------------------------------------    
"""


# Load python modules to be used in the routine
import sys, math
from os.path import join, exists
from StringIO import StringIO
from optparse import OptionParser


# Get header keywords and save in class variables
# ===============================================
def getHeader(image):
    print >> sys.stdout, '\n Getting image header keywords...'
    
    # Input can be a single FITS image, multi-extension image or
    # multi-extension image with particular extension
    if len(image.split( '[', 1 )) > 1:
        ext = image.split('[', 1)[1].replace(']', '')
        image = image.split('[', 1)[0]
    else:
        ext = ''
    
    # Open Header Unit List (HDU) to read header keywords
    try:
        hdulist = pyfits.open(image)
    except:
        print >> sys.stderr, 'Error: Not able to read FITS header. Exiting.'
        sys.exit(-1)

    hdulist.close()

    # Get header parameters - checking number of extensions and using 1st extension
    # in case of multi extension FITS image
    if len( hdulist ) > 1:
        if ext == '':
            hdrdata = hdulist[1].header
        else:
            
            hdrdata = hdulist[int(ext)].header 
    else:
        hdrdata = hdulist[0].header


    # Get CRPIX keyword values
    crpix1 = hdrdata['CRPIX1']
    crpix2 = hdrdata['CRPIX2']

    # Get CRVAL keyword values
    ra0 = hdrdata['CRVAL1']
    dec0 = hdrdata['CRVAL2']

    # Get CD keyword values
    cd11 = hdrdata['CD1_1']
    cd12 = hdrdata['CD1_2']
    cd21 = hdrdata['CD2_1']
    cd22 = hdrdata['CD2_2']


    # Return image header keywords
    return ra0, dec0, crpix1, crpix2, cd11, cd12, cd21, cd22


# Convert RA from degree to hour:min:sec
# ======================================
def degree2hours(degrees):
    hour = int(degrees)
    tmp = (degrees - hour) * 60
    min = int(tmp)
    sec = (tmp - min) * 60 
    
    return '%2d%s%02d%s%2.4f' %(hour, ':', min, ':', sec)


# Translate X,Y image pixel coordinates to sky coordinates RA,DEC
# ===============================================================
def translate(x, y, ra0, dec0, crpix1, crpix2, cd11, cd12, cd21, cd22, degree):
    # Formulas based on IRAF implementation of xy2rd task
    xi = cd11 * (x - crpix1) + cd12 * (y - crpix2)
    eta = cd21 * (x - crpix1) + cd22 * (y - crpix2)
    
    xi = math.radians(xi)
    eta = math.radians(eta)
    ra0 = math.radians(ra0)
    dec0 = math.radians(dec0)
    
    ra = math.atan2(xi, math.cos(dec0) - eta * math.sin(dec0)) + ra0
    dec = math.atan2(eta * math.cos(dec0) + math.sin(dec0), math.sqrt((math.cos(dec0) - eta * math.sin(dec0))**2 + xi**2))
    
    ra = math.degrees(ra)
    dec = math.degrees(dec)
    ra = ra % 360.0
    if ra < 0.0:
        ra = ra + 360.0
        
    if degree == 'no':
        ra = ra / 15.0
        ra = degree2hours(ra)
        dec = degree2hours(dec)
        print >> sys.stdout, 'X = %6.3f%s' %(x, '\t'), '  Y = %6.3f' %y, '  RA = %s' %ra, '  DEC = %s' %dec
    else:    
        print >> sys.stdout, 'X = %6.3f%s' %(x, '\t'), '  Y = %6.3f' %y, '  RA = %3.9f' %ra, '  DEC = %3.9f' %dec
    
    return (x, y, ra, dec)
    
    
# pix2sky routine to proccess x,y pixel pairs
# ===========================================
def pix2sky(image, infile, degree = 'yes', outfile = None):
    # Read image header keywords
    ra0, dec0, crpix1, crpix2, cd11, cd12, cd21, cd22 = getHeader(image)  

    # Open input file with x and y pixel values
    try:
        ifile = open(infile, 'r')
    except:
        print >> sys.stderr, 'Error: Not able to open the input file ', infile, '. Exiting.'
        sys.exit(-1)

    # Set the output file name
    if not outfile:
        if len(infile.rsplit('/')) > 1:
            outfile = join(infile.rsplit('/')[0], 'pix2sky.out')
        else:
            outfile = 'pix2sky.out'
     
    # Open the output file     
    try:
        ofile = open(outfile, 'w')
    except:
        print >> sys.stderr, 'Error: Not able to open the output file ', outfile, ' for writing. Exiting.'
        sys.exit(-1)


    # Write data headers to the output file    
    ofile.write('# ---------------------------------------------------------\n')
    ofile.write('#    X        Y            RA              DEC             \n')
    ofile.write('# ---------------------------------------------------------\n')


    # Write results to the output file
    while 1:
        line = ifile.readline()
        if not line:
            break
        if line[0] != '#':
            res = translate(float(line.split()[0] ), float( line.split()[1] ), ra0, dec0, crpix1, crpix2, cd11, cd12, cd21, cd22, degree)
            ofile.write('%10s%10s%18s%18s%s' %(str(res[0]), str(res[1]), str(res[2]), str(res[3]), '\n'))

    # Close input and output files
    try:        
        ifile.close()
        ofile.close()            
    except:
        print >> sys.stderr, 'Warning: Not able to close the input and the output files.'
        
    print >> sys.stdout, '\n Results written to - ', outfile


# Main function - doing some data validation before calling pix2sky method
# ========================================================================
def main(image, pixelfile, degree = 'yes', outfile = None):
    if not exists(image.split( '[', 1 )[0]):
        print >> sys.stderr, 'Error: Image ', image, ' does not exist. Exiting.'
        sys.exit(-1)

    if not exists(pixelfile):
        print >> sys.stderr, 'Error: Pixel file ', pixelfile, ' does not exist. Exiting.'
        sys.exit(-1)    

    pix2sky(image, pixelfile, degree, outfile)
    
    
# Entry point for PIX2SKY_SERIAL utility
# ======================================
if __name__ == '__main__':
    usage = "Usage: %prog [options] image pixfile"
    description = "Description. Utility to convert X/Y pixel image coordinates to RA/DEC sky coordinates in serial mode."
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

    print >> sys.stdout, '\n Starting processing...'
    
    # Check verbosity
    if not options.verbose:
        output = StringIO()
        old_stdout = sys.stdout
        sys.stdout = output

    # Check if pyraf module is available
    try:
        import pyfits
    except:
        print >> sys.stderr, 'Error: Python module pyfits not found. Exiting.'
        sys.exit(-1)
    
    main(args[0], args[1], options.degree, options.filename)
    
    # Reset verbosity
    if not options.verbose:
        sys.stdout = old_stdout
    
    print >> sys.stdout, '\n Process completed successfully.'   
