#!/usr/bin/env python

"""
    ------------------------------------------------------------------------
    Routine to translate X and Y pixel coordinates to RA, DEC values using
    parallel python module.

    Usage: python pix2sky_pp.py [options] image pixfile

    Input:
        image: Input image with basic header keywords (CRVAL1, CRVAL2, CRPIX1, CRPIX2)
        pixfile: Input file with (x, y) pixel value in column 1 and 2
      
      
    [Options]:
        --help: help
        --version: program version
        --verbose: show result messages
        --quiet: don't show messages
        --filename: output file name (default is pix2sky.dat)
        --degree: (ra, dec) in degrees? (default is yes) else in hh:mm:ss/degree:arcmin:arcsec format
        --ncpus: Number of processors to use [default is maximum number of cores on the machine]
        --scheduler: Type of scheduler [static, guided, dynamic]; default is guided
        
        
    Output:
        pix2sky.dat: Output file with (x, y) and corresponding (ra, dec)
        

    Author:
        Navtej Singh

    Organization:
        Centre for Astronomy, National University of Ireland, Galway, Ireland

    Version:
        24 October 2012     1.0     Parallel python implementation of pix2sky 
    --------------------------------------------------------------------------    
"""


# Load python modules to be used in the routine
import sys, math
from os.path import getsize, join, exists
from StringIO import StringIO
from optparse import OptionParser


# Check if parallel python module is available
try:
    import pp
except:
    print >> sys.stderr, 'Error: Parallel python module not found. Use serial version of this routine. Exiting.'
    sys.exit(-1)


# Create chunks of data to be distributed to multiple processors/cores 
# based on the file size. Arbitrary factor of 20 is choosen    
# ====================================================================
def getchunks(infile, n_cpus, scheduler = 'guided'):
    # Divide input data based on scheduler type
    if scheduler == 'static':
        size = getsize(infile) / n_cpus
    else:
        size = getsize(infile) / (n_cpus * 20)

    # Open input file    
    try:    
        ifile = open(infile)
    except:
        print >> sys.stderr, 'Error: Not able to open ', infile, '. Exiting.'
        sys.exit(-1)

    # Create chunk of data to be distributed to nodes    
    while 1:
        start = ifile.tell()
        ifile.seek(size, 1)
        s = ifile.readline()
        yield start, ifile.tell() - start
        if not s:
            break

    # Close the input file    
    try:    
        ifile.close()
    except:
        print >> sys.stderr, 'Warning: Error closing the file ', ifile
        
    

# Get header keywords of the input image and save in class variables
# ==================================================================
def getHeader(image):
    print >> sys.stdout, '\n Getting image header keywords...'

    # Input can be a single FITS image, multi-extension image or
    # multi-extension image with a particular extension
    if len(image.split( '[', 1 )) > 1:
        ext = image.split( '[', 1 )[1].replace(']', '')
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
    if len(hdulist) > 1:
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

    # Return value to calling program
    return ra0, dec0, crpix1, crpix2, cd11, cd12, cd21, cd22


# Convert RA from degrees to hour:min:sec
# =======================================
def degree2hours(degrees):
    hour = int(degrees)
    tmp = (degrees - hour) * 60
    min = int(tmp)
    sec = (tmp - min) * 60 
    
    return '%2d%s%02d%s%2.4f' %(hour, ':', min, ':', sec)
    

# Translate X,Y pixels to the sky coordinates (RA, DEC)
# =====================================================
def translate(x, y, ra0, dec0, crpix1, crpix2, cd11, cd12, cd21, cd22, degree):
    # Formulas based on the IRAF implementation of xy2rd task
    xi = cd11 * (x - crpix1) + cd12 *  (y - crpix2)
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


# Parallel Python worker function
# ===============================
def worker(indata):
    # Unpack the input python list
    chunk0, chunk1, infile, ra0, dec0, crpix1, crpix2, cd11, cd12, cd21, cd22, hms = indata

    # Open the input file       
    try:
        ifile = open(infile, 'r')
    except:
        print >> sys.stderr, 'Error: Not able to open input file ', infile, '. Exiting.'
        sys.exit(-1)

    # Translate data from pixel to sky coordinates     
    ifile.seek(chunk0)
    result = []
    for line in ifile.read(chunk1).splitlines():
        if line[0] != '#':
            result.append(translate(float(line.split()[0]), float(line.split()[1]), ra0, dec0, crpix1, crpix2, cd11, cd12, cd21, cd22, hms))

    # Return result back to the calling method      
    return result

       
# pix2sky routine to process x,y pixel pairs
# ==========================================
def pix2sky(image, infile, degree = 'yes', outfile = None, ncpus = None, scheduler = 'guided'):
    # Read image headers
    ra0, dec0, crpix1, crpix2, cd11, cd12, cd21, cd22 = getHeader(image)  

    # Set number of server worker processes based either on input
    # parameter or default to number of cores on the machine        
    ppservers = ()
    if ncpus:
        # Creates jobserver with ncpus workers
        job_server = pp.Server(int(ncpus), ppservers = ppservers)
    else:
        # Creates jobserver with automatically detected number of workers
        job_server = pp.Server(ppservers = ppservers)
    
    # Divide data in smaller chunks for better performance based on scheduler type
    chunks = getchunks(infile, job_server.get_ncpus(), scheduler)
    
    # Start the worker processes in parallel
    jobs = []
    for value in chunks:
            indata = (value[0], value[1], infile, ra0, dec0, crpix1, crpix2, cd11, cd12, cd21, cd22, degree)
            jobs.append(job_server.submit(worker, (indata,), (translate,degree2hours,), ("math","sys",)))

    # Append the results   
    results = []
    for job in jobs:
        results.append(job())

    # Set the output file name
    if not outfile:
        if len(infile.rsplit('/')) > 1:
            outfile = join(infile.rsplit('/')[0], 'pix2sky.out')
        else:
            outfile = 'pix2sky.out'

    # Open the output file to write the results   
    try:
        ofile = open(outfile, 'w')
    except:
        print >> sys.stderr, 'Error: Not able to open the outfile file ', outfile, '. Exiting.'
        sys.exit(-1)


    # Format and write records to output file    
    ofile.write('# ---------------------------------------------------------\n')
    ofile.write('#    X     Y       RA             DEC  \n')
    ofile.write('# ---------------------------------------------------------\n')
        
    for result in results:
        for value in result:
            ofile.write('%10s%10s%18s%18s%s' %(str(value[0]), str(value[1]), str(value[2]), str(value[3]), '\n'))

    # Close the output file        
    try:        
        ofile.close()
    except:
        print >> sys.stderr, 'Warning: Not able to close the output file ', outfile
    
    print >> sys.stdout, '\n Results written to file - ', outfile



# Main function - doing some data validation before calling pix2sky method
# ========================================================================
def main(image, pixelfile, degree = 'yes', outfile = None, ncpus = None, scheduler = 'guided'):
    # Check if the image exists 
    if not exists(image.split( '[', 1 )[0]):
        print >> sys.stderr, 'Error: Image ', image, ' does not exist. Exiting.'
        sys.exit(-1)

    # Check if the input pixel file exists 
    if not exists(pixelfile):
        print >> sys.stderr, 'Error: Pixel file ', pixelfile, ' does not exist. Exiting.'
        sys.exit(-1)  

    # Execute the pix2sky task
    pix2sky(image, pixelfile, degree, outfile, ncpus, scheduler)
    

    
# Entry point for PIX2SKY utility
# ===============================
if __name__ == '__main__':
    usage = "Usage: %prog [options] image pixfile"
    description = "Description. Utility to convert X/Y pixel image coordinates to RA/DEC sky coordinates in multiprocessing mode."
    parser = OptionParser( usage = usage, version = "%prog 1.0", description = description )
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
    parser.add_option("-n", "--ncpus", dest = "ncpus", metavar="NCPUS",
                    action="store", help = "number of cpus (cores) for processing [default is maximum cores on the machine]"
                    )
    parser.add_option("-s", "--scheduler", dest = "scheduler", metavar="SCHEDULER",
                    action="store", help = "scheduler for multiprocessing [default is guided]",
                    choices=['guided', 'static'], default = 'guided'
                    )
    (options, args) = parser.parse_args()
    
    # Check for number of input arguments
    if len( args ) != 2:
        parser.error("Error: Incorrect number of arguments.")

    print >> sys.stdout, '\n Starting processing...'
    
    # Check verbosity
    if not options.verbose:
        output = StringIO()
        old_stdout = sys.stdout
        sys.stdout = output

    # Check if pyfits module is available
    try:
        import pyfits
    except:
        print >> sys.stderr, 'Error: Python module pyfits not found. Exiting.'
        sys.exit(-1)
    
    main(args[0], args[1], options.degree, options.filename, options.ncpus, options.scheduler)
    
    # Reset verbosity
    if not options.verbose:
        sys.stdout = old_stdout
    
    print >> sys.stdout, '\n Process completed successfully.'
