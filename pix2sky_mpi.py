'''
    MPI pix2sky routine called from pix2sky.py program.
'''

# Load python modules to be used in the routine
import sys, os
from StringIO import StringIO
from math import cos, sin, atan2, sqrt, radians, degrees


# Check if mpi4py module is available
try:
    from mpi4py import MPI 
except ImportError:
    raise ImportError("Module mpi4py not found.")


# Create chunks of data to be distributed to multiple processors based on file size
def getfilechunks(infile, n_cpus, scheduler = 'guided'):
    if scheduler == 'static':
        size = os.path.getsize(infile) / n_cpus
    else:
        size = os.path.getsize(infile) / (n_cpus * 20)


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
        

# Get header keywords and save in class variables
def getHeader(image):
    print >> sys.stdout, '\n Getting image header keywords...'

    # Check if pyraf module is installed
    try:
        import pyfits
    except:
        print >> sys.stderr, 'Error: Python module pyfits not found. Exiting.'
        sys.exit(-1)
    
    # Input can be a single FITS image, multi-extension image or
    # multi-extension image with particular extension
    if len(image.split('[', 1)) > 1:
        ext = image.split('[', 1 )[1].replace( ']', '')
        image = image.split('[', 1)[0]
    else:
        ext = ''
    
    # Open Header Unit List (HDU) to read header keywords
    try:
        hdulist = pyfits.open(image)
    except:
        print >> sys.stderr, 'Error: Not able to read FITS header. Exiting.'
        exit(-1)

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

    return ra0, dec0, crpix1, crpix2, cd11, cd12, cd21, cd22


# Convert RA from degree to hours:min:sec
def degree2hours(degrees):
    hour = int( degrees )
    tmp = (degrees - hour) * 60
    min = int(tmp)
    sec = (tmp - min) * 60 
    
    return '%2d%s%02d%s%2.4f' %(hour, ':', min, ':', sec)
    

# Translate X,Y pixel coordinates to sky coordinates RA,DEC
def translate(x, y, ra0, dec0, crpix1, crpix2, cd11, cd12, cd21, cd22, degree):
    # Formulas based on IRAF implementation of xy2rd task
    xi = cd11 * (x - crpix1) + cd12 * (y - crpix2)
    eta = cd21 * (x - crpix1) + cd22 * (y - crpix2)
    
    xi = radians(xi)
    eta = radians(eta)
    ra0 = radians(ra0)
    dec0 = radians(dec0)
    
    ra = atan2(xi, cos(dec0) - eta * sin(dec0)) + ra0
    dec = atan2(eta * cos(dec0) + sin(dec0), sqrt((cos(dec0) - eta * sin(dec0))**2 + xi**2))
    
    ra = degrees(ra)
    dec = degrees(dec)
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


# MPI worker function
def worker(indata):
    # Unpack input python list
    chunk0, chunk1, infile, outfile, ra0, dec0, crpix1, crpix2, cd11, cd12, cd21, cd22, degree = indata
    
    # Open input file
    try:
        ifile = open(infile, 'r')
    except:
        print >> sys.stderr, 'Error: Not able to open input file ', infile, '. Exiting.'
        sys.exit(-1)
        
    ifile.seek(chunk0)
    
    # Process data and save in output python list
    result = []
    for line in ifile.read(chunk1).splitlines():
    	if line[0] != '#':
    	    result.append(translate(float(line.split()[0]), float(line.split()[1]), ra0, dec0, crpix1, crpix2, cd11, cd12, cd21, cd22, degree))

    # Return list of calling program        
    return result

       
# pix2sky routine to proccess x,y pixel pairs
def pix2sky(image, infile, degree = 'yes', outfile = None, verbose = False):

    # Start mpi processing
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    n_cpus = comm.Get_size()

    # Divide the data in chunks equal to number of processors/cores
    if rank == 0:
        print >> sys.stdout, '\n Starting processing...'

	if not os.path.exists(image.split( '[', 1 )[0]):
	    print >> sys.stderr, 'Error: Image ', image, ' does not exist. Exiting.'
	    sys.exit(-1)

	if not os.path.exists(infile):
	    print >> sys.stderr, 'Error: Pixel file ', infile, ' does not exist. Exiting.'
	    sys.exit(-1)
	
	# Read image header keywords
	ra0, dec0, crpix1, crpix2, cd11, cd12, cd21, cd22 = getHeader(image)	

	# Divide input data in chunks to be processed by each processor/core
	chunks = getfilechunks(infile, n_cpus, scheduler = 'static')

        data = []
        for chunk in chunks:
            data.append((chunk[0], chunk[1], infile, outfile, ra0, dec0, crpix1, crpix2, cd11, cd12, cd21, cd22, degree))
    else:
        data = None

    
    # Check verbosity
    if not verbose:
        output = StringIO()
        old_stdout = sys.stdout
        sys.stdout = output


    # Scatter the data to available processors (collective communication) and gather the results
    indata = comm.scatter(data, root = 0)
    outdata = worker(indata)
    results = comm.gather(outdata, root = 0)


    if rank == 0:
        # Write headers to the output file
        if outfile:
            outfile = outfile
        else:
            outfile = os.path.join(infile.rsplit('/')[0], 'pix2sky.out')

	try:
	    ofile = open(outfile, 'w')
	except:
		print >> sys.stderr, 'Error: Not able to open file ', outfile, ' for writing. Exiting.'
		sys.exit(-1)

        ofile.write('# ---------------------------------------------------------\n')
        ofile.write('#    X		Y		RA	      	   DEC  \n')
        ofile.write('# ---------------------------------------------------------\n')

        # Write results to output file
        for result in results:
            for data in result:
                ofile.write('%10s%10s%18s%18s%s' %( str(data[0]), str(data[1]), str(data[2]), str(data[3]), '\n'))		    	

        # Close output file        
        try:        
            ofile.close()
        except:
            print >> sys.stderr, 'Warning: Not able to close output file ', outfile

	print >> sys.stdout, '\n Results written to - ', outfile

    # Reset verbosity
    if not verbose:
        sys.stdout = old_stdout
    
    if rank == 0:
        print >> sys.stdout, '\n Process completed successfully.'    