pix2sky
=======

Python routine to convert CCD pixel coordinates to sky coordinates (RA/DEC) on multicore machines.

- Author:       Navtej Singh
- Contact:      n.saini1@nuigalway.ie
- Web Site:     http://astro.nuigalway.ie/staff/navtejs
- Organization: CfA@NUIG <http://astro.nuigalway.ie>

This routine was coded as part of research paper "Parallel astronomical data processing with Python: Recipes for multicore machines", published in Astronomy and Computing. Astro-ph link: http://arxiv.org/abs/1306.0573.

Thank you for downloading pix2sky code. It includes three different versions of the same code - pix2sky_serial.py for running code in serial mode, pix2sky_multi.py for running on multicore/multiprocessor machines and pix2sky_pp.py uses external parallel python library.


- Following requirements should be met to run the code.

    + A Python 2.4/2.5/2.6/2.7/3.0/3.1/3.2 distribution.

    + pyfits python module to handle FITS image files. Download from STScI's
      website (http://www.stsci.edu/resources/software_hardware/pyfits)
 
    + parallel python module in case of pix2sky_pp.py. It can be downloaded
      from www.parallelpython.com.

- Test data is included in data directory


- Execute the following commands to run the code

    + Serial Mode: 
              $python pix2sky_serial.py data/in.fits data/in_xy.cat
              $python pix2sky_serial.py data/in.fits data/in_xy.cat -d no  

    + Multicore Mode: 
              $python pix2sky_multi.py data/in.fits data/in_xy.cat
              $python pix2sky_multi.py data/in.fits data/in_xy.cat -d no

              $python pix2sky_pp.py data/in.fits data/in_xy.cat         
              $python pix2sky_pp.py data/in.fits data/in_xy.cat -d no
    
- For available command line options:

    + $python pix2sky_serial.py --help
    + $python pix2sky_multi.py --help
    + python pix2sky_pp.py --help   
