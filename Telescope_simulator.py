import numpy as np 
import matplotlib.pyplot as plt
import MegaScreen as M
from scipy.fftpack import fft
import imageio
import os
from PIL import Image
from scipy import signal
from scipy import interpolate

def Sum2d(a):
    """"Sum over the rightmost two dimensions of an array of dimension >= 2"""
    return(np.sum(np.sum(a,-1),-1))

def fried_parameter_cm(wavelength,arcseconds_of_seeing_500nm=1.,zenith_angle_deg = 0.):
    """return Fried parameter r0 in cm calculation taken from DFB Practical optical interferometry p58. 
    Find constant of proportionality k using given r0=10cm at 500nm for 1 arcsecond seeing @500nm"""
    r0_500nm_cm = (500e-9/(arcseconds_of_seeing_500nm*(np.pi/(180*3600))))*100
    k = r0_500nm_cm/(500e-9)**(6./5)
    r00 = k*wavelength**(6./5.)
    zenith_angle_rad = np.radians(zenith_angle_deg)
    r0z = r00 * np.cos(zenith_angle_rad)**(3/5.) #p60 DFB POI
    return r0z

def Circ_Aperture_Mask(propagator_size):
    """
    Returns a square 2D array with dimensions propagator_size, and a circular mask with diameter the
    length of the input array
    """

    x = np.arange(0, propagator_size)
    y = np.arange(0, propagator_size)
    arr = np.zeros((y.size, x.size))
    diam = propagator_size
    r= diam/2

    # The two lines below could be merged, but I stored the mask
    # for code clarity.
    mask = (x[np.newaxis,:]-propagator_size/2)**2 + (y[:,np.newaxis]-propagator_size/2)**2 < r**2
    arr[mask] = 1.

    return arr

def quick_complex_pupil(phaseScreen, array_to_propgate_size):
    """
    Takes in the 2D phase screen, converts it to a complex amplitude and returns it
    """

    propagator_size = len(phaseScreen) # gets the dimension of the square phase screen grid 
    telescope_wavefront = np.exp(1j*phaseScreen)
    initial_pupil = Circ_Aperture_Mask(propagator_size)*telescope_wavefront #complex amplitude with A = 1 inside UT aperture, 0 outside. Note, if "telescope_diameter_m" in "temporalPhaseScreens" is the size of the telescope aperture "mask_size_fraction" here should be = 1

    #DM padding 
    pad_left, pad_right = int(np.floor(array_to_propgate_size/2 - initial_pupil.shape[1]/2)), int(np.ceil(array_to_propgate_size/2 - initial_pupil.shape[1]/2))
    pad_up, pad_down = int(np.floor(array_to_propgate_size/2 - initial_pupil.shape[0]/2)), int(np.ceil(array_to_propgate_size/2 - initial_pupil.shape[0]/2))

    thisPupil = np.pad(initial_pupil, pad_width=((pad_left,pad_right),(pad_up, pad_down)), mode='constant')

    return  thisPupil

def Focus_beam(Collimated_Pupil, pad_width = 0):
    """
    Takes the collimated pupil and transforms it to a psf via a Fourier transform 
    """

    Collimated_Pupil_padded = np.pad(Collimated_Pupil,pad_width=int(pad_width),mode='constant') 

    f = np.fft.fft2(Collimated_Pupil_padded) #must be complex amplitude going in here
    fshift = np.fft.fftshift(f)
    intensity_image = (np.abs(fshift))**2
    
    return intensity_image

def angular_to_physical_pixels(angular_pixel_size, focal_length): 
    """
    Takes an angular pixel size in arcseconds/pixel ("/pixel) and converts it to meters/pixel
    or whatever the units of focal length is. 
    """
    plate_scale = 206265/focal_length
    pixel_size_input_image = angular_pixel_size*(1/plate_scale)

    return pixel_size_input_image

if __name__ == '__main__':

    #Contained within the ## below are the parameters you can modify to alter the simulation: 

    ##############################################

    #physical parameters
    input_image_name = "Hubble_Jupiter_image.png" #name of image, in quotation marks including file extension 
    telescope_diameter_m = 203e-3 #in meters
    telescope_focal_length_m = 1.2*3 #in meters
    seeing_arcsec_500nm = 1.5 #in arcseconds
    zenith_angle_deg = 0 #in deg, zero being at the zenith 
    atmosphere = True #True or False, if True simulates atmospheric perturbations. If False simulates purely diffraction effects 
    angular_pixel_size_input_image = 0.05 #arcseconds/pixel (need to calculate this based on angular size of object, number of pixels and scope plate scale) pixel size = "/pixel * 1/plate_scale 
    wavelength = 600e-9 #in meters 
    CCD_pixel_size = 2e-6 #in meters
    CCD_pixel_count = 700 #The pixel width of your simulated CCD

    #simulation parameters
    Num_psfs_to_gen = 10 #number of psfs (and in turn output images) the run will generate 
    pixels_per_ro = 30 #how well you wish to sample your phase screen

    ##############################################

    #open image to convolve 
    im = Image.open(input_image_name)
    im_array = np.asarray(im)

    #makes image square if not square already by cropping
    if np.shape(im_array)[0] != np.shape(im_array)[1]:

        min_dim, max_dim = np.min((np.shape(im_array)[0],np.shape(im_array)[1])), np.max((np.shape(im_array)[0],np.shape(im_array)[1]))

        pixels_to_crop = max_dim - min_dim

        if np.shape(im_array)[0] > np.shape(im_array)[1]:
            left, right = 0, np.shape(im_array)[1]
            top, bottom = np.floor(pixels_to_crop/2), np.shape(im_array)[0] - np.ceil(pixels_to_crop/2)

        if np.shape(im_array)[0] < np.shape(im_array)[1]:
            left, right = np.floor(pixels_to_crop/2), np.shape(im_array)[1] - np.ceil(pixels_to_crop/2)
            top, bottom = 0, np.shape(im_array)[0]

        im = im.crop((left, top, right, bottom)) 

        im_array = np.asarray(im)

    #uncomment here to view the input image
    #plt.imshow(im)
    #plt.show()

    pixel_size_input_image = angular_to_physical_pixels(angular_pixel_size_input_image,telescope_focal_length_m)
    ideal_pixel_size_pupil = (wavelength*telescope_focal_length_m)/(len(im_array)*pixel_size_input_image)
    #print("pixel_size_input_image", pixel_size_input_image)

    r0_cm = fried_parameter_cm(wavelength,arcseconds_of_seeing_500nm=seeing_arcsec_500nm,zenith_angle_deg=zenith_angle_deg)
    telescope_aperture_width_pixels = int(np.ceil((pixels_per_ro/(r0_cm*0.01))*telescope_diameter_m))
    Pixel_size_pupil_plane = telescope_diameter_m/telescope_aperture_width_pixels

    print("telescope_aperture_width_pixels: ", telescope_aperture_width_pixels)

    pixel_size_psf_image_plane = (wavelength*telescope_focal_length_m)/(len(im_array)*Pixel_size_pupil_plane) #term in denominator is the gridwidth in the pupil plane, image plane psf MUST NOT be cropped Verified this is correct!
    #print("pixel_size_psf_image_plane", pixel_size_psf_image_plane)


    phase_screens = []
    if atmosphere == True:

            if int(np.ceil(telescope_aperture_width_pixels)) < 200:
                #print("r0_cm: ",np.round(r0_cm,2))
                ro = telescope_aperture_width_pixels/(telescope_diameter_m/(r0_cm*0.01))  #Based on the assumption that pixelSize = 1. Calculates how many pixels wide a distance of the Fried parameter is based on how many Fried paramters wide the sample box is and how many pixels are present in the sample box
                #print("Pixels/r0: ", np.round(ro,2))
                print("D/r_0: ",np.round(telescope_diameter_m/r0_cm*1e2,2)) #this is the ratio of your telescope size to the Fried parameter 

            else:
                telescope_aperture_width_pixels = 200
                ro = telescope_aperture_width_pixels/(telescope_diameter_m/(r0_cm*0.01))
                print("WARNING: the telescope aperture and seeing conditions do not allow for the sampling of the atm requested.")
                print("Simulation sampling at the maximum resolution ", ro, " Pixels/r0")

            for phaseScreen in M.MegaScreen(numIter=Num_psfs_to_gen, r0=ro, dx=0.5*telescope_aperture_width_pixels, windowShape = (telescope_aperture_width_pixels,telescope_aperture_width_pixels), pixelSize=1): #dx is 1.25x the gird width as this ensures that the tweeter screens of each subsequent phase screen are not correlated
                phase_screens.append(phaseScreen)

    if atmosphere == False:
            for j in range(0,Num_psfs_to_gen):
                phase_screens.append(np.zeros((telescope_aperture_width_pixels, telescope_aperture_width_pixels),dtype=np.complex64))
    print("Phase screens generated")

    #Runs the atmospheric simulations to convert the input phase screens into psfs to be convolved with the image
    image_plane_psfs = []
    for i in range(0,len(phase_screens)):

        complex_amplitude = quick_complex_pupil(phase_screens[i], array_to_propgate_size=len(im_array[0])) #for each loop, complex amplitudes is a array of 2D arrays where each 2D array is the propgated complex amplitude from telescope
        intensity_image = Focus_beam(complex_amplitude)
        image_plane_psfs.append(intensity_image)


    #uncomment this to see the first psf generated
    #plt.imshow(image_plane_psfs[0])
    #plt.show()

    for i in range(0,len(image_plane_psfs)):
        
        x_psf_samples = np.linspace(-pixel_size_psf_image_plane*len(image_plane_psfs[i])/2, pixel_size_psf_image_plane*len(image_plane_psfs[i])/2, len(image_plane_psfs[i]))
        y_psf_samples = np.linspace(-pixel_size_psf_image_plane*len(image_plane_psfs[i])/2, pixel_size_psf_image_plane*len(image_plane_psfs[i])/2, len(image_plane_psfs[i]))

        f = interpolate.interp2d(x_psf_samples, y_psf_samples, image_plane_psfs[i], kind='cubic')

        x_input_image = np.linspace(-pixel_size_input_image*len(im_array)/2, pixel_size_input_image*len(im_array)/2, len(im_array))
        y_input_image = np.linspace(-pixel_size_input_image*len(im_array)/2, pixel_size_input_image*len(im_array)/2, len(im_array))

        resampled_psf = f(x_input_image, y_input_image)
        image_plane_psfs[i] = resampled_psf

    convolved_images = []
    convolved_array_shape = np.shape(signal.convolve(im_array[:,:,0], image_plane_psfs[0]*(1/np.max(image_plane_psfs[0])))) #this line carries out a test convolution to get the shape of the convolved arrays for the variable convolved_array
    for j in range(0,len(image_plane_psfs)):

        convolved_array = np.zeros((convolved_array_shape[0],convolved_array_shape[1],3)) 
        for i in range(0,3):
            temp = signal.convolve(im_array[:,:,i], image_plane_psfs[j]*(1/np.max(image_plane_psfs[j])))
            convolved_array[:,:,i] = temp

        convolved_array = np.uint8((convolved_array)*(254/np.max(convolved_array)))
        convolved_images.append(convolved_array)

    output_images = []
    for j in range(0,len(convolved_images)):

        x_psf_samples = np.linspace(-pixel_size_input_image*len(convolved_images[j])/2, pixel_size_input_image*len(convolved_images[j])/2, len(convolved_images[j]))
        y_psf_samples = np.linspace(-pixel_size_input_image*len(convolved_images[j])/2, pixel_size_input_image*len(convolved_images[j])/2, len(convolved_images[j]))

        x_CCD = np.linspace(-CCD_pixel_size*CCD_pixel_count/2, CCD_pixel_size*CCD_pixel_count/2, CCD_pixel_count)
        y_CCD = np.linspace(-CCD_pixel_size*CCD_pixel_count/2, CCD_pixel_size*CCD_pixel_count/2, CCD_pixel_count)

        output_image = np.zeros((CCD_pixel_count,CCD_pixel_count,3))
        for i in range(0,3):
            f = interpolate.interp2d(x_psf_samples, y_psf_samples, convolved_images[j][:,:,i], kind='cubic')
            temp = f(x_CCD, y_CCD)
            output_image[:,:,i] = temp

        output_image = np.uint8((output_image)*(255/np.max(output_image)))
        output_images.append(output_image)

    current_directory = os.getcwd()
    final_directory = os.path.join(current_directory, r'output_images')
    if not os.path.exists(final_directory):
       os.makedirs(final_directory)

    os.chdir(final_directory)

    imageio.mimsave('Animated_images.gif', output_images, format='.gif') #palettesize = 256 is the max value. The image countours are due to this being too low

    for i in range(0,len(output_images)):
            imageio.imwrite('{}.jpg'.format(i), output_images[i])

    os.chdir(current_directory)