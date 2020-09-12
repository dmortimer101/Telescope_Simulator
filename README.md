# Telescope_Simulator
 A python script for not just simulating the field of view but also atmospheric seeing and diffraction effects of a telescope.


<img src="https://github.com/dmortimer101/Telescope_Simulator/blob/master/Images/Jupiter_Hubble_plus_telescope_3.9m_focal_length_2micron_pixels.gif?raw=true" width="600" height="300" />

Left: is the input, an [image of Jupiter](https://www.nasa.gov/feature/goddard/2019/hubble-new-portrait-of-jupiter) taken by the Hubble space telescope. Right: is the output of the simulation, the image of Jupiter as would be seen by an 8" telescope with a focal length of 1.2 meters, using a 3x Barlow lens and a CCD with 2-micron sized pixels.

Above is an example of the simulation in action. Not only does this script resample the image, it also simulates the effects of diffraction limited resolution and perturbations due to atmospheric seeing. 

The script is written to be easy to read as it is in part used as a teaching tool. There are likely more efficient ways of writing this code but for the sake of clarity they have not been used. 

## Installation: 

This script has the following dependencies 

1. NumPy
2. MegaScreen
3. Matplotlib
4. SciPy
5. imageio
6. Pillow

All of which can be installed by the command pip install "module name". 

## Usage: 

This script is designed to be as simple as possible to get started with. Once the installation described above is complete *Telescope_simulator.py* can be run as is and will by default produce 10 images of Jupiter as seen by an 8 inch telescope with an effective 3.6-meter focal length under average seeing conditions. The parameters which can be reconfigured are under the comments *#physical parameters* (line 90) and *#simulation parameters* (line 102) in the main file, *Telescope_simulator.py*. There are descriptive comments alongside each parameter which are not repeated here however below is a suggested range for each parameter: 

1. *telescope_diameter_m*: of order a few hundreds of mm upwards (100e-3 upwards). Due to the nature of how the atmospheric perturbations are simulated the upper limit on the diameter will depend on the level of seeing, with worse seeing reducing the maximum telescope diameter that can be simulated
2. *telescope_focal_length_m*: of order a few hundreds of mm to a few meters (300e-3 to 5 meters) 
3. *seeing_arcsec_500nm*: between 0.7 and 2 arcseconds is a typical range from excellent to moderately poor seeing 
4. *zenith_angle_deg*: between 0 (at the zenith) and 70 (close to the horizon) degrees 
5. *wavelength*: typically the visible wavelength range, 380-740 nanometers (380e-9 to 740e-9 meters)
6. *CCD_pixel_size*: a few to a few tens of microns (1e-6 to 20e-6 meters) 
7. *CCD_pixel_count*: A few hundred pixels
8. *pixels_per_ro*: leave at 30 by default but always keep above 7 to ensure accuracy in simulating atmospheric perturbations 


## Talk: 
 
If you are interested in diving into the physics behind this simulation here is a talk aimed at the level of the interested amateur: 

[![How to simulate a telescope (in python)](https://img.youtube.com/vi/SuoV_5ACOIo/0.jpg)](https://www.youtube.com/watch?v=SuoV_5ACOIo)

## Use cases: 

### 1. Field of view calculator 

 by setting the following variables:
 
 1. *telescope_focal_length_m* = the effective focal length of the simulated telescope 
 2. *angular_pixel_size_input_image* = the angular width of each pixel (for example if Jupiter is taken to be 50" in diameter on sky, divide this by the number of pixels across it is in the input image) 
 3. *CCD_pixel_size* (in meters) 
 4. *CCD_pixel_count* 
 
the input image is rescaled as if it were sampled by the above described detector attached to the above described telescope.   

 ### 2. Produce diffraction limited images for a given telescope: 

<img src="https://github.com/dmortimer101/Telescope_Simulator/blob/master/Images/Hubble_Jupiter_8_inch_diffraction_limited_Jupiter.png?raw=true" width="600" height="300" />

Left: Hubble Jupiter image. Right: Diffraction limited image for an 8" telescope.

  by setting *atmosphere = False* the image is rescaled as described above and processed to produce a diffraction limited image
  
### 3. Simulate varying levels of atmospheric seeing: 

<img src="https://github.com/dmortimer101/Telescope_Simulator/blob/master/Images/Jupiter_zenith_10_deg_above_horizon.gif?raw=true" width="600" height="300" />

Left: Jupiter observed at Zenith. Right: Jupiter observed 10 degrees above the horizon.

By setting the variable *seeing_arcsec_500nm* you can vary how "good" the seeing is. 2 arcseconds is considered moderately good for an amateur site, 0.7 arcseconds is considered excellent for a professional observatory. The variable *zenith_angle_deg allows* you to simulate how seeing degrades due to observing targets nearer the horizon (hence with an greater airmass). 0 degrees is a target at zenith, 90 degrees a target at the horizon.


### 4. End to end astrophotography simulator 
<img src="https://github.com/dmortimer101/Telescope_Simulator/blob/master/Images/diffraction_limited_Jupiter_atm_perturbed_Jupiter_RegiStax_Jupiter.gif?raw=true" width="600" height="200" />
 
 Left: Diffraction limited image. Middle: individual atmospherically perturbed frames. Right: result of processing middle frames in RegiStax.
 
 The output images can be used as simulated data for astrophotography image processing tools such as RegiStax to test how good an image can be recovered for given seeing conditions and a given number of images. 
  
  
## Future updates: 
 
1. Adaptive optics simulation 
2. Parallelisation of image generation for efficiency 

## Tips for usage: 

1. Ensure *telescope_aperture_width_pixels* (printed in terminal) is above 35 for a given configuration. This can be raised by increasing the variable *pixels_per_ro*.  

2. Objects with a large angular diameter (the width of the full moon or bigger) will return a monochromatic square. This is believed to be due to the fact that the psf is poorly sampled on images where the pixels are so large.
 
