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

All of which can be installed by the command pip install "module name". All that is required to get started after that is to download this repo and run 
Telescope_simulator.py the script is configured to load the image Hubble_Jupiter_image.png by default and simulate the output of a 8" telescope as described above. 
 
All configurable variables are contained in the #### region below the line #Contained within the ## below are the parameters you can modify to alter the simulation:

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
  
  
## Future Updates: 
 
1. Adaptive optics simulation 
2. Parallelisation of image generation for efficiency 

## Tips for usage: 

1. Ensure *telescope_aperture_width_pixels* (printed in terminal) is above 35 for a given configuration. This can be raised by increasing the variable *pixels_per_ro*.  

2. Objects with a large angular diameter (the width of the full moon or bigger) will return a monochromatic square. This is believed to be due to the fact that the psf is poorly sampled on images where the pixels are so large.
 
