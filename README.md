# Telescope_Simulator
 A python script for not just simulating the field of view but also atmospheric seeing and diffraction effects.


<img src="https://github.com/dmortimer101/Telescope_Simulator/blob/master/Images/Jupiter_Hubble_plus_telescope_3.9m_focal_length_2micron_pixels.gif?raw=true" width="600" height="300" />

Left: is the input, a image of Jupiter taken by the Hubble space telescope. Right: is the output of the simulation, the image of Jupiter as would be seen by a 8" telescope with a focal length of 1.2 meters, using a 3x barlow lens and a CCD with 2 micron sized pixels.

Above is an example of the simulation in action. Not only does this script resample the image, it also simulates the effects of diffraction limited resolution and perturbations due to atmospheric seeing. 

## Installation: 

This script has the following dependencies 

1. numpy
2. matplotlib
3. scipy
4. imageio
6. pillow

All of which can be installed by the command pip install <module name>. All that is required to get started after that is to download this repo and run 
Telescope_simulator.py the script is configured to load the image Hubble_Jupiter_image.png by default and simulate the output of a 8" telescope as described above. 
 
All configurable variables are contained in the #### region below the line #Contained within the ## below are the parameters you can modify to alter the simulation:

 
If you are interested in diving into the physics behind this simulation here is a talk about this simulation aimed at the level of the interested amateur. 

## Use cases: 

### 1. Field of view calculator 

 by setting the following variables:
 
 *telescope_focal_length_m* = the effective focal length of the simulated telescope 
 *angular_pixel_size_input_image* = the angular width of each pixel (for example if Jupiter is taken to be 50" in diameter on sky, devide this by the number of pixels accross it is in the input image) 
 *CCD_pixel_size* (in meters) 
 *CCD_pixel_count* 
 
 the input image is rescaled as if it were sampled by the above described detector attached to the above described telescope 
 
 ### 2. Produce diffraction limited images for a given telescope: 

<img src="https://github.com/dmortimer101/Telescope_Simulator/blob/master/Images/Hubble_Jupiter_8_inch_diffraction_limited_Jupiter.png?raw=true" width="600" height="300" />

Left: Hubble Jupiter image. Right: Diffraction limited image for an 8" telescope 

  by setting *atmosphere = False* the image is rescaled as described above and processed to produce a diffraction limited image
  
### 3. Simulate varying levels of atmospheric seeing: 

<img src="https://github.com/dmortimer101/Telescope_Simulator/blob/master/Images/Jupiter_zenith_10_deg_above_horizon.gif?raw=true" width="600" height="300" />


By setting the variable *seeing_arcsec_500nm* you can vary how "good" the seeing is. 2 arcseconds is considered moderatly good for an amature site, 0.7 arcseconds is considered excelent for a professional observatory  
 the variable zenith_angle_deg allows you to simulate how seeing degrades due to observing targets nearer the horizon (hence with an increase airmass). 0 degrees is a target at zenith, 90 degrees a target at the horizon  


### 4. End to end astrophotography simulator 
<img src="https://github.com/dmortimer101/Telescope_Simulator/blob/master/Images/diffraction_limited_Jupiter_atm_perturbed_Jupiter_RegiStax_Jupiter.gif?raw=true" width="600" height="200" />
 
 Left: Diffraction limited image. Middle: individual atmospherically perturbed frames. Right: result of processing middle frames in RegiStax
 
 The output images can be used as simulated data for astrophotography image processing tools such as RegiStax to test how good an image can be recovered for given seeing conditions and a given number of images. 
  
  
## Future Updates: 
 
1. Adaptive optics simulation 
2. parallelisaion of image generation for efficency 

## Tips for usage: 

1. Ensure *telescope_aperture_width_pixels* (printed in terminal) is above 35 for a given configuration. This can be raised by increasing the *variable pixels_per_ro*.  
 
