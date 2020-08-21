# Telescope_Simulator
 A python script for not just simulating the field of view but also atmospheric seeing and diffraction effects.


<img src="https://github.com/dmortimer101/Telescope_Simulator/blob/master/Images/Jupiter_Hubble_plus_telescope_3.9m_focal_length_2micron_pixels.gif?raw=true" width="600" height="300" />


Above is an example of the simulation in action. Left is the input, a image of Jupiter taken by the Hubble space telescope. Right is the output of the simulation, the image of Jupiter as would be seen by a 8" telescope with a focal length of 1.2 meters, using a 3x barlow lens and a CCD with 2 micron sized pixels. Not only does this script resample the image, it also simulates the effects of diffraction limited resolution and perturbations due to atmospheric seeing. 

Installation: 

This script has the following dependencies 

1. numpy
2. matplotlib
3. scipy
4. imageio
6. pillow

All of which can be installed by the command pip install <module name>. All that is required to get started after that is to download this repo and run 
Telescope_simulator.py the script is configured to load the image Hubble_Jupiter_image.png by default and simulate the output of a 8" telescope as described above. 
 
