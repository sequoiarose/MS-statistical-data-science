Files for phase congruency from Kovesi:
- filtergrid
- hysthresh
- lowpassfilter
- nonmaxsup
- perfft2
- phasecongomono
- smoothing

Files I created:
- project_script_phase_II
- project_script_phase_III
*both of these should run as is if the data and code folders are placed in a directory that can run matlab
*image import paths may need to be changed to get it to run. they are currently configured for the set up in Andrade_Project.


project_script_phase_II:

enumerates through the ideal images and parameter values to create animations of phase congruency outputs.

project_script_phase_III:

performs quantitative parameter tuning on both the parasite images and ideal images.
tests two proposed pipelines for including sharpening on the parasite images.

WARNING: the parameter tuning for the parasite images in phase III has a long run time
