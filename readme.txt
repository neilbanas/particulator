Particulator
Neil Banas, 2018


QUICK START

For a simple example of a script that identifies a ROMS run, sets up a particle
experiment, and then tracks the particles, see
examples/example_romsCascadia.m

For more detail on the setup options, see the list of defaults at the start of 
par_release.m, and
examples/example_JdF_verticalOptions.m
examples/example_JdF.m
for further examples.



HOW IT WORKS

par_integrate.m is the particle-tracking code itself. The algorithm consists of reading in one timestep ("frame") of model output at a time, with only two frames held in RAM, and advancing a vector of particles through a series of internal steps between those frames. This approach assumes that reading model output from disk is a slow step, and that it is hard to improve upon the efficiency of the vector operations in matlab.

Output is saved at the same interval as the model input files. Timestepping uses the midpoint method (so that each internal step involves two Euler steps): see takeStep() in par_integrate.m. On each step, particles advance to a new (x,y,z,t) position, and then velocities, tracers, etc. are interpolated at the new location. Various shortcuts are taken to prevent bad values with a minimum of fuss: for example, vertical positions are simply clipped to fall between sigma (z/H) between -1 and 0, as opposed to reflective surface and bottom b.c.'s. When a step is calculated to carry a particle onto land, we simply don't take the step--leave the particle where it is, advance time, and hope that at the next timestep the flow will have changed. (This behaviour can be changed with the "avoidLand" flag in par_release.m.) Some of these shortcuts have been better tested than others.

Note that there are options both to interpolate tracers at the exact (x,y,z) location of the particles, and also to save entire vertical profiles of tracers at the (x,y) locations of the particles. (This is useful when using particle tracks as a basis for the ensemble-Lagrangian NPZ model of Banas et al., JGR, 2016, or other kinds of "plan view" applications.)

modelRun.m is the interface between par_integrate.m and particular hydrodynamic model runs. See modelRun_romsCascadia.m for the details of how ROMS output files are read in and interpolated. To add support for a new model (or a new flavour of ROMS output), DON'T modify par_integrate.m or modelRun.m (if you can at all help it), but rather write a new case, similar to modelRun_romsCascadia.m.

Once a modelRun has been created (let's call it "run"), and a particular frame of model output has been loaded into memory with run.loadFrame(), what was loaded in can be accessed directly as run.F1, bypassing the interp() functions. (This can be useful for debugging.)

It will be a while before I write a proper manual, but please let me know what other questions I should answer in this readme!
