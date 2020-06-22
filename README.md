# Inner_ear_vertex_code
Matlab code for simulations of OoC development

The code simulates the organization of outer hair cells in the cochlea during embryonic development.
Description of the model can be found at:
https://doi.org/10.1101/707422

Random cellular lattices can be found in "random_lattices" folder.
By choosing inner hair cells and pillar cells in a random lattice we create the initial conditions for the described model.
Initial lattices for the simulations are found in "initial_lattices".

To run the simulations open and run the file "runSimulations.m".
The simulations for stage 1 and stage 2 of the model will start to run. 
On a 3.4GHz processor core a single simulation with the default parameters takes roughly 4 minutes. 
A total of 20 simulations (which is the default in "runSimulations.m") will then take 1-1.5 hours.

Each simulation saves snapshots of the cellular lattice during run:
Time points for stage 1 are saved in the folder "stage1_timepoints" in the format "lat(#)_step(#).mat".
Time points for stage 2 are saved in the folder "stage2_timepoints" in the format "lat(#)_step(#).mat".
The first number represents the cellular lattice index and the second number represents the time point.

To plot the lattice at a certain timepoint, load the desired time point ("lat(#)_step(#).mat" file) and use "LatticePresentation(g,0)".
For presentation of the lattice including cells numbers use "LatticePresentation(g,1)".
