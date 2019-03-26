# darkmatter
Dark Matter Project

Included is a simulation that produces a highly-populated recoil spectrum for the LUX experiment (the target) for a given WIMP mass and WIMP-nucleon cross section, and then generates sample datasets for a single run of LUX at each point in a 100x100 mass/cross section phase space. A K.S test is run for the target and each sample, with 50 iterations at each point to produce a mean p-value. These p-values help to quantify the ability of these measurements to identify the mass and cross section.

The velocity distribution on the samples can be varied to explore the systematic error introduced by uncertainty in this. 

'yield.dat' and 'formfactor.dat' represent input files whose raw data is used in the main simulation
- will need to be placed in same folder as the simulation when run

Consult the accompanying report for more information
