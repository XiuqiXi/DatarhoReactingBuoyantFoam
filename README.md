# DatarhoReactingBuoyantFoam
This version of simulation is to calculate the droplet combustion, which is the simplest case, one dimension, including the influence of buoyancy.

`DatarhoReactingBuoyantFoam.C` is modified from `rhoReactingBuoyantFoam.C` including `dataimport.H` `analyticalRegion.H`, `assimilation.H`,

`dataImport.H` is to import the droplet size data from the data file. `Countlines.H` is to count the lines of data, coupled with the data import H file.

`assimilation.H` is the core part to exchange the heat and mass flux between the CFD region and analytical region. `analyticalRegion.H `claims the class, the data and member functions. `analyticalRegionI.H` defines the member functions including the calculation of boundary conditions.
