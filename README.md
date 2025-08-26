# SMLM-SREV

## Overview

This series of scripts is for submitting to the Northwestern Quest Supercomputer Cluster batches of simulations where each simulation is of 10000 fluorophores (dye labels) diffusing within 1 SREV configuration.

In each simulation of fluorophore diffusivity, 10000 mobile dye labels are diffusing about 1 SREV configuration of ellipsoidal nucleosomes that remain frozen/stuck for the whole simulation pipeline. After a certain amount of time, the coordinates of all dye labels are outputted, and a Euclidean distance check is performed for all dye labels at the last outputted timestep to see if the dye label appears to be colocalized with an SREV nucleosome (this is based on the distances between the center of the dye label and the centers of all SREV nucleosomes). If the dye label does appear to be colocalized with a nucleosome, the dye label gets tagged as a dye label that is now "immobilized" to stay bound to that nucleosome. If the dye label does not appear to be colocalized with a nucleosome however, the dye label gets tagged as a dye label that should continue to be mobile in the next simulation iteration. Following this distance check, the dye labels with their updated tags are inputted into a script for continuing the simulation from where it left off. The dye labels that are tagged for immobilization do not move at all; the mobile dye labels continue to diffuse through the SREV configuration. After a certain amount of time, the coordinates of the dye labels are outputted and the distance check is repeated to tag more dye labels for immobilization before the simulation continues. This cycle continues until some user-inputted number of iterations has occurred. A final "wash" step then occurs where all dye labels that had not been tagged for immobilization or colocalized to a nucleosome at that point are "washed" out of the simulated sample, leaving only the dye labels that had bound to the chromatin.

Tweaking the array and working directory parameters will enable each simulation to run in parallel.

## Workflow

The order of operations is the following:

1. Go to the directory containing the SREV configuration (ideally should have 1 SREV configuration per folder).

2. Copy over the fluorophore diffusivity simulation scripts to that directory.

3. Load the LAMMPS Molecular Dynamics and MPI packages.

4. Enter SREV configuration filenames into the fluorophore diffusivity scripts.

5. Initialize the random positions of the dye labels within the simulation box containing the SREV configuration.

5. Run a brief LAMMPS simulation to generate a version of the SREV configuration that lacks overlapping nucleosomes.

6. Run a brief LAMMPS simulation to output the number of SREV nucleosomes that will be used in the simulation. This is vital for properly setting the temperature of the mobile objects at the desired level.

7. Run the first iteration of the simulation. This simulation run will output a trajectory of all 10000 dye labels over the course of this simulation iteration.

For *n* number of user-specified iterations, loop through Steps #8-10 below:

> 8. Input the dye label trajectory file to perform the distance-based colocalization check + appropriate tagging of all dye labels and generate a new dye label input file containing those tags.
> 9. Run a brief simulation to obtain output files that give the number of immobilized dye labels and stuck SREV nucleosomes in the simulation. This is important for recomputing the parameter needed to properly set the temperature of the mobile objects at the desired level.
> 10. Run the next iteration of the simulation with the new dye label input file. This simulation run will output a trajectory of all 10000 dye labels over the course of this simulation iteration.

11. Simulate a "wash" step of the simulated labeled chromatin sample by filtering out all dye labels that have not bound to an SREV nucleosome by that time. Output the coordinates of these dye labels still present in the model chromatin post-wash. This output is your post-wash labeled model chromatin sample. 


## BEFORE YOU BEGIN:

Set up your directory of SREV configurations in the following manner:

- Have a production directory that is for containing all of your SMLM simulations in your set of input SREV configurations. This directory should ONLY contain SREV configurations (configs) that were made with the same input parameters (ex: all the SREV configs in this directory being made with input alpha = 1.15 & phi = 0.12 is good, some SREV configs in this directory being made with input alpha = 1.15 & phi = 0.12 and other SREV configs in this directory being made with input alpha = 1.15 & phi = 0.16 is bad).

- For the *n* number of SREV configs in this production directory, have *n* directories (one for each SREV config). Each SREV config should be assigned its own index that is between 1 to *n*. The SREV config file name should be named as "config-$_i_.dump", where $_i_ is the index of that SREV config. Each of the directories within the production directory should also be assigned its own index that is between 1 to *n*. Each of these directories within the production directory should be named as "Config_$_i_", where $_i_ is the assigned index of that directory. Each of these directories should only contain the SREV config that is of the same index as the directory.

Specify the simulation box dimensions. This should match the SREV config file dimensions. All these scripts are in units of reduced units or "ru", where 1 ru = 9.8 nm.
Go to InitWalkers.f90, ImmobilConversion.f90, and FinalWash.f90, and ensure all boxlength variables are consistent with the box dimensions stipulated in the SREV config file (so 200d0 or 130d0 for 2000 nm long or 1300 nm long, respectively) and that all xlowerbound, ylowerbound, and zlowerbound variables are consistent with the lowerbound dimensions stipulated in the SREV config file (so -65 for -650 nm or -100 for -1000 nm).

## Conceptual Notes

### Why is redefining the temperature parameters important?

Over the course of all the simulation iterations, it is important that you are collecting molecule localizations that occur at the same temperature. When specifying the input target temperature parameters, LAMMPS takes that to specify the temperature of the whole simulation. Since the vast majority of the simulation particles are SREV nucleosome ellipsoids that are kept stuck in place, this complicates the actual energy that gets distributed among the mobile particles to have them move at the specified temperature. Therefore, the input temperature parameters that are given to the LAMMPS simulator program should take this into account so that the mobile particles are actually diffusing at your desired target temperature (in this case the script is set to simulate the objects at a target temperature T = 2.50*).

### What does "ru" mean for how I change out different variables?

Explanations on reduced units used.

## References and Acknowledgements

The simulations scripted here were performed using the LAMMPS Molecular Dynamics Simulator<sup>1</sup>. For more on LAMMPS, consult https://www.lammps.org . Interactions between the particles were modeled using LAMMPS' *pair gayberne* framework for modeling Gay-Berne interaction potentials between ellipsoidal particles<sup>2</sup>.

1. LAMMPS - A flexible simulation tool for particle-based materials modeling at the atomic, meso, and continuum scales, Comp. Phys. Comm. 271, 108171 (2022)

2. Brown, W. M., Petersen, M. K., Plimpton, S. J., & Grest, G. S. (2009). Liquid crystal nanodroplets in solution. *The Journal of Chemical Physics, 130*(4), 044901. https://doi.org/10.1063/1.3058435
