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

<ins>**Specify the simulation box dimensions.**</ins> 

- This should match the SREV config file dimensions. All these scripts are in units of reduced units or "ru", where 1 ru = 10 nm.
Go to InitWalkers.f90, ImmobilConversion.f90, and FinalWash.f90, and ensure all xlength, ylength, and zlength variables are consistent with the box dimensions indicated in the SREV config file and that in the InitWalkers.f90 file all xlowerbound, ylowerbound, and zlowerbound variables are consistent with the lowerbound dimensions stipulated in lines 6-8 of the SREV config file. To identify the xlowerbound, ylowerbound, and zlowerbound information in the SREV config file, understand that line 6 specifies the lower and upper limits of the box in the x direction, line 7 specifies the lower and upper limits of the box in the y direction, and line 8 specifies the lower and upper limits of the box in the z direction. In lines 6-8, the order of the values is lower limit in that dimension and then upper limit in that dimension. Understand that these box dimensions are given in ru and in scientific notation. So, if line 7 says " -1.0000000000000000e+02 1.0000000000000000e+02", that means the box dimension in the y dimension is -100 ru to 100 ru. This means you should go to the lines that define the ylowerbound and ylength variables and set ylowerbound as -100 and ylength as 200 (as y_upper_limit - y_lower_limit = ylength).

<ins>**Update the bash file.**</ins>

Open the execute_sim_in_array.sh file in EdU_Execution_Scripts and in BrdU_Execution_Scripts. Make the following updates:
- In line #3: change the range specified in the definition of array to match the range of SREV configs you wish to simulate SMLM labeling in. So, if you would like to simulate SMLM labeling in 100 SREV configs, change the "_#SBATCH --array=1-10_" line to be "_#SBATCH --array=1-100_".
- In line #13: change the production directory for the location of the directory containing your SREV configs and utilize the SLURM_ARRAY_TASK_ID variable to tell the supercomputer how to locate each SREV config. So, if you have 100 SREV configs in a directory with the name "_/home/my_SMLM_simulations/SREV_configs_made_with_phi_0.12_alpha_1.15/_" and in this directory, each SREV configuration of index _i_ is in a directory with a name of "_SREV_Config_i__" (ex: "SREV_Config1", "SREV_Config2", "SREV_Config3", etc.), change the "_dirname="/projects/p31375/SMLM-SREV-main/BrdU_Phi08/Config_${SLURM_ARRAY_TASK_ID}"_" line to be "_dirname="/home/my_SMLM_simulations/SREV_configs_made_with_phi_0.12_alpha_1.15/SREV_Config_${SLURM_ARRAY_TASK_ID}"_".
- In line 19: change the source directory line for the location of the directory containing the BrdU and EdU SMLM simulation execution scripts. These are the locations of the BrdU_Execution_Scripts folder and EdU_Execution_Scripts folder on your machine/workspace. So, if the BrdU_Execution_Scripts folder and EdU_Execution_Scripts folder are in a directory with the name _"/home/SMLM_Execution_Scripts/BrdU_Execution_Scripts/"_ and _"/home/SMLM_Execution_Scripts/EdU_Execution_Scripts/"_ respectively, in the execute_sim_in_array.sh file located in BrdU_Execution_Scripts, change the "_cp /projects/p31375/SMLM-SREV-main/BrdU_Execution_Scripts/* ._" line to "_cp /home/SMLM_Execution_Scripts/BrdU_Execution_Scripts/* ._" and in the execute_sim_in_array.sh file located in EdU_Execution_Scripts, change the "_cp /projects/p31375/SMLM-SREV-main/EdU_Execution_Scripts/* ._" line to "_cp /home/SMLM_Execution_Scripts/EdU_Execution_Scripts/* ._" . Make sure to preserve the "_* ._" at the end!
- In line 2 and 9: change the line "_#SBATCH --account=p31375_" to specify instead the project you would be running this through on Quest and change the "_#SBATCH --output=log.%j_" to your preferred name for output logs.

## HOW TO EXECUTE:

Go to the directory containing the execution scripts (so EdU_Execution_Scripts or BrdU_Execution_Scripts). After you complete all the steps in this README's "BEFORE YOU BEGIN:" section, type "sbatch execute_sim_in_array.sh". This will submit the SMLM simulation jobs over all your inputted SREV configs to Quest.

## Conceptual Notes

### Why is redefining the temperature parameters important?

Over the course of all the simulation iterations, it is important that you are collecting molecule localizations that occur at the same temperature. When specifying the input target temperature parameters, LAMMPS takes that to specify the temperature of the whole simulation. Since the vast majority of the simulation particles are SREV nucleosome ellipsoids that are kept stuck in place, this complicates the actual energy that gets distributed among the mobile particles to have them move at the specified temperature. Therefore, the input temperature parameters that are given to the LAMMPS simulator program should take this into account so that the mobile particles are actually diffusing at your desired target temperature (in this case the script is set to simulate the objects at a target temperature T = 2.50*).

The equation that computes the temperature to input to produce the desired temperature is the following: number_of_mobile_particles / (number_of_mobile_particles + number_of_immobile_particles) * desired_temp = temp_to_input

## References and Acknowledgements

The simulations scripted here were performed using the LAMMPS Molecular Dynamics Simulator<sup>1</sup>. For more on LAMMPS, consult https://www.lammps.org . Interactions between the particles were modeled using LAMMPS' *pair gayberne* framework for modeling Gay-Berne interaction potentials between ellipsoidal particles<sup>2</sup>.

1. LAMMPS - A flexible simulation tool for particle-based materials modeling at the atomic, meso, and continuum scales, Comp. Phys. Comm. 271, 108171 (2022)

2. Brown, W. M., Petersen, M. K., Plimpton, S. J., & Grest, G. S. (2009). Liquid crystal nanodroplets in solution. *The Journal of Chemical Physics, 130*(4), 044901. https://doi.org/10.1063/1.3058435

We acknowledge the contributions of Dr. Marcelo Carignano in determining the formula for obtaining the desired simulation temperature based on the number of mobile particles and immobilized particles.
