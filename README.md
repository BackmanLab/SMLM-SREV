# SMLM-SREV

This series of scripts is for submitting to Quest batches of simulations where each simulation is of 10000 fluorophores (dye labels) diffusing within 1 SREV configuration.

In each simulation of fluorophore diffusivity, 10000 mobile dye labels are diffusing about 1 SREV configuration of ellipsoidal nucleosomes that remain frozen/stuck for the whole simulation pipeline. After a certain amount of time, the coordinates of all dye labels are outputted, and a Euclidean distance check is performed for all dye labels at the last outputted timestep to see if the dye label appears to be colocalized with an SREV nucleosome (this is based on the distances between the center of the dye label and the centers of all SREV nucleosomes). If the dye label does appear to be colocalized with a nucleosome, the dye label gets tagged as a dye label that is now "immobilized" to stay bound to that nucleosome. Following this distance check, the dye labels with their updated tags are inputted into a script for continuing the simulation from where it left off. The dye labels that are tagged for immobilization do not move at all; the mobile dye labels continue to diffuse through the SREV configuration. After a certain amount of time, the coordinates of the dye labels are outputted and the distance check is repeated to tag more dye labels for immobilization before the simulation continues. This cycle continues until some user-inputted number of iterations has occurred. A final "wash" step then occurs where all dye labels that had not been tagged for immobilization or colocalized to a nucleosome at that point are "washed" out of the simulated sample, leaving only the dye labels that had bound to the chromatin.

Tweaking the array and working directory parameters will enable each simulation to run in parallel.

The order of operations is the following:

1. Go to the directory containing the SREV configuration (ideally should have 1 SREV configuration per folder).

2. Copy over the fluorophore diffusivity simulation scripts to that directory.

3. Load the LAMMPS Molecular Dynamics and MPI packages.

4. Enter SREV configuration filenames into the fluorophore diffusivity scripts.

5. Initialize the random positions of the dye labels within the simulation box containing the SREV configuration.

5. Run a brief LAMMPS simulation to generate a version of the SREV configuration that lacks overlapping nucleosomes.

6. Run a brief LAMMPS simulation to output the number of SREV nucleosomes that will be used in the simulation. This is vital for properly setting the temperature of the mobile objects at the desired level.

7. Run the first iteration of the simulation. This simulation run will output a trajectory of all 10000 dye labels over the course of this simulation iteration.

Then, for n number of user-specified iterations:

	8. Input the dye label trajectory file to perform the distance-based colocalization check + appropriate tagging of all dye labels and generate a new dye label input file containing those tags.

	9. Run a brief simulation to obtain output files that give the number of immobilized dye labels and stuck SREV nucleosomes in the simulation. This is important for recomputing the parameter needed to properly set the temperature of the mobile objects at the desired level.

	10. Run the next iteration of the simulation with the new dye label input file. This simulation run will output a trajectory of all 10000 dye labels over the course of this simulation iteration.


Simulate a "wash" step of the simulated labeled chromatin sample by filtering out all dye labels that have not bound to an SREV nucleosome by that time. Output the coordinates of these dye labels still present in the model chromatin post-wash. This output is your post-wash labeled model chromatin sample. 


BEFORE YOU BEGIN:

Set up your directory of SREV configurations in the following manner:

- Have a directory that is for containing all of your SMLM simulations in your set of input SREV configurations. This directory should ONLY contain SREV configurations that were made with the same input parameters (ex: all the SREV configurations in this directory being made with input alpha = 1.15 & phi = 0.12 is good, some SREV configurations in this directory being made with input alpha = 1.15 & phi = 0.12 and other SREV configurations in this directory being made with input alpha = 1.15 and phi = 0.16 is bad).
  

mkdir /yourprojectdir/set_of_simulations/SREV_configs_made_with_same_SREV_input_params/

for i in {1..number_of_srev_configs}; do
    mkdir /yourprojectdir/set_of_simulations/SREV_configs_made_with_same_SREV_input_params/Config_$i
    cp 



