#!/bin/bash
#SBATCH --account=p31375
#SBATCH --array=1-10
#SBATCH --partition=normal
#SBATCH --time=48:00:00
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=16
#SBATCH --mem=30G
#SBATCH --output=log.%j

#TRADE OUT DIRNAME FOR DIRECTORY TO YOUR SREV CONFIGS

dirname="/projects/p31375/SMLM-SREV-main/BrdU_Phi08/Config_${SLURM_ARRAY_TASK_ID}"

cd $dirname

# TRADE OUT SOURCE DIRECTORY YOU COPY FILES FROM FOR LOCATION OF THESE SCRIPTS

cp /projects/p31375/SMLM-SREV-main/BrdU_Execution_Scripts/* .

# safeguard preventing combination of data from multiple simulations in case you need to rerun anything

rm walkers.xyz
rm walkerinput.in
rm post-immobil*
rm edited-config*
rm log*
rm dump*

# load packages

module purge

module load lammps/20200303-openmpi-4.0.5-intel-19.0.5.281

module load mpi/openmpi-4.1.4-gcc

# automatically trade out config-specific filenames in scripts

lammpsdump="config-relaxed-${SLURM_ARRAY_TASK_ID}.dump"
editedlammpsdump="edited-config-${SLURM_ARRAY_TASK_ID}.dump"
outputwashfile="BrdU-postwash-config-${SLURM_ARRAY_TASK_ID}.xyz"

sed -i "s:^read_dump .*:read_dump ${lammpsdump} 0 x y z add yes box yes:" getconfig_nooverlaps.in
sed -i "s:^dump 1 ellip custom 100 edited-config.*:dump 1 ellip custom 100 ${editedlammpsdump} id type x y z c_0[*]:" getconfig_nooverlaps.in

sed -i "s:^  configfilename = .*:  configfilename = '${lammpsdump}':" InitWalkers.f90
sed -i "s/\r//g" InitWalkers.f90
gfortran -Ofast InitWalkers.f90 -o InitWalkers
./InitWalkers

sed -i "s:^  configfilename = .*:  configfilename = '${editedlammpsdump}':" ImmobilConversion.f90
sed -i "s/\r//g" ImmobilConversion.f90

sed -i "s:^  configfilename = .*:  configfilename = '${editedlammpsdump}':" GetColoc.f90
sed -i "s/\r//g" GetColoc.f90

sed -i "s:^read_dump .*:read_dump ${lammpsdump} 0 x y z add yes box yes:" overlapfilter.in
sed -i "s:^read_dump .*:read_dump ${lammpsdump} 0 x y z add yes box yes:" first_sim.in
sed -i "s:^read_dump .*:read_dump ${lammpsdump} 0 x y z add yes box yes:" continue_sim.in

sed -i "s:^  configfilename = .*:  configfilename = '${editedlammpsdump}':" FinalWash.f90
sed -i "s:^open(unit = 9, file =.*:open(unit = 9, file = '${outputwashfile}'):" FinalWash.f90

# get version of SREV configs that lack all overlapping nucleosomes

mpirun -np 32 lmp -in getconfig_nooverlaps.in

# get needed values for formula for temp to input that produces target mobile temp

mpirun -np 32 lmp -in overlapfilter.in

echo "first ellipsoid particle number"
var=$(sed -n '4,4p' dump.ellipsoid)
echo $var
a=10000
c=$var
d=2.50
e=$(echo "scale=10; $a / ( $a + $c ) * $d" | bc -l)
echo $e
sed -i "s:^fix 1 all nvt temp .*:fix 1 all nvt temp $e $e 1.0:" first_sim.in
rm dump.ellipsoid

# run first iteration of simulation

mpirun -np 32 lmp -in first_sim.in

for i in $(seq 1 10); do
	# below rm command clears nothing when i==1; kept in anyway since it performs rm when needed for rest of i's in loop range
	rm post-immobil-input.in
	# produce version of label input data that marks labels that should be immobilized at the conclusion of that iteration
	gfortran -Ofast ImmobilConversion.f90
	./a.out
	rm walkers.xyz
	rm dump.ellipsoid*
	# get needed values for formula for temp to input that produces target mobile temp
	sed -i "s:^read_dump .*:read_dump ${lammpsdump} 0 x y z add yes box yes:" overlapfilter_continue.in
	mpirun -np 32 lmp -in overlapfilter_continue.in
	var=$(sed -n '4,4p' dump.immob)
	var2=$(sed -n '4,4p' dump.ellipsoid)
	echo "immob particle numbers" # print statements monitoring changes in number of mobile particles
	echo $var
	a=$(echo "10000 - ( $var - $var2 )" | bc -l)
	echo $var2
	echo $a
	c=$var
	d=2.50
	e=$(echo "scale=10; $a / ( $a + $c ) * $d" | bc -l)
	echo $e
	sed -i "s:^fix 1 all nvt temp .*:fix 1 all nvt temp $e $e 1.0:" continue_sim.in
	rm dump.immob
	rm dump.ellipsoid
	# run next iteration of simulation
	mpirun -np 32 lmp -in continue_sim.in
done

gfortran -Ofast FinalWash.f90

./a.out
