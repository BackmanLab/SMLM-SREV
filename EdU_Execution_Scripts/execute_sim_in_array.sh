#!/bin/bash
#SBATCH --account=p31375
#SBATCH --array=1-10
#SBATCH --partition=long
#SBATCH --time=120:00:00
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=16
#SBATCH --mem=30G
#SBATCH --output=log.%j

#cd /projects/p31375/STORM/SphericalSize/MultipleConfigs/FixAndContinue/Phi12/Phi12/Config_7/

dirname="/projects/p31375/STORM/SphericalSize/MultipleConfigs/Set4/Phi08/Config_${SLURM_ARRAY_TASK_ID}"

cd $dirname

cp /projects/p31375/STORM/SphericalSize/MultipleConfigs/FixAndContinue/EdU_Execution_Scripts/* .

rm walkers.xyz
rm walkerinput.in
rm edited-config*

module purge

module load lammps/20200303-openmpi-4.0.5-intel-19.0.5.281

module load mpi/openmpi-4.1.4-gcc

lammpsdump="config-relaxed-${SLURM_ARRAY_TASK_ID}.dump"
editedlammpsdump="edited-config-${SLURM_ARRAY_TASK_ID}.dump"

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

mpirun -np 32 lmp -in getconfig_nooverlaps.in

mpirun -np 32 lmp -in overlapfilter.in

var=$(sed -n '4,4p' dump.ellipsoid)
echo $var
a=10000
c=$var
d=2.50
e=$(echo "scale=10; $a / ( $a + $c ) * $d" | bc -l)
echo $e
sed -i "s:^fix 1 all nvt temp .*:fix 1 all nvt temp $e $e 1.0:" first_sim.in
rm dump.ellipsoid

mpirun -np 32 lmp -in first_sim.in

for i in $(seq 1 10); do
	rm post-immobil-input.in
	gfortran -Ofast ImmobilConversion.f90
	./a.out
	rm walkers.xyz
	rm dump.ellipsoid*
	sed -i "s:^read_dump .*:read_dump ${lammpsdump} 0 x y z add yes box yes:" overlapfilter_continue.in
	mpirun -np 32 lmp -in overlapfilter_continue.in
	var=$(sed -n '4,4p' dump.immob)
	var2=$(sed -n '4,4p' dump.ellipsoid)
	echo $var
	a=$(echo "10000 - ( $var - $var2 )" | bc -l)
	echo $a
	c=$var
	d=2.50
	e=$(echo "scale=10; $a / ( $a + $c ) * $d" | bc -l)
	echo $e
	sed -i "s:^fix 1 all nvt temp .*:fix 1 all nvt temp $e $e 1.0:" continue_sim.in
	rm dump.immob
	rm dump.ellipsoid
	mpirun -np 32 lmp -in continue_sim.in
done
