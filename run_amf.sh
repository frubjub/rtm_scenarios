#!/bin/bash
# Copyright Stephen Broccardo, 2018

# read in function definitions
. ./build_control_inp
. ./build_control_ac_inp
. ./build_man_aer_inp
. ./build_control_geom_inp

# append the AOT:
AOT=01

# some admin:
# change SCIATRAN_DIR to be the correct path for your setup!
SCENARIO_DIR=$HOME/doas/sciatran/steve_scenarios
SCIATRAN_DIR=$HOME/doas/sciatran/Execute-2.2.2
cd $SCIATRAN_DIR
ulimit -s unlimited

# the paths of the control files
CONTROL_INP=$SCIATRAN_DIR/control.inp
CONTROL_AC_INP=$SCIATRAN_DIR/control_ac.inp
MAN_AER_INP=$SCIATRAN_DIR/man_aer.inp
CONTROL_GEOM_INP=$SCIATRAN_DIR/control_geom.inp

# the paths of the output files
OUTPUT_FILES="$SCIATRAN_DIR/DATA_OUT/AEROSOL.OUT $SCIATRAN_DIR/DATA_OUT/vert_col.dat $SCIATRAN_DIR/DATA_OUT/output_map.inf $SCIATRAN_DIR/DATA_OUT/amf.dat $SCIATRAN_DIR/DATA_OUT/SCIATRAN_SCENARIO.OUT $SCIATRAN_DIR/errors.log"

# remove the old control files and output files
# there is another script (janitor) for removing old files from the output directories
rm -f $CONTROL_INP $CONTROL_AC_INP $MAN_AER_INP $CONTROL_GEOM_INP
rm -f $OUTPUT_FILES


#for scenario in scenario185 scenario19 scenario20 scenario21 scenario22 scenario23
#for scenario in scenario01 scenario02 scenario03 scenario04 scenario05 scenario06 scenario07 scenario08 scenario09 scenario10 scenario11 scenario12 scenario125 scenario13 scenario14 scenario15 scenario16 scenario17 scenario18 scenario185 scenario19 scenario20 scenario21 scenario22 scenario23 scenario24 scenario25 scenario26 scenario27 scenario28 scenario29 scenario30 scenario31 scenario32 scenario33 scenario34 scenario35 scenario36
for scenario in scenario_s10 
do
	for albedo in albedo002 albedo005 albedo008 albedo011
	do
		for ssa in ssa082 ssa090 ssa098 # ssa100
		do 
			build_control_inp $scenario $albedo $CONTROL_INP amf
			build_control_ac_inp $scenario $CONTROL_AC_INP $ssa
			build_man_aer_inp $scenario $MAN_AER_INP $ssa
			build_control_geom_inp $scenario $CONTROL_GEOM_INP amf

			# run the model, then move the output files to the relevant
			# directory, and name them according to who is running the 
			# model, so we can compare Steve's and Alex's results:

			./SCIA_INTEL.exe && 
			for file in ./DATA_OUT/*.OUT;
	       		do mv $file $SCENARIO_DIR/$scenario/amf/output/$albedo/$ssa/`basename $file .OUT`_$AOT.out;
		       	done && 
		       	for file in ./DATA_OUT/*.dat; 
			do mv $file $SCENARIO_DIR/$scenario/amf/output/$albedo/$ssa/`basename $file .dat`_$AOT.dat; 
			done && 
			mv ./errors.log $SCENARIO_DIR/$scenario/amf/output/$albedo/$ssa/errors_$AOT.log 
			mv ./DATA_OUT/output_map.inf $SCENARIO_DIR/$scenario/amf/output/$albedo/$ssa/output_map_$AOT.inf
			mv control.inp control_ac.inp man_aer.inp control_geom.inp \
					$SCENARIO_DIR/$scenario/amf/output/$albedo/$ssa/
			
		done

	done
done

