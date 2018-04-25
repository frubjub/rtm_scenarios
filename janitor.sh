#!/bin/bash
# Copyright (C) 2018 Stephen Broccardo

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

# cleans out all the results from previous runs of sciatran

SCENARIO_DIR=/home/steve/doas/sciatran/steve_scenarios

for scenario in 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36
do
	for albedo in 002 005 007 008 009 011
	do
		for ssa in 082 090 098 100
		do
			rm -vf $SCENARIO_DIR/scenario$scenario/amf/output/albedo$albedo/ssa$ssa/*
			rm -vf $SCENARIO_DIR/scenario$scenario/block_amf/output/albedo$albedo/ssa$ssa/*
		done
	done
done

