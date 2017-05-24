#!/bin/bash

PROG="../test/parsing_xml"

INPUT="../test/input_files/fail_parsing_1.xml  ../test/input_files/solar_flux.dat ../test/input_files/CH4_hv_cs.dat"

$PROG $INPUT

