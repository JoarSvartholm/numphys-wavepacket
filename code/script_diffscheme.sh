#!/bin/bash
#
# Sample script to run the wavepacket program repeatedly
# modifying the parameter file at each execution
#
# This script must be an executable file, the command
# ls -l script_sb.sh
# must show something like
# -rwxr-xr-x 1 cldi0001 adusers 484 mar 10  2016 script_sb.sh
#    *  *  * with x's
# If this is not the case, then do
# chmod +x script_sb.sh
#
# Note:  The parameters for the Gaussian wave packet are read by my function
# initialize_wf from a file.  Therefore, I execute my program using
# ./single_barrier single_barrier.in gauss.in
# where single_barrier.in is the wavepacket program parameter file and
# gauss.in is the file read by initialize_wf.  Modifications to this script
# may be necessary depending on how you set the value of k0 in your program.
#
# Claude Dion, March 2017
#

# Create an empty file that will contain the results as a function of k0
# *** will erase the file if it already exists ***
echo -n > ../data/diff_RT_all.dat

# Set initial value of k0
K=150

# Loop over values of k0
while [ $K -le 350 ]; do

    # Copy template file containing initialize_wf parameters except k0
    cp input_files/diff_scheme_template.in input_files/diff_$K.in

    # Add results_file with k0-dependent file name
    echo results_file = ../data/diff_k0_$K.dat >> input_files/diff_$K.in
    echo wf_output_text = ../data/diff_barrier_wave_k0_$K.dat >> input_files/diff_$K.in

    # Print on screen current value of k0
    echo k0 = $K
    # Execute program with timing information, sending all
    # standard output to a  k0-dependent file
    time ./diffScheme input_files/diff_$K.in $K > output_files/diff_k0_$K.out

    # Add value of k0 to file single_barrier_RT_all.dat
    echo -n $K" " >> ../data/diff_RT_all.dat
    # Add last time point of output from user_observe
    tail -1 ../data/diff_barrier_RT_k0_$K.dat >> ../data/diff_RT_all.dat

    # Increment k0 before looping
    let K=K+10
done
