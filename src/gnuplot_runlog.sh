#!/bin/bash
#Script by Angela Shipman, v1 Sep 18 2023

#Pretty Formatting
#COLORS
red=$(tput setaf 1)
green=$(tput setaf 2)
yellow=$(tput setaf 3)
blue=$(tput setaf 4)
magenta=$(tput setaf 5)
cyan=$(tput setaf 6)
bold=$(tput bold)
reset=$(tput sgr0)

#SCRIPT BEGINS
echo -e "${green}${bold}What would you like to plot on the X axis?\n\n${reset}${green}1) Total Energy of the system [a.u.]\n2) Step Number\n3) Conserved energy [hartree]\n4) CPU time per MD step [s] \n5) Energy drift per atom [K]\n6) Potential energy [hartree]\n7) Kinetic energy [hartree]\n8) Instan. Temperature [K]\n9) Avg. Temperature [K]\n10) Pressure [bar]"
echo -e "Your answer (#):${reset}"
read x_axis

echo -e "${magenta}${bold}What would you like to plot on the Y axis?${reset}${magenta}\n\n1) Total Energy of the system [a.u.]\n2) Step Number\n3) Conserved energy [hartree]\n4) CPU time per MD step [s] \n5) Energy drift per atom [K]\n6) Potential energy [hartree]\n7) Kinetic energy [hartree]\n8) Instan. Temperature [K]\n9) Avg. Temperature [K]\n10) Pressure [bar]"
echo -e "Your answer (#):${reset}"
read y_axis

#Grab Information
#x axis
if [ $x_axis == 1 ] ;
then
	grep "Total FORCE_EVAL" runlog.log | awk '{printf $9 "\n"}' > x_axis.txt
	x_label="Total Energy [a.u.]"
fi

if [ $x_axis == 2 ] ;
then
	grep "Step number" runlog.log | awk '{printf $4 "\n"}' > x_axis.txt
	x_label="Step Number"
fi

if [ $x_axis == 3 ] ;
then
        grep "Conserved quantity" runlog.log | awk '{printf $5 "\n"}' > x_axis.txt
	x_label="Conserved Energy [hartree]"
fi

if [ $x_axis == 4 ] ;
then
        grep "CPU time per MD step" runlog.log | awk '{printf $8 "\n"}' > x_axis.txt
	x_label="CPU time per MD step [s]"
fi

if [ $x_axis == 5 ] ;
then
        grep "Energy drift per atom" runlog.log | awk '{printf $7 "\n"}' > x_axis.txt
	x_label="Energy drift per atom [K]"
fi

if [ $x_axis == 6 ] ;
then
        grep "Potential energy" runlog.log | awk '{printf $5 "\n"}' > x_axis.txt
	x_label="Potential Energy [hartree]"
fi

if [ $x_axis == 7 ] ;
then
        grep "Kinetic energy" runlog.log | awk '{printf $5 "\n"}' > x_axis.txt
	x_label="Kinetic Energy [hartree]"
fi

if [ $x_axis == 8 ] ;
then
	grep "Temperature" runlog.log | awk '{printf $4 "\n"}' > x_axis.txt
	x_label="Instantaneous Temperature [K]"
fi

if [ $x_axis == 9 ] ;
then
	grep "Temperature" runlog.log | awk '{printf $5 "\n"}' > x_axis.txt
	x_label="Average Temperature [K]"
fi

if [ $x_axis == 10 ] ;
then
        grep "Pressure" runlog.log | awk '{printf $4 "\n"}' > x_axis.txt
	x_label="Pressure [bar]"
fi

#y_axis
if [ $y_axis == 1 ] ;
then
        grep "Total FORCE_EVAL" runlog.log | awk '{printf $9 "\n"}' > y_axis.txt
        y_label="Total Energy [a.u.]"
fi

if [ $y_axis == 2 ] ;
then
        grep "Step number" runlog.log | awk '{printf $4 "\n"}' > y_axis.txt
	y_label="Step Number"
fi

if [ $y_axis == 3 ] ;
then
        grep "Conserved quantity" runlog.log | awk '{printf $5 "\n"}' > y_axis.txt
        y_label="Conserved Energy [hartree]"
fi

if [ $y_axis == 4 ] ;
then
        grep "CPU time per MD step" runlog.log | awk '{printf $8 "\n"}' > y_axis.txt
        y_label="CPU time per MD step [s]"
fi

if [ $y_axis == 5 ] ;
then
        grep "Energy drift per atom" runlog.log | awk '{printf $7 "\n"}' > y_axis.txt
        y_label="Energy drift per atom [K]"
fi

if [ $y_axis == 6 ] ;
then
        grep "Potential energy" runlog.log | awk '{printf $5 "\n"}' > y_axis.txt
        y_label="Potential Energy [hartree]"
fi

if [ $y_axis == 7 ] ;
then
        grep "Kinetic energy" runlog.log | awk '{printf $5 "\n"}' > y_axis.txt
        y_label="Kinetic Energy [hartree]"
fi

if [ $y_axis == 8 ] ;
then
	grep "Temperature" runlog.log | awk '{printf $4 "\n"}' > y_axis.txt
	y_label="Instantaneous Temperature [K]"
fi

if [ $y_axis == 9 ] ;
then
	grep "Temperature" runlog.log | awk '{printf $5 "\n"}' > y_axis.txt
	y_label="Average Temperature [K]"
fi

if [ $y_axis == 10 ] ;
then
        grep "Pressure" runlog.log | awk '{printf $4 "\n"}' > y_axis.txt
        y_label="Pressure [bar]"
fi

# Combine the two files into a single file with two columns
paste x_axis.txt y_axis.txt > data.txt

# Clean up the temporary files
rm x_axis.txt y_axis.txt

gnuplot <<EOF
set title "${y_label} vs ${x_label}"
set xlabel "${x_label}"
set ylabel "${y_label}"
plot "data.txt" using 1:2 with lines title "${y_label}"
pause mouse
EOF
