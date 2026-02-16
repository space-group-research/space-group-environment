#!/bin/bash
#Script by Angela Shipman, v3 Feb 15 2024

###HOW TO USE
#1. Type 'plot'
#2. Script searches for a runlog.log, stdout.*, *energy.dat file, or allows you to type in a custom name.
#3. Select the file. Script tries to detect if the file is from CP2K or MPMC.
#4(a). If CP2K is detected, script searches for the following job types: Energy Force, Geo Opt, Cell Opt, or MD
#4(b). If MPMC is detected, script discerns if it is an averaged core file (stdout.*) or if it contains individual core information (*.energy.dat). If it is the latter, it will prompt you to choose which core to plot and parse the file accordingly.
#5. Plots your data!

#**Requires Gnuplot

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
#Check for runlog.log or stdout.#
echo -e "${cyan}Searching for ${bold}runlog.log${reset}${cyan} or ${bold}stdout.*${reset}${cyan} files..."

if [ ! -f runlog.log ] ;
then
	if [ ! -f stdout.* ] &> /dev/null ;
	then
		echo -e "${red}No ${bold}runlog.log${reset}${red} or ${bold}stdout.*${reset}${red} file detected. Please ensure there is an output file in this directory.${reset}"
		exit
	fi
fi

#set output_file

list=`ls | grep "stdout./*\|runlog.log\|/*.energy.dat/*"`
list+=`echo -e "\nCUSTOM"`

num_opts=`echo $list | wc -w`
norm=$(( $num_opts - 1 ))

echo -e "${cyan}${bold}Detected the following output files. Please select the one you'd like to plot, or type in a custom filename if you don't see the desired file below.${reset}${cyan}"

select file in $list
do
	if (( $REPLY <= $norm ));
	then
		output_file="$file"
                echo -e "\n${green}${bold}Selected $file.${reset}"
	elif (( $REPLY == $num_opts ));
	then
		others=`ls`
		echo -e "\n${cyan}${bold}All other files in this directory are listed below. Type the entire file name you would like to plot below.${reset}${cyan}"
		for i in $others;
		do
			echo $i
		done
		echo -e "${cyan}${bold}Custom File Name:${reset}"
		read custom_file
		output_file="$custom_file"
		echo -e "\n${green}${bold}Selected $custom_file.${reset}"
	fi
	break
done

#Check for Program
echo -e "\n${magenta}Detecting program...${reset}"
check=`head -n 1 $output_file | awk '{print $1}'`

if [ "$check" = "DBCSR|" ] ;
then
	echo -e "${green}${bold}CP2K${reset}${green} detected.${reset}"
	program=1 #CP2K
elif [ "$check" = "MPMC" ] ;
then
	echo -e "${green}${bold}MPMC${reset}${green} detected.${reset}"
	program=2 #MPMC
elif [ "$check" = "#step" ] ;
then
	echo -e "${green}${bold}MPMC (individual processor data file)${reset}${green} detected.${reset}"
	program=3 #MPMC_individualcores
else
	echo -e "${yellow}${bold}Cannot detect program from default options. Which program is the runlog from?${reset}${yellow}\n\n1) CP2K\n2) MPMC\n3) MPMC (Individual Core File)\nYour answer (#):"
read program
fi


if [ $program == 1 ] ;
then
	if [ -f cutoff.sh ] ; #cut off exists
	then
		echo -e "${cyan}${bold}A CUTOFF script has been detected. Would you like to:${reset}${cyan}\n1) Plot the ${bold}FIRST ${reset}${cyan}energy calculated as a function of CUTOFF Ry\n2) Plot the ${bold}LAST${reset}${cyan} energy calculated as a function of CUTOFF Ry ${bold}(COMPLETED JOBS ONLY)\n${reset}${cyan}3) Plot the ${bold}DIFFERENCE IN FIRST ENERGY${reset}${cyan} as a function of CUTOFF Ry\n4) Plot the ${bold}DIFFERENCE IN LAST ENERGY${reset}${cyan} as a function of CUTOFF Ry ${bold}(COMPLETED JOBS ONLY)${reset}"
		read a1

		if [ $a1 == 1 ] ; #plot first energy
		then
			output_filepaths=`find -name "runlog.log" | sort`
        		echo $output_filepaths | grep -o "[0-9]*" > x_axis.txt
			x_label="Cutoff Energy (Ry)"
        		grep "     1 OT" $output_filepaths | awk '{print $8}' > y_axis.txt
			y_label="First Energy [a.u.]"

		elif [ $a1 == 2 ] ; #plot last energy
		then
			output_filepaths=`find -name "runlog.log" | sort`
			echo $output_filepaths | grep -o "[0-9]*" > x_axis.txt
			x_label="Cutoff Energy (Ry)"
			grep -B 2 "*** SCF run converged in" $output_filepaths | awk '{print $8}' | awk 'NR%4==1 || NR%4==2' | awk 'NR%2==1' > y_axis.txt
			y_label="Last Energy [a.u.]"

		elif [ $a1 == 3 ] ;
		then
			output_filepaths=`find -name "runlog.log" | sort`
                        echo $output_filepaths | grep -o "[0-9]*" > x_axis.txt
                        x_label="Cutoff Energy (Ry)"
                        grep "     1 OT" $output_filepaths | awk '{print $8}' > y_axis0.txt
			awk 's{print ((s>$0)?s-$0:$0-s)}{s=$0}' y_axis0.txt > y_axis.txt
			rm y_axis0.txt
                        y_label="Difference in First Energy [a.u.]"

		elif [ $a1 == 4 ] ;
		then
			output_filepaths=`find -name "runlog.log" | sort`
                        echo $output_filepaths | grep -o "[0-9]*" > x_axis.txt
                        x_label="Cutoff Energy (Ry)"
                        grep -B 2 "*** SCF run converged in" $output_filepaths | awk '{print $8}' | awk 'NR%4==1 || NR%4==2' | awk 'NR%2==1' > y_axis0.txt
			awk 's{print ((s>$0)?s-$0:$0-s)}{s=$0}' y_axis0.txt > y_axis.txt
			rm y_axis0.txt
                        y_label="Difference in Last Energy [a.u.]"

		else
			echo "${red}Invalid answer. Enter a number between 1-4."
		fi
	elif [ -f rel_cutoff.sh ] ; #rel cut off exists
	then
                echo -e "${cyan}${bold}A REL CUTOFF script has been detected. Would you like to:${reset}${cyan}\n1) Plot the ${bold}FIRST${reset}${cyan} energy calculated as a function of REL CUTOFF Ry\n2) Plot the ${bold}LAST${reset}${cyan} energy calculated as a function of REL CUTOFF Ry ${bold}(COMPLETED JOBS ONLY)\n${reset}${cyan}3) Plot the ${bold}DIFFERENCE IN FIRST ENERGY${reset}${cyan} as a function of REL CUTOFF Ry\n4) Plot the ${bold}DIFFERENCE IN LAST ENERGY${reset}${cyan} as a function of REL CUTOFF Ry ${bold}(COMPLETED JOBS ONLY)${reset}"
		read a1

                if [ $a1 == 1 ] ; #plot first energy
                then
                        output_filepaths=`find -name "runlog.log" | sort`
			organize=`echo "$output_filepaths" | sort -n -t / -k 2,2 -k 3,3`
                        echo $organize | grep -o "[0-9]*" | sort -n > x_axis.txt
                        x_label="Rel Cutoff Energy (Ry)"
                        grep "     1 OT" $organize | awk '{print $8}' > y_axis.txt
                        y_label="First Energy [a.u.]"

                elif [ $a1 == 2 ] ; #plot last energy
                then
                        output_filepaths=`find -name "runlog.log" | sort`
			organize=`echo "$output_filepaths" | sort -n -t / -k 2,2 -k 3,3`
                        echo $organize | grep -o "[0-9]*" | sort -n > x_axis.txt
                        x_label="Rel Cutoff Energy (Ry)"
                        grep -B 2 "*** SCF run converged in" $organize | awk '{print $8}' | awk 'NR%4==1 || NR%4==2' | awk 'NR%2==1' > y_axis.txt
                        y_label="Last Energy [a.u.]"
		elif [ $a1 == 3 ] ;
		then
			output_filepaths=`find -name "runlog.log" | sort`
			organize=`echo "$output_filepaths" | sort -n -t / -k 2,2 -k 3,3`
                        echo $organize | grep -o "[0-9]*" | sort -n > x_axis.txt
                        x_label="Rel Cutoff Energy (Ry)"
                        grep "     1 OT" $organize | awk '{print $8}' > y_axis0.txt
			awk 's{print ((s>$0)?s-$0:$0-s)}{s=$0}' y_axis0.txt > y_axis.txt
                        rm y_axis0.txt
                        y_label="Difference in First Energy [a.u.]"

		elif [ $a1 == 4 ] ;
		then
			output_filepaths=`find -name "runlog.log" | sort`
			organize=`echo "$output_filepaths" | sort -n -t / -k 2,2 -k 3,3`
                        echo $organize | grep -o "[0-9]*" | sort -n > x_axis.txt
                        x_label="Rel Cutoff Energy (Ry)"
                        grep -B 2 "*** SCF run converged in" $organize | awk '{print $8}' | awk 'NR%4==1 || NR%4==2' | awk 'NR%2==1' > y_axis0.txt
                        awk 's{print ((s>$0)?s-$0:$0-s)}{s=$0}' y_axis0.txt > y_axis.txt
                        rm y_axis0.txt
			y_label="Difference in Last Energy [a.u.]"
	       	
		else
                        echo "${red}Invalid answer. Enter a number between 1-4."
                fi

	else #detect job type
		job_type=`grep "Run type" runlog.log | tail -1 | awk '{print $4}'`
		
		if [ "$job_type" == "ENERGY_FORCE" ] ;
		then
			echo "${bold}This is an ENERGY FORCE job on CP2K.${reset}"

			echo -e "${green}${bold}What would you like to plot on the X axis?\n\n${reset}${green}1) Total Energy of the system [a.u.]\n2) Step Number\n3) Time of OT Step [sec]\n4) Convergence Energy [a.u.] \n5) Change in Energy per OT Step"
			echo -e "Your answer (#):${reset}"
			read x_axis

			echo -e "${yellow}${bold}What would you like to plot on the Y axis?${reset}${yellow}\n\n1) Total Energy of the system [a.u.]\n2) Step Number\n3) Time of OT Step [sec]\n4) Convergence Energy [a.u.] \n5) Change in Energy per OT Step"
			echo -e "Your answer (#):${reset}"
			read y_axis

			#x axis
			if [ $x_axis == 1 ] ;
			then
				cp runlog.log tmp00.txt
				sed -i "s/                    /     0.00000000     /g" tmp00.txt
				grep " OT " tmp00.txt | awk '{print $7}' | tail -n+3 > x_axis.txt
				x_label="Total Energy [a.u.]"
				rm tmp00.txt

			elif [ $x_axis == 2 ] ;
			then
				grep " OT " runlog.log | awk '{print $1}' | tail -n+3 > x_axis.txt
				x_label="Step Number"
			

			elif [ $x_axis == 3 ] ;
			then
				grep " OT " runlog.log | awk '{print $5}' | tail -n+3 > x_axis.txt
				x_label="Time of OT Step [sec]"
			

			elif [ $x_axis == 4 ] ;
			then
				cp runlog.log tmp00.txt
				sed -i "s/                    /     0.00000000     /g" tmp00.txt
				grep " OT " tmp00.txt | awk '{print $6}' | tail -n+3 > x_axis.txt
				rm tmp00.txt
				x_label="Convergence Energy [a.u.]"
			

			elif [ $x_axis == 5 ] ;
			then
				cp runlog.log tmp00.txt
				sed -i "s/                    /     0.00000000     /g" tmp00.txt
				grep " OT " tmp00.txt | awk '{print $8}' | tail -n+3 > x_axis.txt
				rm tmp00.txt
				x_label="Change in Energy [a.u.]"
			else
				echo "${red}Invalid entry. Type a number between 1-5.${reset}"
				break
			fi
			
			#y axis

			if [ $y_axis == 1 ] ;
			then
				cp runlog.log tmp00.txt
				sed -i "s/                    /     0.00000000     /g" tmp00.txt
				grep " OT " tmp00.txt | awk '{print $7}' | tail -n+3 > y_axis.txt
				y_label="Total Energy [a.u.]"
				rm tmp00.txt

			elif [ $y_axis == 2 ] ;
			then
				grep " OT " runlog.log | awk '{print $1}' | tail -n+3 > y_axis.txt
				y_label="Step Number"
			

			elif [ $y_axis == 3 ] ;
			then
				grep " OT " runlog.log | awk '{print $5}' | tail -n+3> y_axis.txt
				y_label="Time of OT Step [sec]"
			

			elif [ $y_axis == 4 ] ;
			then
				cp runlog.log tmp00.txt
				sed -i "s/                    /     0.00000000     /g" tmp00.txt
				grep " OT " tmp00.txt | awk '{print $6}' | tail -n+3 > y_axis.txt
				sed -i "s/0.00000000/ /g" y_axis.txt
				rm tmp00.txt
				y_label="Convergence Energy [a.u.]"


			elif [ $y_axis == 5 ] ;
			then
				cp runlog.log tmp00.txt
				sed -i "s/                    /     0.00000000     /g" tmp00.txt
				grep " OT " tmp00.txt | awk '{print $8}' | tail -n+3 > y_axis.txt
				rm tmp00.txt
				y_label="Change in Energy [a.u.]"
			else
				echo "${red}Invalid entry. Enter a number between 1-5.${reset}"
			fi
		
		elif [ "$job_type" == "GEO_OPT" ] ;
		then
			echo "${bold}This is a GEO OPT job on CP2K.${reset}"

			echo -e "${green}${bold}What would you like to plot on the X axis?\n\n${reset}${green}1) Total Energy (from FORCE EVAL) [a.u.]\n2) Step Number\n3) Used Time [sec]\n4) Real Energy Change [a.u.]"
			echo -e "Your answer (#):${reset}"
			read x_axis

			echo -e "${yellow}${bold}What would you like to plot on the Y axis?${reset}${yellow}\n\n1) Total Energy (from FORCE EVAL) [a.u.]\n2) Step Number\n3) Used Time [sec]\n4) Real Energy Change [a.u.]"
			echo -e "Your answer (#):${reset}"
			read y_axis

			#x axis
			if [ $x_axis == 1 ] ; #total energy
			then
				grep "ENERGY| Total FORCE_EVAL" runlog.log | awk '{print $9}' > x_axis.txt
				x_label="Total Energy [a.u.]"	
			elif [ $x_axis == 2 ] ; #step number
			then
				grep "Informations at step" runlog.log | awk '{print $6}' > x_axis.txt
				x_label="Step Number"
			elif [ $x_axis == 3 ] ; #used time
			then
				grep "Used time" runlog.log | awk '{print $4}' > x_axis.txt
				x_label="Used Time [sec]"
			elif [ $x_axis == 4 ] ; #real energy change
			then
				grep "Real energy change" runlog.log | awk '{print $5}' > x_axis.txt
				x_label="Real Energy Change [a.u.]"
			fi

			#y axis
			if [ $y_axis == 1 ] ; #total energy
			then
				grep "ENERGY| Total FORCE_EVAL" runlog.log | awk '{print $9}' > y_axis.txt
				y_label="Total Energy [a.u.]"
			elif [ $y_axis == 2 ] ; #step number
			then
				grep "Informations at step" runlog.log | awk '{print $6}' > y_axis.txt
				y_label="Step Number"
			elif [ $y_axis == 3 ] ; #used time
			then
				grep "Used time" runlog.log | awk '{print $4}' > y_axis.txt
				y_label="Used Time [sec]"
			elif [ $y_axis == 4 ] ; #real energy change
			then
				grep "Real energy change" runlog.log | awk '{print $5}' > y_axis.txt
				y_label="Real Energy Change [a.u.]"
			fi

		elif [ "$job_type" == "CELL_OPT" ] ;
		then
			echo "${bold}This is a CELL OPT job on CP2K.${reset}"

			#Plot: Step number, total energy FORCE_EVAL [a.u.], Lattice lengths (a, b, c), lattice angles (alp, bet, gam), volume, STEP FOR MORE UNCOMMON THINGS, Internal Pressure [bar], Real energy change, Used Time, 

			echo -e "${green}${bold}What would you like to plot on the X axis?\n\n${reset}${green}1) Total Energy (from FORCE EVAL) [a.u.]\n2) Step Number\n3) Lattice Length a [Angstroms] \n4) Lattice Length b [Angstroms]\n5) Lattice Length c [Angstroms]\n6) Lattice Angle Alpha [degrees]\n7) Lattice Angle Beta [degrees]\n8) Lattice Angle Gamma [degrees]\n9) Volume [Angstroms^3]\n10) Used Time [sec]\n11) Real Energy Change [a.u.]\n12) Internal Pressure [bar]"
			echo -e "Your answer (#):${reset}"
			read x_axis

			echo -e "${yellow}${bold}What would you like to plot on the Y axis?${reset}${yellow}\n\n1) Total Energy (from FORCE EVAL) [a.u.]\n2) Step Number\n3) Lattice Length a [Angstroms] \n4) Lattice Length b [Angstroms]\n5) Lattice Length c [Angstroms]\n6) Lattice Angle Alpha [degrees]\n7) Lattice Angle Beta [degrees]\n8) Lattice Angle Gamma [degrees]\n9) Volume [Angstroms^3]\n10) Used Time [sec]\n11) Real Energy Change [a.u.]\n12) Internal Pressure [bar]"
			echo -e "Your answer (#):${reset}"
			read y_axis

			#x axis
			if [ $x_axis == 1 ] ; #total energy
			then
				grep "Total Energy               =" runlog.log | awk '{print $4}' > x_axis.txt
				x_label="Total Energy [a.u.]"	
			elif [ $x_axis == 2 ] ; #step number
			then
				grep "Informations at step" runlog.log | awk '{print $6}' > x_axis.txt
				x_label="Step Number"
			elif [ $x_axis == 3 ] ; #len a
			then
				grep "Vector a" runlog.log | awk '{print $10}' | head -n -1 > x_axis.txt
				x_label="Lattice Length a [Angstroms]"
			elif [ $x_axis == 4 ] ; #len b
			then
				grep "Vector b" runlog.log | awk '{print $10}' | head -n -1 > x_axis.txt
				x_label="Lattice Length b [Angstroms]"
			elif [ $x_axis == 5 ] ; #len c
			then
				grep "Vector c" runlog.log | awk '{print $10}' | head -n -1 > x_axis.txt
				x_label="Lattice Length c [Angstroms]"
			elif [ $x_axis == 6 ] ; #alpha
			then
				grep "CELL| Angle (b,c), alpha" runlog.log | awk '{print $6}' | head -n -1 > x_axis.txt
				x_label="Lattice Angle Alpha [degrees]"
			elif [ $x_axis == 7 ] ; #beta
			then
				grep "CELL| Angle (a,c), beta" runlog.log | awk '{print $6}' | head -n -1 > x_axis.txt
				x_label="Lattice Angle Beta [degrees]"
			elif [ $x_axis == 8 ] ; #gamma
			then
				grep "CELL| Angle (a,b), gamma" runlog.log | awk '{print $6}' | head -n -1 > x_axis.txt
				x_label="Lattice Angle Gamma [degrees]"
			elif [ $x_axis == 9 ] ; #volume
			then
				grep "Volume" runlog.log | awk '{print $4}' | head -n -1 > x_axis.txt
				x_label="Volume [Angstroms^3]"
			elif [ $x_axis == 10 ] ; #used time
			then
				grep "Used time" runlog.log | awk '{print $4}' > x_axis.txt
				x_label="Used Time [sec]"
			elif [ $x_axis == 11 ] ; #real energy change
			then
				grep "Real energy change" runlog.log | awk '{print $5}' > x_axis.txt
				x_label="Real Energy Change [a.u.]"
			elif [ $x_axis == 12 ] ; #internal pressure
			then
				grep "Internal Pressure" runlog.log | awk '{print $5}' > x_axis.txt
				x_label="Internal Pressure [bar]"
			fi

			#y axis
			if [ $y_axis == 1 ] ; #total energy
			then
				grep "Total Energy               =" runlog.log | awk '{print $4}' > y_axis.txt
				y_label="Total Energy [a.u.]"
			elif [ $y_axis == 2 ] ; #step number
			then
				grep "Informations at step" runlog.log | awk '{print $6}' > y_axis.txt
				y_label="Step Number"
			elif [ $y_axis == 3 ] ; #len a
			then
				grep "Vector a" runlog.log | awk '{print $10}' | head -n -1 > y_axis.txt
				y_label="Lattice Length a [Angstroms]"
			elif [ $y_axis == 4 ] ; #len b
			then
				grep "Vector b" runlog.log | awk '{print $10}' | head -n -1 > y_axis.txt
				y_label="Lattice Length b [Angstroms]"
			elif [ $y_axis == 5 ] ; #len c
			then
				grep "Vector c" runlog.log | awk '{print $10}' | head -n -1 > y_axis.txt
				y_label="Lattice Length c [Angstroms]"
			elif [ $y_axis == 6 ] ; #alpha
			then
				grep "CELL| Angle (b,c), alpha" runlog.log | awk '{print $6}' | head -n -1 > y_axis.txt
				y_label="Lattice Angle Alpha [degrees]"
			elif [ $y_axis == 7 ] ; #beta
			then
				grep "CELL| Angle (a,c), beta" runlog.log | awk '{print $6}' | head -n -1 > y_axis.txt
				y_label="Lattice Angle Beta [degrees]"
			elif [ $y_axis == 8 ] ; #gamma
			then
				grep "CELL| Angle (a,b), gamma" runlog.log | awk '{print $6}' | head -n -1 > y_axis.txt
				y_label="Lattice Angle Gamma [degrees]"
			elif [ $y_axis == 9 ] ; #volume
			then
				grep "Volume" runlog.log | awk '{print $4}' | head -n -1 > y_axis.txt
				y_label="Volume [Angstroms^3]"
			elif [ $y_axis == 10 ] ; #used time
			then
				grep "Used time" runlog.log | awk '{print $4}' > y_axis.txt
				y_label="Used Time [sec]"
			elif [ $y_axis == 11 ] ; #real energy change
			then
				grep "Real energy change" runlog.log | awk '{print $5}' > y_axis.txt
				y_label="Real Energy Change [a.u.]"
			elif [ $y_axis == 12 ] ; #internal pressure
			then
				grep "Internal Pressure" runlog.log | awk '{print $5}' > y_axis.txt
				y_label="Internal Pressure [bar]"
			fi
		
		elif [ "$job_type" == "MD" ] ;
		then
			echo "${bold}This is an MD job on CP2K."

			echo -e "${green}${bold}What would you like to plot on the X axis?\n\n${reset}${green}1) Total Energy of the system [a.u.]\n2) Step Number\n3) Conserved energy [hartree]\n4) CPU time per MD step [s] \n5) Energy drift per atom [K]\n6) Potential energy [hartree]\n7) Kinetic energy [hartree]\n8) Instan. Temperature [K]\n9) Avg. Temperature [K]\n10) Pressure [bar]"
			echo -e "Your answer (#):${reset}"
			read x_axis

			echo -e "${yellow}${bold}What would you like to plot on the Y axis?${reset}${yellow}\n\n1) Total Energy of the system [a.u.]\n2) Step Number\n3) Conserved energy [hartree]\n4) CPU time per MD step [s] \n5) Energy drift per atom [K]\n6) Potential energy [hartree]\n7) Kinetic energy [hartree]\n8) Instan. Temperature [K]\n9) Avg. Temperature [K]\n10) Pressure [bar]"
			echo -e "Your answer (#):${reset}"
			read y_axis

			#x axis
			if [ $x_axis == 1 ] ;
			then
				grep -B 3 "Step number" runlog.log | grep "ENERGY| Total FORCE_EVAL" | awk '{printf $9 "\n"}' > x_axis.txt
				x_label="Total Energy [a.u.]"

			elif [ $x_axis == 2 ] ;
			then
				grep "Step number" runlog.log | awk '{printf $4 "\n"}' > x_axis.txt
				x_label="Step Number"

			elif [ $x_axis == 3 ] ;
			then
				grep "Conserved quantity" runlog.log | awk '{printf $5 "\n"}' > x_axis.txt
				x_label="Conserved Energy [hartree]"

			elif [ $x_axis == 4 ] ;
			then
				grep "CPU time per MD step" runlog.log | awk '{printf $8 "\n"}' > x_axis.txt
				x_label="CPU time per MD step [s]"

			elif [ $x_axis == 5 ] ;
			then
				grep "Energy drift per atom" runlog.log | awk '{printf $7 "\n"}' > x_axis.txt
				x_label="Energy drift per atom [K]"

			elif [ $x_axis == 6 ] ;
			then
				grep "Potential energy" runlog.log | grep -v "MD_INI" | grep -v "MD_PAR" | awk '{printf $5 "\n"}' > x_axis.txt
				x_label="Potential Energy [hartree]"

			elif [ $x_axis == 7 ] ;
			then
				grep "Kinetic energy" runlog.log | grep -v "MD_INI" | grep -v "MD_PAR" | awk '{printf $5 "\n"}' > x_axis.txt
				x_label="Kinetic Energy [hartree]"

			elif [ $x_axis == 8 ] ;
			then
				grep "Temperature" runlog.log | grep -v "MD_INI" | grep -v "MD_PAR" | awk '{printf $4 "\n"}' > x_axis.txt
				x_label="Instantaneous Temperature [K]"

			elif [ $x_axis == 9 ] ;
			then
				grep "Temperature" runlog.log | grep -v "MD_INI" | grep -v "MD_PAR" | awk '{printf $5 "\n"}' > x_axis.txt
				x_label="Average Temperature [K]"

			elif [ $x_axis == 10 ] ;
			then
				grep "Pressure" runlog.log | grep -v "MD_INI" | grep -v "MD_PAR" | awk '{printf $4 "\n"}' > x_axis.txt
				x_label="Pressure [bar]"
			fi #close x axis generation

			#y_axis
			if [ $y_axis == 1 ] ;
			then
				grep -B 3 "Step number" runlog.log | grep "ENERGY| Total FORCE_EVAL" | awk '{printf $9 "\n"}' > y_axis.txt
				y_label="Total Energy [a.u.]"

			elif [ $y_axis == 2 ] ;
			then
				grep "Step number" runlog.log | awk '{printf $4 "\n"}' > y_axis.txt
				y_label="Step Number"

			elif [ $y_axis == 3 ] ;
			then
				grep "Conserved quantity" runlog.log | awk '{printf $5 "\n"}' > y_axis.txt
				y_label="Conserved Energy [hartree]"

			elif [ $y_axis == 4 ] ;
			then
				grep "CPU time per MD step" runlog.log | awk '{printf $8 "\n"}' > y_axis.txt
				y_label="CPU time per MD step [s]"

			elif [ $y_axis == 5 ] ;
			then
				grep "Energy drift per atom" runlog.log | awk '{printf $7 "\n"}' > y_axis.txt
				y_label="Energy drift per atom [K]"

			elif [ $y_axis == 6 ] ;
			then
				grep "Potential energy" runlog.log | grep -v "MD_INI" | grep -v "MD_PAR" | awk '{printf $5 "\n"}' > y_axis.txt
				y_label="Potential Energy [hartree]"

			elif [ $y_axis == 7 ] ;
			then
				grep "Kinetic energy" runlog.log | grep -v "MD_INI" | grep -v "MD_PAR" | awk '{printf $5 "\n"}' > y_axis.txt
				y_label="Kinetic Energy [hartree]"

			elif [ $y_axis == 8 ] ;
			then
				grep "Temperature" runlog.log | grep -v "MD_INI" | grep -v "MD_PAR" | awk '{printf $4 "\n"}' > y_axis.txt
				y_label="Instantaneous Temperature [K]"

			elif [ $y_axis == 9 ] ;
			then
				grep "Temperature" runlog.log | grep -v "MD_INI" | grep -v "MD_PAR" | awk '{printf $5 "\n"}' > y_axis.txt
				y_label="Average Temperature [K]"

			elif [ $y_axis == 10 ] ;
			then
				grep "Pressure" runlog.log | grep -v "MD_INI" | grep -v "MD_PAR" | awk '{printf $4 "\n"}' > y_axis.txt
				y_label="Pressure [bar]"
			fi #close y axis generation
		else
			echo "Unsupported CP2K job type."
		fi #close job type actions
	fi  #close cp2k
#---------------------------------------------------------------------MPMC
elif [ $program == 2 ] ; #MPMC
then

	echo -e "${green}${bold}What would you like to plot on the X axis?\n\n${reset}${green}1) Step Number\n2) Boltzmann Factor (BF)\n3) Acceptance Ratio (AR)\n4) Potential Energy [K]\n5) Electrostatic Energy [K]\n6) Repulsion/Dispersion Energy [K]\n7) N (# sorbed molecules in system)\n8) Density of sorbed gas in material [g/cm^3]\n9) wt %\n10) wt % (ME)\n11) qst [kJ/mol]\n12) Energy/particle [kJ/mol]\n13) Compressibility [atm^-1]\n14) Bulk Modulus [GPa]\n15) Heat Capacity [kJ/mol K]\n16) Para Spin Ratio [%]"
        echo -e "Your answer (#):${reset}"
        read x_axis

        echo -e "${yellow}${bold}What would you like to plot on the Y axis?${reset}${yellow}\n\n1) Step Number\n2) Boltzmann Factor (BF)\n3) Acceptance Ratio (AR)\n4) Potential Energy [K]\n5) Electrostatic Energy [K]\n6) Repulsion/Dispersion Energy [K]\n7) N (# sorbed molecules in system)\n8) Density of sorbed gas in material [g/cm^3]\n9) wt %\n10) wt % (ME)\n11) qst [kJ/mol]\n12) Energy/particle [kJ/mol]\n13) Compressibility [atm^-1]\n14) Bulk Modulus [GPa]\n15) Heat Capacity [kJ/mol K]\n16) Para Spin Ratio [%]"
        echo -e "Your answer (#):${reset}"
        read y_axis

	#x-axis
	if [ $x_axis == 1 ] ;
                        then
                                grep "Completed step" $output_file | awk '{print $4}' | awk -F"/" '{print $1}' > x_axis.txt
				x_label="Step Number"

                        elif [ $x_axis == 2 ] ;
                        then
                                grep "BF" $output_file | tail -n +2 | awk '{print $4}' > x_axis.txt
				x_label="Boltzmann Factor"

                        elif [ $x_axis == 3 ] ;
                        then
                                grep "AR" $output_file | tail -n +2 | awk '{print $4}' > x_axis.txt
				x_label="Acceptance Ratio"

                        elif [ $x_axis == 4 ] ;
                        then
                                grep "potential energy" $output_file | tail -n +2 | awk '{print $5}' > x_axis.txt
				x_label="Potential Energy [K]"

                        elif [ $x_axis == 5 ] ;
                        then
                                grep "electrostatic energy" $output_file | tail -n +2 | awk '{print $5}' > x_axis.txt
				x_label="Electrostatic Energy [K]"

                        elif [ $x_axis == 6 ] ;
                        then
                                grep "repulsion/dispersion energy" $output_file | tail -n +3 | awk '{print $5}' > x_axis.txt
				x_label="Repulsion/Dispersion Energy [K]"

                        elif [ $x_axis == 7 ] ;
                        then
                                grep "N =" $output_file | tail -n +3 | awk '{print $4}' > x_axis.txt
				x_label="N [# sorbed molecules in system]"

                        elif [ $x_axis == 8 ] ;
                        then
				grep "density" $output_file | tail -n +3 | awk '{print $4}' > x_axis.txt
                                x_label="Density of Sorbed Molecules [g/cm^3]"

                        elif [ $x_axis == 9 ] ;
                        then
				grep "wt % =" $output_file | tail -n +3 | awk '{print $5}' > x_axis.txt
                                x_label="wt %"

                        elif [ $x_axis == 10 ] ;
                        then
				grep "wt % (ME) =" $output_file | tail -n +3 | awk '{print $6}' > x_axis.txt
				x_label="wt % (ME)"

			elif [ $x_axis == 11 ] ;
                        then
				grep "qst" $output_file | tail -n +2 | awk '{print $4}' > x_axis.txt
				x_label="qst [kJ/mol]"

			elif [ $x_axis == 12 ] ;
                        then
                                grep "Energy/particle" $output_file | tail -n +2 | awk '{print $4}' > x_axis.txt
				x_label="Energy/particle [kJ/mol]"

			elif [ $x_axis == 13 ] ;
                        then
                                grep "compressibility" $output_file | tail -n +2 | awk '{print $4}' > x_axis.txt
				x_label="Compressibility [atm^(-1)]"

			elif [ $x_axis == 14 ] ;
                        then
                                grep "bulk modulus" $output_file | tail -n +2 | awk '{print $5}' > x_axis.txt
				x_label="Bulk Modulus [GPa]"

			elif [ $x_axis == 15 ] ;
                        then
                                grep "heat capacity" $output_file | tail -n +2 | awk '{print $5}' > x_axis.txt
				x_label="Heat Capacity [kJ/mol K]"

			elif [ $x_axis == 16 ] ;
                        then
				grep "para spin ratio" $output_file | tail -n +3 | awk '{print $6}' > x_axis.txt
                                x_label="Para Spin Ratio [%]"

                        fi #close x axis generation


        #y-axis
        if [ $y_axis == 1 ] ;
                        then
                                grep "Completed step" $output_file | awk '{print $4}' | awk -F"/" '{print $1}' > y_axis.txt
                                y_label="Step Number"

                        elif [ $y_axis == 2 ] ;
                        then
                                grep "BF" $output_file | tail -n +2 | awk '{print $4}' > y_axis.txt
                                y_label="Boltzmann Factor"

                        elif [ $y_axis == 3 ] ;
                        then
                                grep "AR" $output_file | tail -n +2 | awk '{print $4}' > y_axis.txt
                                y_label="Acceptance Ratio"

                        elif [ $y_axis == 4 ] ;
                        then
                                grep "potential energy" $output_file | tail -n +2 | awk '{print $5}' > y_axis.txt
                                y_label="Potential Energy [K]"

                        elif [ $y_axis == 5 ] ;
                        then
                                grep "electrostatic energy" $output_file | tail -n +2 | awk '{print $5}' > y_axis.txt
                                y_label="Electrostatic Energy [K]"

                        elif [ $y_axis == 6 ] ;
                        then
                                grep "repulsion/dispersion energy" $output_file | tail -n +3 | awk '{print $5}' > y_axis.txt
                                y_label="Repulsion/Dispersion Energy [K]"

                        elif [ $y_axis == 7 ] ;
                        then
                                grep "N =" $output_file | tail -n +3 | awk '{print $4}' > y_axis.txt
                                y_label="N [# sorbed molecules in system]"

                        elif [ $y_axis == 8 ] ;
                        then
                                grep "density" $output_file | tail -n +3 | awk '{print $4}' > y_axis.txt
                                y_label="Density of Sorbed Molecules [g/cm^3]"

                        elif [ $y_axis == 9 ] ;
                        then
                                grep "wt % =" $output_file | tail -n +3 | awk '{print $5}' > y_axis.txt
                                y_label="wt %"

                        elif [ $y_axis == 10 ] ;
                        then
                                grep "wt % (ME) =" $output_file | tail -n +3 | awk '{print $6}' > y_axis.txt
                                y_label="wt % (ME)"

                        elif [ $y_axis == 11 ] ;
                        then
                                grep "qst" $output_file | tail -n +2 | awk '{print $4}' > y_axis.txt
                                y_label="qst [kJ/mol]"

                        elif [ $y_axis == 12 ] ;
                        then
                                grep "Energy/particle" $output_file | tail -n +2 | awk '{print $4}' > y_axis.txt
                                y_label="Energy/particle [kJ/mol]"

                        elif [ $y_axis == 13 ] ;
                        then
                                grep "compressibility" $output_file | tail -n +2 | awk '{print $4}' > y_axis.txt
                                y_label="Compressibility [atm^(-1)]"

                        elif [ $y_axis == 14 ] ;
                        then
                                grep "bulk modulus" $output_file | tail -n +2 | awk '{print $5}' > y_axis.txt
                                y_label="Bulk Modulus [GPa]"

                        elif [ $y_axis == 15 ] ;
                        then
                                grep "heat capacity" $output_file | tail -n +2 | awk '{print $5}' > y_axis.txt
				y_label="Heat Capacity [kJ/mol K]"

                        elif [ $y_axis == 16 ] ;
                        then
                                grep "para spin ratio" $output_file | tail -n +3 | awk '{print $6}' > y_axis.txt
                                y_label="Para Spin Ratio [%]"

                        fi #close x axis generation

			

#---------------------------------------------------------------------MPMC (INDIVIDUAL CORE FILE)
elif [ $program == 3 ] ; #MPMC INDIVIDUAL CORE FILE
then
	#determine number of cores
	num_cores=`cat $output_file | awk '{print $1}' | tr " " "\n" | sort | uniq -c | sort -k2nr | awk '{printf("%s\t%s\n",$2,$1)}END{print}' | head -n -3 | awk '{print $2}' | uniq -c | awk '{print $2}'`

	if (( $num_cores > 1 ));
	then
		core_list=$(for i in $(seq 1 $num_cores);
		do
			echo "Core-#$i" ; 
		done)
	
		echo -e "${cyan}Which core would you like to plot data from?${reset}"
		select core in $core_list
		do
			core_sel=$REPLY
			echo "You have selected Core #$core_sel."
			break
		done
		
		subtitle="Core #${core_sel}"
		
		#skip every BLAH number of lines for specific core
		touch tmp1_ams2024.txt
		space=$(( $core_sel + 2 ))
		remaining=`tail -n +${space} $output_file >> tmp1_ams2024.txt`
	
		skip=$num_cores	
		parsed=`awk -v NUM=$skip 'NR % NUM == 1' tmp1_ams2024.txt >> tmp2_ams2024.txt`
		output_file=tmp2_ams2024.txt
	else
		echo -e "Only 1 core detected."
		remaining=`tail -n +2 $output_file >> tmp1_ams2024.txt`
		output_file=tmp1_ams2024.txt	
	fi #end core check
	#generate x y data
	echo -e "${green}${bold}What would you like to plot on the X axis?\n\n${reset}${green}1) Step Number\n2) Total Energy [K]\n3) Coulombic Energy [K]\n4) RD\n5) Polarization Contribution [K]\n6) VDW Contribution [K]\n7) Kinetic Energy [K]\n8) Kinetic Temperature [K]\n9) N [# molecules adsorbed] \n10) Spin Ratio [%]\n11) Volume [Angstroms^3]\n12) Core Temperature [K]"
        echo -e "Your answer (#):${reset}"
        read x_axis

        echo -e "${yellow}${bold}What would you like to plot on the Y axis?${reset}${yellow}\n\n1) Step Number\n2) Total Energy [K]\n3) Coulombic Energy [K]\n4) RD\n5) Polarization Contribution [K]\n6) VDW Contribution [K]\n7) Kinetic Energy [K]\n8) Kinetic Temperature [K]\n9) N [# molecules adsorbed] \n10) Spin Ratio [%]\n11) Volume [Angstroms^3]\n12) Core Temperature [K]"
        echo -e "Your answer (#):${reset}"
        read y_axis

	
	if [ $x_axis == 1 ] ;
                        then
                                cat $output_file | awk '{print $1}' > x_axis.txt
				x_label="Step Number"

                        elif [ $x_axis == 2 ] ;
                        then
                                cat $output_file | awk '{print $2}' > x_axis.txt
				x_label="Energy [K]"

                        elif [ $x_axis == 3 ] ;
                        then
                                cat $output_file | awk '{print $3}' > x_axis.txt
				x_label="Coulombic Energy [K]"

                        elif [ $x_axis == 4 ] ;
                        then
                                cat $output_file | awk '{print $4}' > x_axis.txt
				x_label="RD"

                        elif [ $x_axis == 5 ] ;
                        then
                                cat $output_file | awk '{print $5}' > x_axis.txt
				x_label="Polarization Contribution [K]"

                        elif [ $x_axis == 6 ] ;
                        then
                                cat $output_file | awk '{print $6}' > x_axis.txt
				x_label="VDW Contribution [K]"

                        elif [ $x_axis == 7 ] ;
                        then
                                cat $output_file | awk '{print $7}' > x_axis.txt
				x_label="Kinetic Energy [K]"

                        elif [ $x_axis == 8 ] ;
                        then
                                cat $output_file | awk '{print $8}' > x_axis.txt
				x_label="Kinetic Temperature [K]"

                        elif [ $x_axis == 9 ] ;
                        then
                                cat $output_file | awk '{print $9}' > x_axis.txt
				x_label="N [# molecules adsorbed]"

                        elif [ $x_axis == 10 ] ;
                        then
                                cat $output_file | awk '{print $10}' > x_axis.txt
				x_label="Spin Ratio [%]"

                        elif [ $x_axis == 11 ] ;
                        then
                                cat $output_file | awk '{print $11}' > x_axis.txt
				x_label="Volume [Angstroms^3]"

                        elif [ $x_axis == 12 ] ;
                        then
                                cat $output_file | awk '{print $12}' > x_axis.txt
				x_label="Core Temperature [K]"

	fi #end x-axis

	if [ $y_axis == 1 ] ;
                        then
                                cat $output_file | awk '{print $1}' > y_axis.txt
                                y_label="Step Number"

                        elif [ $y_axis == 2 ] ;
                        then
                                cat $output_file | awk '{print $2}' > y_axis.txt
                                y_label="Energy [K]"

                        elif [ $y_axis == 3 ] ;
                        then
                                cat $output_file | awk '{print $3}' > y_axis.txt
                                y_label="Coulombic Energy [K]"

                        elif [ $y_axis == 4 ] ;
                        then
                                cat $output_file | awk '{print $4}' > y_axis.txt
                                y_label="RD"

                        elif [ $y_axis == 5 ] ;
                        then
                                cat $output_file | awk '{print $5}' > y_axis.txt
                                y_label="Polarization Contribution [K]"

                        elif [ $y_axis == 6 ] ;
                        then
                                cat $output_file | awk '{print $6}' > y_axis.txt
                                y_label="VDW Contribution [K]"

                        elif [ $y_axis == 7 ] ;
                        then
                                cat $output_file | awk '{print $7}' > y_axis.txt
                                y_label="Kinetic Energy [K]"

                        elif [ $y_axis == 8 ] ;
                        then
                                cat $output_file | awk '{print $8}' > y_axis.txt
                                y_label="Kinetic Temperature [K]"

                        elif [ $y_axis == 9 ] ;
                        then
                                cat $output_file | awk '{print $9}' > y_axis.txt
                                y_label="N [# molecules adsorbed]"

                        elif [ $y_axis == 10 ] ;
                        then
                                cat $output_file | awk '{print $10}' > y_axis.txt
                                y_label="Spin Ratio [%]"

                        elif [ $y_axis == 11 ] ;
                        then
                                cat $output_file | awk '{print $11}' > y_axis.txt
                                y_label="Volume [Angstroms^3]"

                        elif [ $y_axis == 12 ] ;
                        then
                                cat $output_file | awk '{print $12}' > y_axis.txt
                                y_label="Core Temperature [K]"

	fi #end y-axis
	
	#clear tmp files
	if [ -f tmp1_ams2024.txt ];
	then
		rm tmp1_ams2024.txt
	fi
	if [ -f tmp2_ams2024.txt ];
	then
		rm tmp2_ams2024.txt
	fi

else
	echo "${red}Invalid response. Please type a number between 1-3.${reset}"
	exit
fi

#---------------------------------------------------------------------Plot things
# Combine the two files into a single file with two columns

echo -e "${bold}Plotting...${reset}"

if [ ! -f x_axis.txt ] || [ ! -f y_axis.txt ] ;
then
        echo -e "${red}Error: No x axis or no y axis data was generated.${reset}"
        exit
fi

paste x_axis.txt y_axis.txt > data.txt

if [ "$program" == "1" ] && [ "$job_type" == "ENERGY_FORCE" ] && [ "$x_axis" == "4" ] || [ "$y_axis" == "4" ] || [ "$x_axis" == "5" ] || [ "$y_axis" == "5" ] ;
then
	sed '2~3d;3~3d' data.txt > tmp11.txt
	mv tmp11.txt data.txt
fi


#Clean up the temporary files
rm x_axis.txt y_axis.txt

echo -e "${green}${bold}Done!${reset}"
echo -e "${cyan}${reset}"

echo -e "${cyan}${bold}Generating the following data.csv...${reset}"
echo -e "${cyan}$x_label    $y_label"
cat data.txt
echo -e "${reset}"

echo "$x_label,$y_label" > data.csv
sed 's/ \+/,/g' data.txt >> data.csv

(cat<<EOF
set title "${y_label} vs ${x_label}\n ${subtitle}"
set xlabel "${x_label}"
set ylabel "${y_label}"
plot "data.txt" using 1:2 with lines title "${y_label}"
pause -1 "Press Ctrl-D to exit."
EOF
cat /dev/tty) | gnuplot

