#!/bin/bash
#Script by Angela Shipman, v2 Oct 9 2023

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
echo -e "${yellow}${bold}Which program is the runlog from?${reset}${yellow}\n\n1) CP2K\n2) MPMC\n3) ORCA\n\nYour answer (#):"
read program

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
				grep -B 3 "Step number" runlog.log | grep "ENERGY| Total FORCE_EVAL" | awk '{printf $9 "\n"}' | sort -u > x_axis.txt
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
				grep "ENERGY| Total FORCE_EVAL" runlog.log | awk '{printf $9 "\n"}' > y_axis.txt
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

	echo "blah"

#---------------------------------------------------------------------ORCA
elif [ $program == 3 ] ; #ORCA
then
	echo "blah"


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
set title "${y_label} vs ${x_label}"
set xlabel "${x_label}"
set ylabel "${y_label}"
plot "data.txt" using 1:2 with lines title "${y_label}"
pause -1 "Press Ctrl-D to exit."
EOF
cat /dev/tty) | gnuplot

