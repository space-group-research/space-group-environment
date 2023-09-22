#UPDATE ENTRY (UE)
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

#Check for necessary files
echo -e "@--------------------------------------------------@"
echo "${bold}Updating entry...${reset}"
#Check for submit script and program...
if [ ! -f entry.txt ] ;
then
    echo -e "${bold}${red}You have not made an entry yet, so there is not one to update. Use command 'ce' to create an entry.${reset}"
else
program=`grep "Program |" entry.txt | awk '{print $4}'`


#-------------------------------------------------------------------------------------------------------------BUILD VARIABLES BASED ON PROGRAM
#@---------------------------------------------------------CP2K--------------------------------------------------------@

if [ "$program" == "CP2K" ] ;
then    
#BUILD VARIABLES - CP2K

    #Check which computer the job is for
    if grep -q "SBATCH" *.sh
    then
            echo "${bold}This is a job on Bridges. Searching for *.sh, runlog.log, and slurm.*.out files.${reset}"
            #Bridges
            job_id=`find "slurm"* | grep -o '[0-9]\+'`
            #submit_time=`sacct -X -j $job_id --format Start | grep "T" | awk '{printf $0 " "}'`
    fi
    if grep -q "BSUB" *.sh
    then
            echo "${bold}This is a job on Hazel. Searching for *.sh, runlog.log, and stdout.* files.${reset}"
            #Hazel
            job_id=`grep "PROGRAM PROCESS ID" runlog.log | head -1 | awk '{printf $9 " "}'`
            submit_time=`grep "was submitted from host" stdout.* | awk '{printf $15" "$16" "$17" "$18" "$19" "}'`
    fi

    #Generic CP2K variables:
    pickup_time=`grep "PROGRAM STARTED AT" runlog.log | awk '{printf $8" "$9", "}'`
    number_run=`grep -o ',' <<<"$pickup_time" | grep -c .`
    final_time=`grep "PROGRAM ENDED AT" runlog.log | awk '{printf $8" "$9", "}'`
    final_energy_au=`grep "ENERGY|" runlog.log | tail -1 | awk '{printf $9 " "}'`
    if [ ! -z "$final_energy_au" ];
    then
        final_energy_eV=`echo "$final_energy_au * 27.211324570273" | bc -l`
    fi

    #Determine CP2K Job Type & Grab new info
    if grep -q "Run Type:  ENERGY_FORCE" entry.txt ;
    then
        job_type="ENERGY_FORCE"
        echo "${bold}This job is an ENERGY_FORCE job.${reset}"
        total_steps=`grep " OT " runlog.log | awk '{print $1}' | tail -1`
        avg_time=`grep " OT " runlog.log | awk '{sum+=$5} END {print sum/NR,"seconds"}'`
    fi

    if grep -q "Run Type:  GEO_OPT" entry.txt ;
    then
        job_type="GEO_OPT"
        echo "${bold}This job is a GEO_OPT job.${reset}"
                total_steps=`grep "Informations at step" runlog.log | tail -1 | grep -o '[0-9]\+'`
                avg_time=`grep "Used time" runlog.log | awk '{sum+=$4} END {print sum/NR,"seconds"}'`
                if [[ $(grep "GEOMETRY OPTIMIZATION COMPLETED" runlog.log) ]] ;
                then
                        #finished="Job completed!"
                        final_cell_lengths=`grep -A 11 "GEOMETRY OPTIMIZATION COMPLETED" runlog.log | tail -7 | grep Vector| awk '{printf $10 " "}'`
                        final_cell_angles=`grep -A 11 "GEOMETRY OPTIMIZATION COMPLETED" runlog.log | tail -7 | grep Angle | awk '{printf $6 " "}'`
                        final_cell_volume=`grep -A 11 "GEOMETRY OPTIMIZATION COMPLETED" runlog.log | tail -7 | grep Volume | awk '{printf $4 " "}'`

                else
                        echo -e "${red}${bold}Job did not complete! Grabbing info from last recorded step... ${reset}"
                        #finished="Job did not reach completion."
                        final_cell_lengths=`grep -A 12 "OPTIMIZATION STEP" runlog.log | grep "CELL|" | tail -9 | grep Vector | awk '{printf $10 " "}'`
                        final_cell_angles=`grep -A 12 "OPTIMIZATION STEP" runlog.log | grep "CELL|" | tail -9 | grep Angle | awk '{printf $6 " "}'`
                        final_cell_volume=`grep -A 12 "OPTIMIZATION STEP" runlog.log | grep "CELL|" | tail -9 | grep Volume | awk '{printf $4 " "}'`
                fi
    fi

    if grep -q "Run Type:  CELL_OPT" entry.txt ;
    then
        job_type="CELL_OPT"
        echo "${bold}This job is a CELL_OPT job.${reset}"
                total_steps=`grep "Informations at step" runlog.log | tail -1 | grep -o '[0-9]\+'`
                avg_time=`grep "Used time" runlog.log | awk '{sum+=$4} END {print sum/NR,"seconds"}'`
                if [[ $(grep "GEOMETRY OPTIMIZATION COMPLETED" runlog.log) ]] ;
                then
                        #finished="Job completed!"
                        final_cell_lengths=`grep -A 11 "GEOMETRY OPTIMIZATION COMPLETED" runlog.log | tail -7 | grep Vector| awk '{printf $10 " "}'`
                        final_cell_angles=`grep -A 11 "GEOMETRY OPTIMIZATION COMPLETED" runlog.log | tail -7 | grep Angle | awk '{printf $6 " "}'`
                        final_cell_volume=`grep -A 11 "GEOMETRY OPTIMIZATION COMPLETED" runlog.log | tail -7 | grep Volume | awk '{printf $4 " "}'`

                else
                        echo -e "${red}${bold}Job did not complete! Grabbing info from last recorded step... ${reset}"
                        #finished="Job did not reach completion."
                        final_cell_lengths=`grep -A 12 "OPTIMIZATION STEP" runlog.log | grep "CELL|" | tail -9 | grep Vector | awk '{printf $10 " "}'`
                        final_cell_angles=`grep -A 12 "OPTIMIZATION STEP" runlog.log | grep "CELL|" | tail -9 | grep Angle | awk '{printf $6 " "}'`
                        final_cell_volume=`grep -A 12 "OPTIMIZATION STEP" runlog.log | grep "CELL|" | tail -9 | grep Volume | awk '{printf $4 " "}'`
                fi
    fi

    if grep -q "Run Type:  MD" entry.txt ;
    then
        job_type="MD"
            echo "${bold}This is an MD job.${reset}"
            avg_time=`grep "CPU time per MD step" runlog.log | awk '{print $9}' | tail -1`
        total_steps=`grep "Step number" runlog.log | tail -1 | grep -o '[0-9]\+'`
        step_size=`grep "TIMESTEP " *.inp | awk '{print $2}'`
        final_sim_time=$(echo "scale=4; $step_size*$total_steps" | bc)
        final_avg_temp=`grep "Temperature" runlog.log | awk '{print $5}' | tail -1`
        final_avg_pres=`grep "Pressure" runlog.log | awk '{print $5}' | tail -1`
        final_avg_ke=`grep "Kinetic energy" runlog.log | awk '{print $6}' | tail -1`
        final_avg_pe=`grep "Potential energy" runlog.log | awk '{print $6}' | tail -1`
        final_avg_drift=`grep "Energy drift per atom" runlog.log | awk '{print $8}' | tail -1`
    fi
#UPDATE ENTRY - CP2K
#Check for entry to update / if there is an entry, update it with variables. If not, throw error.
    if test -z entry.txt;
    then
            echo "${bold}${red}You have not made an entry yet, so there is not one to update. Use command 'ce' to create an entry.${reset}"
    else
            if test -z runlog.log ;
            then
                echo "${bold}${red}The job has not started yet. Wait for a runlog to be generated.${reset}"
            else
            #Check if there is already a results section - replace lines if so.
            if [[ $(grep "R E S U L T S" entry.txt) ]] ;
                    then
            
            #CP2K Update Results Section
                sed -i "s/This job has been submitted a total of.*/This job has been submitted a total of $number_run time(s)./g" entry.txt
                sed -i "s/Date and Time of Job Submission:.*/Date and Time of Job Submission: $submit_time/g" entry.txt
                sed -i "s/Date and Time of Job Pick Up:.*/Date and Time of Job Pick Up: $pickup_time/g" entry.txt
                sed -i "s/Date and Time of Job Completion.*/Date and Time of Job Completion: $final_time/g" entry.txt
                sed -i "s/Total Number of Steps.*/Total Number of Steps: $total_steps/g" entry.txt
                sed -i "s/Average Time Per Step.*/Average Time Per Step [s]: $avg_time/g" entry.txt
                sed -i "s/Final Energy.*/Final Energy: $final_energy_au a.u. or $final_energy_eV eV/g" entry.txt
                if [ "$job_type" == "GEO_OPT" ] ;
                then
                    sed -i "s/Final Cell Lengths.*/Final Cell Lengths (Angstrom): $final_cell_lengths/g" entry.txt
                    sed -i "s/Final Cell Angles.*/Final Cell Angles (alpha beta gamma): $final_cell_angles/g" entry.txt
                    sed -i "s/Final Cell Volume.*/Final Cell Volume (Angstrom^3): $final_cell_volume/g" entry.txt
                fi                
                if [ "$job_type" == "CELL_OPT" ] ;
                then
                    sed -i "s/Final Cell Lengths.*/Final Cell Lengths (Angstrom): $final_cell_lengths/g" entry.txt
                    sed -i "s/Final Cell Angles.*/Final Cell Angles (alpha beta gamma): $final_cell_angles/g" entry.txt
                    sed -i "s/Final Cell Volume.*/Final Cell Volume (Angstrom^3): $final_cell_volume/g" entry.txt
                fi
                if [ "$job_type" == "MD" ] ;
                then
                    sed -i "s/Length of simulation.*/Length of simulation [fs]: $final_sim_time/g" entry.txt
                    sed -i "s/Final Average Temperature.*/Final Average Temperature [K]: $final_avg_temp/g" entry.txt
                    sed -i "s/Final Average Pressure.*/Final Average Pressure [bar]: $final_avg_pres/g" entry.txt
                    sed -i "s/Final Average Kinetic Energy.*/Final Average Kinetic Energy [hartree]: $final_avg_ke/g" entry.txt
                    sed -i "s/Final Average Potential Energy.*/Final Average Potential Energy [hartree]: $final_avg_pe/g" entry.txt
                    sed -i "s/Final Average Energy Drift per Atom.*/Final Average Energy Drift per Atom [K]: $final_avg_drift/g" entry.txt
                fi
                    else
                #If there is no results section, write one.
                            echo -e "@------------------------------------@\n             R E S U L T S\n@------------------------------------@" >> entry.txt
                            echo -e "*If fields are blank, ensure the job has finished.\n" >> entry.txt
            #CP2K Write Results Section
                echo -e "This job has been submitted a total of $number_run time(s)." >> entry.txt
                echo -e "Date and Time of Job Submission: " $submit_time >> entry.txt
                echo -e "Date and Time of Job Pick Up: " $pickup_time >> entry.txt
                echo -e "Date and Time of Job Completion: " $final_time >> entry.txt

                echo -e "\nTotal Number of Steps: " $total_steps >> entry.txt
                echo -e "Average Time Per Step [s]: " $avg_time >> entry.txt

                echo -e "\nInformation for last completed step (Step $total_steps)" >> entry.txt
                echo -e "Final Energy:" $final_energy_au "a.u. / " $final_energy_eV "eV" >> entry.txt
                if [ "$job_type" == "GEO_OPT" ] ;
                then
                    echo -e "Final Cell Lengths (a b c): " $final_cell_lengths >> entry.txt
                    echo -e "Final Cell Angles (alpha beta gamma): " $final_cell_angles >> entry.txt
                    echo -e "Final Cell Volume (Angstrom^3): " $final_cell_volume "\n" >> entry.txt
                fi
                if [ "$job_type" == "CELL_OPT" ] ;
                then
                    echo -e "Final Cell Lengths (a b c): " $final_cell_lengths >> entry.txt
                    echo -e "Final Cell Angles (alpha beta gamma): " $final_cell_angles >> entry.txt
                    echo -e "Final Cell Volume (Angstrom^3): " $final_cell_volume "\n" >> entry.txt
                fi
                if [ "$job_type" == "MD" ] ;
                then
                    echo -e "Length of simulation [fs]:" $final_sim_time >> entry.txt
                    echo -e "Final Average Temperature [K]: " $final_avg_temp >> entry.txt
                    echo -e "Final Average Pressure [bar]:" $final_avg_pres >> entry.txt
                    echo -e "Final Average Kinetic Energy [hartree]:" $final_avg_ke >> entry.txt
                    echo -e "Final Average Potential Energy [hartree]:" $final_avg_pe >> entry.txt
                    echo -e "Final Average Energy Drift per Atom [K]:" $final_avg_drift >> entry.txt
                fi
            fi
        fi
    fi
fi

#@---------------------------------------------------------MPMC--------------------------------------------------------@

if [ "$program" == "MPMC" ] ;
then
#BUILD VARIABLES - MPMC
    if [ ! -f runlog.log ] ;
        then
            echo "${bold}${red}The job has not started yet. Wait for a runlog to be generated.${reset}"
        else
        version=`grep "MPMC" runlog.log | awk '{print $6}'`
        pickup_time=`grep "processes started on" runlog.log | awk '{print $8" "$9}'`
        final_time=`grep "OUTPUT" runlog.log | tail -18 | grep "Root collecting statistics at" | awk '{print $6" "$7" "$8" "$9" "$10}'`
        e_per_p=`grep "OUTPUT" runlog.log | tail -18 | grep "Energy" | awk '{print $4}'`
        total_steps=`grep "Completed step" runlog.log | tail -1 | awk '{print $4}' | grep -o '[0-9]\+' | head -1`
        avg_time=`grep "OUTPUT" runlog.log | tail -18 | grep "sec/step" | awk '{print $2}'`
        qst=`grep "OUTPUT" runlog.log | tail -18 | grep "qst" | awk '{print $4}'`
        n=`grep "OUTPUT" runlog.log | tail -18 | grep "N" | awk '{print $4" "$5" "$6" "$7}'`
        wt=`grep "OUTPUT" runlog.log | tail -18 | grep "wt %" | grep -v "(ME)" | awk '{print $5" "$6" "$7" "$8}'`
        wt_ME=`grep "OUTPUT" runlog.log | tail -18 | grep "(ME)" | awk '{print $6" "$7" "$8" "$9}'`
        if [ "$(grep "Completed step" runlog.log | awk '{print $5}' | grep -o '[0-9]\+'.'[0-9]\+' | tail -1)" == "100.000" ] ; 
        then 
            echo -e "${green}${bold}Job completed!${reset}"
            finish="COMPLETED"
        else
            echo -e "${red}${bold}Job did not complete!${reset}"
            finish="INCOMPLETE"
        fi
#WRITE ENTRY - MPMC
#Check for entry to update / if there is an entry, update it with variables. If not, throw error.
        #Check if there is already a results section - replace lines if so.
        if [[ $(grep "R E S U L T S" entry.txt) ]] ;
        then    
        #MPMC Update Results Section
            sed -i "s/Program | Version.*/Program | Version: MPMC | $version/g" entry.txt
            sed -i "s/Date and Time of Job Pick Up.*/Date and Time of Job Pick Up: $pickup_time/g" entry.txt
            sed -i "s/Date and Time of Job End.*/Date and Time of Job End: $final_time/g" entry.txt
            sed -i "s/Completion Status.*/Completion Status: $finish/g" entry.txt

            sed -i "s/Completed Steps.*/Completed Steps: $total_steps/g" entry.txt
            sed -i "s@Time Per Step.*@Time Per Step [sec/step]: $avg_time@g" entry.txt

            sed -i "s/Information for last completed step.*/Information for last completed step (Step $total_steps)/g" entry.txt
            sed -i "s@Energy per particle.*@Energy per particle [kJ/mol]: $e_per_p@g" entry.txt
            sed -i "s@qst .*@qst [kJ/mol]: $qst@g" entry.txt
            sed -i "s/Molecules per unit cell.*/Molecules per unit cell [N]: $n/g" entry.txt
            sed -i "s/Correct Weight Percent.*/Correct Weight Percent: $wt/g" entry.txt
            sed -i "s/ME Weight Percent.*/ME Weight Percent: $wt_ME/g" entry.txt
        else
        #If there is no results section, write one.
            echo -e "@------------------------------------@\n             R E S U L T S\n@------------------------------------@" >> entry.txt
            echo -e "*If fields are blank, ensure the job has finished.\n" >> entry.txt
        #MPMC Write Results Section
            sed -i "s/Program | Version.*/Program | Version: MPMC | $version/g" entry.txt
            echo -e "Date and Time of Job Pick Up: " $pickup_time >> entry.txt
            echo -e "Date and Time of Job End: " $final_time >> entry.txt
            echo -e "Completion Status:" $finish >> entry.txt

            echo -e "\nCompleted Steps: " $total_steps >> entry.txt
            echo -e "Time Per Step [sec/step]: " $avg_time >> entry.txt

            echo -e "\nInformation for last completed step (Step $total_steps)" >> entry.txt
            echo -e "Energy per particle [kJ/mol]: " $e_per_p >> entry.txt
            echo -e "qst [kJ/mol]: " $qst >> entry.txt
            echo -e "Molecules per unit cell [N]:" $n >> entry.txt
            echo -e "Correct Weight Percent: " $wt >> entry.txt
            echo -e "ME Weight Percent: " $wt_ME >> entry.txt
        fi
    fi
fi


#@-------------------------------------------------------OPEN--MM------------------------------------------------------@
#Customize and Build your OpenMM variables here! - CUSTOMIZE ME!

if [ "$program" == "OPENMM" ] ;
then
    #if [ ! -f output.file ] ;
    #then
        #echo "${bold}${red}The job has not started yet. Wait for a runlog to be generated.${reset}"
    #else
    
#BUILD VARIABLES - OPEN MM (CUSTOMIZE ME!)
    #Warn user ue.sh hasn't been customized yet.
    echo -e "${bold}${red}You haven't customized your variables in the ${green}'ue.sh'${red} script in $entry_path yet! Look for designated OpenMM sections with ${magenta}'CUSTOMIZE ME!' ${red}strings. ${reset}"
    #Uncomment line below and start grabbing your variables!
    #your_custom_variables=``

#WRITE ENTRY - OPEN MM (CUSTOMIZE ME!)

         #Check if there is already a results section - replace lines if so.
        if [[ $(grep "R E S U L T S" entry.txt) ]] ;
        then    
        #OpenMM Update Results Section - CUSTOMIZE ME!
            echo -e "${red}${bold}Update the ${magenta}'OpenMM Update Results Section'${red} in the ${green}'ue.sh'${red} script (located in ${magenta}$entry_path ${red})!${reset}"
            #Uncomment lines below! This section updates an already-existing results section. This is so the entry can be updated on the fly.
            #sed -i "s/Your Original Lines from Writing Results section below/Your updated lines (see other program sections for examples of how sed works)/g" entry.txt
        else
            #If there is no results section, write one.
            echo -e "@------------------------------------@\n             R E S U L T S\n@------------------------------------@" >> entry.txt
            echo -e "*If fields are blank, ensure the job has finished.\n" >> entry.txt
        
            #OpenMM Write Results Section - CUSTOMIZE ME!
            echo -e "${red}${bold}Update the ${magenta}'OpenMM Write Results Section'${red} in the ${green}'ue.sh'${red} script (located in ${magenta}$entry_path ${red})!${reset}"
            #Uncomment the lines below! This section writes a RESULTS section if one hasn't already been written.
            #echo -e "Your Custom Information Here with: $your_custom_variables" >> entry.txt
        fi
        #fi - include if want a check for output files.
fi
#@---------------------------------------------------------ORCA--------------------------------------------------------@
if [ "$program" == "ORCA" ] ;
then
    if [ ! -f runlog.log ] ;
        then
        echo "${bold}${red}The job has not started yet. Wait for a runlog to be generated.${reset}"
        else
        #read job type
        job_type=`grep "Run Type:" entry.txt | awk '{print $3}'`
        version=`grep "Program Version" runlog.log | awk '{print $3}'`
    #BUILD VARIABLES - ORCA
        #Runlog variables
        final_time=`grep "TOTAL RUN TIME:" runlog.log | awk '{print $4" "$5" "$6" "$7" "$8" "$9" "$10" "$11" "$12" "$13}'`
        if [ "$job_type" == "Geo" ] ;
        then
            total_steps=`grep "ORCA GEOMETRY RELAXATION STEP" runlog.log | wc -l`
            avg_time=`grep "The final length of the internal step" runlog.log | awk '{sum+=$9} END {print sum/NR,"seconds"}'`
        fi
        if grep "****ORCA TERMINATED NORMALLY****" runlog.log &> /dev/null ;
        then
            finish="COMPLETE"
        else
            finish="INCOMPLETE"
        fi
        final_sng_pt_e=`grep "FINAL SINGLE POINT ENERGY" runlog.log | tail -1 | awk '{print $5}'`
        total_dipole_moment=`grep "Total Dipole Moment" runlog.log | awk '{print $1" "$2" "$3" "$4" "$5" (X) "$6" (Y) "$7" (Z) "}'`
        dipole_magnitude=`grep "Magnitude (a.u.)" runlog.log | awk '{print $4}'`
        if [ ! -f orca.xyz ] ;
        then
            echo -e "${bold}${red}No orca.xyz file detected.${reset}"
        else
            final_coord=`grep -A 1000 "Coordinates from ORCA-job" orca.xyz | awk 'NR > 1 { print }'`
        fi

#WRITE ENTRY - ORCA
#Check for entry to update / if there is an entry, update it with variables. If not, throw error.
        echo "${bold}${red}The job has not started yet. Wait for a runlog to be generated.${reset}"
        #Check if there is already a results section - replace lines if so.
        if [[ $(grep "R E S U L T S" entry.txt) ]] ;
        then
        
        #ORCA Update Results Section
            sed -i "s/Program | Version.*/Program | Version: ORCA | $version/g" entry.txt    
            sed -i "s/Total Time of Job.*/Total Time of Job: $final_time/g" entry.txt
            sed -i "s/Total Number of Steps.*/Total Number of Steps: $total_steps/g" entry.txt
            sed -i "s/Average Time Per Step.*/Average Time Per Step [s]: $avg_time/g" entry.txt
            sed -i "s/Final Single Point Energy.*/Final Single Point Energy: $final_sng_pt_e Eh/g" entry.txt
            sed -i "s/Total Dipole Moment.*/Total Dipole Moment (a.u.): $total_dipole_moment >> entry.txt/g" entry.txt
            sed -i "s/Dipole Magnitude.*/Dipole Magnitude (a.u.): $dipole_magnitude/g" entry.txt
        else    
        #If there is no results section, write one.
            echo -e "@------------------------------------@\n             R E S U L T S\n@------------------------------------@" >> entry.txt
            echo -e "*If fields are blank, ensure the job has finished.\n" >> entry.txt
        #ORCA Write Results Section
            sed -i "s/Program | Version.*/Program | Version: ORCA | $version/g" entry.txt
            echo -e "Total Time of Job: " $final_time >> entry.txt

            echo -e "\nTotal Number of Steps: " $total_steps >> entry.txt
            echo -e "Average Time Per Step [s]: " $avg_time >> entry.txt

            echo -e "\nInformation for last completed step (Step $total_steps)" >> entry.txt
            echo -e "Final Single Point Energy:" $final_sng_pt_e "Eh" >> entry.txt
            echo -e "Total Dipole Moment (a.u.):" $total_dipole_moment >> entry.txt
            echo -e "Dipole Magnitude (a.u.):" $dipole_magnitude >> entry.txt
        fi
    fi
fi


echo "${bold}Entry has been updated.${reset}"
echo "${magenta}If some fields are blank, ensure the job has finished.${reset}"
fi
echo -e "@--------------------------------------------------@"
