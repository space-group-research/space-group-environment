#CREATE ENTRY (CE)

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

source $entry_path/var_storage.sh
echo -e "${bold}Creating entry...${reset}"
echo -e "${bold}Searching for ${green}submit script${reset} ${bold}and ${green}input script${reset}${bold} files...${reset}"

#Check if Entry already exists
if test -f entry*.txt;
then
        echo "${red}${bold}WARNING: A starting entry already exists for this job. Are you sure you want to overwrite it?${reset}"
        select yn in "Yes" "No";
        do
                case $yn in
                        Yes )
                                echo -e "\n${bold}Overwriting previous entry...${reset}\n"
                break;;
                        No )
                                echo -e "\n${bold}Exiting entry creation.\n${reset}"
                exit;;
                esac
        done
fi

#Get User Input
        if test -z "$default_material" ;
        then
                echo -e ${green}"\nWhat material are you simulating?\n${red}${bold}DO NOT USE ' IN YOUR ANSWER.\n${reset}${green}Your answer:${reset}"
                read material
                echo -e ${green}"\nWould you like to set this as your default material? (You can change this at any time.)${bold} [Select '1' or '2']${reset}\n"
                select yn in "Yes" "No";
                do
                        case $yn in
                                Yes )
                                        echo -e "${bold} New default material set.${reset}"
                    export default_material="$material"
                                        break;;
                                No )
                                       break;;
                        esac
                done
        else
                echo -e "\n${green}What material are you simulating? ${bold}${red}\nDO NOT USE ' IN YOUR ANSWER.\n${reset}${green}${bold}If you want to use your default material, type '1'.${reset}\n${green}Your default material is currently: ${magenta}$default_material \n${green}Your answer:"${reset}
                read material
                if [[ "$material" != "1" ]] ;
                then
                        echo -e "\n${bold}$material has been set as the material for this job.${reset}\n"
                        echo -e ${green}"Would you like to set this material as your new default material? ${bold} [Select '1' or '2'] ${reset}"
                        select yn in "Yes" "No";
                        do
                                case $yn in
                                        Yes )
                                                echo -e "${bold} New default material set.${reset}"
                        export default_material="$material"
                                                break;;
                                        No )
                                                break;;
                                esac
                        done
                fi
        fi
        echo -e "${green}\nIs it loaded with gas? ${bold}[Answer '1' or '2']${reset}\n"
        select yn in "Yes" "No";
        do
                case $yn in
                        Yes )
                                echo -e "\n${green}What is the gas? (e.g. water would be 'H2O')\n${red}Do not use ' in your answer.${green}\nYour answer: ${reset}"
                                read gas
                                echo -e "\n${green}What is the concentration of gas? (e.g. 5 M, 4 per uc, etc.)\n${red}Do not use ' in your answer. ${green}\nYour answer: ${reset}"
                                read concentration
                                break;;
                        No )
                                unset gas
                                unset concentration
                                break;;
                esac
        done
        echo -e "\n${green}What are you hoping to learn from this job (if anything unique from routine)?\n${red}${bold}DO NOT USE ' IN YOUR ANSWER.\n${reset}${bold}${green}Type '1' if you'd like to use your previous job's description."${reset}
        if test -z "$previous_description" ;
        then
                echo -e "\n${red}No descriptions stored on this login yet.${green}\nYour answer: ${reset}"
        else
                echo -e "\n${green}${bold}Last Recorded Description:${reset}${magenta} $previous_description${green}\nYour answer: ${reset}"
        fi
        read description
        if [[ "$description" != "1" ]] ;
        then
                echo "${bold}Updated description.${reset}"
                export previous_description="$description"
        fi
echo -e "@--------------------------------------------------@"



#Build universal variables (regardless of program)
        user=`last | grep -o $USER | head -1`
        dir=`pwd`
        dt=$(date '+%d/%m/%Y %H:%M:%S')
    increment_entry_num=$(($entry_number+1))

#Check What Program Is Used
if [ -f *.inp ] ;
then
check=`head -n 1 *.inp | awk '{print $1}'`

        if grep -q "&GLOBAL" *.inp || grep -q "&GLOBAL" *-1.restart &> /dev/null ; #assuming all cp2k *.inp files start with &GLOBAL
        then
                echo -e "${bold}You are using CP2K${reset}"
                program="CP2K"

        elif grep  "job_name" *.inp &> /dev/null ; #assuming all mpmc *.inp files contains 'job_name'
        then
                echo -e "${bold}You are using MPMC.${reset}"
                program="MPMC"

        elif grep "import openmm" *.py &> /dev/null || grep "import openff" *.py &> /dev/null ; #assuming all open mm jobs have a *.py file with 'import openmm'
        then
                echo -e "${bold}You are using OPENMM.${reset}"
                program="OPENMM"

        elif [ "$check" = "!" ] ; #assuming all ORCA *.inp files start with '!'
        then
                echo -e "${bold}You are using ORCA.${reset}"
                program="ORCA"
        else
                echo -e "${bold}${magenta}*.inp file found, but cannot detect a program in list of default programs. Writing basic entry.${reset}"
                program="No program detected."
        fi
else
        echo -e "${bold}${cyan}No *.inp file found. Writing basic entry.${reset}"
        program="No program detected."
fi

#--------------------------------------------------------------------------------------------------------------------BUILD PROGRAM-UNIQUE ENTRIES
if [ "$program" == "No program detected." ] ;
then

#WRITE ENTRY - No Program Detected
    #Introduction
    #Titles are 36 -, four spaces between words.
    echo -e "\n\n\n\n" > entry.txt
    echo -e "________________________________________________________________\nENTRY # $entry_number$index\n" >> entry.txt
    echo -e "@------------------------------------@\n      E N T R Y    D E T A I L S\n@------------------------------------@" >> entry.txt
    echo -e "By: " $user >> entry.txt
    echo -e "Date of Entry (DD/MM/YYYY HR:MIN:SEC): " $dt >> entry.txt
    echo -e "Working Directory: "$dir >> entry.txt
    echo -e "Program | Version: " $program >> entry.txt
    if [[ $description == "1" ]]
    then
        echo -e "Description: " $previous_description >> entry.txt
    else
        echo -e "Description: " $description >> entry.txt
    fi
    echo -e "\n" >> entry.txt

    #System Details
    echo -e "@------------------------------------@\n   M A T E R I A L    D E T A I L S\n@------------------------------------@" >> entry.txt

    if [[ $material == "1" ]]
    then
        echo "Material: " $default_material >> entry.txt
    else
        echo -e "Material: " $material >> entry.txt
    fi

    if test -z $gas;
    then
        echo "Gas Loading: None." >> entry.txt
    else
        echo "Gas Loading: " $gas ',' $concentration >> entry.txt
    fi
    echo -e "\n" >> entry.txt

fi

#@---------------------------------------------------------CP2K--------------------------------------------------------@

if [ "$program" == "CP2K" ] ;
then

#BUILD VARIABLES - CP2K
    version=`grep "source" *.sh | grep -o '[0-9]\+.[0-9]\+'`    
    
    #Check for restart file
    if [ -f *-1.restart ] ;
        then
                echo -e "${bold}${magenta}*-1.restart file detected.${reset}"
                if [ ! -f *.inp ] ;
                then
                        echo "${red}${bold}No input file (*.inp) detected. Grabbing information from *-1.restart file instead. Please note that some information may be missing and you may need to refer to the original input script.${reset}"
                fi
                lattice_a=`grep -A 3 "&CELL" *-1.restart | grep "A"`
                lattice_b=`grep -A 3 "&CELL" *-1.restart | grep "B"`
                lattice_c=`grep -A 3 "&CELL" *-1.restart | grep "C "`
                lattice_len="[Restart file cell matrix]\n$lattice_a\n$lattice_b\n$lattice_c"
                lattice_ang="\nRestart files do not contain cell angles, and Bash isn't designed for trigonometry. Please refer to original input script, calculated angles in runlog.log, or calculate them by hand from the cell matrix above."
                multiplicity=`grep "MULTIPLICITY " *-1.restart | awk '{printf $2 "\n"}'`
                job_type=`grep  "RUN_TYPE " *-1.restart | awk '{printf $2 "\n"}'`
                functional1=`grep -A 2 "&XC_FUNCTIONAL " *-1.restart | grep "&END " | awk '{printf $2 "\n"}'`
                functional2=`grep -A 8 "&VDW_" *-1.restart | grep "^ *TYPE " | awk '{printf $2 "\n"}'`
                atoms=`grep "&KIND" *-1.restart | awk '{printf $2 "\n"}'`
                basis_sets=`grep -A 2 "&KIND" *-1.restart | grep "BASIS_SET " | awk '{printf $2 "\n"}'`
        else
                echo "No *-1.restart file detected."
        fi

    #Check for input file, overwrite data gathered from restart file. Otherwise, keep data from restart file.
        if [ -f *.inp ] ;
        then
                echo -e "${bold}${green}*.inp file detected.${reset}"
                lattice_len=`grep "ABC " *.inp | awk '{printf $2" "$3" "$4 "\n"}'`
                lattice_ang=`grep "ALPHA_BETA_GAMMA " *.inp | awk '{printf $2" "$3" "$4 "\n"}'`
                multiplicity=`grep "MULTIPLICITY " *.inp | awk '{printf $2 "\n"}'`
                job_type=`grep "RUN_TYPE " *.inp | awk '{printf $2 "\n"}'`
                functional1=`grep "&XC_FUNCTIONAL " *.inp | awk '{printf $2 "\n"}'`
                functional2=`grep -A 8 "&VDW_" *.inp | grep "^ *TYPE " | awk '{printf $2 "\n"}'`
                atoms=`grep "&KIND" *.inp | awk '{printf $2 "\n"}'`
                basis_sets=`grep -A 2 "&KIND" *.inp | grep "BASIS_SET " | awk '{printf $2 "\n"}'`
        else
                echo "No *.inp file detected."
        fi

        #Check which computer the job is for
        if grep -q "SBATCH" *.sh ;
        then
                echo "${bold}This is a job on Bridges. Modifying search for variables...${reset}"
                cpu_request=`grep "mpirun -np " *.sh | awk '{printf $3 "\n"}'`
                requested_time=`grep "SBATCH -t" *.sh | awk '{printf $3 "\n"}'`
        fi
        if grep -q "BSUB" *.sh
        then
                echo "${bold}This is a job on Hazel. Modifying search for variables...${reset}"
                cpu_request=`grep "#BSUB -n" *.sh | awk '{printf $3 "\n"}'`
                requested_time=`grep "#BSUB -W " *.sh | awk '{printf $3 "\n"}'`
        fi

        #Check the job type
        if [ "$job_type" == "ENERGY_FORCE" ] ;
        then
                echo -e "${bold}This is an ENERGY_FORCE job.${reset}"
        fi

        if [ "$job_type" == "GEO_OPT" ] ;
        then
                echo -e "${bold}This is a GEO_OPT job.${reset}"
        fi

        if [ "$job_type" == "CELL_OPT" ] ;
        then
                echo -e "${bold}This is a CELL_OPT job.${reset}"
        fi

        if [ "$job_type" == "MD" ] ;
        then
                echo -e "${bold}This is an MD job.${reset}"
                ensemble=`grep "ENSEMBLE" *.inp | awk '{print $2}'`
                timestep=`grep "TIMESTEP " *.inp | awk '{print $2}'`
                thermostat=`grep -A 1 "&THERMOSTAT" *.inp | grep "TYPE" | awk '{print $2}'`
                temperature=`grep "TEMPERATURE" *.inp | awk '{print $2}'`
        fi

#WRITE ENTRY
    #Introduction
    #Titles are 36 -, four spaces between words.
    echo -e "\n\n\n\n" > entry.txt
    echo -e "________________________________________________________________\nENTRY # $entry_number$index\n" >> entry.txt
    echo -e "@------------------------------------@\n      E N T R Y    D E T A I L S\n@------------------------------------@" >> entry.txt
    echo -e "By: " $user >> entry.txt
    echo -e "Date of Entry (DD/MM/YYYY HR:MIN:SEC): " $dt >> entry.txt
    echo -e "Working Directory: "$dir >> entry.txt
    echo -e "Program | Version: " $program " | " $version >> entry.txt
    if [[ $description == "1" ]]
    then
            echo -e "Description: " $previous_description >> entry.txt
    else    
            echo -e "Description: " $description >> entry.txt
    fi
    echo -e "\n" >> entry.txt

    #System Details
    echo -e "@------------------------------------@\n   M A T E R I A L    D E T A I L S\n@------------------------------------@" >> entry.txt

    if [[ $material == "1" ]]
    then
            echo "Material: " $default_material >> entry.txt
    else
            echo -e "Material: " $material >> entry.txt
    fi
    echo -e "Initial Cell Lengths (a b c): " $lattice_len >> entry.txt
    echo -e "Initial Cell Angles (alpha beta gamma): " $lattice_ang >> entry.txt
    echo -e "Multiplcity: " $multiplicity >> entry.txt

    if test -z $gas;
    then
        echo "Gas Loading: None." >> entry.txt
    else
        echo "Gas Loading: " $gas ',' $concentration >> entry.txt
    fi
    echo -e "\n" >> entry.txt

    #Job Details
    echo -e "@------------------------------------@\n      I N P U T    D E T A I L S\n@------------------------------------@" >> entry.txt
    echo -e "Run Type: " $job_type >> entry.txt
    echo -e "Functional: " $functional1 $functional2 >> entry.txt
    if [ "$job_type" == "MD" ]
    then
        echo -e "Ensemble: " $ensemble >> entry.txt
        echo -e "Timestep [fs]: " $timestep >> entry.txt
        echo -e "Thermostat Type: " $thermostat >> entry.txt
        echo -e "Temperature [K]: " $temperature >> entry.txt
    fi
    echo -e "Basis Sets: " >> entry.txt
    printf '%s\n' "$atoms" "$basis_sets" | pr -2 -T >> entry.txt
    echo -e "\n" >> entry.txt

    #Submit Details
    echo -e "@------------------------------------@\n     S U B M I T    D E T A I L S\n@------------------------------------@" >> entry.txt
    echo -e "Requested CPUs: " $cpu_request >> entry.txt
    echo -e "Set Time Limit: $requested_time hours\n" >> entry.txt

echo -e "${bold}Entry created! When${green} runlog.log ${reset}${bold}and ${green}stdout${reset}${bold} files are generated, use ${magenta}ue${reset} ${bold}command to update the entry with more information.${reset}"

fi








#@---------------------------------------------------------MPMC--------------------------------------------------------@
if [ "$program" == "MPMC" ] ;
then

#BUILD VARIABLES - MPMC
    version="(Use 'ue' command to grab program version.)"
           
    if [ ! -f *.inp ] ; #Check if *.inp does not exist
    then
        echo "${red}${bold}No input file (*.inp) detected. Please note that some information may be missing and you may need to refer to the original input script.${reset}"
    
    else #if *.inp does exist
            echo -e "${bold}${green}*.inp file detected.${reset}"
        
        if grep -q "abcbasis" *.inp &> /dev/null ;
        then
            matrix="false"
            lattice_len=`grep "abcbasis" *.inp | awk '{print $2" "$3" "$4}'`
            lattice_ang=`grep "abcbasis" *.inp | awk '{print $5" "$6" "$7}'`
        elif grep "CRYST1" *.pqr | grep -v *.traj.pqr | grep -v *.restart.pqr ;
        then
            matrix="false"
            lattice_len=`grep "CRYST1" *.pqr | grep -v *.traj.pqr | grep -v *.restart.pqr | awk '{print $2" "$3" "$4}'`
            lattice_ang=`grep "CRYST1" *.pqr | grep -v *.traj.pqr | grep -v *.restart.pqr | awk '{print $5" "$6" "$7}'`
        elif grep "REMARK BOX" *pqr | grep -v *.traj.pqr | grep -v *.restart.pqr ;
        then
            matrix="true"
            lattice_len=`grep "REMARK BOX" *pqr | grep -v *.traj.pqr | grep -v *.restart.pqr | awk '{print $5" "$6" "$7}'`
        else
            echo -e "${bold}No lattice lengths or angles found in input scripts.${reset}"
        fi

            ensemble=`grep "ensemble" *.inp | awk '{print $2}'`
        temperature=`grep "temperature" *.inp | awk '{print $2}'`
        pressure=`grep "pressure" *.inp | awk '{print $2}'`
        numsteps=`grep "numsteps" *.inp | awk '{print $2}'`
        corrtime=`grep "corrtime" *.inp | awk '{print $2}'`
    fi

    #Check which computer the job is for
    if grep -q "SBATCH" *.sh &> /dev/null ;
    then
        echo "${bold}This is a job on Bridges. Modifying search for variables...${reset}"
        cpu_request=`grep "mpirun -np " *.sh | awk '{printf $3 "\n"}'`
        requested_time=`grep "SBATCH -t" *.sh | awk '{printf $3 "\n"}'`
    elif grep -q "BSUB" *.sh &> /dev/null ;
    then
        echo "${bold}This is a job on Hazel. Modifying search for variables...${reset}"
        cpu_request=`grep "#BSUB -n" *.sh | awk '{printf $3 "\n"}'`
        requested_time=`grep "#BSUB -W " *.sh | awk '{printf $3 "\n"}'`
    else
        echo -e "${bold}${red}No submit script (*.sh) found.${reset}"
    fi


#WRITE ENTRY - MPMC
    #Introduction
    #Titles are 36 -, four spaces between words.
    echo -e "\n\n\n\n" > entry.txt
    echo -e "________________________________________________________________\nENTRY # $entry_number$index\n" >> entry.txt
    echo -e "@------------------------------------@\n      E N T R Y    D E T A I L S\n@------------------------------------@" >> entry.txt
    echo -e "By: " $user >> entry.txt
    echo -e "Date of Entry (DD/MM/YYYY HR:MIN:SEC): " $dt >> entry.txt
    echo -e "Working Directory: "$dir >> entry.txt
    echo -e "Program | Version: " $program " | " $version >> entry.txt
    if [[ $description == "1" ]]
    then
        echo -e "Description: " $previous_description >> entry.txt
    else
        echo -e "Description: " $description >> entry.txt
    fi
    echo -e "\n" >> entry.txt

    #System Details
    echo -e "@------------------------------------@\n   M A T E R I A L    D E T A I L S\n@------------------------------------@" >> entry.txt

    if [[ $material == "1" ]]
    then
        echo "Material: " $default_material >> entry.txt
    else
        echo -e "Material: " $material >> entry.txt
    fi
    echo -e "Initial Cell Lengths (a b c): " $lattice_len >> entry.txt
    echo -e "Initial Cell Angles (alpha beta gamma): " $lattice_ang >> entry.txt

    if test -z $gas;
    then
        echo "Gas Loading: None." >> entry.txt
    else
        echo "Gas Loading: " $gas ',' $concentration >> entry.txt
    fi
    echo -e "\n" >> entry.txt

    #Input file Details
    echo -e "@------------------------------------@\n      I N P U T    D E T A I L S\n@------------------------------------@" >> entry.txt

    echo -e "Ensemble: " $ensemble >> entry.txt
    echo -e "Designated Number of Steps: " $numsteps >> entry.txt
    echo -e "Correlation Time [# steps to average]:" $corrtime  >> entry.txt
    echo -e "Temperature [K]: " $temperature >> entry.txt
    echo -e "Pressure [atm]: " $pressure >> entry.txt
    echo -e "\n" >> entry.txt

    #Submit Details
    echo -e "@------------------------------------@\n     S U B M I T    D E T A I L S\n@------------------------------------@" >> entry.txt
    echo -e "Requested CPUs: " $cpu_request >> entry.txt
    echo -e "Set Time Limit: $requested_time hours\n" >> entry.txt

#Reminder user to update file.
        echo -e "${bold}Entry created! When${green} runlog.log, .pqr ${reset}${bold}and ${green}.dat${reset}${bold} files are generated, use ${magenta}ue${reset} ${bold}command to update the entry with more information.${reset}"

fi


#@-------------------------------------------------------OPEN--MM------------------------------------------------------@
if [ "$program" == "OPENMM" ] ;
then

#BUILD VARIABLES - OPENMM
#CUSTOMIZE ME!
    #Warn user that this script hasn't been customized yet.
    echo -e "${magenta}${bold}Since OpenMM lacks a standard format, entries will largely be left blank except for your notes. If you have a personal format, you may write your own variables and customize entries into designated OpenMM sections in the ${green}ce.sh${magenta} and ${green}ue.sh${magenta} scripts in your ${reset}${bold}$entry_path ${magenta}directory. Search for ${red}CUSTOMIZE ME!${magenta} strings.${reset}"

    write_your_custom_variables_here=``
    #Common Variables Below:
    version="Unknown"
    lattice_len=``
        lattice_ang=``
        cpu_request=``
        requested_time=``


#WRITE ENTRY - OPENMM
    #Introduction - CUSTOMIZE ME!
    #Titles are 36 -, four spaces between words.
    echo -e "\n\n\n\n" > entry.txt
    echo -e "________________________________________________________________\nENTRY # $entry_number$index\n" >> entry.txt
    echo -e "@------------------------------------@\n      E N T R Y    D E T A I L S\n@------------------------------------@" >> entry.txt
    echo -e "By: " $user >> entry.txt
    echo -e "Date of Entry (DD/MM/YYYY HR:MIN:SEC): " $dt >> entry.txt
    echo -e "Working Directory: "$dir >> entry.txt
    echo -e "Program | Version: " $program " | " $version >> entry.txt
    if [[ $description == "1" ]]
    then
        echo -e "Description: " $previous_description >> entry.txt
    else
        echo -e "Description: " $description >> entry.txt
    fi
    echo -e "\n" >> entry.txt
    echo -e "Since OpenMM lacks a standard format, entries will largely be left blank except for your notes. If you have a personal format, you may write your own variables and customize entries into designated OpenMM sections in the 'ce.sh' and 'ue.sh' scripts in your $entry_path directory. Search for 'CUSTOMIZE ME!' strings." >> entry.txt
    echo -e "\n" >> entry.txt

    echo -e "@------------------------------------@\n   M A T E R I A L    D E T A I L S\n@------------------------------------@" >> entry.txt
    #THIS WRITES THE MATERIAL DETAILS SECTION - CUSTOMIZE ME!
    #IF YOU WANT TO CUSTOMIZE FOR OPENMM, uncomment the lines below:
    #echo -e "Your text" $your_custom_variables >> entry.txt

    if [[ $material == "1" ]]
    then
        echo "Material: " $default_material >> entry.txt
    else
        echo -e "Material: " $material >> entry.txt
    fi
    #echo -e "Initial Cell Lengths (a b c): " $lattice_len >> entry.txt
    #echo -e "Initial Cell Angles (alpha beta gamma): " $lattice_ang >> entry.txt

    if test -z $gas;
    then
        echo "Gas Loading: None." >> entry.txt
    else
        echo "Gas Loading: " $gas ',' $concentration >> entry.txt
    fi
    echo -e "\n" >> entry.txt

    echo -e "@------------------------------------@\n      I N P U T    D E T A I L S\n@------------------------------------@" >> entry.txt

    #This writes the INPUT DETAILS SECTION - CUSTOMIZE ME!
    #IF YOU WANT TO CUSTOMIZE FOR OPENMM, uncomment the line below and edit:
    #echo -e "Your custom lines here with $your_custom_variables" >> entry.txt
    echo -e "\n" >> entry.txt

    echo -e "@------------------------------------@\n     S U B M I T    D E T A I L S\n@------------------------------------@" >> entry.txt
    #This writes the  SUBMIT SCRIPT DETAILS - CUSTOMIZE ME!
    #echo -e "Your custom submit script details here with $your_custom_variables." >> entry.txt
    #echo -e "Requested CPUs: " $cpu_request >> entry.txt
    #echo -e "Set Time Limit: $requested_time hours\n" >> entry.txt
    echo -e "\n" >> entry.txt

    #Reminder user to update file.
        echo -e "${bold}Entry created! When${green} OUTPUT ${reset}${bold} files are generated, use ${magenta}ue${reset} ${bold}command to update the entry with more information.${reset}"

fi
#@---------------------------------------------------------ORCA--------------------------------------------------------@
if [ "$program" == "ORCA" ] ;
then
#BUILD VARIABLES - ORCA
    #Inp variables
    basis_set=`head -1 *.inp | awk '{print $3}'`
    functional1=`head -1 *.inp | awk '{print $2}'`
    cpu_request=`grep -i "nprocs" *.inp | awk '{print $2}'`
    requested_mem=`grep -i "%MaxCore" *.inp | awk '{print $2}'`
    #Grab job type
    if head -1 *.inp | grep -i "opt" &> /dev/null ; 
    then
        job_type="Geo Opt"
    elif head -1 *.inp | grep -i "MD" &> /dev/null ;
        then
                job_type="MD"
        timestep=`grep "timestep" *.inp | awk -F'[ _]' '{print $2}'`
        temperature=`grep "initvel" *.inp | awk -F'[ _]' '{print $2}'`
        thermostat=`grep "thermostat" *.inp | awk -F'[ _]' '{print $2}'`
        ensemble="NVT"
    elif head -1 *.inp | grep -i "EnGrad" &> /dev/null ;
    then
        job_type="Energy + Gradient"
    elif head -1 *.inp | grep -i "Freq" &> /dev/null ;
    then
        job_type="Vibrational Frequencies"
    fi

    
    #runlog.log variables
    version="(Once runlog.log is generated, update this entry (ue) for the program version.)"

    #xyz variables
    if [ -f orca.xyz ] || [ -f orca_trj.xyz ] ;
    then
        xyz_file=`ls *.xyz | grep -v orca_trj.xyz | grep -v orca.xyz`
        num_atoms=`cat $xyz_file | head -1`
                atoms=`cat $xyz_file | awk 'NR > 2 { print }' | awk '{print $1}' | sort -u`
        
    else
        if [ ! -f *.xyz ] ;
        then
            echo -e "${red}${bold}No *.xyz file detected.${reset}"
        else
            num_atoms=`cat *.xyz | head -1`
            atoms=`cat *.xyz | awk 'NR > 2 { print }' | awk '{print $1}' | sort -u`
        fi
    fi

    #Submit variables
    cpu_request=`grep -i "nprocs" *.inp | awk '{print $2}'`
    requested_mem=`grep -i "%MaxCore" *.inp | awk '{print $2}'`

#WRITE ENTRY - ORCA
    #Titles are 36 -, four spaces between words.
    echo -e "\n\n\n\n" > entry.txt
    echo -e "________________________________________________________________\nENTRY # $entry_number$index\n" >> entry.txt
    echo -e "@------------------------------------@\n      E N T R Y    D E T A I L S\n@------------------------------------@" >> entry.txt
    echo -e "By: " $user >> entry.txt
    echo -e "Date of Entry (DD/MM/YYYY HR:MIN:SEC): " $dt >> entry.txt
    echo -e "Working Directory: "$dir >> entry.txt
    echo -e "Program | Version: " $program " | " $version >> entry.txt
    if [[ $description == "1" ]]
    then
        echo -e "Description: " $previous_description >> entry.txt
    else
        echo -e "Description: " $description >> entry.txt
    fi
    echo -e "\n" >> entry.txt

    #System Details
    echo -e "@------------------------------------@\n   M A T E R I A L    D E T A I L S\n@------------------------------------@" >> entry.txt

    if [[ $material == "1" ]]
    then
        echo "Material: " $default_material >> entry.txt
    else
        echo -e "Material: " $material >> entry.txt
    fi


    if test -z $gas;
    then
        echo "Gas Loading: None." >> entry.txt
    else
        echo "Gas Loading: " $gas ',' $concentration >> entry.txt
    fi
    echo -e "\n" >> entry.txt

    #Job Details
    echo -e "@------------------------------------@\n      I N P U T    D E T A I L S\n@------------------------------------@" >> entry.txt
    echo -e "Run Type: " $job_type >> entry.txt
        echo -e "Functional: " $functional1 $functional2 >> entry.txt
    echo -e "Basis Set: " $basis_set >> entry.txt
        if [ "$job_type" == "MD" ]
        then
                echo -e "Ensemble: " $ensemble >> entry.txt
                echo -e "Timestep [fs]: " $timestep >> entry.txt
                echo -e "Thermostat Type: " $thermostat >> entry.txt
                echo -e "Temperature [K]: " $temperature >> entry.txt
        fi
        echo -e "Number of Atoms:" $num_atoms >> entry.txt
    echo -e "Atoms: \n$atoms" >> entry.txt
        echo -e "\n" >> entry.txt

    #Submit Details
    echo -e "@------------------------------------@\n     S U B M I T    D E T A I L S\n@------------------------------------@" >> entry.txt

    echo -e "Requested CPUs: " $cpu_request >> entry.txt
    echo -e "Set Memory Limit [MB]: $requested_mem \n" >> entry.txt

    #Reminder user to update file.
    echo -e "${bold}Entry created! When${green} runlog.log ${reset}${bold}and ${green}orca.xyz${reset}${bold} files are generated, use ${magenta}ue${reset} ${bold}command to update the entry with more information.${reset}"
fi

#Final Export Defaults - NO TOUCHY
sed -i "s/default_material.*/default_material='$default_material'/g" $entry_path/var_storage.sh
sed -i "s/previous_description.*/previous_description='$previous_description'/g" $entry_path/var_storage.sh
sed -i "s/entry_number.*/entry_number=$increment_entry_num /g" $entry_path/var_storage.sh

#End
echo -e "@--------------------------------------------------@"
