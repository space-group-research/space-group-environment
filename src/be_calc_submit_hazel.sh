#!/bin/bash
#Script by Angela Shipman, v2-HAZEL Aug 21, 2026
#Now supports both LSF and Slurm!

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

#=====================================================================
# SCHEDULER HELPERS (LSF/bsub or Slurm/sbatch)
#=====================================================================

#Detect whether a submit script is written for LSF (#BSUB) or Slurm (#SBATCH)
detect_scheduler() {
        local script="$1"
        if grep -q "#SBATCH" "$script" 2>/dev/null; then
                echo "slurm"
        elif grep -q "#BSUB" "$script" 2>/dev/null; then
                echo "lsf"
        else
                echo "unknown"
        fi
}

#Sets the global $scheduler / $scheduler_name and pulls cores/hrs/queue/memory
#out of $submit_script, using the right directives for whichever scheduler
#was detected.
parse_submit_script() {
        scheduler=$(detect_scheduler "$submit_script")

        if [ "$scheduler" = "slurm" ]; then
                scheduler_name="Slurm"
                cores=`grep "#SBATCH --ntasks=" $submit_script | sed 's/.*--ntasks=//'`
                hrs=`grep "#SBATCH --time=" $submit_script | sed 's/.*--time=//'`
                queue=`grep "#SBATCH --partition=" $submit_script | sed 's/.*--partition=//'`
                memory=`grep "#SBATCH --mem=" $submit_script | sed 's/.*--mem=//' | sed 's/[^0-9]*$//'`
        elif [ "$scheduler" = "lsf" ]; then
                scheduler_name="LSF"
                cores=`grep "#BSUB -n" $submit_script | awk '{print $3}'`
                hrs=`grep "#BSUB -W" $submit_script | awk '{print $3}'`
                queue=`grep "#BSUB -q" $submit_script | awk '{print $3}'`
                memory=`grep "#BSUB -M" $submit_script | awk '{print substr($3, 1, length($3)-3)}'`
        else
                echo -e "${red}${bold}ERROR:${reset}${red} Could not detect a supported scheduler in ${bold}$submit_script${reset}${red}.\nMake sure it contains either ${bold}#BSUB${reset}${red} (LSF) or ${bold}#SBATCH${reset}${red} (Slurm) directives.${reset}"
                exit
        fi
}

#Writes a job name into a submit script using the right directive for the
#detected scheduler
set_job_name() {
        local script="$1"
        local name="$2"
        if [ "$scheduler" = "slurm" ]; then
                sed -i "s/#SBATCH --job-name=.*/#SBATCH --job-name=${name}/g" "$script"
        else
                sed -i "s/#BSUB -J.*/#BSUB -J ${name}/g" "$script"
        fi
}

#Submits a job using the right command for the detected scheduler
submit_job() {
        local script="$1"
        if [ "$scheduler" = "slurm" ]; then
                sbatch "$script"
        else
                bsub < "$script"
        fi
}

#Rewrites cores/hours/queue/memory in guest.sh using the right directives
#for the detected scheduler. Assumes guest.sh already exists (copied from
#$submit_script) before this is called.
update_guest_script() {
        if [ "$scheduler" = "slurm" ]; then
                sed -i "s/#SBATCH --time=.*/#SBATCH --time=$guest_hours/g" guest.sh
                sed -i "s/#SBATCH --ntasks=.*/#SBATCH --ntasks=$guest_cores/g" guest.sh
                sed -i "s/#SBATCH --partition=.*/#SBATCH --partition=$guest_queue/g" guest.sh
                sed -i "s/#SBATCH --mem=.*/#SBATCH --mem=${guest_memory}G/g" guest.sh
        else
                sed -i "s/#BSUB -W.*/#BSUB -W $guest_hours/g" guest.sh
                sed -i "s/#BSUB -n.*/#BSUB -n $guest_cores/g" guest.sh
                sed -i "s/#BSUB -q.*/#BSUB -q $guest_queue/g" guest.sh
                sed -i "s/#BSUB -M.*/#BSUB -M ${guest_memory}GB!/g" guest.sh
                sed -i "s/\[mem\=.*GB\/task\]/\[mem\=${guest_memory}GB\/task\]/g" guest.sh
        fi
}

#=====================================================================

#Get user info
echo -e "${bold}${yellow}NOTICE:${reset}${yellow} Make sure you've run your *.xyz files through PDB Wizard's SORT option (available in version 0.3.1) so that all guest molecules are organized at the bottom of the page.\n${reset}"

if [  -f *.cp ] ;
then
        cp_file=`ls *.cp`
        submit_script=`awk 'f;/submit_script/{f=1}' $cp_file | head -1`

        job_type=`grep "RUN_TYPE" *.inp | awk '{print $2}'`
        host_multiplicity=`grep "MULTIPLICITY" *.inp | awk '{print $2}'`
        parse_submit_script

        atoms1_mof=`awk 'f;/mof_elements/{f=1}' $cp_file | head -1`
        guest_cores=`awk 'f;/guest_cores/{f=1}' $cp_file | head -1`
        guest_hours=`awk 'f;/guest_hours/{f=1}' $cp_file | head -1`
        guest_queue=`awk 'f;/guest_queue/{f=1}' $cp_file | head -1`
        guest_memory=`awk 'f;/guest_memory/{f=1}' $cp_file | head -1`
        guest=`awk 'f;/guest_name/{f=1}' $cp_file | head -1`
        num_guest_atoms=`awk 'f;/num_guest_atoms/{f=1}' $cp_file | head -1`
        num_guests=`awk 'f;/num_guests/{f=1}' $cp_file | head -1`
        guest_elements=`awk 'f;/guest_elements/{f=1}' $cp_file | head -1`
        guest_multiplicity=`awk 'f;/guest_multiplicity/{f=1}' $cp_file | head -1`

        echo -e "${green}Detected *.cp file. Does the following all sound correct?${reset}\n\n\n${bold}SUBMIT SCRIPT: ${reset}$submit_script\n${bold}SCHEDULER DETECTED:${reset} $scheduler_name\n${bold}JOB TYPE:${reset} $job_type\n\n${bold}${cyan}=====HOST DETAILS=====${reset}"
        echo -e "${bold}MULTIPLICITY: ${reset}$host_multiplicity \n${bold}JOB CORES:${reset} $cores \n${bold}JOB HOURS: ${reset}$hrs \n${bold}JOB QUEUE: ${reset}$queue\n${bold}MEMORY (GB) PER CORE: ${reset}$memory\n${bold}ELEMENTS: ${reset}$atoms1_mof\n"
        echo -e "\n${cyan}${bold}=====GUEST DETAILS=====${reset}"
        echo -e "${bold}NAME:${reset} $guest\n${bold}MULTIPLICITY: ${reset}$guest_multiplicity \n${bold}JOB CORES:${reset} $guest_cores \n${bold}JOB HOURS: ${reset}$guest_hours \n${bold}JOB QUEUE: ${reset}$guest_queue\n${bold}MEMORY (GB) PER CORE: ${reset}$guest_memory\n${bold}ELEMENTS: ${reset}$guest_elements\n${bold}NUM ATOMS:${reset} $num_guest_atoms\n${bold}NUM MOLECULES:${reset} $num_guests"
        echo -e "\n\n${bold}${yellow}PROCEED?${reset}"
        select yn in "Yes" "No"; do
                case $yn in
                        Yes) break;;
                        No) exit;;
                esac
        done
else
        echo -e "No *.cp file"
        echo -e "${green}Please type the ${bold}full filename of your submit script${reset}${green}. e.g. Type: cp2k.sh ${reset}"
        read submit_script
        echo -e "\n\n${bold}${red}WARNING:${reset}${red} Make sure you have the following:\n1) Your *.xyz file that's been run through PDB wizard's SORT option.\n2) Your *.inp file with the host's multiplicity\n3) Your $submit_script script (LSF ${bold}#BSUB${reset}${red} or Slurm ${bold}#SBATCH${reset}${red} format) which has the desired cores, queue, and run times for your HOST jobs. You will be prompted for GUEST job details.\n\n${reset}${bold}DO YOU HAVE ALL OF THESE ITEMS?${reset}"
        select yn in "Yes" "No"; do
                        case $yn in
                        Yes) break;;
                        No) exit;;
                esac
        done

        if [ ! -f $submit_script ] ;
        then
                echo "${red}You are missing a submit script. Please create a submit script with the desired cores, queue, and and run times for your HOST jobs that will be generated. (You will be prompted for GUEST JOB details.)${reset}"
                exit
        fi
        if [ ! -f *.xyz ] ;
        then
                echo "${red}You are missing your start-pos.xyz file. Please make sure you run it through PDB Wizard's SORT option before using it here.${reset}"
                exit
        fi
        if [ ! -f *.inp ] ;
        then
                echo "${red}You are missing your *.inp file. Be sure to include one here with your HOST'S multiplicity.${reset}"
                exit
        fi

        job_type=`grep "RUN_TYPE" *.inp | awk '{print $2}'`
        host_multiplicity=`grep "MULTIPLICITY" *.inp | awk '{print $2}'`
        parse_submit_script

        echo -e "\n${yellow}According to your current files, the 2 host jobs will be submitted with the following:\n${reset}${bold}SCHEDULER DETECTED:${reset} $scheduler_name\n${bold}JOB TYPE:${reset} $job_type \n${bold}HOST MULTIPLICITY: ${reset}$host_multiplicity \n${bold}HOST JOB CORES:${reset} $cores \n${bold}HOST JOB HOURS: ${reset}$hrs \n${bold}HOST JOB QUEUE: ${reset}$queue\n${bold}HOST JOB MEMORY PER CORE: ${reset}${memory} GB\n${yellow}${bold}Proceed?${reset}"
        select yn in "Yes" "No"; do
                case $yn in
                        Yes) break;;
                        No) exit;;
                esac
        done

        echo -e "\n${green}${bold}How many cores ${reset}${green}would you like to allocate for the ${bold}guest-only jobs?${reset}"
        read guest_cores

        if [ "$scheduler" = "slurm" ]; then
                time_hint="Use the format HH:MM:SS. e.g. 2 hrs would be: 02:00:00"
        else
                time_hint="Use the format HR:MM. e.g. 2 hrs would be: 02:00"
        fi
        echo -e "\n${yellow}${bold}How many hours ${reset}${yellow}would you like to allocate for the ${bold}guest-only jobs?${reset}${yellow} $time_hint"
        read guest_hours

        if [ "$scheduler" = "slurm" ]; then
                queue_hint="This will be used as the Slurm --partition."
        else
                queue_hint="Check spelling!"
        fi
        echo -e "\n${yellow}${bold}What queue ${reset}${yellow}would you like to submit the ${bold}guest-only jobs to?${reset}${yellow} $queue_hint"
        read guest_queue

        echo -e "\n${green}${bold}How many GB of memory ${reset}${green}would you like to allocate to the ${bold}guest-only jobs${reset}${green}? Type an integer."
        read guest_memory

        echo -e "\n${green}What is your ${bold}guest ${reset}${green}molecule? ${bold}Type a name.${reset}${green}\n(This will be used to name the guest molecule directories.)"
        read guest

        echo -e "\n${yellow}What ${bold}elements${reset}${yellow} are in ${bold}$guest${reset}${yellow}? \n${bold}Type in elements separated with spaces.${reset}${yellow} Please match the case sensitivity as the elements appear in your input script.\nEx #1. For H2O, type: H O\nEx #2. For N2, type: N"
        read guest_elements

        echo -e "\n${green}How many ${bold}total atoms${reset}${green} are in an individual ${bold}$guest${reset}${green}? ${bold}Type an integer.${reset}${green}\nEx #1. For H2O, type: 3\nEx #2. For N2, type: 2"
        read num_guest_atoms

        echo -e "\n${yellow}How many ${bold}total molecules${reset}${yellow} of ${bold}$guest${reset}${yellow} are in the system that you're trying to calculate binding energies for?${bold} Type an integer.${reset}${yellow}\nEx #1. If there are 12 waters in your system, type: 12\nEx #2. If there are 423 N2 in your system, type: 423"
        read num_guests

        echo -e "\n${green}For a lone ${bold}$guest${reset}${green} molecule, what is it's ${bold}multiplicity? Type an integer.${reset}${green} \ne.g. H2O has multiplicity of 1."
        read guest_multiplicity

        echo -e "\n${yellow}What ${bold}elements${reset}${yellow} is your host made of (negate ${bold}$guest${reset}${yellow} atoms)? ${bold}Type in elements separated by spaces,${reset}${yellow} and match the case sensitivity of the elements listed in your input script. e.g. H O C N Cu"
        read atoms1_mof

fi

total_jobs=$(( $num_guests + 2 ))

echo -e "${cyan}${bold}\n\n\n\n\nFINAL CHECK! ${reset}\n${cyan}Does this all sound right?\n\nAccording to your input, a total of ${bold}$total_jobs jobs ${reset}${cyan}will be submitted using ${bold}$scheduler_name${reset}${cyan}. \n\n${bold}HOST JOBS:${reset}${cyan} Host jobs will be submitted to the ${bold}$queue queue${reset}${cyan}, all of which will run for ${bold}$hrs hrs${reset}${cyan} with ${bold}$cores cores${reset}${cyan} (${memory}GB per core${reset}${cyan}). The host system is comprised of ${bold}$atoms1_mof elements ${reset}${cyan}for a total ${bold}host multiplicity of $host_multiplicity. \n\n${bold}GUEST JOBS:${reset}${cyan} There are ${bold}$num_guests $guest molecules${reset}${cyan} in your host system. The ${bold}$guest ${reset}${cyan}molecule is made of ${bold}$num_guest_atoms atoms${reset}${cyan}, which comprise of ${bold}$guest_elements elements${reset}${cyan}, for a total ${bold}guest multiplicity of $guest_multiplicity${reset}${cyan}. Guest jobs will be submitted to the ${bold}$guest_queue queue${reset}${cyan} with ${bold}$guest_cores cores${reset}${cyan} (${guest_memory}GB per core for LSF, ${guest_memory}GB total for Slurm) for ${bold}$guest_hours hours. \n\n${yellow}SUBMIT JOBS?${reset}"
select yn in "Submit!" "Cancel"; do
        case $yn in
                Submit!) break;;
                Cancel) exit;;
        esac
done

#turn guest_elements into array
guest_elements2="$guest_elements"
set -- $guest_elements2

#turn atoms1_mof into an array
atoms2_mof="$atoms1_mof"
set -- $atoms2_mof

#Math for finding lines where waters are.
total_lines=`cat *.xyz | wc -l`
lines_guest=$(( $num_guests * $num_guest_atoms ))
mof_lines=$(( $total_lines - $lines_guest - 2 ))

#create mof and everything_real directories with corresponding ghosts
mkdir all_real
cp *.inp *.xyz $submit_script all_real
cd all_real
set_job_name $submit_script all_real
submit_job $submit_script
cd ..
echo "Submitting job for all real atoms..."

mkdir host_only
cp *.inp *.xyz $submit_script host_only
cd host_only
base=$(( 2 + $mof_lines ))
host=$(( base + 1 ))
for ((i=$host; i<=$total_lines; i++));
do
        for elem in $guest_elements2;
        do
                sed -i "${i}s/${elem} /${elem}_ghost /g" *.xyz
        done
done
set_job_name $submit_script host_real
submit_job $submit_script
cd ..
echo "Submitting job for host-only..."
#make a version of xyz file that is all ghost.
cp *.xyz all_ghost.xyz

for i in $atoms2_mof;
do
        sed -i "s/$i /${i}_ghost /g" all_ghost.xyz
done

#Update input script for guests with multiplicity:
cp *.inp guest.inp;
sed -i "s/MULTIPLICITY.*/MULTIPLICITY $guest_multiplicity/g" guest.inp;

#update guest script for guest hrs/cores/queue/memory (scheduler-specific):
cp $submit_script guest.sh
update_guest_script

#create guest directories with corresponding ghosts
for ((i=1;i<=$num_guests;i++));
do
        skip=$(( $i * $num_guest_atoms ))
        base_skip=$(( $base + $skip ))

        mkdir ${i}${guest}
        cp guest.sh guest.inp all_ghost.xyz ${i}${guest}
        cd ${i}${guest}

        for ((x=1;x<=$num_guest_atoms;x++));
        do
                line_guest=$(( $base_skip - $x + 1 ))
                for elem in $guest_elements2;
                do
                        sed -i "${line_guest}s/${elem}_ghost /${elem} /g" all_ghost.xyz
                done
        done

        mv all_ghost.xyz start-pos.xyz
        mv guest.sh $submit_script
        set_job_name $submit_script ${i}${guest}
        submit_job $submit_script
        cd ..
        echo -e "Submitting job for $guest calculation #$i..."
done
rm guest.inp all_ghost.xyz guest.sh

echo "${bold}${green}Complete! All energy jobs running!${reset}"
