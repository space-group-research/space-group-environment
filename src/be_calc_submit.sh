#!/bin/bash
#Script by Angela Shipman, v1 Nov 20, 2023

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

#Get user info
echo -e "${bold}${yellow}NOTICE:${reset}${yellow} Make sure you've run your *.xyz files through PDB Wizard's SORT option (available in version 0.3.1) so that all guest molecules are organized at the bottom of the page.\n${reset}"

echo -e "${green}Please type the ${bold}full filename of your submit script${reset}${green}. e.g. Type: cp2k.sh ${reset}"
read submit_script

echo -e "${bold}${red}WARNING:${reset}${red} Make sure you have the following:\n1) Your *.xyz file that's been run through PDB wizard's SORT option.\n2) Your *.inp file with the host's multiplicity\n3) Your $submit_script script which has the desired cores, queue, and run times for ALL of the jobs.\n\n${reset}${bold}DO YOU HAVE ALL OF THESE ITEMS?${reset}"
select yn in "Yes" "No"; do
                case $yn in
                Yes) break;;
                No) exit;;
        esac
done

if [ ! -f $submit_script ] ;
then
        echo "${red}You are missing a submit script. Please create a submit script with the desired cores, queue, and and run times for ALL of the jobs that will be generated.${reset}"
        exit
fi
if [ ! -f *.xyz ] ;
then
        echo "${red}You are missing your start-pos.xyz file. Please make sure you run it through PDB Wizard's SORT option before using it here.${reset}"
        exit
fi
if [ ! -f *.inp ] ;
then
        echo "${red}You are missing your *.inp file. Be sure to include one here with your MOF's multiplicity.${reset}"
        exit
fi

job_type=`grep "RUN_TYPE" *.inp | awk '{print $2}'`
host_multiplicity=`grep "MULTIPLICITY" *.inp | awk '{print $2}'`
cores=`grep "#BSUB -n" $submit_script | awk '{print $3}'
hrs=`grep "#BSUB -W" $submit_script | awk '{print $3}'`
queue=`grep "#BSUB -q" $submit_script | awk '{print $3}'`

echo -e "\n${yellow}According to your current files, all jobs will be submitted with the following:\n${reset}${bold}JOB TYPE:${reset} $job_type \n${bold}HOST MULTIPLICITY: ${reset}$host_multiplicity \n${bold}CORES:${reset} $cores \n${bold}HOURS: ${reset}$hrs \n${bold}QUEUE: ${reset}$queue \n\n${yellow}${bold}Proceed?${reset}"
select yn in "Yes" "No"; do
        case $yn in
                Yes) break;;
                No) exit;;
        esac
done

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

total_jobs=$(( $num_guests + 2 ))

echo -e "${cyan}${bold}\n\n\n\n\nFINAL CHECK! ${reset}\n${cyan}Does this all sound right?\n\nAccording to your input, a total of ${bold}$total_jobs jobs ${reset}${cyan}will be submitted to the ${bold}$queue queue${reset}${cyan}, all of which will run for ${bold}$hrs hrs${reset}${cyan} with ${bold}$cores cores.${reset}${cyan}  \n\nThere are ${bold}$num_guests $guest molecules${reset}${cyan} in your host system. The ${bold}$guest ${reset}${cyan}molecule is made of ${bold}$num_guest_atoms atoms${reset}${cyan}, which comprise of ${bold}$guest_elements elements${reset}${cyan}, for a total ${bold}guest multiplicity of $guest_multiplicity${reset}${cyan}. The host system is comprised of ${bold}$atoms1_mof elements ${reset}${cyan}for a total ${bold}host multiplicity of $host_multiplicity. \n\nSUBMIT JOBS?${reset}"
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
lines_guest=$(( $num_guests*$num_guest_atoms ))
mof_lines=$(( $total_lines - $lines_guest - 2 ))

#create mof and everything_real directories with corresponding ghosts
mkdir all_real
cp *.inp *.xyz $submit_script all_real
cd all_real
sed -i "s/#BSUB -J.*/#BSUB -J all_real/g" $submit_script
#bsub < $submit_script
cd ..
echo "Submitting job for all real atoms..."
mkdir host_only
cp *.inp *.xyz $submit_script host_only
cd host_only
base=$(( 2 + $mof_lines ))
for ((i=$base; i<=$total_lines; i++));
do
        for elem in $guest_elements2;
        do
                sed -i "${i}s/ ${elem} / ${elem}_ghost /g" *.xyz
        done
done
sed -i "s/#BSUB -J.*/#BSUB -J host_real/g" $submit_script
#bsub < $submit_script
cd ..
echo "Submitting job for host-only..."
#make a version of xyz file that is all ghost.
cp *.xyz all_ghost.xyz

for i in $atoms2_mof;
do
        sed -i "s/ $i / ${i}_ghost /g" all_ghost.xyz
done

#Update input script for guests with multiplicity:
cp *.inp guest.inp;
sed -i "s/MULTIPLICITY.*/MULTIPLICITY $guest_multiplicity/g" guest.inp;

#create guest directories with corresponding ghosts
for ((i=1;i<=$num_guests;i++));
do
        skip=$(( $i * $num_guest_atoms ))
        base_skip=$(( $base + $skip ))

        mkdir ${i}${guest}
        cp $submit_script guest.inp all_ghost.xyz ${i}${guest}
        cd ${i}${guest}

        for ((x=1;x<=$num_guest_atoms;x++));
        do
                line_guest=$(( $base_skip - $x + 1 ))
                for elem in $guest_elements2;
                do
                        sed -i "${line_guest}s/ ${elem}_ghost / ${elem} /g" all_ghost.xyz
                done
        done

        mv all_ghost.xyz start-pos.xyz
        sed -i "s/#BSUB -J.*/#BSUB -J ${i}${guest}/g" $submit_script
        #bsub < $submit_script
        cd ..
        echo -e "Submitting job for $guest calculation #$i..."
done
rm guest.inp all_ghost.xyz

echo "${bold}${green}Complete! All energy jobs running!${reset}"
