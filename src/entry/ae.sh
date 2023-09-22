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


user=`last | grep -o $USER | head -1`
source $entry_path/var_storage.sh

echo -e "@--------------------------------------------------@"
echo -e "${bold}Compiling entries into one file...${reset}"
echo -e "COMPLETE ENTRY LIST BY: $user" > $entry_path/all_entries.txt
echo -e "Compiled on: $(date '+%d/%m/%Y %H:%M:%S')" >> $entry_path/all_entries.txt


find $user_directory -name entry.txt | xargs cat > $entry_path/tmp.txt
grab=`gawk -v RS="\n\n\n\n" -v ORS="\n\n\n\n\n\0" 1 $entry_path/tmp.txt |sort -z |tr -d 'q\0'`
echo -e "$grab" >> $entry_path/all_entries.txt
rm $entry_path/tmp.txt

echo -e "${bold}All entries have been compiled in file ${green}all_entries.txt${reset}${bold} in ${magenta}$entry_path${reset}.\n${green}Remember you can access this folder with ${bold}cdentry${reset}${green}.${reset}"
echo -e "@--------------------------------------------------@"
