#NOTATE ENTRY (ne)

#Fancy Formatting
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
echo -e "@--------------------------------------------------@"
echo "${green}If applicable, which ${bold}lab notebook${reset}${green} and which ${bold}page of this lab notebook${reset}${green} is this job discussed on?${reset}"
    read pg_number
echo "${bold}${green}Enter any notes, conclusions or thoughts to the end of this entry.${reset}"
    read conclusions

if grep -q "N O T E S  &  C O N C L U S I O N S" entry.txt
then
    echo -e "NOTE: $conclusions" >> entry.txt
    echo -e "Corresponding Lab Notebook and Page # (if any): $pg_number" >> entry.txt
    echo -e "Notes signed (DD/MM/YYYY): $(date '+%d/%m/%Y %H:%M:%S') by $user\n" >> entry.txt

else
    echo -e "\n@-------------------------------------@\n  N O T E S  &  C O N C L U S I O N S\n@-------------------------------------@" >> entry.txt
    echo -e "NOTE: $conclusions" >> entry.txt
    echo -e "Corresponding Lab Notebook and Page # (if any): $pg_number" >> entry.txt
    echo -e "Notes signed (DD/MM/YYYY): $(date '+%d/%m/%Y %H:%M:%S') by $user\n" >> entry.txt
fi

echo -e "${bold}${green}Entry updated with notes!${reset}"
echo -e "@--------------------------------------------------@"
