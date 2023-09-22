#USER DEFINED VARIABLES

export user_directory='/update/me/with/your/path/'

if [ $hostsystem = 'chem' ]
then
    export user_directory="/home/$USER"
elif [ $hostsystem = 'bridges2' ]
then
    export user_directory="/ocean/projects/che220043p/$USER"
elif [ $hostsystem = 'hazel' ]
then
    export user_directory="/share/ssp/$USER"
fi

#This sets the directory that ae.sh will search for ALL entries you've ever made on the machine. On your local computer, this can be '~' to search the entire computer. For an HPC, you'll likely want to set this to your user directory.
#Common User Directories:
#Hazel: '/share/ssp/YourUnityID'
#Bridges: '/ocean/projects/che220043p/YourUsername'

export index="-L"
#The index gives unique identifiers to the names of your entries. For example, I use "L" for entries created on my local computer, "H" on Hazel, and "B" on Bridges. You may want to use one based on project, material, etc... do what suits you! Leave it blank if you don't want an index at all.

#AUTO-GENERATED (No need to touch!)
export default_material=''
export entry_number=1 
export previous_description=''
