Entry Script written by Angela Shipman.
V.1. compiled 09/12/2023.

@------------------------------------@
   W H A T   A R E   E N T R I E S ?
@------------------------------------@

Entries will allow you to keep notes on your jobs, all while grabbing information about the job from your input files (*.inp) (if it cannot find an *.inp file, it will grab information from the restart files), submit scripts (cp2k*.sh, *.sh), runlogs (runlog.log), and other files depending on the program you're using.

Compatible Programs:      Required Files for ce                   Required Files for ue (All require entry.txt)
    - CP2K              - *.inp or *-1.restart, and *.sh        - runlog.log, stdout.*
    - MPMC              - *.inp, *.pqr, *.sh                    - runlog.log
    - OpenMM            - *.py                                  - customize!
    - Orca              - *.inp, *.xyz                          - runlog.log, orca.xyz

If you do not see your program on this list, feel free to make your own with the OpenMM template in the ce.sh and ue.sh scripts.

@------------------------------------@
     H O W    T O    S E T    U P
@------------------------------------@

Step 1. Place the `entry` directory where you want it. You will want the working directory for this in step 3. The entry directory will contain the following:
        entry
        |
    |---READ_ME.txt (you are here)
        |---all_entries.txt
        |---ae.sh
        |---ce.sh
        |---ue.sh
        |---ne.sh
        |---var_storage.sh

Step 2. Copy the below script into your bashrc:

        #ENTRY COMMANDS
        entry_path=/insert/your/path/to/entry/folders/here/entry
        alias cdentry='cd insert/your/path/to/entry/folders/here/entry'

        ce() { #Create Entry
                export entry_path
                sh $entry_path/var_storage.sh
                sh $entry_path/ce.sh
        }

        ue() { #Update Entry
                export entry_path
                sh $entry_path/var_storage.sh
                sh $entry_path/ue.sh
        }

        ne() { #Notate Entry
                export entry_path
                sh $entry_path/var_storage.sh
                sh $entry_path/ne.sh
        }

        ae() { #Compile all entries into one file
                export entry_path
                sh $entry_path/var_storage.sh
                sh $entry_path/ae.sh
        }

Step 3. Within your bashrc, edit the `entry_path` and `alias cdentry` variables with the working directory where you placed the `entry` directory.
E.g. 
entry_path=/home/amshipma/Programs/entry
alias cdentry='cd /home/amshipma/Programs/entry'

Step 4. Within the var_storage.sh file, customize the following user settings:

USER DIRECTORY: Edit the `user_directory` variable with a directory such that when you want to compile ALL ENTRIES (see `ae` command below), the computer can access any entry you've written on the computer.
For example, I will always write entries on Hazel within my own user directory, so I set:
user_directory=/share/ssp/amshipma

INDEX: Edit the `index` variable with your preferred naming scheme. This will display after the entry number.
For example, I want a different naming scheme depending on which computer I write the entry on. So on my local computer, the index is '-L' and on Hazel, the index is '-H'. The result is that on my local computer, the entries will be titled Entry #418-L, #419-L, etc... and my Hazel entries will be titled Entry #219-H, #220-H, etc. You may want a different organization, perhaps based on project or something, so set this index as you please or leave it blank.



@------------------------------------@
        H O W    T O    U S E
@------------------------------------@

Go to the entry directory - 'cdentry'
        This is an alias that has been set up so you can easily access your entry directory.

Create Entry - 'ce'
        Create an entry using the 'ce' command. You will be prompted with questions that fill the fields of the entry.

Update Entry - 'ue'
        Once the job begins, you may use the 'ue' command to populate information throughout the course of the job. Note that some fields will be empty if the job hasn't finished.

Notate Entry - 'ne'
        At any point after entry creation, you can use the 'ne' command to add any further custom notes, especially conclusions about the completed job. You will also be prompted about any corresponding pages in your lab notebook if that is helpful to you. If not, press enter at this prompt and move on.

All Entries - 'ae'
        It may be useful to compile all entries on the computer into one file to search for items. Use the 'ae' command to overwrite the all_entries.txt file with an updated list of all of your entries.

NOTE: If you do not have a unique index set up, entry numbers may duplicate if you use entries across multiple HPCs (e.g. Bridges AND Hazel). You can update your entry number and index in the var_storage.sh file at any time as needed.


You're now ready to use the entries system! If you have any trouble, e-mail Angela Shipman: amshipma@ncsu.edu
Happy entry writing!
