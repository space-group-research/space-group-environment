#!/usr/bin/bash


# --- Exit early if not interactive ---
[[ $- == *i* ]] || return

# --- Determine Host System ---

hostname=`hostname`

if echo $hostname | grep --quiet "chem.ncsu.edu"
then
    hostsystem='chem'
elif echo $hostname | grep --quiet "bridges2.psc.edu"
then
    hostsystem='bridges2'
elif [ $hostname = 'login01' ] || [ $hostname = 'login02' ] || [ $hostname = 'login03' ] || [ $hostname = 'login04' ]
then
    hostsystem='hazel'
else
    hostsystem='unknown'
fi

export hostsystem

# --- Aliases ---

# pdbwizard and other programs

export PATH=$PATH:~/bin

alias pdbwizard="python3 ~/bin/pdb_wizard.py"
alias pdb_wizard="python3 ~/bin/pdb_wizard.py"
alias plot='bash ~/bin/gnuplot_runlog.sh'
alias bashrc='vim ~/.bashrc'
alias bashup='source ~/.bashrc'

# Easier navigation: .., ..., and ~
alias ..="cd .."
alias ...="cd ../.."
alias ~="cd ~"

# List dir contents aliases
# ref: http://ss64.com/osx/ls.html
# Long form no user group, color
alias ls="ls --color=auto"
alias l="ls -oG"
# Order by last modified, long form no user group, color
alias lt="ls -toG"
# List all except . and ..., color, mark file types, long form no user group, file size
alias la="ls -AGFoh"
# List all except . and ..., color, mark file types, long form no use group, order by last modified, file size
alias lat="ls -AGFoth"

# --- Bash Exports ---

# Make vim the default editor
export EDITOR="vim"

# Ignore duplicate commands in the history
export HISTCONTROL=ignoredups

# Increase the maximum number of lines contained in the history file
# (default is 500)
export HISTFILESIZE=10000

# Increase the maximum number of commands to remember
# (default is 500)
export HISTSIZE=10000

# Don't clear the screen after quitting a manual page
export MANPAGER="less -X"

# Make new shells get the history lines from all previous
# shells instead of the default "last window closed" history
export PROMPT_COMMAND="history -a; $PROMPT_COMMAND"

#ENTRY COMMANDS
entry_path=~/bin/entry
alias cdentry='cd ~/bin/entry'

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

# --- File Management ---

deletecorefiles(){
        find -type f -name 'core.*'
        echo -e "THE ABOVE FILES WILL BE DELETED. CONTINUE? ( Y or N )?"
        read answer1
        case $answer1 in
                Yes | Y )
                        echo -e "ARE YOU SURE? ಠ_ಠ ( Y or N )?"
                        read answer2
                        case $answer2 in
                                Yes | Y )
                                        find -type f -name 'core.*' -delete
                                        echo "CORE FILES DELETED."
                                ;;
                                No | N )
                                ;;
                        esac
                        ;;
                No | N )
                        ;;
        esac
}


# --- System Specific Exports ---

if [ $hostsystem = 'chem' ]
then

    if [ `alias | grep "alias hazel=" | wc -l` = 0 ]
    then
        alias hazel="ssh login.hpc.ncsu.edu"
    fi
    if [ `alias | grep "alias bridges=" | wc -l` = 0 ]
    then
        alias bridges="ssh bridges2.psc.edu"
    fi
    
    export CUDA_VISIBLE_DEVICES=0
    export OMP_NUM_THREADS=12
    
    ulimit -s unlimited
    
    play() {
        xdg-open $1
    }
    
    open() {
        xdg-open $1
    }

    fromhazel(){
        scp -r $USER@login.hpc.ncsu.edu:$1 ./new
        mv new/* ./
        rm -r new/
    }

    tohazel(){
        echo "Enter path: "
        read destination
        scp -r $1 $USER@login.hpc.ncsu.edu:$destination/$1
    }

    tobridges2(){
        echo "Enter path: "
        read destination
        scp -r $1 $USER@bridges2.psc.edu:$destination/$1
    }

    update_hazel(){
        printf "Enter path on Hazel to sync to CURRENT directory:\n%s\n" `pwd`
        read source_path
        rsync -avh $USER@login.hpc.ncsu.edu:$source_path/* .
    }

    update_bridges(){
        printf "Enter path on Bridges2 to sync to CURRENT directory:\n%s\n" `pwd`
        read source_path
        rsync -avh $USER@bridges2.psc.edu:$source_path/ .
    }
    
elif [ $hostsystem = 'bridges2' ]
then

    export WORK="/ocean/projects/che220043p"
    alias cdw='cd $WORK'
    alias be_calc='sh ~/bin/be_calc_submit_bridges.sh'
    alias binding_energy='sh ~/bin/be_calc_submit_bridges.sh'
    
    # Print job ID and directory
    bjobsdir() {
        all_jobs=`squeue -u $USER | wc -l`
        all_jobs=`echo "$all_jobs - 1" | bc`

        for i in `squeue -u $USER | awk '{print $1}'| tail -$all_jobs`
        do
                dir=`scontrol show job $i| grep "WorkDir="`
                echo $i $dir
        done
    }
    
elif [ $hostsystem = 'hazel' ]
then

    export WORK="/share/ssp/$USER"
    export LSB_BJOBS_FORMAT="id user queue:6 stat:4 exec_host:8 job_name:30 submit_time start_time"
    alias cdw='cd $WORK'
    alias wjobs='bjobs -w'
    alias memory='bjobs -r -X -o "jobid queue cpu_used run_time avg_mem max_mem slots "'
    alias space='qstat -a space'
    alias job_count='qstat -a | grep -c $USER'
    alias touch=~/bin/touch
    alias be_calc='sh ~/bin/be_calc_submit_hazel.sh'
    alias binding_energy='sh ~/bin/be_calc_submit_hazel.sh'
    
    # Print job ID and directory
    bjobsdir() {
        no_jobs=`bjobs | wc -l`
        all_jobs=`echo "$no_jobs - 1" | bc`
        printf "jobID\tSTATUS\tDIRECTORY\n"
        for i in `bjobs | awk '{print $1}'| tail -$all_jobs`
        do
                stat=`bjobs $i | awk 'NR > 1 {print $3}'`
                big=`bjobs -l $i | grep -A 3 "CWD"`
                small=${big#*CWD }
                smaller=`echo $small|sed "s/ *//g"`
                smallest=${smaller%%>*}
                final=${smallest##*<}
                printf "%s\t%s\t%s\n" $i $stat $final
                #echo $i $stat $final
        done
            }
            
    # Kill jobs located in current directory
    bkill_here()
    {
        base=`pwd`
        job_id=`bjobsdir | grep $base | awk {'print $1'}`

        for i in $job_id; do
            bkill $i
        done
    }

    # Check jobs running in current directory
    jobs_here()
    {
        base=`pwd`
        job_id=`bjobsdir | grep $base | awk {'print $1'} | sort`

        for i in $job_id; do
            bjobs -w $i | grep $i | awk {'print $1, $4, $7'}
        done
    }
    
fi

source ~/.bash_prompt

