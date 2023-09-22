#!/usr/bin/bash

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

# Easier navigation: .., ..., and ~
alias ..="cd .."
alias ...="cd ../.."
alias ~="cd ~"

# List dir contents aliases
# ref: http://ss64.com/osx/ls.html
# Long form no user group, color
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


# --- System Specific Exports ---

if [ $hostsystem = 'chem' ]
then

    alias hazel="ssh login.hpc.ncsu.edu"
    alias bridges="ssh bridges2.psc.edu"
    
    export CUDA_VISIBLE_DEVICES=0
    export OMP_NUM_THREADS=12
    
    ulimit -s unlimited
    
elif [ $hostsystem = 'bridges2' ]
then

    export WORK="/ocean/projects/che220043p"
    alias cdw='cd $WORK'
    
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
    
    # Print job ID and directory
    bjobsdir()
    {
        job_id=`qstat -a | grep $USER | awk {'print $3'}`

        for i in $job_id; do
            job_dir=`echo $i | xargs bjobs -l | grep -A 1 "CWD" | head -n 2 | paste -d " " - - | grep -o "CWD <*[/a-zA-Z0-9.-\_]* *[/a-zA-Z0-9.-\_]*>" | sed "s/ //g" | sed "s/CWD//g" | sed "s/^<//g" | sed "s/>$//g"`
            echo $i $job_dir
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

