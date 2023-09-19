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

# --- Aliases ---

# pdbwizard
alias pdbwizard="python3 ~/pdb_wizard.py"
alias pdb_wizard="python3 ~/pdb_wizard.py"

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
    
elif [ $hostsystem = 'hazel' ]
then

    export WORK="/share/ssp/$USER"
    
fi

source ~/.bash_prompt

