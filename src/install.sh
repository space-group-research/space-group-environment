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

echo "installing on $hostsystem ..."

mkdir -p ~/bin
mkdir -p ~/bin/entry

if [ $hostsystem = 'hazel' ]
then
    echo "~/bin/touch"
    gcc src/touch.c -o ~/bin/touch

    echo "~/be_calc_submit_hazel.sh"
    cp src/be_calc_submit_hazel.sh ~/bin
fi

if [ $hostsystem = 'bridges2' ]
then
    echo "~/be_calc_submit_bridges.sh"
    cp src/be_calc_submit_bridges.sh ~/bin
fi

echo "~/pdb_wizard.py"
cp src/pdb_wizard.py ~/bin

echo "gnuplot_runlog.sh"
cp src/gnuplot_runlog.sh ~/bin

echo "~/.spacebashrc"
cp src/spacebashrc.sh ~/.spacebashrc

echo "entry"
cp src/entry/ae.sh ~/bin/entry
cp src/entry/ce.sh ~/bin/entry
cp src/entry/ne.sh ~/bin/entry
cp src/entry/ue.sh ~/bin/entry
cp src/entry/READ_ME.txt ~/bin/entry
if ! [ -f ~/bin/entry/var_storage.sh ]
then
    cp src/entry/var_storage.sh ~/bin/entry
fi

if ! grep --quiet "source ~/.spacebashrc" ~/.bashrc
then
    echo "adding hook in ~/.bashrc"
    echo "" >> ~/.bashrc
    echo "source ~/.spacebashrc" >> ~/.bashrc
    echo "" >> ~/.bashrc
fi

if ! [ -f ~/.bash_prompt ]
then
    echo "~/.bash_prompt"
    cp src/bash_prompt.sh ~/.bash_prompt
fi

if ! [ -f ~/.vmdrc ]
then
    echo "~/.vmdrc"
    cp src/vmdrc ~/.vmdrc
fi

echo "done!"
echo "you can reload your environment by running 'source ~/.bashrc' or 'bashup'"
echo "or opening a new terminal"


