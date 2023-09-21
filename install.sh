#!/usr/bin/bash

echo "installing ..."

mkdir -p ~/bin

echo "~/pdb_wizard.py"
cp src/pdb_wizard.py ~/bin

echo "gnuplot_runlog.sh"
cp src/gnuplot_runlog.sh ~/bin

echo "~/.spacebashrc"
cp src/spacebashrc.sh ~/.spacebashrc

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

