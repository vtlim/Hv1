
# Purpose:  Combine NAMD fepout files for bar4fep.py processing.
# Usage:    combineExtended.sh orig.fepout extended.fepout combined.fepout 

mkdir ../results-extended

head -n-1 $1 > $3   # take orig file without the "Free energy change" last line
tail -n+4 $2 >> $3  # take new file without the "NEW FEP WINDOW" header

cp $3 ../results-extended

