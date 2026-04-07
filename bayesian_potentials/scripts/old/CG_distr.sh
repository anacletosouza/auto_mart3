#!/bin/bash

# ============================
# ARGUMENTS
# ============================
BONDS_NDX="$1"
ANGLES_NDX="$2"
DIHEDRALS_NDX="$3"
XTC_FILE="$4"
TPR_FILE="$5"
OUT_PREFIX="$6"

# ============================
# Bonds and constraints
# ============================
DIR="${OUT_PREFIX}_bonds"
rm -rf $DIR
mkdir $DIR

grep -E '^[[:space:]]*[0-9]+[[:space:]]+[0-9]+' $BONDS_NDX > bonds_list.txt
NBONDS=$(wc -l < bonds_list.txt)
echo "Found $NBONDS bonds"

REPORT="$DIR/report_bonds.txt"
rm -f $REPORT
touch $REPORT

IBOND=0
while [ $IBOND -lt $NBONDS ]; do
    bead_pair=$(sed -n "$((IBOND+1))p" bonds_list.txt)
    echo "Processing bond $IBOND: $bead_pair"
    
    echo "[ bond_$IBOND ]" > temp_bond_$IBOND.ndx
    echo $bead_pair >> temp_bond_$IBOND.ndx
    
    echo 0 | gmx distance -f $XTC_FILE \
                          -n temp_bond_$IBOND.ndx \
                          -s $TPR_FILE \
                          -oall $DIR/bond_$IBOND.xvg \
                          -xvg none &> $DIR/bond_$IBOND.log
    
    if [ $? -ne 0 ] || [ ! -f $DIR/bond_$IBOND.xvg ]; then
        echo "ERROR on bond $IBOND" >> $DIR/errors.log
        cat $DIR/bond_$IBOND.log >> $DIR/errors.log
    else
        echo "---- bond $IBOND (${bead_pair}) ----" >> $DIR/data_bonds.txt
        awk '/Average distance/ {print $3} /Standard deviation/ {print $3}' $DIR/bond_$IBOND.log >> $DIR/data_bonds.txt
        gmx analyze -f $DIR/bond_$IBOND.xvg -dist $DIR/distr_bond_$IBOND.xvg -xvg none -bw 0.001 

        echo "$IBOND: $bead_pair" >> $REPORT
    fi
    
    rm -f temp_bond_$IBOND.ndx
    let IBOND=$IBOND+1
done
rm -f bonds_list.txt


# ============================
# Angles
# ============================
DIR="${OUT_PREFIX}_angles"
rm -rf $DIR
mkdir $DIR

grep -E '^[[:space:]]*[0-9]+[[:space:]]+[0-9]+[[:space:]]+[0-9]+' $ANGLES_NDX > angles_list.txt
NANGLES=$(wc -l < angles_list.txt)
echo "Found $NANGLES angles"

REPORT="$DIR/report_angles.txt"
rm -f $REPORT
touch $REPORT

IANG=0
while [ $IANG -lt $NANGLES ]; do
    bead_trio=$(sed -n "$((IANG+1))p" angles_list.txt)
    echo "Processing angle $IANG: $bead_trio"
    
    echo "[ angle_$IANG ]" > temp_angle_$IANG.ndx
    echo $bead_trio >> temp_angle_$IANG.ndx
    
    echo 0 | gmx angle -f $XTC_FILE \
                       -n temp_angle_$IANG.ndx \
                       -ov $DIR/ang_$IANG.xvg &> $DIR/ang_$IANG.log
    
    if [ $? -ne 0 ] || [ ! -f $DIR/ang_$IANG.xvg ]; then
        echo "ERROR on angle $IANG" >> $DIR/errors.log
        cat $DIR/ang_$IANG.log >> $DIR/errors.log
    else
        echo "---- ang $IANG (${bead_trio}) ----" >> $DIR/data_angles.txt
        awk '/< angle >/ {print $5} /Std. Dev./ {print $4}' $DIR/ang_$IANG.log >> $DIR/data_angles.txt
        gmx analyze -f $DIR/ang_$IANG.xvg -dist $DIR/distr_ang_$IANG.xvg -xvg none -bw 1.0 

        echo "$IANG: $bead_trio" >> $REPORT
    fi
    
    rm -f temp_angle_$IANG.ndx
    let IANG=$IANG+1
done
rm -f angles_list.txt


# ============================
# Dihedrals
# ============================
DIR="${OUT_PREFIX}_dihedrals"
rm -rf $DIR
mkdir $DIR

grep -E '^[[:space:]]*[0-9]+[[:space:]]+[0-9]+[[:space:]]+[0-9]+[[:space:]]+[0-9]+' $DIHEDRALS_NDX > dihedrals_list.txt
NDIHEDRALS=$(wc -l < dihedrals_list.txt)
echo "Found $NDIHEDRALS dihedrals"

REPORT="$DIR/report_dihedrals.txt"
rm -f $REPORT
touch $REPORT

IDIH=0
while [ $IDIH -lt $NDIHEDRALS ]; do
    bead_quartet=$(sed -n "$((IDIH+1))p" dihedrals_list.txt)
    echo "Processing dihedral $IDIH: $bead_quartet"
    
    valid=1
    for num in $bead_quartet; do
        if [ $num -lt 1 ] || [ $num -gt 30 ]; then
            valid=0
            echo "ERROR: bead $num out of range (1-30)" >> $DIR/errors.log
        fi
    done
    
    if [ $valid -eq 1 ]; then
        echo "[ dihedral_$IDIH ]" > temp_dih_$IDIH.ndx
        echo $bead_quartet >> temp_dih_$IDIH.ndx
        
        echo 0 | gmx angle -type dihedral \
                           -f $XTC_FILE \
                           -n temp_dih_$IDIH.ndx \
                           -ov $DIR/dih_$IDIH.xvg &> $DIR/dih_$IDIH.log
        
        if [ $? -ne 0 ] || [ ! -f $DIR/dih_$IDIH.xvg ]; then
            echo "ERROR on dihedral $IDIH" >> $DIR/errors.log
            cat $DIR/dih_$IDIH.log >> $DIR/errors.log
        else
            echo "---- dih $IDIH (${bead_quartet}) ----" >> $DIR/data_dihedrals.txt
            awk '/< angle >/ {print $5} /Std. Dev./ {print $4}' $DIR/dih_$IDIH.log >> $DIR/data_dihedrals.txt
            gmx analyze -f $DIR/dih_$IDIH.xvg -dist $DIR/distr_dih_$IDIH.xvg -xvg none -bw 1.0 

            echo "$IDIH: $bead_quartet" >> $REPORT
        fi
    fi
    
    rm -f temp_dih_$IDIH.ndx
    let IDIH=$IDIH+1
done
rm -f dihedrals_list.txt

echo "Processing completed!"
echo "Results in: ${OUT_PREFIX}_bonds/, ${OUT_PREFIX}_angles/, ${OUT_PREFIX}_dihedrals/"
echo "Check errors.log in each directory if needed."
