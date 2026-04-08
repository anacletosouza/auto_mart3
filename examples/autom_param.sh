#!/bin/bash

set -e
set -u
set -o pipefail

# ----------------------------------------
# Configurações principais
# ----------------------------------------
rm -rf automatization_of_potentials

ITER=20
BASEDIR=/grain/anacleto/projects/project_1_virion_simulation/CG_simulations/glycoprotein/ASN_FA2/6_CG_MD_simulation
WORKDIR=$BASEDIR/automatization_of_potentials

NDX_BONDS=$BASEDIR/NDX/bonds.ndx
NDX_ANGLES=$BASEDIR/NDX/angles.ndx
NDX_DIHEDRALS=$BASEDIR/NDX/dihedrals.ndx

SCRIPTS_PATH=$BASEDIR/scripts
MDP=$BASEDIR/mdp

export OMP_NUM_THREADS=12
export GMX_GPU_DD_COMMS=true
export GMX_GPU_PME_PP_COMMS=true
export GMX_FORCE_UPDATE_DEFAULT_GPU=true
export GMX_ENABLE_DIRECT_GPU_COMM=TRUE
export GMX_DISABLE_GPU_TIMING=TRUE

mkdir -p $WORKDIR

# ----------------------------------------
# Criação do diretório inicial iter_0
# ----------------------------------------
if [ ! -d "$WORKDIR/iter_0" ]; then
    mkdir -p $WORKDIR/iter_0

    cp $BASEDIR/MDRUN/topol.top $WORKDIR/iter_0/
    cp $BASEDIR/MDRUN/molecule.gro $WORKDIR/iter_0/
    cp $BASEDIR/MDRUN/solv_ions_CG.gro $WORKDIR/iter_0/
    cp $BASEDIR/MDRUN/ASN_FA2.ndx $WORKDIR/iter_0/
    cp $BASEDIR/MDRUN/posre_ASN_FA2.itp $WORKDIR/iter_0/
    cp $BASEDIR/MDRUN/*.ndx $WORKDIR/iter_0/ 2>/dev/null
    cp $BASEDIR/MDRUN/*.tpr $WORKDIR/iter_0/ 2>/dev/null
    cp $BASEDIR/MDRUN/FA2_final.itp $WORKDIR/iter_0/
    cp -r $BASEDIR/ff_files $WORKDIR/iter_0/
fi

# ----------------------------------------
# Loop principal
# ----------------------------------------
for ((i=0; i<$ITER; i++)); do
    echo "========================================="
    echo "Starting iteration $i"
    echo "========================================="

    CURRENT_ITER=iter_$i
    NEXT_ITER=iter_$((i+1))
    CURDIR=$WORKDIR/$CURRENT_ITER

    cd $CURDIR
    echo "PWD: $(pwd)"

    mkdir -p stat figures

    # -------------------------------
    # MD Simulation
    # -------------------------------
    gmx grompp -f $MDP/minimization.mdp -p topol.top -c solv_ions_CG.gro -o em.tpr -r solv_ions_CG.gro -maxwarn 2
    gmx mdrun -deffnm em -v -pin on -ntmpi 1 -ntomp 12

    gmx grompp -f $MDP/CG_nvt_1000.mdp -p topol.top -c em.gro -o nvt.tpr -n index.ndx -r em.gro -maxwarn 2
    gmx mdrun -deffnm nvt -ntmpi 1 -ntomp 12 -v

    gmx grompp -f $MDP/CG_npt_1000.mdp -p topol.top -c nvt.gro -o npt.tpr -n index.ndx -r nvt.gro -t nvt.cpt -maxwarn 2
    gmx mdrun -deffnm npt -ntmpi 1 -ntomp 12 -v

    gmx grompp -f $MDP/CG_md.mdp -p topol.top -c npt.gro -o md.tpr -n index.ndx -r npt.gro -t npt.cpt -maxwarn 2
    gmx mdrun -deffnm md -ntmpi 1 -ntomp 12 -v

    # -------------------------------
    # Processamento da trajetória
    # -------------------------------
    echo "0" | gmx trjconv -s solv_ions.tpr -f md.xtc -pbc whole -o md_whole.xtc
    echo "1 0" | gmx trjconv -s solv_ions.tpr -f md_whole.xtc -o md_center_fit.xtc -fit rot+trans -n ASN_FA2.ndx

    # -------------------------------
    # Análise de distribuições
    # -------------------------------
    bash $SCRIPTS_PATH/CG_distr.sh \
        $NDX_BONDS $NDX_ANGLES $NDX_DIHEDRALS \
        md_center_fit.xtc solv_ions.tpr
    
    rm -fr *.xtc *.trr nvt* npt* md* *.cpt *.edr *.log \#* solv_ions.tpr topol.top ff_files em* *.xvg *.ndx solv_ions_CG.gro posre_ASN_FA2.itp molecule.gro ASN_FA2.tpr

    echo "Checking distributions..."
    ls -lah _bonds/ _angles/ _dihedrals/ || echo "WARNING: missing distribution dirs"

    # -------------------------------
    # Plot + statistiscs
    # -------------------------------
    python3 $SCRIPTS_PATH/plot_distributions.py \
        --bonds_dir _bonds \
        --angles_dir _angles \
        --dihedrals_dir _dihedrals \
        --figures_dir figures \
        --bond_out stat/bond_${i}.tsv \
        --angle_out stat/angle_${i}.tsv \
        --dihedral_out stat/dihedral_${i}.tsv

    echo "Checking TSV..."
    ls -lah stat/

    if [ ! -f stat/bond_${i}.tsv ]; then
        echo "ERROR: bond_${i}.tsv missing"
        exit 1
    fi

    # -------------------------------
    # Update of the potentials if this is not ultimate
    # -------------------------------
    if [ $i -lt $((ITER-1)) ]; then
        NEXTDIR=$WORKDIR/$NEXT_ITER
        mkdir -p $NEXTDIR

        python3 $SCRIPTS_PATH/generate_itp_update.py \
            --bonds_ref $BASEDIR/AA_REF/bond_ref.tsv \
            --angles_ref $BASEDIR/AA_REF/angle_ref.tsv \
            --dihedrals_ref $BASEDIR/AA_REF/dihedral_ref.tsv \
            --bonds_sim stat/bond_${i}.tsv \
            --angles_sim stat/angle_${i}.tsv \
            --dihedrals_sim stat/dihedral_${i}.tsv \
            --atoms_json $JSON_TOPOLOGY_MARTINI3 \
            --ndx_bounds $NDX_BONDS \
            --ndx_angles $NDX_ANGLES \
            --ndx_dihedrals $NDX_DIHEDRALS \
            --molecule_name FA2 \
            --dihedrals_target \
            --itp_out $NEXTDIR/FA2_final.itp

        # copy to next iteration
        cp $BASEDIR/MDRUN/topol.top $NEXTDIR/
        cp $BASEDIR/MDRUN/molecule.gro $NEXTDIR/
        cp $BASEDIR/MDRUN/solv_ions_CG.gro $NEXTDIR/
        cp $BASEDIR/MDRUN/*.ndx $NEXTDIR/ 2>/dev/null
        cp $BASEDIR/MDRUN/posre_ASN_FA2.itp $NEXTDIR/
        cp -r $BASEDIR/ff_files $NEXTDIR/
        cp $BASEDIR/MDRUN/*.tpr $NEXTDIR/
    fi

done

echo "========================================="
echo "DONE"
echo "Final result: $WORKDIR/iter_$((ITER-1))/carb_optimized.itp"
