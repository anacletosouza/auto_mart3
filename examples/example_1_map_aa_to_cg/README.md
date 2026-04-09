Examples of usage

In pathway bayesian_potentials/examples/example_1_map_aa_to_cg:

(1) create one folder "results"
mkdir -p results

(2) enter into results
cd results

(3) if you want to make step-by-step, we suggest you apply the module map followed by gen-top. 
    Otherwise, you may run all steps with the module pipeline.

### `bayesian-potentials map`

Maps an atomistic trajectory to coarse-grained coordinates.

```bash
bayesian-potentials map \
                        --index_cg       ../ndx/cg.ndx   \
                        --aa_tpr         ../aa_md_data/md.tpr \
                        --aa_xtc         ../aa_md_data/md.xtc \
                        --output_mapped  mapped.xtc      \
                        --output_cg_gro  cg.gro          \
                        --remove_pbc                     \
                        --verbose
```

### `bayesian-potentials gen-top`

Generates a CG topology file (ITP) for GROMACS.

```bash
bayesian-potentials gen-top \
                            --path_ff           ../ff_files                    \
                            --ff                martini_v3.0.0.itp             \
                            --ions              martini_v3.0.0_ions_v1.itp     \
                            --solvent           martini_v3.0.0_solvents_v1.itp \
                            --itp_ligand        cg.itp                         \
                            --name_molecule     "molecule"                     \
                            --number_molecule   1                              \
                            --output_topol      topol.top
```

### `bayesian-potentials pipeline`

**Complete pipeline** - runs mapping, ITP generation, topology creation, and `grompp` in one command.

```bash
bayesian-potentials pipeline \
                             --cg_ndx             ../ndx/cg.ndx \
                             --aa_tpr             ../aa_md_data/md.tpr \
                             --aa_xtc             ../aa_md_data/md.xtc \
                             --output_mapped      mapped.xtc \
                             --output_cg_gro      cg.gro \
                             --aa_itp             ../aa_md_data/carb.itp \
                             --output_cg_itp      cg.itp \
                             --input_mdp          ../mdp/minimization.mdp \
                             --output_dir         ../results/ \
                             --output_topol       topol.top \
                             --path_ff            ../ff_files \
                             --name_molecule      "molecule" \
                             --number_molecule    1
