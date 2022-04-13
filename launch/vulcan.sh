#!/usr/bin/env sh
# Author: Thiago Trevizam Dorini
######################################
#-------Parameters------#
PROJECT="PDIN"
namework="surf_energy_final_pdin"
machine="explor"
partition_explor="mysky" # On Explor --"std", "sky", "mysky", "freestdepyc", "mystdcasijl", "freestdcas"
partition_occigen="HSW24" # On Occigen -- "BDW28", "HSW24"

vasp="True"
qe="False" # Por enquanto so funciona com a funcao "relax". Deixar todas as outras em "False"
#------Calculation------#
new="True"
ncore_test="False"
dftu="False"

slab="False"
bader="False"

#########
md="False"
steps_production="50"
steps_final="200"
temperature_final="300"
##
gamma="False"
relax="True"
scf="False"
hse="False"
dos="False"
##
cohp="False"
nbands="500"
##
stm="False"
##
adsorption="False"
plane_of_separation="24" # In the z direction
calculate_chgdiff="False"
##
enthalpy="False"
coh="False"
##
surf_energy="False"
count_relation="1"
numb_of_layers="8"
##
relax_far="False"
##
#-------Parameters------#
encut="500"
xc="pbe"
######RESULTS######
results="False" # Sometimes it does not work, use it carefully.
######################################
here=$(pwd)
export here
if [[ $machine == "jeanzay" ]]; then
    work="$SCRATCH/$PROJECT/${namework}"
elif [[ $machine == "idefix" ]]; then
    work="/users/trevizam1/$PROJECT/${namework}"
else
    work="$SCRATCHDIR/$PROJECT/${namework}"
fi
database="${namework}.db"

python vulcan.py
. ./structure_size
source ~/.bashrc

if [[ $enthalpy == "True" ]]; then
    echo -e "from vulcan import *\nx = Configure('$here' ,'$work', '$database', '$xc', $encut)\nx.get_bulks()" > bulks.py
    python bulks.py
    rm bulks.py
fi

for structure in $( seq 0 $STRUCTURE_SIZE )
do
    export structure

    #----------Creating the vasp.slurm----------#
    if [[ $machine == "explor" ]]; then
        echo -e "from vulcan import *\nx = Calculo('$here' ,'$work', '$database', '$xc', $encut)\nx.launch_explor($structure, partition='$partition_explor')" > launch_$structure.py

    elif [[ $machine == "occigen" ]]; then
        echo -e "from vulcan import *\nx = Calculo('$here' ,'$work', '$database', '$xc', $encut)\nx.launch_occigen($structure, partition='$partition_occigen')" > launch_$structure.py

    elif [[ $machine == "jeanzay" ]]; then
        echo -e "from vulcan import *\nx = Calculo('$here' ,'$work', '$database', '$xc', $encut)\nx.launch_jeanzay($structure)" > launch_$structure.py

    elif [[ $machine == "cobalt" ]]; then
        echo -e "from vulcan import *\nx = Calculo('$here' ,'$work', '$database', '$xc', $encut)\nx.launch_cobalt($structure)" > launch_$structure.py
    fi
    #----------Creating the vasp.slurm----------#

    if [[ $machine == "idefix" ]]; then
        echo -e "from vulcan import *\nx = Calculo('$here' ,'$work', '$database', '$xc', $encut)\nx.run_step_relax($structure, vasp=$vasp, qe=$qe, slab=$slab, md=$md, steps_production=float('$steps_production'), steps_final=float('$steps_final'), temperature_final=float('$temperature_final'), gamma=$gamma, relax=$relax, scf=$scf, hse=$hse,  dos=$dos, stm=$stm, cohp=$cohp, nbands=float('$nbands'),
        adsorption=$adsorption, plane_of_separation=float('$plane_of_separation'), enthalpy=$enthalpy, coh=$coh, surf_energy=$surf_energy, count_relation=int('$count_relation'), numb_of_layers=int('$numb_of_layers'), calculate_chgdiff=$calculate_chgdiff, bader=$bader, idefix=True, relax_far=$relax_far, ncore_test=$ncore_test, new=$new, dftu=$dftu, results=$results)" > run_$structure.py

    else
        echo -e "from vulcan import *\nx = Calculo('$here' ,'$work', '$database', '$xc', $encut)\nx.run_step_relax($structure, vasp=$vasp, qe=$qe, slab=$slab, md=$md, steps_production=float('$steps_production'), steps_final=float('$steps_final'), temperature_final=float('$temperature_final'), gamma=$gamma, relax=$relax, scf=$scf, hse=$hse,  dos=$dos, stm=$stm, cohp=$cohp, nbands=float('$nbands'),
        adsorption=$adsorption, plane_of_separation=float('$plane_of_separation'), enthalpy=$enthalpy, coh=$coh, surf_energy=$surf_energy, count_relation=int('$count_relation'), numb_of_layers=int('$numb_of_layers'), calculate_chgdiff=$calculate_chgdiff, relax_far=$relax_far, bader=$bader, ncore_test=$ncore_test, new=$new, dftu=$dftu, results=$results)" > run_$structure.py
    fi

    mkdir -p $work; mv launch_$structure.py run_$structure.py $work; cp *.py *.vasp *.cif *.db $work 2>/dev/null; cd $work
    echo -e "$here" > homepath

    #------CANCEL SCRIPT-----#

    #cat > verify.sh < EOF
    ##!/usr/bin/env sh
    #total_time=20 #hours -- total time

    #while [ total_time -gt 2 ]
    #do
    #sleep 3600
    #let "total_time=total_time-1"
    #done
    #EOF

    #------LAUNCH------#
    if [[ $machine != "idefix" ]]; then
        python launch_$structure.py
    fi

    if [[ $partition_explor == "mysky" && $machine == "explor" ]]; then
        echo "explor!"
        sbatch -A uwp34 run_$structure.slurm > output_sbatch; awk '{ print $4 }' output_sbatch > $here/jobid_$structure

    elif [[ $machine == "cobalt" ]]; then
        echo "submitting in cobalt!"
        ccc_msub run_$structure.slurm > output_sbatch; awk '{ print $4 }' output_sbatch > $here/jobid_$structure

    elif [[ $machine == "idefix" ]]; then
        echo "submitting in idefix!"
        python run_$structure.py

    else
        echo "Running on Jeanzay or Explor!"
        sbatch run_$structure.slurm > output_sbatch; awk '{ print $4 }' output_sbatch > $here/jobid_$structure
    fi
    #------LAUNCH------#

    cat run_$structure.py >> resume_run_and_launch
    cat launch_$structure.py >> resume_run_and_launch
    cat run_$structure.slurm >> resume_run_and_launch
    rm launch_$structure.py run_$structure.slurm

    cd $here

    echo -e "$work" > workpath

done
