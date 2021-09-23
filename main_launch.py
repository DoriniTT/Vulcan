#!/usr/bin/env python

# Author: Thiago Trevizam Dorini

import os
import numpy as np
import time
import json
import math as m
import subprocess as sb
from subprocess import check_output
from ase.spacegroup import get_spacegroup
from ase.visualize import view
from ase.build import surface, bulk, make_supercell, cut, stack, add_vacuum, molecule, rotate, add_adsorbate
from ase import Atoms
from ase.io import read, write
from ase.calculators.vasp import Vasp
from ase.calculators.espresso import Espresso
from ase.db import connect
from ase.constraints import FixAtoms
from ase.io.vasp import read_vasp_xdatcar

###############################---Functions to modify---###############################

###############################################################################
# Do you want to restart from the database?
r_from_data = False
# Give a list of structures id's if you do not want to restart for all.
list_of_structures_id_to_restart = []
###############################################################################

def structures(restart_from_database=r_from_data, visualize=False) -> list:

    global struc

    if visualize == True:
        cwd = '.'
    else:
        cwd = os.environ.get("here")

    try:
        sb.run(["mkdir "+cwd+"/structures"], shell=True)
    except:
        print("Getting structures from the structures folder")

    #################### RESTARTING FROM THE DATABASE ########################

    if restart_from_database == True:

        db = connect('')
        c = db.count()
        if list_of_structures_id_to_restart == True:
            struc = [ db.get_atoms(id=i) for i in range(c+1) ]
        else:
            struc = [ db.get_atoms(id=i) for i in list_of_structures_id_to_restart ]
        #print(struc[0].get_potential_energy())

    ##########################################################################

    else:
        if len(os.listdir(cwd+"/structures")) == 0:

            #################### DEFINE THE STRUCTURES HERE ########################

            struc = [0]

            #######################################################################

            for n, s in enumerate(struc):
                write(cwd+"/structures/"+str(s.get_chemical_formula(mode='hill'))+"_"+str(n)+".vasp", s)
        else:
            import glob
            structures = glob.glob(cwd+'/structures/*.vasp')
            struc = [ read(i, format="vasp") for i in structures ]

    if visualize == True:
        cwd = os.getcwd()
        view(struc)

    mycontent = '''STRUCTURE_SIZE="'''+str(len(struc) - 1)+'''"
'''
    with open("structure_size", "w") as f:
        f.write(mycontent)

    #a1 = bulk('Cu', 'fcc', a=3.6)
    #struc = [a1]
    return struc

def calculators(vasp=False, qe=False):

    if vasp == True:

        abinit = Vasp(istart=0, algo = 'Normal', nelm=40,setups='recommended',\
            ediff=1e-5, ediffg=-0.03, ismear=0, xc = 'pbe',\
            prec = 'Accurate', lreal = 'Auto', kpts=(1,1,1), gamma=True,\
            ibrion=1, isif=2, nsw=500, laechg = True, lorbit = 11, lvtot=True, lwave=False,\
            ncore=1, encut=400)

    elif qe == True:

        pseudopotentials = {'Na': 'Na.pbe-spn-kjpaw_psl.1.0.0.UPF',\
                'Cl': 'Cl.pbe-n-kjpaw_psl.1.0.0.UPF', 'Al': 'Al.pbe-n-kjpaw_psl.1.0.0.UPF',\
                    'Co': 'Co.pbe-spn-kjpaw_psl.0.3.1.UPF', 'O': 'O.pbe-n-kjpaw_psl.1.0.0.UPF',\
                        'H': 'H.pbe-kjpaw_psl.1.0.0.UPF'}

        abinit = Espresso(pseudopotentials=pseudopotentials,\
            tstress=True, tprnfor=True, calculation='relax', etot_conv_thr=1e-5, \
                    forc_conv_thr=0.03, nstep=500, ecutwfc=80, ecutrho=800, occupations='smearing',\
                        degauss=0.002, smearing = 'marzari-vanderbilt', mixing_beta = 0.7, \
                            ion_dynamics = 'damp', electron_maxstep = 40, scf_must_converge=False,\
                                adaptive_thr=True)

    return abinit

###############################---Functions to modify---###############################

class Configure():

    def __init__(self, homedir, workdir, database_name, xc, ecut):

        self.homedir = homedir
        self.workdir = workdir
        self.database_name = database_name
        self.xc = xc
        self.ecut = ecut

    def get_structures(self) -> list:

        struc = structures()
        return struc

    def get_calculators(self, vasp=False, qe=False):

        calc = calculators(vasp=vasp, qe=qe)
        return calc

    def set_env(self, n, vasp=False, qe=False, k=None):

        global calc, cell, namefile, spgroup, unique_identifier

        calc = self.get_calculators(vasp=vasp, qe=qe)
        pwd = os.getcwd()

        # setting variables
        if k != None:
            cell = structures[n][k].cell.cellpar()
            spgroup = get_spacegroup(structures[n][k]).symbol.replace(" ","").replace("/", "_")
            namefile = structures[n][k].get_chemical_formula(mode='hill')
        else:
            cell = structures[n].cell.cellpar()
            spgroup = get_spacegroup(structures[n]).symbol.replace(" ","").replace("/", "_")
            namefile = structures[n].get_chemical_formula(mode='hill')

        unique_identifier = str(namefile)+"_"+str(spgroup)+"_"+self.xc+"_"+str(self.ecut)

        if not os.path.isdir(unique_identifier):
            sb.run(["mkdir "+unique_identifier], shell=True)
        os.chdir(unique_identifier)

    def database(self, structure, namefile):

        global db, ID
        db = connect(self.homedir + "/" + self.database_name)
        try:
            ID = db.get(unique=unique_identifier).id
        except:
            ID = None
        if ID == None:
            db.write(structure, name=str(namefile), xc=self.xc, unique=unique_identifier) # Write in the database.
        else:
            db.update(id=ID, atoms=structure)

    def set_kpts(self, slab=False, precision=6) -> tuple:

        # Defining kpts
        if slab == False:
            k = (max(1,int(precision*2*m.pi/cell[0])), \
                  max(1,int(precision*2*m.pi/cell[1])), \
                  max(1,int(precision*2*m.pi/cell[2])))
        else:
            k = (max(1,int(precision*2*m.pi/cell[0])), \
                  max(1,int(precision*2*m.pi/cell[1])), \
                  1)
        return k

    def recuperar(self, namefile="test", md=False, cwd=None, vasp=False, final=False):

        if not os.path.isdir("backup_"+namefile):
            sb.run(["mkdir backup_"+namefile], shell=True)

        import fnmatch

        if vasp == True:
            inputs = ["XDATCAR", "CONTCAR", "POSCAR", "OUTCAR", "OSZICAR", "INCAR", "KPOINTS"]

            numb = len(fnmatch.filter(os.listdir("backup_"+namefile), 'INCAR*'))
            for i in inputs:
                sb.run(["cp "+i+" backup_"+namefile+"/"+i+""+str(numb)], shell=True)

        if final == True:
            sb.run(["cp CHGCAR backup_"+namefile+"/CHGCAR_final"], shell=True)
            try:
                sb.run(["cp AECCAR0 backup_"+namefile+"/AECCAR0_final"], shell=True)
                sb.run(["cp AECCAR1 backup_"+namefile+"/AECCAR1_final"], shell=True)
                sb.run(["cp AECCAR2 backup_"+namefile+"/AECCAR2_final"], shell=True)
            except:
                pass

        if cwd != None:
            sb.run(["cp -r backup_* " + cwd], shell=True)

class Calculo(Configure):

    def __init__(self, homedir, workdir, database_name, xc, ecut):
        super().__init__(homedir, workdir, database_name, xc, ecut)

    #---TESTING--#
    def check_calculation(self, n, get_energy_function):

        with open(self.homedir+"/jobid_"+str(n), "r") as f:
            jobid = f.read().replace('\n', '')

        done=False
        tries=0
        while done == False and tries < 2:

            get_energy_function

            output = str(check_output(u'sacct -a --jobs {}'.format(jobid), shell=True))
            if u'vasp_gam' in output:
                output = output.split('vasp_gam')[1]
            elif u'vasp_std' in output:
                output = output.split('vasp_std')[1]

            if u'COMPLETED' in output:
                done=True
            elif u'CANCELLED' in output or u'FAILED' in output or u'OUT_OF_ME+' in output:
                tries=tries+1
    #---TESTING---#

    def calculation(self, n, vasp=False, qe=False, adsorption=False, new=False, restart=True, dftu=False, idefix=False, check_calc=False) -> None:

        global energy, s

        if dftu == True:
            ldau_luj = {'In': {'L': 2, 'U': 7.0, 'J': 0}, 'Cr': {'L': 2, 'U': 3.0, 'J': 0}, 'Mn': {'L': 2, 'U': 4.3, 'J': 0}, 'Fe': {'L': 2, 'U': 3.0, 'J': 0}, 'Co': {'L': 2, 'U': 4.3, 'J': 0}, 'Ni': {'L': 2, 'U': 4.3, 'J': 0}}
            calc.set(ldau_luj=ldau_luj, ldautype=2)

        # Manual parameter to change:
        if vasp == True:
            calc.set(xc=self.xc)

        # Launch calculation
        if os.path.exists("CONTCAR") and os.path.getsize("CONTCAR") > 0 and \
              restart == True and new == False or (os.path.exists('result_espresso.out') \
                    and os.path.getsize("result_espresso.out") > 0):
            print("========================Restarting!=========================")
            if vasp == True and adsorption == False:
                print("Here!")
                contcar = read("CONTCAR", format="vasp")
            elif vasp == True and adsorption == True:
                if y == 0:
                    contcar = read("CONTCAR_structure", format="vasp")
                elif y == 1:
                    contcar = read("CONTCAR_substrate", format="vasp")
                elif y == 2:
                    contcar = read("CONTCAR_oxide", format="vasp")
            elif qe == True:
                contcar = read("result_espresso.out", format="espresso-out")
            contcar.calc = calc

            if vasp == True and idefix == False and check_calc == True:
                self.check_calculation(n, contcar.get_potential_energy()) # TESTING!
                if adsorption == False:
                    calc.write_json(str(namefile) + "_state.json")
                    #self.recuperar(namefile=str(unique_identifier), final=True)
            else:
                energy = contcar.get_potential_energy()
                if adsorption == False:
                    calc.write_json(str(namefile) + "_state.json")
                    #self.recuperar(namefile=str(unique_identifier), final=True)
            try:
                self.database(contcar, str(namefile))
            except:
                print("There is an error in the database!")
            s = contcar
        else:
            structures[n].calc = calc
            if vasp == True and idefix==False and check_calc == True:
                self.check_calculation(n, structures[n].get_potential_energy()) # TESTING!
                #self.recuperar(namefile=str(unique_identifier))
            else:
                energy = structures[n].get_potential_energy()
                #self.recuperar(namefile=str(unique_identifier))
                calc.write_json(str(namefile) + "_state.json")
            #energy = structures[n].get_potential_energy()
            try:
                if adsorption == True and y == 0:
                    self.database(structures[n], str(namefile))
                elif adsorption == False:
                    self.database(structures[n], str(namefile))
                else:
                    pass
            except:
                print("There is an error in the database!")
            s = structures[n]

        return s

    def md(self, n, steps_production=5, steps_final=10, temperature_final=5000, vasp=False, qe=False, new=False, idefix=False) -> None:

        if idefix:
            calc.set(command='mpirun -np 12 vasp_gam')
        elif kpts == (1, 1, 1):
            calc.set(command='srun vasp_gam')
        else:
            calc.set(command='srun vasp_std')

        md_dir = self.homedir+"/md_files/"
        sb.run(["mkdir "+md_dir], shell=True)

        #Production phase
        calc.set(nsw=steps_production)
        calc.set(mdalgo=2, ibrion=0, smass=0, tebeg=0, teend=temperature_final,\
                 potim=0.05, algo='Fast', ispin=1, isym=0, nwrite=0, lwave=False,\
                 prec='Accurate')
        try:
            q = read_vasp_xdatcar(file='XDATCAR_production', index=slice(None))
        except:
            q = [0]
        if int(len(q)) + 2 < steps_production:
            new_steps = int(steps_production) - int(len(q))
            calc.set(nsw=new_steps)

            self.calculation(n, vasp=vasp, qe=qe, new=False, idefix=idefix)
            sb.run(["cat XDATCAR >> XDATCAR_production"], shell=True)
            sb.run(["awk <PCDAT >>"+md_dir+"PCDAT_prod."+str(steps_production)+"_"+str(unique_identifier)+"_"+str(n)+"fs ' NR==8 {pcskal=$1} NR==9 {pcfein=$1} NR>=13 {line=line+1; print (line-0.5)*pcfein/pcskal,$1} '"], shell=True)
            sb.run(["grep 'free  energy' OUTCAR | awk ' {print $5}' >> "+md_dir+"prod_energy_"+str(unique_identifier)+"_"+str(n)+".dat"], shell=True)
        else:
            sb.run(["touch "+md_dir+"finished_production_"+str(n)], shell=True)

        #Final phase
        calc.set(nsw=steps_final)
        calc.set(mdalgo=2, ibrion=0, smass=0, tebeg=temperature_final, teend=temperature_final,\
                 potim=0.25, algo='Fast', ispin=1, isym=0, nwrite=0, lwave=False,\
                 prec='Accurate')

        try:
            q = read_vasp_xdatcar(file='XDATCAR_final', index=slice(None))
        except:
            q = [0]
        if int(len(q))+2 < steps_final:
            new_steps = int(steps_final) - int(len(q))
            calc.set(nsw=new_steps)

            self.calculation(n, vasp=vasp, qe=qe, new=False, idefix=idefix)
            sb.run(["cat XDATCAR >> XDATCAR_final"], shell=True)
            sb.run(["awk <PCDAT >>"+md_dir+"PCDAT_final."+str(steps_final)+"_"+str(unique_identifier)+"_"+str(n)+"fs ' NR==8 {pcskal=$1} NR==9 {pcfein=$1} NR>=13 {line=line+1; print (line-0.5)*pcfein/pcskal,$1} '"], shell=True)
            sb.run(["grep 'free  energy' OUTCAR | awk ' {print $5}' >> "+md_dir+"final_energy_"+str(unique_identifier)+"_"+str(n)+".dat"], shell=True)
        else:
            sb.run(["touch "+md_dir+"finished_final_"+str(n)], shell=True)

    def run_step_relax(self, n, vasp=False, qe=False, slab=False, new=False, dftu=False, \
                       md=False, steps_production=50, steps_final=200, temperature_final=300,\
                       gamma=True, relax=True, dos=False, idefix=False, ncore_test=False,\
                       stm=False, eint=[-1], \
                       adsorption=False, plane_of_separation=0, calculate_chgdiff=False, bader=False,\
                       cohp=False, nbands=1000, \
                       results=True) -> None :

        global structures, kpts, y
        structures = self.get_structures()
        self.set_env(n, vasp=vasp, qe=qe)

        kpts = (1, 1, 1)
        if vasp == True:
            inputs = ["INCAR", "KPOINTS", "POTCAR", "IBZKPT", "vasprun.xml", "OUTCAR", "OSZICAR", "CONTCAR", "POSCAR", "XDATCAR"]

        #----------MD---------#
        if md == True:
            if vasp == True:
                self.md(n, vasp=vasp, qe=qe, new=new, idefix=idefix,\
                        steps_production=steps_production, steps_final=steps_final, temperature_final=temperature_final)

                inputs.append("PCDAT")
                for i in inputs:
                    sb.run(["cp "+i+" "+i+"_md"], shell=True)
        #----------MD---------#

        calc.set(nsw=1000, isym=2)
        #--------GAMMA-------#
        if vasp == True:
            calc.set(mdalgo=None, ibrion=1, smass=None, tebeg=None, teend=None, potim=0.5)
        if gamma == True:
            if idefix:
                calc.set(command='mpirun -np 12 vasp_gam')
            else:
                calc.set(command='srun vasp_gam')
            try:
                x = str(check_output('grep "reached required" OUTCAR_gam', shell=True))
            except:
                x = 1
            if isinstance(x, int):
                self.calculation(n, vasp=vasp, qe=qe, new=False, idefix=idefix)
                for i in inputs:
                    sb.run(["cp "+i+" "+i+"_gam"], shell=True)
                if idefix:
                    sb.run(["touch ../finished_calculation"], shell=True)
        #--------GAMMA-------#

        #--------RELAX-------#
        if relax == True:
            kpts = self.set_kpts(slab=slab)

            if vasp == True:
                if slab==True:
                    calc.set(idipol=3, ldipol=True) # Define this!
                else:
                    calc.set(idipol=None, ldipol=None) # Define this!

                if idefix:
                    calc.set(command='mpirun -np 12 vasp_std')
                else:
                    calc.set(command='srun vasp_std')
                calc.set(kpts=kpts)
                try:
                    x = str(check_output('grep "reached required" OUTCAR_relax', shell=True))
                except:
                    x = 1
                if isinstance(x, int):
                    self.calculation(n, vasp=vasp, qe=qe, new=False, idefix=idefix)
                    if idefix:
                        sb.run(["touch ../finished_calculation"], shell=True)
                    for i in inputs:
                        sb.run(["cp "+i+" "+i+"_relax"], shell=True)

            elif qe == True:
                try:
                    x = str(check_output('grep "!" result_expresso.out', shell=True))
                except:
                    x = 1
                if isinstance(x, int):
                    self.calculation(n, vasp=vasp, qe=qe, new=False, idefix=idefix)
        #--------RELAX-------#

        #--------DOS-------#
        if dos == True:
            if idefix:
                calc.set(command='mpirun -np 12 vasp_std')
            else:
                calc.set(command='srun vasp_std')
            kpts = self.set_kpts(slab=slab)
            kpts = (2*kpts[0], 2*kpts[1], 2*kpts[2])
            calc.set(icharg=11, kpts=kpts, ismear=-5, prec='Accurate', nsw=0)
            self.calculation(n, vasp=vasp, qe=qe, new=False)
            if idefix:
                sb.run(["touch ../finished_calculation"], shell=True)
            if vasp == True:
                inputs.append("PROCAR")
                for i in inputs:
                    sb.run(["cp "+i+" "+i+"_dos"], shell=True)

            # Getting DOS results
            if results == True:
                calc_load = Vasp(restart=True, directory=".")
                energies, dos = calc_load.get_dos()

                import matplotlib.pyplot as plt
                fig, ax = plt.subplots()
                ax.plot(energies, dos)
                ax.set(xlabel="Energy (eV)", ylabel="DOS")
                ax.grid()

                dir_dos = self.homedir+"/dos_results/"
                if not os.path.isdir(dir_dos):
                    sb.run(["mkdir "+dir_dos], shell=True)
                fig.savefig(dir_dos+unique_identifier+"_"+str(n)+".pdf")
        #--------DOS-------#

        #--------STM-------#
        if stm == True:
            if idefix:
                calc.set(command='mpirun -np 12 vasp_std')
            else:
                calc.set(command='srun vasp_std')
            kpts = self.set_kpts(slab=slab)
            try:
                x = str(check_output('grep "Total CPU time used" OUTCAR_scf_stm', shell=True))
            except:
                x = 1
            if isinstance(x, int):
                calc.set(nsw=0, nelm=120, kpts=kpts, prec='Accurate')
                self.calculation(n, vasp=vasp, qe=qe, new=False)
                if vasp == True:
                    for i in inputs:
                        sb.run(["cp "+i+" "+i+"_scf_stm"], shell=True)

            try:
                x = str(check_output('grep "Total CPU time used" OUTCAR_par_stm', shell=True))
            except:
                x = 1
            if isinstance(x, int):
                calc.set(nsw=0, kpts=kpts, prec='Accurate', lsepk=False, lsepb=False, lpard=True, nbmod=-3, eint=eint)
                calc.initialize(structures[n])
                calc.write_input(structures[n])
                sb.run(["srun vasp_std"], shell=True)
                #self.calculation(n, vasp=vasp, qe=qe, new=False)

                sb.run(["cp OUTCAR OUTCAR_par_stm"], shell=True)

                if vasp == True:
                    for i in inputs:
                        sb.run(["cp "+i+" "+i+"_par_stm"], shell=True)

            if idefix:
                sb.run(["touch ../finished_calculation"], shell=True)
        #--------STM-------#

        #--------COHP-------#
        if cohp == True:
            if idefix:
                calc.set(command='mpirun -np 12 vasp_std')
            else:
                calc.set(command='srun vasp_std')
            kpts = self.set_kpts(slab=slab)
            kpts = (2*kpts[0], 2*kpts[1], 2*kpts[2])
            calc.set(lwave=True, isym=0, nbands=nbands, prec='Accurate', nsw=0, kpts=kpts, ismear=-5)
            try:
                x = str(check_output('grep "Total CPU time used" OUTCAR_cohp', shell=True))
            except:
                x = 1
            if isinstance(x, int):
                self.calculation(n, vasp=vasp, qe=qe, new=False)

            if vasp == True:
                inputs.append("PROCAR")
                for i in inputs:
                    sb.run(["cp "+i+" "+i+"_cohp"], shell=True)

            lobsterin = '''COHPstartEnergy -10
COHPendEnergy 5
userecommendedbasisfunctions
            '''
            with open("lobsterin", "w") as f:
                f.write(lobsterin)

            sb.run(["lobster"], shell=True)

            # Getting DOS results
            if results == True:
                calc_load = Vasp(restart=True, directory=".")
                energies, dos = calc_load.get_dos()

                import matplotlib.pyplot as plt
                fig, ax = plt.subplots()
                ax.plot(energies, dos)
                ax.set(xlabel="Energy (eV)", ylabel="DOS")
                ax.grid()

                dir_dos = self.homedir+"/dos_results/"
                if not os.path.isdir(dir_dos):
                    sb.run(["mkdir "+dir_dos], shell=True)
                fig.savefig(dir_dos+unique_identifier+"_"+str(n)+".pdf")

            if idefix:
                sb.run(["touch ../finished_calculation"], shell=True)
        #--------COHP-------#

        #--------ADSORPTION-------#
        #Does not work if new=False !
        if adsorption == True:

            inputs.append("CHGCAR")
            inputs.append("AECCAR0")
            inputs.append("AECCAR2")
            kpts = self.set_kpts(slab=slab)
            calc.set(nsw=0, nelm=150)

            if vasp == True:
                if idefix:
                    calc.set(command='mpirun -np 12 vasp_std')
                else:
                    calc.set(command='srun vasp_std')
                calc.set(kpts=kpts, lvtot=True, ldipol=True, idipol=3)

                ## GETTING THE THREE STRUCTURES
                ## IMPORTANT: Define the plane (z) that separates the substrate from the oxide!
                # For the sixton:
                substrate = structures[n].copy()
                del substrate[[atom.index for atom in substrate if atom.position[2]>plane_of_separation]]
                oxide = structures[n].copy()
                del oxide[[atom.index for atom in oxide if atom.position[2]<plane_of_separation]]

                try:
                    y = str(check_output('grep "Total CPU time used" OUTCAR_structure', shell=True))
                except:
                    y = 0
                if isinstance(y, int):
                    try:
                        self.calculation(n, vasp=vasp, qe=qe, adsorption=True, new=True)
                    except:
                        pass
                    for i in inputs:
                        sb.run(["cp "+i+" "+i+"_structure"], shell=True)
                        sb.run(["mv LOCPOT LOCPOT_structure"], shell=True)

                structures[n] = substrate
                calc.set(lvtot=False, ldipol=None, idipol=None)
                try:
                    y = str(check_output('grep "Total CPU time used" OUTCAR_substrate', shell=True))
                except:
                    y = 1
                if isinstance(y, int):
                    try:
                        self.calculation(n, vasp=vasp, qe=qe, adsorption=True, new=True)
                    except:
                        pass
                    for i in inputs:
                        sb.run(["cp "+i+" "+i+"_substrate"], shell=True)

                structures[n] = oxide
                try:
                    y = str(check_output('grep "Total CPU time used" OUTCAR_oxide', shell=True))
                except:
                    y = 2
                if isinstance(y, int):
                    try:
                        self.calculation(n, vasp=vasp, qe=qe, adsorption=True, new=True)
                    except:
                        pass
                    for i in inputs:
                        sb.run(["cp "+i+" "+i+"_oxide"], shell=True)

                    if idefix:
                        sb.run(["touch ../finished_calculation"], shell=True)

                result_substrate = read("OUTCAR_substrate", format="vasp-out")
                result_oxide = read("OUTCAR_oxide", format="vasp-out")
                result_structure = read("OUTCAR_structure", format="vasp-out")
                e_substrate = result_substrate.get_potential_energy()
                e_oxide = result_oxide.get_potential_energy()
                e_structure = result_structure.get_potential_energy()

                e_adsorption = (e_structure - e_oxide - e_substrate)/(cell[0]*cell[1])
                db.update(id=ID, energy_substrate=e_substrate, energy_oxide=e_oxide, adhesive_energy=e_adsorption)

                with open(self.homedir+"/adh_energy.txt", "w") as f:
                    f.write(str(e_adsorption)+' eV/A\n')

                if calculate_chgdiff == True:
                    sb.run(["chgsum.pl CHGCAR_structure CHGCAR_oxide 1 -1; mv CHGCAR_sum chgdiff; chgsum.pl chgdiff CHGCAR_substrate 1 -1; mv CHGCAR_sum chgdiff"], shell=True)

                if bader == True:
                    sb.run(["chgsum.pl AECCAR0_substrate AECCAR2_substrate; mv CHGCAR_sum CHGCAR_diff; bader CHGCAR_substrate -ref CHGCAR_diff; mv ACF.dat ACF_substrate.dat; mv BCF.dat BCF_substrate.dat"], shell=True)
                    sb.run(["chgsum.pl AECCAR0_oxide AECCAR2_oxide; mv CHGCAR_sum CHGCAR_diff; bader CHGCAR_oxide -ref CHGCAR_diff; mv ACF.dat ACF_oxide.dat; mv BCF.dat BCF_oxide.dat"], shell=True)
                    sb.run(["chgsum.pl AECCAR0_structure AECCAR2_structure ; mv CHGCAR_sum CHGCAR_diff; bader CHGCAR_structure -ref CHGCAR_diff; mv ACF.dat ACF_structure.dat; mv BCF.dat BCF_structure.dat"], shell=True)

                savedir = self.homedir+"/CALC_STATES/"+str(unique_identifier)+"_"+str(n)
                if not os.path.isdir(savedir):
                    sb.run(["mkdir -p "+savedir], shell=True)
                sb.run(["cp OUTCAR_substrate OUTCAR_oxide OUTCAR_structure chgdiff LOCPOT_structure ACF_structure.dat BCF_structure.dat ACF_oxide.dat BCF_oxide.dat BCF_substrate.dat ACF_substrate.dat "+savedir], shell=True)
        #--------ADSORPTION-------#

        #--------NCORE TEST-------#
        if ncore_test == True:
            if idefix:
                calc.set(command='mpirun -np 12 vasp_std')
            else:
                calc.set(command='srun vasp_gam')
            ncore_list = [1, 2, 4, 8, 16]
            kpts = (1, 1, 1)
            times = {}
            for i in ncore_list:
                sb.run(["mkdir ncore"+str(i)], shell=True)
                os.chdir("ncore"+str(i))
                calc.set(kpts=kpts, nsw=3, ncore=i)
                time_i = time.time()
                try:
                    self.calculation(n, vasp=vasp, qe=qe, new=False)
                except:
                    print("There was a problem with ncore_"+str(i))
                time_f = time.time()
                times["ncore"+str(i)] = str(round(time_f-time_i,2))+"s"
                os.chdir("..")

            with open(self.homedir + "/" + namefile + "_" + str(n) + "_" +"time.txt", "w") as r:
                r.write(json.dumps(times)+'\n')
        #--------NCORE TEST-------#

        #--------FINAL VERIFICATION--------#
        if md == False or adsorption == False:
            try:
                x = str(check_output('grep "reached required" OUTCAR', shell=True))
            except:
                x = 1
            if isinstance(x, str):
                sb.run(["cp OUTCAR OUTCAR_FINAL"], shell=True)

        # Writing the calculator's state
        savedir = self.homedir+"/CALC_STATES/"+str(unique_identifier)+"_"+str(n)
        if not os.path.isdir(savedir):
            sb.run(["mkdir -p "+savedir], shell=True)
        calc.write_json(savedir+"/"+str(unique_identifier)+"_"+str(n)+"_state.json")

        analyse_script = '''from ase.calculators.vasp import Vasp
from ase.io import read, write

x = Vasp(restart=True, directory='.')
x.read_results()
print(x)
        '''
        with open(savedir+"/analyse_run.py", "w") as f:
            f.write(analyse_script)

        try:
            for i in inputs:
                if relax == True:
                    sb.run(["cp "+i+"_relax "+savedir], shell=True)
                if gamma == True:
                    sb.run(["cp "+i+"_gam "+savedir], shell=True)
                if md == True:
                    sb.run(["cp "+i+"_md "+savedir], shell=True)
                if dos == True:
                    sb.run(["cp "+i+"_dos "+savedir], shell=True)
                if cohp == True:
                    sb.run(["cp "+i+"_cohp "+savedir], shell=True)
                if stm == True:
                    sb.run(["cp "+i+"_par_stm "+savedir], shell=True)
        except:
            sb.run(["cp INCAR KPOINTS POTCAR CONTCAR IBZKPT OUTCAR vasprun.xml "+savedir], shell=True)

        if cohp == True:
            import glob
            outputs_lobster = glob.glob('*.lobster')
            outputs_lobster.append("lobsterout")
            for o in outputs_lobster:
                sb.run(["cp "+o+" "+savedir], shell=True)

        if bader == True and adsorption == False:
            sb.run(["chgsum.pl AECCAR0 AECCAR2; bader CHGCAR -ref CHGCAR_sum; cp ACF.dat BCF.dat "+savedir], shell=True)

###########################  JOB SCRIPTS #################################

    def launch_jeanzay(self, n):

        structures = self.get_structures()
        processors = len(structures[n].get_chemical_symbols())
        if processors > 40:
            while int(processors) % 40 != 0:
                processors -= 1
        else:
            while int(processors) % 40 != 0:
                processors += 1
        #processors = processors*3
        name = structures[n].get_chemical_formula(mode='hill')

        processors=40
        myscript_content='''#!/bin/sh
#SBATCH --nodes='''+str(int(processors/40))+'''               # 1 node reserved
#SBATCH --ntasks-per-node=40    # 40 MPI tasks
##SBATCH --mem=0       # 1 OpenMP thread
#SBATCH --hint=nomultithread    # Disable hyperthreading
#SBATCH -J '''+name+'''
#SBATCH -o '''+name+'''.out
#SBATCH -e '''+name+'''.out
#SBATCH --time=20:00:00
#SBATCH -A fcv@cpu        # Account to specify for multiproject user
##SBATCH --qos=qos_cpu-dev   # Uncomment for runs under 2 hours

module purge
ve
module load vasp/5.4.4-mpi-vtst

## relax
python run_$structure.py
'''
        with open("run_"+str(n)+".slurm", "w") as f:
            f.write(myscript_content)

    def launch_occigen(self, n, partition='HSW24'):

        structures = self.get_structures()
        processors = len(structures[n].get_chemical_symbols())
        if processors > 24:
            while int(processors) % 24 != 0:
                processors -= 1
        else:
            while int(processors) % 24 != 0:
                processors += 1
        #processors = processors*3
        name = structures[n].get_chemical_formula(mode='hill')

        myscript_content='''#!/bin/sh
#SBATCH -J '''+name+'''
#SBATCH --nodes='''+str(int(processors/24))+'''
#SBATCH --tasks='''+str(processors)+'''
#SBATCH --ntasks-per-node=24
#SBATCH --time=23:00:00
#SBATCH --output vasp.output
##SBATCH --constraint=BDW28
#SBATCH --constraint='''+partition+'''
##SBATCH --mem=118000MB

module purge
module load intel/17.2 python/3.6.3

## relax
python run_$structure.py
'''
        with open("run_"+str(n)+".slurm", "w") as f:
            f.write(myscript_content)

    def launch_explor(self, n, partition="std"):

        structures = self.get_structures()
        processors = len(structures[n].get_chemical_symbols())
        if processors > 32:
            while int(processors) % 32 != 0:
                processors -= 1
        else:
            while int(processors) % 32 != 0:
                processors += 1
        name = structures[n].get_chemical_formula(mode='hill')

        myscript_content='''#!/bin/sh
#SBATCH --nodes='''+str(int(processors/32))+'''
##SBATCH -A uwp34
#SBATCH -p '''+partition+'''
#SBATCH -J '''+name+'''
#SBATCH -o '''+name+'''.out
#SBATCH -e '''+name+'''.out
## Number of tasks requested
## Number of tasks requested per node
#SBATCH --ntasks-per-node=32
#SBATCH -t 2-23:00:00
#SBATCH --output output

module purge
module load quantum-espresso/6.6 vasp/5.4.4/std

python run_$structure.py >> output_python
cp '''+name+'''.out '''+self.homedir+'''
'''
        with open("run_"+str(n)+".slurm", "w") as f:
            f.write(myscript_content)

    def launch_cobalt(self, n):

        structures = self.get_structures()
        processors = len(structures[n].get_chemical_symbols())
        if processors > 40:
            while int(processors) % 40 != 0:
                processors += 1
        else:
            while int(processors) % 40 != 0:
                processors += 1
        name = structures[n].get_chemical_formula(mode='hill')

        myscript_content='''#!/bin/bash
#MSUB -q skylake
#MSUB -T 28800
#MSUB -N '''+str(int(processors/40))+'''
#MSUB -n '''+str(int(processors))+'''
#MSUB -r qetest
#MSUB -o log.output
#MSUB -e log.error

module purge
module load espresso/6.4.1
module load python3/3.7.5

python3 run_$structure.py
'''
        with open("run_"+str(n)+".slurm", "w") as f:
            f.write(myscript_content)

if __name__ == "__main__":
    print("Initializing")
    structures(visualize=False)
