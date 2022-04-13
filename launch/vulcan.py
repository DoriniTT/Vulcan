#!/usr/bin/env python

# Author: Thiago Trevizam Dorini

import os
import numpy as np
import time
import json
import math as m
import pint
import glob
import subprocess as sb
import datetime
from launch.vulcan_helper import * # My functions!
from launch.vulcan_structures import * # My functions!
from launch.vulcan_custom_calculators import * # My functions!
from collections import defaultdict
from subprocess import check_output
from ase.spacegroup import get_spacegroup
from ase.visualize import view
from ase.build import surface, bulk, make_supercell, cut, stack, add_vacuum, molecule, rotate, add_adsorbate, fcc111
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
name_of_database = ''
###############################################################################

def structures(restart_from_database=r_from_data) -> list:

    global struc

    #################### RESTARTING FROM THE DATABASE ########################

    if restart_from_database == True:

        db = connect(name_of_database)
        c = db.count()
        if list_of_structures_id_to_restart == True:
            struc = [ db.get_atoms(id=i) for i in range(c+1) ]
        else:
            struc = [ db.get_atoms(id=i) for i in list_of_structures_id_to_restart ]

    ##########################################################################

    #################### DEFINE THE STRUCTURES HERE ########################

    s = read("pdin.vasp")
    s1 = s.copy()
    s1 = supercell(s1, matrix=[1, 1, 8])
    s1 = cut(s1, nlayers=13)
    add_vacuum(s1, 10)
#     view(s1)

    s2 = s.copy()
    s2 = cut(s2, origo=[0.5, 0.5, 0.5])
    s2 = supercell(s2, matrix=[1, 1, 8])
    s2 = cut(s2, nlayers=13)
    add_vacuum(s2, 10)
#     view(s2)

    struc = [s1, s2]

    #######################################################################

    mycontent = '''STRUCTURE_SIZE="'''+str(len(struc) - 1)+'''"
'''
    with open("structure_size", "w") as f:
        f.write(mycontent)

    return struc

###############################---Functions to modify---###############################

class Configure():

    def __init__(self, homedir, workdir, database_name, xc, encut):

        self.homedir = homedir
        self.workdir = workdir
        self.database_name = database_name
        self.xc = xc
        self.encut = encut

    def get_structures(self) -> list:

        struc = structures()
        return struc

    def get_calculators(self, vasp=False, qe=False):

        calc = general_calculator(vasp=vasp, qe=qe, encut=self.encut)
        return calc

    def verify_outcar(self, outcar_name, reached=False, time=False):

        try:
            if reached == True:
                x = str(check_output('grep "reached required" '+outcar_name, shell=True))

            elif time == True:
                x = str(check_output('grep "Total CPU time used" '+outcar_name, shell=True))
        except:
            x = 1

        return x

    def set_env(self, n, vasp=False, qe=False, k=None):

        global cell, namefile, spgroup, unique_identifier

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

        unique_identifier = str(namefile)+"_"+str(spgroup)+"_"+self.xc+"_"+str(self.encut)+"_"+str(n)

        if not os.path.isdir(unique_identifier):
            sb.run(["mkdir "+unique_identifier], shell=True)
        os.chdir(unique_identifier)

    def get_bulks(self):

        # Script by Florian Brix.
        from pymatgen.ext.matproj import MPRester

        if not os.path.isdir(self.homedir+"/bulk_for_enthalpies"):
            sb.run(["mkdir "+self.homedir+"/bulk_for_enthalpies"], shell=True)

        if len(os.listdir(self.homedir+'/bulk_for_enthalpies') ) == 0:

            structures = self.get_structures()
            for structure in structures:
                namefile = structure.get_chemical_formula(mode='hill')
                sb.run(["mkdir "+self.homedir+"/bulk_for_enthalpies/"+str(namefile)], shell=True)
                formule=structure.get_chemical_symbols()
                try:
                  for j in range(len(formule)):
                     if j>0 and formule[j]==formule[j-1]:
                        pass
                     else:
                         bulk=MPRester(api_key='jMQJq6g41qTTTTqqqnoI')
                         sortie=bulk.get_data(str(formule[j]), prop="e_above_hull")
                         for k in range(len(sortie)):
                             if sortie[k].get('e_above_hull')==0:
                                 meilleurbulkid=(sortie[k].get('material_id'))
                         sortie=bulk.get_data(str(formule[j]),data_type='vasp', prop="cif")
                         for k in range(len(sortie)):
                             if sortie[k].get('material_id')==meilleurbulkid:
                                 meilleurbulkcif=sortie[k].get('cif')
                         filename=str(formule[j]+'.cif')
                         with open(self.homedir+"/bulk_for_enthalpies/"\
                                   +str(namefile)+"/"+filename, 'w') as g:
                             g.write((meilleurbulkcif))
                except:
                    print('le code doit avoir acces a internet pour acceder  aux donnees bulk l execution une fois hors noeuds de calcul peut suffire')
                    pass

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

    def set_kpts(self, slab=False, precision=7) -> tuple:

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

class Calculo(Configure):

    def __init__(self, homedir, workdir, database_name, xc, encut):
        super().__init__(homedir, workdir, database_name, xc, encut)

    #---TESTING--#
    #def check_calculation(self, n, get_energy_function):
    def check_calculation(self, n):

        with open(self.homedir+"/jobid_"+str(n), "r") as f:
            jobid = f.read().replace('\n', '')

        done=False
        tries=0
        while done == False and tries <= 2:

            get_energy_function

            output = str(check_output(u'sacct -a --jobs {}'.format(jobid), shell=True))
            output = output.split('\\n')[-2]

            if u'COMPLETED' in output:
                done=True
            else:
                tries = tries+1
                if tries == 2:
                    sb.run(["cd "+self.homedir+"/; sh vulcan.sh"], shell=True)
                    sb.run(["scancel "+str(jobid)], shell=True)

    #---TESTING---#

    def calculation(self, n, vasp=False, qe=False, adsorption=False, new=False, restart=True, dftu=False, idefix=False, check_calc=False) -> None:

        global energy, s

        if dftu == True:
            ldau_luj = {'Ti': {'L': 2, 'U': 1.0, 'J': 0}, 'In': {'L': 2, 'U': 7.0, 'J': 0}, 'Cr': {'L': 2, 'U': 3.0, 'J': 0}, 'Mn': {'L': 2, 'U': 4.3, 'J': 0}, 'Fe': {'L': 2, 'U': 3.0, 'J': 0}, 'Co': {'L': 2, 'U': 4.3, 'J': 0}, 'Ni': {'L': 2, 'U': 4.3, 'J': 0}}
            calc.set(ldau_luj=ldau_luj, ldautype=2)

        # Manual parameter to change:
        if vasp == True:
            calc.set(xc=self.xc)

        # Launch calculation
        if os.path.exists("CONTCAR") and os.path.getsize("CONTCAR") > 0 and \
              restart == True or (os.path.exists('result_espresso.out') \
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
            else:
                energy = contcar.get_potential_energy()
                if adsorption == False:
                    calc.write_json(str(namefile) + "_state.json")
            try:
                self.database(contcar, str(namefile))
            except:
                print("There is an error in the database!")
            s = contcar

        elif new:
            structures[n].calc = calc
            if vasp == True and idefix==False and check_calc == True:
                self.check_calculation(n, structures[n].get_potential_energy()) # TESTING!
            else:
                energy = structures[n].get_potential_energy()
                calc.write_json(str(namefile) + "_state.json")
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

        return energy

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
                 potim=1, algo='Normal', ispin=1, isym=0, nwrite=0, lwave=False,\
                 prec='Accurate')
        try:
            q = read_vasp_xdatcar(file='XDATCAR_production', index=slice(None))
        except:
            q = [0]
        if int(len(q)) + 2 < steps_production:
            new_steps = int(steps_production) - int(len(q))
            calc.set(nsw=new_steps)

            self.calculation(n, vasp=vasp, qe=qe, new=new, idefix=idefix)
            sb.run(["cat XDATCAR >> XDATCAR_production"], shell=True)
            sb.run(["awk <PCDAT >>"+md_dir+"PCDAT_prod."+str(new_steps)+"_"+str(unique_identifier)+"_"+str(n)+"fs ' NR==8 {pcskal=$1} NR==9 {pcfein=$1} NR>=13 {line=line+1; print (line-0.5)*pcfein/pcskal,$1} '"], shell=True)
            sb.run(["grep 'free  energy' OUTCAR | awk ' {print $5}' >> "+md_dir+"prod_energy_"+str(unique_identifier)+"_"+str(n)+".dat"], shell=True)
        else:
            sb.run(["touch "+md_dir+"finished_production_"+str(n)], shell=True)

        #Final phase
        calc.set(nsw=steps_final)
        calc.set(mdalgo=2, ibrion=0, smass=0, tebeg=temperature_final, teend=temperature_final,\
                 potim=2, algo='Fast', ispin=1, isym=0, nwrite=0, lwave=False,\
                 prec='Accurate')

        try:
            q = read_vasp_xdatcar(file='XDATCAR_final', index=slice(None))
        except:
            q = [0]
        if int(len(q))+2 < steps_final:
            new_steps = int(steps_final) - int(len(q))
            calc.set(nsw=new_steps)

            self.calculation(n, vasp=vasp, qe=qe, new=new, idefix=idefix)
            sb.run(["cat XDATCAR >> XDATCAR_final"], shell=True)
            sb.run(["awk <PCDAT >>"+md_dir+"PCDAT_final."+str(new_steps)+"_"+str(unique_identifier)+"_"+str(n)+"fs ' NR==8 {pcskal=$1} NR==9 {pcfein=$1} NR>=13 {line=line+1; print (line-0.5)*pcfein/pcskal,$1} '"], shell=True)
            sb.run(["grep 'free  energy' OUTCAR | awk ' {print $5}' >> "+md_dir+"final_energy_"+str(unique_identifier)+"_"+str(n)+".dat"], shell=True)
        else:
            sb.run(["touch "+md_dir+"finished_final_"+str(n)], shell=True)

    def run_step_relax(self, n, vasp=False, qe=False, slab=False, new=False, dftu=False, \
                       md=False, steps_production=50, steps_final=200, temperature_final=300,\
                       gamma=True, relax=True, scf=True, dos=False, idefix=False, ncore_test=False,\
                       stm=False, eint=[-1], \
                       hse=False, \
                       enthalpy=False, \
                       coh=False, \
                       adsorption=False, plane_of_separation=0, calculate_chgdiff=False, bader=False,\
                       cohp=False, nbands=1000, \
                       surf_energy=False, count_relation=1, numb_of_layers=8, \
                       relax_far=False, \
                       results=True) -> None :

        global structures, kpts, y, calc
        structures = self.get_structures()
        calc = self.get_calculators(vasp=vasp, qe=qe)
        if bader == True:
            calc.set(laechg=True)

        calc_copy = [calc].copy()
        self.set_env(n, vasp=vasp, qe=qe)
        db = connect(self.homedir + "/" + self.database_name)

        kpts = (1, 1, 1)
        if vasp == True:
            inputs = ["INCAR", "KPOINTS", "POTCAR", "IBZKPT", "vasprun.xml", "OUTCAR", "OSZICAR", "CONTCAR", "POSCAR", "XDATCAR", "CHGCAR"]
        if slab==True:
            calc.set(idipol=3, ldipol=True) # Define this!
        else:
            calc.set(idipol=None, ldipol=None) # Define this!

        #----------MD---------#
        if md == True:
            if vasp == True:
                self.md(n, vasp=vasp, qe=qe, new=new, idefix=idefix,\
                        steps_production=steps_production, steps_final=steps_final, temperature_final=temperature_final)

                inputs.append("PCDAT")
                for i in inputs:
                    sb.run(["cp "+i+" "+i+"_md"], shell=True)
        #----------MD---------#

        calc = calc_copy.copy()
        calc = calc[0]

        #--------GAMMA-------#
        if gamma == True:
            if idefix:
                calc.set(command='mpirun -np 12 vasp_gam')
            else:
                calc.set(command='srun vasp_gam')

            x = self.verify_outcar("OUTCAR_gam", reached=True)
            if isinstance(x, int):
                self.calculation(n, vasp=vasp, qe=qe, new=new, idefix=idefix, dftu=dftu)
                for i in inputs:
                    sb.run(["cp "+i+" "+i+"_gam"], shell=True)
                if idefix:
                    sb.run(["touch ../finished_calculation"], shell=True)
        #--------GAMMA-------#

        calc = calc_copy.copy()
        calc = calc[0]

        #--------RELAX-------#
        if relax == True:
            kpts = self.set_kpts(slab=slab)

            if vasp == True:
                if idefix:
                    calc.set(command='mpirun -np 12 vasp_std')
                elif kpts == (1, 1, 1):
                    calc.set(command='srun vasp_gam')
                else:
                    calc.set(command='srun vasp_std')
                calc.set(kpts=kpts)

                x = self.verify_outcar("OUTCAR_relax", reached=True)
                if isinstance(x, int):

                    if os.path.exists("mystate_relax.json") and os.path.getsize("mystate_relax.json") > 0:
                        structures[n].calc = calc
                        energy = structures[n].get_potential_energy()
                        calc.write_json('mystate_relax.json')
                    else:
                        self.calculation(n, vasp=vasp, qe=qe, new=new, idefix=idefix, dftu=dftu)
                    if idefix:
                        sb.run(["touch ../finished_calculation"], shell=True)
                    for i in inputs:
                        sb.run(["cp "+i+" "+i+"_relax"], shell=True)
                        sb.run(["cp CHGCAR CHGCAR_relax"], shell=True)

            elif qe == True:
                try:
                    x = str(check_output('grep "!" result_expresso.out', shell=True))
                except:
                    x = 1
                if isinstance(x, int):
                    self.calculation(n, vasp=vasp, qe=qe, new=new, idefix=idefix, dftu=dftu)
        #--------RELAX-------#

        x = self.verify_outcar("OUTCAR_relax", reached=True)
        if isinstance(x, int):
            try:
                sb.run(["cp CONTCAR_relax CONTCAR"], shell=True)
            except:
                print("File CONTCAR_relax does not exist!")

        calc = calc_copy.copy()
        calc = calc[0]
        #--------SCF-------#
        if scf == True:
            kpts = self.set_kpts(slab=slab)

            if vasp == True:
                if slab==True:
                    calc.set(idipol=3, ldipol=True) # Define this!
                else:
                    calc.set(idipol=None, ldipol=None) # Define this!

                if idefix:
                    calc.set(command='mpirun -np 12 vasp_std')
                elif kpts == (1, 1, 1):
                    calc.set(command='srun vasp_gam')
                else:
                    calc.set(command='srun vasp_std')

                calc.set(kpts=kpts, nsw=0, nelm=150, lreal=False)

                x = self.verify_outcar("OUTCAR_scf", time=True)
                if isinstance(x, int):
                    self.calculation(n, vasp=vasp, qe=qe, new=new, idefix=idefix, dftu=dftu)
                    if idefix:
                        sb.run(["touch ../finished_calculation"], shell=True)
                    for i in inputs:
                        sb.run(["cp "+i+" "+i+"_scf"], shell=True)

            elif qe == True:
                try:
                    x = str(check_output('grep "!" result_expresso.out', shell=True))
                except:
                    x = 1
                if isinstance(x, int):
                    self.calculation(n, vasp=vasp, qe=qe, new=new, idefix=idefix, dftu=dftu)
        #--------SCF-------#

        calc = calc_copy.copy()
        calc = calc[0]
        #--------DOS-------#
        if dos == True:
            if idefix:
                calc.set(command='mpirun -np 12 vasp_std')
            elif kpts == (1, 1, 1):
                calc.set(command='srun vasp_gam')
            else:
                calc.set(command='srun vasp_std')
            kpts = self.set_kpts(slab=slab)
            kpts = [2*kpts[0], 2*kpts[1], 2*kpts[2]]
            if slab == True:
                kpts[2] = 1
            calc.set(icharg=11, kpts=kpts, ismear=-5, prec='Accurate', nsw=0, nedos=5000)
            self.calculation(n, vasp=vasp, qe=qe, new=new, dftu=dftu)
            if idefix:
                sb.run(["touch ../finished_calculation"], shell=True)
            if vasp == True:
                inputs.append("PROCAR")
                inputs.append("DOSCAR")
                for i in inputs:
                    sb.run(["cp "+i+" "+i+"_dos"], shell=True)

            # Getting DOS results
            if results == True:

                try:
                    elements = set(structures[n].get_chemical_symbols())

                    with open("total_dos.vp", "w") as f:
                        f.write("11\n")
                        f.write("111\n")
                    sb.run(["vaspkit <total_dos.vp"], shell=True)
                    with open("ldos_element.vp", "w") as f:
                        f.write("11\n")
                        f.write("116\n")
                    sb.run(["vaspkit <ldos_element.vp"], shell=True)

                    sb.run(["sed 's/\#//' LDOS_ELEMENTS.dat > ldos_elements.dat"], shell=True)
                    df = pd.read_csv('TDOS.dat', delim_whitespace=True)
                    df2 = pd.read_csv("ldos_elements.dat", delim_whitespace=True)
                    energy_tot = df["#Energy"]
                    dos_tot = np.array(df["TDOS"])/(120+4+12+18)

                    number_of_atoms = [120, 4, 12, 18]
                    dos_elements = [ df2[i] for i in elements ]
                    dos_elements = [ np.array(df2[i])/n for i, n in zip(elements, number_of_atoms) ]
                    #dos_al = df2["Al"]

                    fig, ax = plt.subplots()
                    ax.plot(energy_tot, dos_tot, c="k", label="TDOS")

                    if len(dos_elements) == 2:
                        colors = ["red", "blue"]
                    elif len(dos_elements) == 3:
                        colors = ["red", "blue", "green"]
                    elif len(dos_elements) == 4:
                        colors = ["red", "blue", "green", "yellow"]
                    else:
                        print("number of elements bigger than 4!")

                    for i, j, c in zip(dos_elements, elements, colors):
                        ax.plot(energy_tot, i, color=c, label="LDOS-"+str(j))
                    ax.set(xlabel=r"E - E$_F$ (eV)", ylabel=r"DOS")
                    ax.set_xlim(-10,5)
                    ax.legend()
                    ax.grid()
                    plt.axvline(x=0, c="k", linestyle='--')

                    dir_dos = self.homedir+"/dos_results/"
                    if not os.path.isdir(dir_dos):
                        sb.run(["mkdir "+dir_dos], shell=True)
                    fig.set_size_inches(13.5, 10.5)
                    fig.savefig(dir_dos+unique_identifier+"_"+str(n)+".pdf", dpi=200)
                    #fig.savefig(fname="dos_images/total_and_ldos.pdf")
                    #plt.show()
                except:
                    print("Could not get the DOS results, do it manually!")

        #--------DOS-------#

        calc = calc_copy.copy()
        calc = calc[0]
        #--------COHP-------#
        if cohp == True:
            kpts = self.set_kpts(slab=slab)
            kpts = [2*kpts[0], 2*kpts[1], 2*kpts[2]]

            if idefix:
                calc.set(command='mpirun -np 12 vasp_std')
            elif kpts == (1, 1, 1):
                calc.set(command='srun vasp_gam')
            else:
                calc.set(command='srun vasp_std')
            if slab == True:
                kpts[2] = 1
            calc.set(lwave=True, isym=0, nbands=nbands, prec='Accurate', nsw=0, kpts=kpts, ismear=-5, nelm=150, nedos=3000)

            x = self.verify_outcar("OUTCAR_cohp", time=True)
            if isinstance(x, int):
                self.calculation(n, vasp=vasp, qe=qe, new=new, dftu=dftu)

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

            if idefix:
                sb.run(["touch ../finished_calculation"], shell=True)
        #--------COHP-------#

        calc = calc_copy.copy()
        calc = calc[0]
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
                elif kpts == (1, 1, 1):
                    calc.set(command='srun vasp_gam')
                else:
                    calc.set(command='srun vasp_std')
                calc.set(kpts=kpts, lvtot=True, ldipol=True, idipol=3)

                ## GETTING THE THREE STRUCTURES
                ## IMPORTANT: Define the plane (z) that separates the substrate from the oxide!
                # For the sixton:
                slab = structures[n].copy()
                substrate = structures[n].copy()
                del substrate[[atom.index for atom in substrate if atom.position[2]>plane_of_separation]]
                oxide = structures[n].copy()
                del oxide[[atom.index for atom in oxide if atom.position[2]<plane_of_separation]]

                structures[n] = slab
                y = self.verify_outcar("OUTCAR_structure", time=True)
                if y == 1:
                    y -= 1

                if (os.path.exists("CONTCAR_relax") and os.path.getsize("CONTCAR_relax")) > 0:
                    for i in inputs:
                        sb.run(["cp "+i+"_relax "+i+"_structure"], shell=True)
                    y = 'Already calculated'
                elif (os.path.exists("CONTCAR_scf") and os.path.getsize("CONTCAR_scf")) > 0:
                    for i in inputs:
                        sb.run(["cp "+i+"_scf"+i+"_structure"], shell=True)
                    y = 'Already calculated'

                if isinstance(y, int):
                    try:
                        self.calculation(n, vasp=vasp, qe=qe, adsorption=True, new=new, dftu=dftu)
                    except:
                        pass
                    for i in inputs:
                        sb.run(["cp "+i+" "+i+"_structure"], shell=True)
                        sb.run(["mv LOCPOT LOCPOT_structure"], shell=True)

                sb.run(["rm CONTCAR"], shell=True)
                calc.set(lvtot=None, ldipol=None, idipol=None)
                structures[n] = oxide
                y = self.verify_outcar("OUTCAR_oxide", time=True)
                if y == 1:
                    y += 1
                if isinstance(y, int):
                    try:
                        self.calculation(n, vasp=vasp, qe=qe, adsorption=True, new=new, dftu=dftu)
                    except:
                        pass
                    for i in inputs:
                        sb.run(["cp "+i+" "+i+"_oxide"], shell=True)

                    if idefix:
                        sb.run(["touch ../finished_calculation"], shell=True)

                sb.run(["rm CONTCAR"], shell=True)
                structures[n] = substrate
                calc.set(lvtot=False, ldipol=None, idipol=None)

                y = self.verify_outcar("OUTCAR_substrate", time=True)
                if isinstance(y, int):
                    try:
                        self.calculation(n, vasp=vasp, qe=qe, adsorption=True, new=new, dftu=dftu)
                    except:
                        pass
                    for i in inputs:
                        sb.run(["cp "+i+" "+i+"_substrate"], shell=True)

                sb.run(["cp CONTCAR_structure CONTCAR"], shell=True)

                result_substrate = read("OUTCAR_substrate", format="vasp-out")
                result_oxide = read("OUTCAR_oxide", format="vasp-out")
                result_structure = read("OUTCAR_structure", format="vasp-out")
                e_substrate = result_substrate.get_potential_energy()
                e_oxide = result_oxide.get_potential_energy()
                e_structure = result_structure.get_potential_energy()

                e_adsorption = (e_structure - e_oxide - e_substrate)/(cell[0]*cell[1]*m.sin(m.radians(cell[5])))

                self.database(result_structure, str(namefile))
                try:
                    db.update(id=ID, energy_substrate=e_substrate, energy_oxide=e_oxide, adhesive_energy=e_adsorption)
                except:
                    print("Could not put in the database")

                with open(self.homedir+"/adh_energy_"+str(n)+".txt", "w") as f:
                    f.write('Adhesive energy: '+str(e_adsorption)+' eV/A2\n')
                    f.write('Total energy '+str(e_structure)+' eV\n')

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

        calc = calc_copy.copy()
        calc = calc[0]
        #--------FORMATION ENTHALPY-------#
        # Given a specific structure, calculate its formation energy.

        if enthalpy == True:

            bulks = sorted(glob.glob(self.homedir+"/bulk_for_enthalpies/"+str(namefile)+"/*.cif"))
            bulks_final = [ read(i, format="cif") for i in bulks ]

            if vasp == True:

                kpts = self.set_kpts(slab=slab)
                calc.set(kpts=kpts, nsw=500, isif=3, ediffg=-0.02)

                if idefix:
                    calc.set(command='mpirun -np 12 vasp_std')
                elif kpts == (1, 1, 1):
                    calc.set(command='srun vasp_gam')
                else:
                    calc.set(command='srun vasp_std')

                energies_per_bulk, cohesive_energies = defaultdict(dict), defaultdict(dict)

                for xcpot in ["pbe", "optpbe-vdw"]:

                    sb.run(["mkdir "+str(xcpot)],shell=True)
                    os.chdir(str(xcpot))
                    calc.set(xc=xcpot)

                    formula_struc = structures[n].get_chemical_formula(mode='hill')
                    print(formula_struc)
                    x = self.verify_outcar("OUTCAR_enthalpy", reached=True)
                    if isinstance(x, int):
                        energ = self.calculation(n, vasp=vasp, qe=qe, new=new, idefix=idefix, dftu=dftu)
                        db.write(structures[n], type_calc="enthalpy")
                        if idefix:
                            sb.run(["touch ../finished_calculation"], shell=True)
                        for i in inputs:
                            sb.run(["cp "+i+" "+i+"_enthalpy"], shell=True)
                    else:
                        energ = read("OUTCAR_enthalpy", format="vasp-out").get_potential_energy()
                    #energies_per_bulk.update = {formula_struc: {xcpot: energ}}
                    energies_per_bulk[formula_struc][xcpot] = energ

                    # To the final formula:
                    # To do:
                    e_tot = read("OUTCAR_enthalpy", format="vasp-out").get_potential_energy()
                    tot_number_of_atoms = len(structures[n].get_chemical_symbols())
                    elements_in_structure = []
                    [elements_in_structure.append(x) for x in structures[n].get_chemical_symbols()if x not in elements_in_structure]
                    #elements_in_structure = list(set(s.get_chemical_symbols()))
                    n_atoms_per_element = {}

                    for el in elements_in_structure:
                        number = structures[n].get_chemical_symbols().count(el)
                        n_atoms_per_element[el] = number

                    # Number of each atom in the structure:

                    ####################

                    for bulk in bulks_final:
                        precision=6
                        k = (max(1,int(precision*2*m.pi/bulk.get_cell()[0][0])), \
                            max(1,int(precision*2*m.pi/bulk.get_cell()[1][1])), \
                            max(1,int(precision*2*m.pi/bulk.get_cell()[2][2])))
                        calc.set(isif=3, ibrion=2, ediffg=-0.02, kpts=k)
                        bulk.calc = calc
                        n_of_atoms_bulk = len(bulk.get_chemical_symbols())
                        element = bulk.get_chemical_symbols()[0]

                        sb.run(["mkdir "+str(element)],shell=True)
                        os.chdir(str(element))

                        x = self.verify_outcar("OUTCAR_"+str(bulk.get_chemical_formula(mode='hill')), reached=True)
                        if isinstance(x, int):
                            e_bulk = bulk.get_potential_energy()
                            for i in inputs:
                                sb.run(["cp "+i+" "+i+"_"+str(bulk.get_chemical_formula(mode='hill'))], shell=True)
                            db.write(bulk, type_calc="bulks_enthalpy")
                        else:
                            e_bulk = read("OUTCAR_"+str(bulk.get_chemical_formula(mode='hill')), format="vasp-out").get_potential_energy()
                        energies_per_bulk[element][xcpot] = e_bulk

                        if coh == True:

                            # Cohesion energy
                            kpts = (1, 1, 1)
                            at = bulk.get_chemical_symbols()[0]
                            atomo = Atoms(at, positions=[(0, 0, 0)], cell=[13, 14, 15], pbc=True)
                            calc.set(isif=None, ibrion=2, nsw=0, nelm=1000, kpts=kpts)
                            atomo.calc = calc
                            x = self.verify_outcar("OUTCAR_at_"+str(atomo.get_chemical_formula(mode='hill')), time=True)
                            if isinstance(x, int):
                                e_atomo = atomo.get_potential_energy()
                                for i in inputs:
                                    sb.run(["cp "+i+" "+i+"_at_"+str(atomo.get_chemical_formula(mode='hill'))], shell=True)
                                db.write(atomo, type_calc="atomo_cohesion")
                            else:
                                e_atomo = read("OUTCAR_at_"+str(atomo.get_chemical_formula(mode='hill'))).get_potential_energy()

                            cohesive = (e_bulk - n_of_atoms_bulk * e_atomo)/n_of_atoms_bulk
                            cohesive_energies[element][xcpot] = cohesive
                            with open(self.homedir+"/cohesive_energy_"+str(bulk.get_chemical_formula(mode='hill')), "w") as f:
                                f.write("Cohesive energy of bulk "+str(atomo.get_chemical_formula(mode='hill')+" is: "+str(cohesive))+"\n")
                                f.write("Number of atoms bulk: "+str(n_of_atoms_bulk)+"\n")
                                f.write("E_bulk: "+str(e_bulk)+"\n")
                                f.write("E_atomo: "+str(e_atomo)+"\n")

                            sb.run(["rm CHG* vasprun*"], shell=True)

                        os.chdir("..")

                    os.chdir("..")
                sb.run(["cp CONTCAR_enthalpy CONTCAR"], shell=True)

                with open(self.homedir+'/results.json', 'w') as convert_file:
                    convert_file.write(json.dumps(energies_per_bulk))
                with open(self.homedir+'/cohesive_energies.json', 'w') as convert_file:
                    convert_file.write(json.dumps(cohesive_energies))
                # Calculating enthalpy of formation
                #hf = ( e_tot - (  ) )/tot_number_of_atoms

        #--------FORMATION ENTHALPY-------#

        calc = calc_copy.copy()
        calc = calc[0]
        #--------STM-------#
        if stm == True:

            kpts = self.set_kpts(slab=slab)
            if idefix:
                calc.set(command='mpirun -np 12 vasp_std')
            elif kpts == (1, 1, 1):
                calc.set(command='srun vasp_gam')
            else:
                calc.set(command='srun vasp_std')

            x = self.verify_outcar("OUTCAR_scf_stm", time=True)
            if isinstance(x, int):
                calc.set(nsw=0, nelm=120, kpts=kpts, prec='Accurate', ismear=1, lwave=True, algo="Fast", istart=0, ibrion=2)
                self.calculation(n, vasp=vasp, qe=qe, new=new, dftu=dftu)
                if vasp == True:
                    for i in inputs:
                        sb.run(["cp "+i+" "+i+"_scf_stm"], shell=True)

            x = self.verify_outcar("OUTCAR_par_stm", time=True)
            if isinstance(x, int):
                for voltage in [-1, -0.5, 0.5, 1]:
                    calc.set(nsw=1, kpts=kpts, prec='Accurate', lsepk=False, lsepb=False, lpard=True, nbmod=-3, eint=[voltage], ismear=2, lplane=True, nsim=4, npar=8, lscalu=False, istart=2, ibrion=1)
                    calc.initialize(structures[n])
                    calc.write_input(structures[n])
                    sb.run(["cp CONTCAR POSCAR"], shell=True)
                    sb.run(["srun vasp_std"], shell=True)

                    sb.run(["cp PARCHG PARCHG_"+str(voltage)], shell=True)
                    sb.run(["cp OUTCAR OUTCAR_par_stm_"+str(voltage)], shell=True)
                    #self.calculation(n, vasp=vasp, qe=qe, new=False)

                if vasp == True:
                    for i in inputs:
                        sb.run(["cp "+i+" "+i+"_par_stm"], shell=True)

            if idefix:
                sb.run(["touch ../finished_calculation"], shell=True)
        #--------STM-------#

        calc = calc_copy.copy()
        calc = calc[0]
        #--------HSE-------#
        if hse == True:
            if idefix:
                calc.set(command='mpirun -np 12 vasp_std')
            else:
                calc.set(command='srun vasp_std')
            kpts = self.set_kpts(slab=slab)
            try:
                x = str(check_output('grep "Total CPU time used" OUTCAR_scf', shell=True))
            except:
                x = 1
            if isinstance(x, int):
                calc.set(nsw=0, nelm=120, kpts=kpts, prec='Accurate', ismear=1, lwave=True, algo="Fast", istart=0, ibrion=2)
                self.calculation(n, vasp=vasp, qe=qe, new=new, dftu=dftu)
                if vasp == True:
                    for i in inputs:
                        sb.run(["cp "+i+" "+i+"_scf_stm"], shell=True)

            try:
                x = str(check_output('grep "Total CPU time used" OUTCAR_hse', shell=True))
            except:
                x = 1
            if isinstance(x, int):
                calc.set(istart=2, xc='hse03')
                self.calculation(n, vasp=vasp, qe=qe, new=new, dftu=dftu)

                if vasp == True:
                    for i in inputs:
                        sb.run(["cp "+i+" "+i+"_par_stm"], shell=True)

            if idefix:
                sb.run(["touch ../finished_calculation"], shell=True)
        #--------HSE-------#

        calc = calc_copy.copy()
        calc = calc[0]
        #--------SURFACE ENERGY-------#
        if surf_energy == True:

            import pint
            ureg = pint.UnitRegistry()

            kpts = self.set_kpts(slab=slab)

            # For the bulk ------------------------------------

            if idefix:
                calc.set(command='mpirun -np 12 vasp_std')
            elif kpts == (1, 1, 1):
                calc.set(command='srun vasp_gam')
            else:
                calc.set(command='srun vasp_std')

            calc.set(isif=3, ibrion=2, ediffg=-0.02, kpts=kpts, nsw=200)

            x = self.verify_outcar("OUTCAR_relax", reached=True)
            if isinstance(x, int):
                self.calculation(n, vasp=vasp, qe=qe, new=new, idefix=idefix, dftu=dftu)
                sb.run(["cp OUTCAR OUTCAR_relax; cp CONTCAR CONTCAR_relax"], shell=True)
            energy_bulk = read("OUTCAR_relax", format="vasp-out").get_potential_energy()

            structures[n] = read("CONTCAR_relax", format="vasp")
            c = FixAtoms(indices=[atom.index for atom in structures[n] if atom.position[2] <= 0.3])
            structures[n].set_constraint(c)

            # For the surfaces ----------------------------------
            area = structures[n].get_cell_lengths_and_angles()[0]*structures[n].get_cell_lengths_and_angles()[1]

            energies = []
            # Need to find the value of "count_relation" according to the number of layers in the bulk!!
            layers = range(2, numb_of_layers)
            for sup_cells in layers:
                s = structures[n].copy()
                s = make_supercell(s, [[1, 0, 0], [0, 1, 0], [0, 0, sup_cells]])
                s = cut(s, nlayers=sup_cells+count_relation)
                s.set_cell([s.get_cell()[0][0], s.get_cell()[1][1], s.get_cell()[2][2]+15])
                n_in = s.get_chemical_symbols().count('In')
                n_pd = s.get_chemical_symbols().count('Pd')

                count_relation += 1

                precision=6
                kpts = (max(1,int(precision*2*m.pi/s.get_cell()[0][0])), \
                    max(1,int(precision*2*m.pi/s.get_cell()[1][1])), \
                    1)
                calc.set(isif=2, nsw=200, nelm=60, kpts=kpts, prec='Accurate', ismear=1, lwave=False, algo="Fast", istart=0, ibrion=2)
                s.calc = calc

                if idefix:
                    calc.set(command='mpirun -np 12 vasp_std')
                elif kpts == (1, 1, 1):
                    calc.set(command='srun vasp_gam')
                else:
                    calc.set(command='srun vasp_std')

                x = self.verify_outcar("OUTCAR_surf"+str(sup_cells), time=True)
                if isinstance(x, int):
                    energy = s.get_potential_energy() * ureg.eV
                    sb.run(["cp OUTCAR OUTCAR_surf"+str(sup_cells)], shell=True)
                else:
                    energy = read("OUTCAR_surf"+str(sup_cells), format="vasp-out").get_potential_energy() * ureg.eV
                #energies.append(energy.magnitude/n_atoms)
                energies.append(energy.magnitude)

                sb.run(["rm CHG* vasprun*"], shell=True)

            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt
            from sklearn.linear_model import LinearRegression

            energies.sort(reverse=True)
            energies = np.array(energies)
            nlayers = np.array(layers).reshape((-1, 1))
            print(nlayers)
            print(energies)
            model = LinearRegression()
            model.fit(nlayers, energies)
            I = model.intercept_
            print('intercept:', model.intercept_)
            print('slope:', model.coef_)
            gamma_eV = I/(2*area) * ureg["eV/angstrom**2"]
            gamma = gamma_eV.to("J/meter**2")
            print(gamma_eV)
            print(gamma)

            fig, ax = plt.subplots()
            ax.scatter(nlayers, energies)
            plt.text(0.7,0.95,r'$\gamma$: '+str(round(gamma,4))+r' J/m$^2$',horizontalalignment='left', verticalalignment='baseline', transform = ax.transAxes)
            plt.text(0.7,0.9,r'$\gamma$: '+str(round(gamma_eV,4))+r' eV/\AA$^2$',horizontalalignment='left', verticalalignment='baseline', transform = ax.transAxes)
            plt.text(0.7,0.85,'intercept: '+str(round(model.intercept_, 4)),horizontalalignment='left', verticalalignment='baseline', transform = ax.transAxes)
            plt.text(0.7,0.8,'slope: '+str(round(model.coef_[0], 4)),horizontalalignment='left', verticalalignment='baseline', transform = ax.transAxes)
            ax.grid()
            ax.set(xlabel="Number of layers", ylabel=r"E$_{\rm tot}$ (eV)")
            ax.plot(nlayers, model.predict(nlayers))
            fig.savefig(self.homedir+"/surf_energy.pdf")
            #sb.run(["evince "+self.homedir+"/surf_energy.pdf"], shell=True)
            #plt.show()

        #--------SURFACE ENERGY-------#

        calc = calc_copy.copy()
        calc = calc[0]
        #--------RELAX FAR STRUCTURES-------#
        if relax_far == True:
            kpts = self.set_kpts(slab=slab)

            if vasp == True:
                if slab==True:
                    calc.set(idipol=3, ldipol=True) # Define this!
                else:
                    calc.set(idipol=None, ldipol=None) # Define this!

                if idefix:
                    calc.set(command='mpirun -np 12 vasp_std')
                elif kpts == (1, 1, 1):
                    calc.set(command='srun vasp_gam')
                else:
                    calc.set(command='srun vasp_std')
                calc.set(kpts=kpts)

                x = self.verify_outcar("OUTCAR_relax1", reached=True)
                if isinstance(x, int):
                    calc.set(isif=7, nsw=20)
                    self.calculation(n, vasp=vasp, qe=qe, new=new, idefix=idefix, dftu=dftu)
                    for i in inputs:
                        sb.run(["cp OUTCAR OUTCAR_relax1"], shell=True)

                x = self.verify_outcar("OUTCAR_relax2", reached=True)
                if isinstance(x, int):
                    calc.set(isif=2, nsw=20)
                    self.calculation(n, vasp=vasp, qe=qe, new=new, idefix=idefix, dftu=dftu)
                    for i in inputs:
                        sb.run(["cp OUTCAR OUTCAR_relax2"], shell=True)

                x = self.verify_outcar("OUTCAR_relax3", reached=True)
                if isinstance(x, int):
                    calc.set(isif=7, nsw=20)
                    self.calculation(n, vasp=vasp, qe=qe, new=new, idefix=idefix, dftu=dftu)
                    for i in inputs:
                        sb.run(["cp OUTCAR OUTCAR_relax3"], shell=True)

                x = self.verify_outcar("OUTCAR_relax4", reached=True)
                if isinstance(x, int):
                    calc.set(isif=2, nsw=20)
                    self.calculation(n, vasp=vasp, qe=qe, new=new, idefix=idefix, dftu=dftu)
                    for i in inputs:
                        sb.run(["cp OUTCAR OUTCAR_relax4"], shell=True)
                        sb.run(["cp OUTCAR_relax"], shell=True)

        #--------RELAX FAR STRUCTURE-------#

        x = self.verify_outcar("OUTCAR_relax", reached=True)
        if isinstance(x, int):
            try:
                sb.run(["cp CONTCAR_relax CONTCAR"], shell=True)
            except:
                print("File CONTCAR_relax does not exist!")

        calc = calc_copy.copy()
        calc = calc[0]
        #--------NCORE TEST-------#
        if ncore_test == True:
            if idefix:
                calc.set(command='mpirun -np 12 vasp_std')
            else:
                calc.set(command='srun vasp_gam')
            ncore_list = [2, 4, 8, 16]
            kpts = (1, 1, 1)
            times = {}
            for i in ncore_list:
                sb.run(["mkdir ncore"+str(i)], shell=True)
                os.chdir("ncore"+str(i))
                calc.set(kpts=kpts, nsw=3, ncore=i)
                time_i = time.time()
                try:
                    self.calculation(n, vasp=vasp, qe=qe, new=new)
                except:
                    print("There was a problem with ncore_"+str(i))
                time_f = time.time()
                times["ncore"+str(i)] = str(round(time_f-time_i,2))+"s"
                os.chdir("..")

            with open(self.homedir + "/" + namefile + "_" + str(n) + "_" +"time.txt", "w") as r:
                r.write(json.dumps(times)+'\n')
        #--------NCORE TEST-------#

        #--------FINAL VERIFICATION--------#
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

        # Saving everything
        sb.run(["cp *_* "+savedir], shell=True)

        if cohp == True:
            sb.run(["cp *.lobster lobsterout "+savedir], shell=True)

        if bader == True and adsorption == False:
            sb.run(["chgsum.pl AECCAR0 AECCAR2; bader CHGCAR -ref CHGCAR_sum; cp ACF.dat BCF.dat "+savedir], shell=True)

        if enthalpy == True:
            sb.run(['rsync -nrv --include="*/" --include="OUTCAR*" --include="CONTCAR*" --exclude="*" . ' + savedir], shell=True)

###########################  JOB SCRIPTS #################################

    def launch_jeanzay(self, n, p=40):

        p = 40
        time = '20:00:00'
        structures = self.get_structures()
        processors = len(structures[n].get_chemical_symbols())
        if processors > p:
            while int(processors) % p != 0:
                processors -= 1
        else:
            while int(processors) % p != 0:
                processors += 1
        #processors = processors*3
        name = structures[n].get_chemical_formula(mode='hill')

        #processors=2*p
        myscript_content='''#!/bin/sh
#SBATCH --nodes='''+str(int(processors/p))+'''               # 1 node reserved
#SBATCH --ntasks-per-node=40    # 40 MPI tasks
##SBATCH --mem=0       # 1 OpenMP thread
#SBATCH --hint=nomultithread    # Disable hyperthreading
#SBATCH -J '''+name+'''
#SBATCH -o '''+name+'''.out
#SBATCH -e '''+name+'''.out
#SBATCH --time='''+time+'''
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

        p = 24
        structures = self.get_structures()
        processors = len(structures[n].get_chemical_symbols())
        if processors > p:
            while int(processors) % p != 0:
                processors -= 1
        else:
            while int(processors) % p != 0:
                processors += 1
        #processors = processors*3
        name = structures[n].get_chemical_formula(mode='hill')

        #processors=2*p
        myscript_content='''#!/bin/sh
#SBATCH -J '''+name+'''
#SBATCH --nodes='''+str(int(processors/p))+'''
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

        p = 32
        structures = self.get_structures()
        processors = len(structures[n].get_chemical_symbols())
        if processors > p:
            while int(processors) % p != 0:
                processors -= 1
        else:
            while int(processors) % p != 0:
                processors += 1
        name = structures[n].get_chemical_formula(mode='hill')

        processors = 2*p
        myscript_content='''#!/bin/sh
#SBATCH --nodes='''+str(int(processors/p))+'''
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

        p = 40
        structures = self.get_structures()
        processors = len(structures[n].get_chemical_symbols())
        if processors > p:
            while int(processors) % p != 0:
                processors += 1
        else:
            while int(processors) % p != 0:
                processors += 1
        name = structures[n].get_chemical_formula(mode='hill')

        processors = 2*p
        myscript_content='''#!/bin/bash
#MSUB -q skylake
#MSUB -T 28800
#MSUB -N '''+str(int(processors/p))+'''
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
    structures()
