#!/usr/bin/env python

# Author: Thiago Trevizam Dorini

from __future__ import division, print_function, unicode_literals, \
    absolute_import
import os
import numpy as np
import glob
import subprocess as sb
import itertools
from launch.vulcan_helper import *
from ase.visualize import view
from ase.build import surface, bulk, make_supercell, cut, stack, add_vacuum, molecule, rotate, add_adsorbate, fcc111
from ase import Atoms
from ase.io import read, write
from ase.calculators.vasp import Vasp
from ase.calculators.espresso import Espresso
from ase.db import connect
from ase.constraints import FixAtoms
from ase.io.vasp import read_vasp_xdatcar
from pymatgen.io.ase import AseAtomsAdaptor as ase
from pymatgen.analysis.adsorption import AdsorbateSiteFinder as asf
from pymatgen.core.structure import Structure, Molecule
import subprocess as sb
from mpinterfaces.calibrate import CalibrateSlab
from mpinterfaces.interface import Interface
from mpinterfaces.transformations import *
from mpinterfaces.utils import *
from pymatgen.core import Structure
from pymatgen.io.ase import AseAtomsAdaptor as ase
from pymatgen.core.surface import SlabGenerator, ReconstructionGenerator, generate_all_slabs
from ase.build import make_supercell, surface
from pymatgen.io.ase import AseAtomsAdaptor as ase

class MyStructures_Approx():

    def __init__(self, structure_file=None):

        self.structure_file = structure_file

    def get_poscars_glob(self):

        files = sorted(glob.glob("POSCAR*"))
        s = [ read(f, format="vasp") for f in files ]
        return s

class MyStructures_PdIn():

    def __init__(self, structure_file, hkl_oxide=[0, 0, 1], molecule='HCOOH', thickness=5):

        self.structure_file = structure_file
        self.hkl_oxide = hkl_oxide
        self.molecule = molecule
        self.thickness = thickness

    def get_or_view_all(self, view_struc=True):

        here = os.getcwd()
        path_stack = "matched_poscars/"+str(self.hkl_oxide[0])+str(self.hkl_oxide[1])+str(self.hkl_oxide[2])+"/stack_"+str(self.thickness)+"A"
        os.chdir(path_stack)

        import glob
        files = sorted(glob.glob("stacked*"))
        struc_ase = [ read(i) for i in files ]
        if view_struc:
            view(struc_ase)

        os.chdir(here)

    def mpinterf_in2o3(self, oxide_file, hkl_oxide=[0, 0, 1], hkl_subst=[0, 0, 1], max_area=400, max_mismatch=0.05, max_angle=1, r1r2_tol=0.01, poscars=True, symmetrize=True):

        # import mpinterfaces
        s = read(self.structure_file)
        oxide = read(oxide_file)
        s = ase.get_structure(s)
        oxide = ase.get_structure(oxide)

        def slabs(atoms, hkl, thick=self.thickness, vac=5, oxide=True):

            sa_sub = SpacegroupAnalyzer(atoms)
            s = sa_sub.get_conventional_standard_structure()
            # slab_f = Interface(s, hkl=hkl, min_thick=thick, min_vac=vac, to_unit_cell=True, primitive=False, from_ase=False, force_normalize=True)
            if (hkl_oxide == [1,1,1]) and (oxide==True):
                slab_f = SlabGenerator(s, hkl, thick, vac, lll_reduce=True, in_unit_planes=False, center_slab=False)
            else:
                slab_f = SlabGenerator(s, hkl, thick, vac, lll_reduce=False, in_unit_planes=False)

            if oxide == False:
                slabgen = slab_f.get_slabs(symmetrize=True, repair=True)
            else:
                slabgen = slab_f.get_slabs(bonds={("In3+", "O2-"): 3}, symmetrize=symmetrize, repair=False)

            return slabgen

        sub = slabs(s, hkl_subst, thick=12, vac=5, oxide=False)

        if symmetrize:
            thickness = 15
        else:
            thickness = self.thickness # We can get thiner oxides here!

        if hkl_oxide == [1, 1, 0]:
            ox = slabs(oxide, hkl_oxide, thick=thickness, vac=5)
        elif hkl_oxide == [0, 0, 1]:
            ox = slabs(oxide, hkl_oxide, thick=thickness, vac=5)
        elif hkl_oxide == [1, 1, 1]:
            ox = slabs(oxide, hkl_oxide, thick=thickness, vac=10)

        # test = ase.get_atoms(ox[0])
        # view(test)
        # sys.exit()
        # Visualize structures
        # list_ox, list_sub = [], []
        # for i, j in zip(sub, ox):
            # j.get_orthogonal_c_slab()
            # sub2 = ase.get_atoms(i)
            # ox2 = ase.get_atoms(j)
            # list_ox.append(ox2)
            # list_sub.append(sub2)
        # view(list_ox)
        # view(list_sub)

        sub = sub[0] # Substrate terminated on Pd.
        sub_ase = ase.get_atoms(sub)
        # view(sub_ase)

        # For the database
        path_stack = "matched_poscars/"+str(hkl_oxide[0])+str(hkl_oxide[1])+str(hkl_oxide[2])+"/stack_"+str(thickness)+"A"
        sb.run(["rm "+path_stack+"/pdin_"+str(thickness)+"A.json"], shell=True)
        ### Main function
        good_structures = []
        for n, oxide_structure in enumerate(ox):

            substrate_slab_aligned, mat2d_slab_aligned = get_aligned_lattices(sub, oxide_structure, max_area=max_area, max_mismatch=max_mismatch, max_angle_diff=max_angle, r1r2_tol=r1r2_tol)
            # Checking whether it maintains the In2O3 stoichiometry:
            mat2d_ase = ase.get_atoms(mat2d_slab_aligned)
            ind = mat2d_ase.get_chemical_symbols().count("In")
            o = mat2d_ase.get_chemical_symbols().count("O")

            if (ind/o == 2/3) or (hkl_oxide == [0, 0, 1]) or (hkl_oxide == [1, 0, 0]) or (hkl_oxide == [0, 1, 0]):

                print("Good stoichiometry! (Or working on the oxide surface (001))")
                sb.run(["mkdir -p matched_poscars/"+str(hkl_oxide[0])+str(hkl_oxide[1])+str(hkl_oxide[2])], shell=True)
                os.chdir("matched_poscars/"+str(hkl_oxide[0])+str(hkl_oxide[1])+str(hkl_oxide[2]))
                substrate_slab_aligned.to(fmt='"poscar', filename='POSCAR'+str(n)+'_substrate_aligned.vasp')
                mat2d_slab_aligned.to(fmt='poscar', filename='POSCAR'+str(n)+'_mat2d_aligned.vasp')
                os.chdir("../..")

                # Stacking structures

                path = "matched_poscars/"+str(hkl_oxide[0])+str(hkl_oxide[1])+str(hkl_oxide[2])
                subst = read(path+'/POSCAR'+str(n)+'_substrate_aligned.vasp')
                ox = read(path+'/POSCAR'+str(n)+'_mat2d_aligned.vasp')
                k = stack(subst, ox, maxstrain=0.5, distance=2.8)
                k.set_cell(subst.get_cell())
                max_position = max([ i.position[2] for i in k ])
                k.cell[2][2] = max_position+10
                k.wrap()
                # view(k)
                # sys.exit()

                sb.run(["mkdir "+path_stack], shell=True)
                write(path_stack+"/stacked"+str(n)+".vasp", k)
                db = connect(path_stack+"/pdin_"+str(thickness)+"A.json")
                db.write(k, orientation=int(str(hkl_oxide[0])+str(hkl_oxide[1])+str(hkl_oxide[2])), thickness_oxide=thickness)
                good_structures.append(k)

            else:
                print("Pass")
        # view(good_structures)

        # orientations = [[1, 1, 0], [0, 0, 1]]
        # for ori in orientations:
            # x = MyStructures_PdIn("pdin.cif").mpinterf_in2o3("in2o3.cif", hkl_oxide=ori, hkl_subst=[0, 0, 1], max_area=1000)
        ####
        # sb.run(["vis -f matched_poscars/POSCAR1_substrate_aligned.vasp; vis -f matched_poscars/POSCAR1_mat2d_aligned.vasp"], shell=True)

        # Creating a Slab from the structure:
        #slab1 = ase.get_structure(mat2d_slab_aligned)

        # slab1 = SlabGenerator(mat2d_slab_aligned, hkl, min_slab_size=0.01, min_vacuum_size=0.001).get_slabs(symmetrize=False, repair=False)
        # for n, i in enumerate(slab1):
            # i.to(fmt='poscar', filename='POSCAR_test')
            # #if n < 3:
            # #    sb.run(["vis -f POSCAR_test"], shell=True)
            # #substrate_slab_aligned.to(fmt='poscar', filename='POSCAR_substrate_aligned.vasp')

        # print(slab1)
        # return substrate_slab_aligned, mat2d_slab_aligned

        # return s

    def create_supercell(self):

        s = read(self.structure_file)
        struc = []
        for i in [2, 3, 4]:
            st = s.copy()
            st = supercell(st, matrix=[1,i,1])
            struc.append(st)

        for j in [2]:
            st = s.copy()
            st = supercell(st, matrix=[j,j,1])
            struc.append(st)

        st = s.copy()
        st = supercell(st, matrix=[2,3,1])
        struc.append(st)

        struc_pym = [ ase.get_structure(i) for i in struc ] # Get pymatgen structures as well

        return struc, struc_pym

    def get_adsorbed_structures(self, distance=3):

        list_of_structures_ase, list_of_structures_pymat = self.create_supercell()
        adsorbate = molecule(self.molecule)
        for atom in adsorbate:
            if atom.symbol == "C":
                c = atom.position
                print(type(c[0]))

        def delete_H(ads):

            for atom in ads:
                if atom.symbol == "H":
                    del ads[atom.index]
                    break

        ads_struc_final = []
        for s in list_of_structures_pymat:

            for p in [80]:
                for t in [0, 90]:

                    afind = asf(s)
                    adsorbate = molecule(self.molecule)

                    delete_H(adsorbate)

                    adsorbate.rotate(p, 'z')
                    adsorbate.rotate(t, 'y')
#                     adsorbate.euler_rotate(phi=p, theta=t, psi=t, center='COU')


    #                 adsorbate = read(self.molecule, format="vasp")
    #                 adsorbate.set_cell(None)
                    adsorbate = ase.get_molecule(adsorbate)

                    ads_structs = afind.generate_adsorption_structures(molecule=adsorbate, find_args={"distance": distance})
                    ads_structs_ase_hcooh = [ ase.get_atoms(i) for i in ads_structs ]
                    c = constraint(ads_structs_ase_hcooh[0], less_position=17) # Constraint
                    ads_structs_ase_hcooh[0].set_constraint(c) # Add constraint
                    ads_struc_final.append(ads_structs_ase_hcooh[0]) # Append to final list

        return ads_struc_final

    def all_adsorbed_structures(self, supercell=[1, 3, 1], distance=2):

        s = Structure.from_file(self.structure_file)
        s.make_supercell(supercell)

        afind = asf(s) # Find adsorption sites

        final_struc = []

        if self.molecule == "H":

            adsorbate = molecule(self.molecule)
            adsorbate = ase.get_molecule(adsorbate)
            ads_structs = afind.generate_adsorption_structures(molecule=adsorbate, find_args={"distance": distance}) # Generate structures with selected adsorption sites
            ads_structs_ase_h = [ ase.get_atoms(i) for i in ads_structs ]

            c = constraint(ads_structs_ase_h[0], less_position=17) # Constraint
            for s in ads_structs_ase_h:
                s.set_constraint(c)

            final_struc = final_struc + ads_structs_ase_h

        elif self.molecule == "HCOOH":

            for angles_z in [30, 60, 90]:
                for angles_x in [0, 45, 90]:

                    adsorbate = molecule(self.molecule)

                    for atom in adsorbate:
                        if atom.symbol == "H":
                            del adsorbate[atom.index]
                            break
#                     adsorbate.euler_rotate()
#                     adsorbate.rotate(angles_z, 'z')
#                     adsorbate.rotate(angles_x, 'y')

                    adsorbate = ase.get_molecule(adsorbate)
                    ads_structs = afind.generate_adsorption_structures(molecule=adsorbate) # Generate structures with selected adsorption sites
                    ads_structs_ase_hcooh = [ ase.get_atoms(i) for i in ads_structs ]
                    c = constraint(ads_structs_ase_hcooh[0], less_position=17) # Constraint
                    for s in ads_structs_ase_hcooh:
                        s.set_constraint(c)

                    final_struc = final_struc + ads_structs_ase_hcooh

        return final_struc

def run_pdin_in2o3_interfaces():

    orientations = [[1, 1, 0], [0, 0, 1]]
    # orientations = [[0, 0, 1], [1, 1, 0]]

    ########### For visualizing only:
    # for ori in orientations:

        # for th in [1, 15]:
            # x = MyStructures_PdIn("pdin.cif", hkl_oxide=ori, thickness=th).get_or_view_all()
    # sys.exit()
    #################################

    for ori in orientations:
        # Thin oxides
        for th in [1]:
        # for th in [8, 15, 20]:
            x = MyStructures_PdIn("pdin.cif", thickness=th).mpinterf_in2o3("in2o3.cif", hkl_oxide=ori, hkl_subst=[0, 0, 1], max_area=1000, symmetrize=False)

        # Thick oxides
        x = MyStructures_PdIn("pdin.cif").mpinterf_in2o3("in2o3.cif", hkl_oxide=ori, hkl_subst=[0, 0, 1], max_area=1000, symmetrize=True)

        # IMPORTANT:
        # If we use "symmetrized = False", we can reduce the oxide thickness, otherwise we need to put at least 15A.

if __name__ == "__main__":

    run_pdin_in2o3_interfaces()
