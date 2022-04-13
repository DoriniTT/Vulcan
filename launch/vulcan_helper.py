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
import pandas as pd
import matplotlib.pyplot as plt
from ase.spacegroup import get_spacegroup
from ase.visualize import view
from ase.build import surface, bulk, make_supercell, cut, stack, add_vacuum, molecule, rotate, add_adsorbate, fcc111
from ase import Atoms
from ase.io import read, write
from ase.io.vasp import write_vasp_xdatcar
from ase.io.animation import write_gif
from ase.calculators.vasp import Vasp
from ase.calculators.espresso import Espresso
from ase.db import connect
from ase.constraints import FixAtoms
from ase.io.vasp import read_vasp_xdatcar

def get_structs_database(database):

    db = connect(database)
    c = db.count()
    strucs = [ db.get_atoms(id=i) for i in range(1, c+1) ]
    return strucs

def stack_structures(substrate, support, distance=1.4, vacuum=15, fix="substrate") -> object:

    for atom1, atom2 in zip(substrate, support):
        atom1.tag = 0
        atom2.tag = 1

    a_support, b_support = support.get_cell()[0][0], support.get_cell()[1][1]
    a_substrate, b_substrate = substrate.get_cell()[0][0], substrate.get_cell()[1][1]
    area_support = a_support * b_support
    area_substrate = a_substrate * b_substrate

    adiff = 100*(a_substrate - a_support)/a_substrate
    bdiff = 100*(b_substrate - b_support)/b_substrate
    area_diff = 100*(area_substrate - area_support)/area_substrate

    print("Lattice par 'a' diff (in %): ", adiff)
    print("Lattice par 'b' diff (in %): ", bdiff)
    print("Area diff (in %): ", area_diff)

    if fix == "substrate":
        support.set_cell([substrate.get_cell()[0][0], substrate.get_cell()[0][0], support.get_cell()[2][2]])
    elif fix == "support":
        substrate.set_cell([support.get_cell()[0][0], support.get_cell()[0][0], substrate.get_cell()[2][2]])

    zmin = support.positions[:, 2].min()
    zmax = substrate.positions[:, 2].max()
    support.positions += (0, 0, zmax - zmin + distance)
    c = support.positions[:, 2].max()
    interface = substrate + support
    interface.set_cell([interface.get_cell()[0][0], interface.get_cell()[1][1], c])
    add_vacuum(interface, vacuum)
    return interface

def constraint(atoms, less_position=None, symbol=None):

    if less_position != None:
        c = FixAtoms(indices=[atom.index for atom in atoms if atom.position[2] <= less_position])
    elif symbol != None:
        c = FixAtoms(indices=[atom.index for atom in atoms if atom.symbol == symbol])
    #atoms.set_constraint(c)
    return c

def change_basis(rotate_2d=True, at1=None, at2=None, at3=None) -> list:
    # at1 is a list with the coordinates of the atom at the new origin.
    # at2 is the atom at the "a" vector.
    # at3 is the atom at the "b" vector.
    # Ex: [[0, 1, 0], [1, 0, 0]]

    if rotate_2d == True:

        # Terminar isso!
        P = [[at1[0]-at2[0], at1[1]-at2[1], at1[2]-at2[2]], \
                  [at1[0]-at3[0], at1[1]-at3[1], at1[2]-at3[2]], \
                                                     [0, 0, 1]]

    v1, v2 = [1, 0, 0], [0, 1, 0]

    c1, c2 = np.linalg.solve(P, v1), np.linalg.solve(P, v2)
    return [c1, c2, [0, 0, 10]]

def delete_atoms(atoms, sym=None, height=None):

    if sym != None:
        del atoms[[atom.index for atom in atoms if atom.symbol == sym]]
    elif height != None:
        del atoms[[atom.index for atom in atoms if atom.positions[2] <= heigth]]

def supercell(atoms, matrix=[1, 1, 1]):

    atoms = make_supercell(atoms, [[matrix[0], 0, 0], [0, matrix[1], 0], [0, 0, matrix[2]]])
    return atoms

def create_substrate(symbol="Pt", supercell=(1,2,1), a=3.99, orthogonal=True):

    pt = fcc111("Pt", supercell, a=3.99, vacuum=10, orthogonal=True, periodic=True)
    return pt

def from_files(path, format_file="vasp"):

    files = sorted(glob.glob(path))
    structures = [ read(i, format=format_file) for i in files ]
    return structures

def open_json(jsonfile):

    # Opening JSON file
    with open(jsonfile) as json_file:
        data = json.load(json_file)
    return data

def rsync_files(file_to_include, towhere, origin="."):

    sb.run(['rsync -nrv --include="*/" --include="'+str(file_to_include)+'" --exclude="*" ' + str(origin) + ' ' + str(towhere)], shell=True)

def make_movie(fd, atoms, gif=False):

    write_vasp_xdatcar(fd, atoms)
    if gif:
        x = read_vasp_xdatcar(file=fd, index=slice(None))
        write_gif("movie.gif", x, interval=1)

def write_image_png(atoms, supercell=(1,1,1), filetype="png", see=False):

    sb.run(["mkdir figures"], shell=True)
    if isinstance(atoms, list):
        for n, a in enumerate(atoms):
            namefig = a.get_chemical_formula(mode="hill")
            write("figures/"+namefig+"_"+str(n)+"."+filetype, a*supercell, rotation='-20z,-60x')
            sb.run(["display figures/"+namefig+"_"+str(n)+"."+filetype+" &"], shell=True)

    else:
        namefig = atoms.get_chemical_formula(mode="hill")
        write("figures/"+namefig+"."+filetype, atoms*supercell, rotation='20z,-80x')

    if see:
        view(atoms)

def oxidation_in_o(database, criteria=2, latex_table=False):

    strucs = get_structs_database(database)

    info_dict = {}
    info_list = []
    for n, s in enumerate(strucs):

        name = s.get_chemical_formula(mode='hill')
        pd_last_layer_height = max([ i.position[2] for i in s if i.symbol == "Pd" ])

        ind_p = [ (i.position[0], i.position[1], i.position[2]) for i in s if i.symbol == "In" \
                                                            and i.position[2] > pd_last_layer_height - 1 ]

        o_p = [ (i.position[0], i.position[1], i.position[2]) for i in s if i.symbol == "O" ]

        ox = []
        for xyz_ind in ind_p:

            count = 0
            for xyz_o in o_p:

                d = m.sqrt( ( xyz_ind[0] - xyz_o[0] )**2 + ( xyz_ind[1] - xyz_o[1] )**2 + ( xyz_ind[2] - xyz_o[2] )**2 )
                if d < criteria:
                    count += 1

            ox.append(count)

        oxidation_in = round(sum(ox)/len(ox), 1)
        info_dict[n] = oxidation_in
        info_list.append([n, name, oxidation_in])

        if latex_table:

            headers = ["Structure number", "Composition", "$\gamma$ (J/m$^2$)"]

    return info_dict, info_list

class Lowdin():

    def __init__(self, pathposcar, path_charge_file):

        self.pathposcar = pathposcar
        self.path_charge_file = path_charge_file

    def analyse_poscar(self):

        poscar = read(self.pathposcar, format="vasp")
        list_elements = sorted(set(poscar.get_chemical_symbols()))

        return list_elements

    def read_charge_lobster(self):

        df = pd.read_csv(self.path_charge_file, delimiter=r"\s\s+", header=0, skiprows=1, engine='python')
        df.drop(df.tail(1).index,inplace=True) # Delete last row

        elements = self.analyse_poscar()
        lowedin_dict, mulliken_dict = {}, {}
        for e in elements:

            df_elements = df.loc[df['atom'] == e]
            lowedin_dict[str(e)] = round(df_elements['Loewdin charge'].mean(axis=0), 2)
            mulliken_dict[str(e)] = round(df_elements['Mulliken charge'].mean(axis=0), 2)

        return lowedin_dict, mulliken_dict

class Bader():

    def __init__(self, pathposcar, path_acf_file):

        self.pathposcar = pathposcar
        self.path_acf_file = path_acf_file

    def analyse_poscar(self):

        # global pd_last_layer_height, poscar
        global poscar

        poscar = read(self.pathposcar, format="vasp")

        # pd_last_layer_height = max([ i.position[2] for i in poscar if i.symbol == "Pd" ])

        list_elements = sorted(set(poscar.get_chemical_symbols()), key=poscar.get_chemical_symbols().index)
        numb_elements = []
        for e in list_elements:
            numb_elements.append(poscar.get_chemical_symbols().count(e))

        return list_elements, numb_elements

    def read_acf(self, only_oxide=False):

        df = pd.read_csv(self.path_acf_file, delimiter=r"\s\s+", header=0, skiprows=0, names=["#", "X", "Y", "Z", "Charge", "Min_dist", "atomic_vol"], engine='python')
        df.drop(df.head(1).index,inplace=True) # Delete last row
        df.drop(df.tail(4).index,inplace=True) # Delete last row

        elements, numb_of_each_atom = self.analyse_poscar()
        baders = {}
        for i, (e, n) in enumerate(zip(elements, numb_of_each_atom)):

            if i == 0:
                df_element = df.loc[df.index <= n]
                if only_oxide:
                    df_element = df.loc[(df.index <= n) & (df["Z"] >= pd_last_layer_height - 1)]
                else:
                    df_element = df.loc[(df.index <= n)]
                baders[str(e)] = round(df_element['Charge'].mean(axis=0), 2)
                k = n
            else:
                j = k+n

        return baders
        # Rumpling example:
        # ind_p = np.array([ i.position[2] for i in poscar if i.symbol == "In" \
                                                            # and i.position[2] > pd_last_layer_height + 0.1 ])

        # o_p = np.array([ i.position[2] for i in poscar if i.symbol == "O" ])
        # ind_mean = np.mean(ind_p)
        # o_mean = np.mean(o_p)
        # rumpling = ind_mean - o_mean

        # return baders, rumpling

# if __name__ == "__main__":
