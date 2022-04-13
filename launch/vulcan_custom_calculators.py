#!/usr/bin/env python

# Author: Thiago Trevizam Dorini

from ase.calculators.vasp import Vasp

def general_calculator(vasp=False, qe=False, encut=400):

    if vasp == True:

        abinit = Vasp(istart=0, algo = 'Fast', nelm=80,setups='recommended',\
            ediff=1e-5, ispin=2, ediffg=-0.01, ismear=1, xc = 'pbe',\
            prec = 'Accurate', lreal = 'Auto', kpts=(1,1,1), gamma=True,\
            ibrion=2, isif=2, nsw=0, laechg = False, lorbit = 11, lvtot=False, lwave=False,\
            ncore=8, encut=encut)

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

def vibration(encut=500):

    abinit = Vasp(istart=0, algo = 'Fast', nelm=60,setups='recommended',\
        ediff=1e-8, ediffg=-0.01, ismear=2, xc = 'pbe',\
        prec = 'Accurate', lreal = False, kpts=(1,1,1), gamma=True,\
        ibrion=5, nfree=2, isif=2, nsw=500, ncore=4, encut=encut)
    return abinit

# if __name__ == "__main__":
