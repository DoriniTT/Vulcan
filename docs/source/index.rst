Welcome to Vulcan's documentation!
===================================

**Vulcan** is a Python code to make life easier when working with VASP. It optimized all steps in the calculation process, from creating the structure to be studied to making a backup of the results on a SQL database. It was developed to work exclusively on VASP (for now) and in a remote supercomputer cluster with a queue system, and makes heavy use of the Atomic Simulation Environment (ASE).

.... figure:: vulcan_flow.png
   ..:scale: 50 %

   ..This is the caption of the figure (a simple paragraph).

..The flowchart show a schematic representation of how the code works. The three initial steps need to be defined as follows: (1) The code takes as initial input a list containing the ASE Atoms objects that one intends to study. The user can use any available method to construct its structure, but at the end it must be converted to the Atoms object. (2) Configure the VASP calculator according to the input structures and the type of calculation that it is intended to be performed. Keywords that are dependent on the type of calculation, such as number of ionic steps, construction of the charge density, projection of atomic orbitals, etc., are automatically changed depending on what calculation you want to perform, and only those parameters that are system dependent, such as how to determine the partial occupancies, size of the plane wave basis set, dipol corrections, etc., are kept fixed. The $k$-point grid is determined automatically depending on the size of the systems, but one can also set it manually. (3) Configure the job script function to work with the queue system of your cluster. All these three steps are detailed in the manual.

..After defining the inputs, the user can define all the calculations intended, which will run in parallel for each structure and in sequence for each type of calculation. If all calculations finished correctly, the code will make a backup of some selected important files in the \$HOME folder for postprocessing. If one has identified that some of the calculations have not finished well, simply running the code again will restart those calculations.

..At the end, with all important files gathered, there are several implemented functions to analyse the outputs, such as the Bader charge, charge density difference, DOS, COHP, etc., and can be easily adapted depending on each analysis intended.

.. note::

   This project is under active development.

Contents
--------

.. toctree::

   usage
   configuration
   useful_scripts
   calculations
