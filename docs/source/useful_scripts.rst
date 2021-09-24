.. _useful_scripts:

Useful scripts
==============

I developed a few scripts that can help in the result's analisys of this code, specially when it comes to the final results on the SQLite3 database, so we'll begin by this script.
It is a good idea to put the folder where you git cloned these files in your path. So in your ``~/.bashrc`` you should add:

::

    export PATH=${PATH}:~/path_to_folder_with_scripts

This way you can use this script anywhere.

.. _search_db:

search_db.py
------------

This script is used to get some simple information on the structures that a .db file has.
A basic usage is with the command:

.. code-block:: console

    $ search_db -f {database_name}.db

For a sample database where I have the results of the H2O molecule, the output is:

::

    There are 1 structures in this database!
    ID | Formula | T. Energy (eV) | Max force (eV/A) | a (A) | b (A) | c (A)
    1   H2O   -14.22402066   0.015  7.94 7.94 7.94


    Only valuable for structures with the same composition:
    The structure 1 is the most stable with the total energy of -14.22402066 eV

Another interesting thing that can be done is select and visualize the structures that you want by ID. Here I'll select ID = 1 (which is the only structure on the database):

::

    search_db -f h2o.db -ids 1 -v True

It will display not only the information regarding this structure but also will use the ``view`` function of ASE to plot its structure.

.. _forces:

forces.py
---------

This is a simple script to output the forces according to the OUTCAR file. It is useful to have an idea of how far the structure is from convergence when doing a structure relaxation.
You only need to specify the -c flag, in which the program will search for all atoms with forces higher than this (units in eV/Angs.).

::

    forces -c 0.02
