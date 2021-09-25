.. _usage:

Usage
=====

.. _requirements:

Requirements
------------

- Python 3.x;
- Numpy;
- Matplotlib;
- ASE (Atomic Simulation Environment).

.. _installation:

Installation
------------

It is only necessary to get two files to use this program for now: vulcan.py and vulcan.sh. They can be found at my GitHub page:
https://github.com/DoriniTT/Vulcan

Initialization
--------------

Since we will be using the ASE calculators, you need to set your environment so that the program knows where to find the right executables and files.
For this, you can go to:

https://wiki.fysik.dtu.dk/ase/ase/calculators/vasp.html?highlight=vasp#ase.calculators.vasp.Vasp

.. note::

    The Vulcan program only supports VASP for now, but I'm working on adding more ab initio codes.

At its core, this program is really simple. It has the vulcan.sh file to manage the environment for the vulcan.py file. After the configuration step, all you have to do to launch the code is:

.. code-block:: console

   $ sh vulcan.sh

But first, we need to know how to configure these files for the type of calculation that we intend to do, which is what we will see in the next section.
