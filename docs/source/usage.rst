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

The code is divided into two essential files ("vulcan.py" and "vulcan.sh") and three complementary files ("vulcan_custom_calculators.py", "vulcan_helper.py", and "vulcan_structures.py"). These files can be obtained by cloning my github repertory:

.. code-block:: console

   $ git clone https://github.com/DoriniTT/Vulcan

I use the complementary files as **modules** for the essential files, so you'll need to add where to find these files in your ``PYTHONPATH``. Essentially, these files are helpful for two things: (1) Creating the structures you desire with the "vulcan_helper.py" and "vulcan_structures.py" and (2) Defining the calculator with the "vulcan_custom_calculators.py". There are a few examples inside these files that I used, so fell free to use them as a template for your case. 

By default, I store these files in ``python/launch/``. Therefore, I suggest you put them somewhere you'll easily find:

.. code-block:: console

   $ mkdir -p <path-the-git-repo>/python/launch
   $ touch <path-the-git-repo>/python/launch/__init__.py

And put the following in your ~/.bashrc:

.. code-block:: console

   $ export PYTHONPATH=${PYTHONPATH}:<path-to-git-repo>/python

This way these complementary files are acessible to the vulcan.py file and you are ready to go.

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
