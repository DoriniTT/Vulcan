.. _calculation:

Types of calculations
=====================

There are already a few available routines that you can perform using this code. It is relativelly simple to implement new routines based on the ones already available.

Let's explain what each routine is doing.

TO BE DONE

.. _gamma:

Gamma-point
-----------

Considering that the structure is far from fully relaxed, the idea of this function is to get a rough relaxation with the :math:`{\\gamma}`-point. 
Before each calculation, the only thing that it does is verify if there is the line ``reached required`` on the OUTCAR_gam file. If it does not, the calculation will begin either from the CONTCAR file (if it exists), or from a new POSCAR.

.. _relax:

Relaxation
----------

The same idea as on the :math:`{\\gamma}`-point calculation but now using a higher number of *k*-points. This time it will verify the ``reached required`` phrase on the OUTCAR_relax file. 
The number of *k*-points is automatically set by using the following formula:

.. math::

   k = (precision*2\pi/a, precision*2\pi/b, precision*2\pi/c)

where ``precision = 6`` by default. You can set this higher if needed. If the flag ``slab`` is set to True in the vulcan.sh file, ``k[2] = 1`` is used.

At the end of the calculation, if everything goes well, it will save the files defined in the ``input`` variable as ``*_relax``. You can add or remove files in this variable as you please.

.. _dos:

Density of States
-----------------

It uses the same idea as the `relax` function but setting different inputs for the INCAR. It will double the number of *k*-points and set the following inputs: ``icharg=11, ismear=-5, prec='Accurate', nsw=0``.
If the flag ``results`` is True on the vulcan.sh file, it will also plot the total Density of States to a file called ``dos_results`` in the launch directory.

.. _cohp:

COHP
----

Similar idea to the Density of States calculation, but with different flags following a COHP calculation. As can be seen in the Lobster manual http://schmeling.ac.rwth-aachen.de/cohp/index.php?menuID=6, one very important parameter is the number of bands, which you need to set manually on the vulcan.sh file.  

After a successful calculation, this routine will also calculate the COHP using the ``lobster`` program, so be sure to have the binary file on your path. (You can get the script for free on their website).
I set a very simple ``lobsterin`` file in the program as using only the recommended basis functions from the POTCAR file, feel free to set your own file according to what you need. You only need to modify the `lobsterin` string on the vulcan.py script.

.. _stm:

STM
---

.. _adsorption:

Adsorption
----------

.. _md:

Molecular dynamics
------------------
