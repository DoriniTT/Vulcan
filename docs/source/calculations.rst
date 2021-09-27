.. _calculation:

Types of calculations
=====================

There are already a few available routines that you can perform using this code. It is relativelly simple to implement new routines based on the ones already available.

Let's explain what each routine is doing.

TO BE DONE

.. _gamma:

Gamma-point
-----------

Considering that the structure is far from fully relaxed, the idea of this function is to get a rough relaxation with the :math:`{\\Gamma}`-point. 
Before each calculation, the only thing that it does is verify if there is the line ``reached required`` on the OUTCAR_gam file. If it does not, the calculation will begin either from the CONTCAR file (if it exists), or from a new POSCAR.

.. _relax:

Relaxation
----------

The same idea as on the :math:`{\\Gamma}`-point calculation but now using a higher number of *k*-points. This time it will verify the ``reached required`` phrase on the OUTCAR_relax file. 
The number of *k*-points is automatically set by using the following formula:

.. math::

   k = (precision*2\pi/a, precision*2\pi/b, precision*2\pi/c)

where ``precision = 6`` by default. You can set this higher if needed. If the flag ``slab`` is set to True in the vulcan.sh file, ``k[2] = 1`` is used.

.. _md:

Molecular dynamics
------------------

.. _dos:

Density of States
-----------------

.. _cohp:

COHP
----

.. _adsorption:

Adsorption
----------

.. _stm:

STM
---
