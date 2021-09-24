.. _configuration1:

Configuration - vulcan.py
=========================

This is the most important part. Let's do a step-by-step general configuration with a example structure.

I'll assume that you have access to a supercomputer to launch your calculations as a job (i.e. Slurm). I'll show later how to launch your calculations in your desktop machine later. 
To begin with, create a new directory in your $HOME (~) for this project and put both the vulcan.py and vulcan.sh there.

.. code-block:: console

   $ mkdir ~/tutorial-vulcan
   $ cd ~/tutorial-vulcan

Assuming that these files are already in your home:

.. code-block:: console

   $ cp ~/vulcan.* ~/tutorial-vulcan

After this step, we can start configuring these files.

.. _structure:

Structure generation
--------------------

In general, the hole python script is centered around the ``structures`` function, which returns a list with all the structures that we want to study.
All you have to do is create an ``Atoms`` ASE object and put it in the list.
In the vulcan.py script, this needs to be done under the **DEFINE THE STRUCTURES HERE** commentary. 

.. note::
    
    If you are not familiar with the Atomic Simulation Environment, please take a look at their documentation on how to create an **Atoms** object.

We will begin creating an Cu-FCC structure:

.. code:: python

    cu = bulk('Cu', 'fcc', a=3.6)
    struc = [cu]

Here we will put only one structure to make things easier, but the program is made so that you can put as much structures as you want.
That's all you need for this part. It becames more complex according to the structures that you want to create, but this part is entirely up to you.

.. _calculation:

Calculation setup
-----------------

Now you need to focus on the ``calculators`` function. I have already put some default parameters for a DFT calculation with VASP, but you need to modify it according to what you need. For the functions implemented on this code, these default parameters usually works just fine.
You can look at the VASP wiki for a description of all input parameters for VASP: https://www.vasp.at/wiki/index.php/The_VASP_Manual

.. _launch:

Launch script setup
-------------------

The last thing that you need to create is a ``launch_{your-machine}`` function (it is located right after the ``run_step_relax`` function). I have 4 functions created depending on the machine that I'm using. You need to create your own function based on a JOB script that you want.
In my case, I automitized the number of processors that the JOB script will launch based on each structure's number of atoms. You can take advantage on the code's defined variables to adjust how you automitize your JOB.

**Important:** The only line that you cannot touch in the JOB script is the: 

.. code:: python

    python run_$structure.py 

Please leave this line as the last one in your job, otherwise the code will not work.

.. _optional:

Optional setup
--------------

Although this is a optional, I highly recommend you to take a look on the ``run_step_relax`` function on the vulcan.py script. Depending on the type of calculation you intend to do, this is the function that you need to modify. 
Again, it has default parameters that works well for some cases.

.. note::

    The **run_step_relax** function is the most versitile function in the code. Please try to understand its format if you want to add a new type of calculation.

.. _configuration2:

Configuration - vulcan.sh
=========================


