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

.. code:: python3

    #################### DEFINE THE STRUCTURES HERE ########################

    cu = bulk('Cu', 'fcc', a=3.6)
    struc = [cu]

    #######################################################################

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

.. code-block:: console

    $ python run_$structure.py 

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

Now, this file's configuration will be devided into two parts: A **permanent** configuration, in which you need to configure only once for your machine, and the defining of the variables.

For the **permanent** configuration, the first thing you need is to define the ``machine`` bash variable with the same name as your supercomputer.
Let's say you are in a cluster called **samba**:

::

    machine="samba"

Now you also need to add a ``elif`` under the ``Creating the vasp.slurm`` commentary. 

.. note::
    
    Remember that you needed to create a function for the JOB in the vulcan.py script. In our example, this function should be called ``launch_samba``.

Therefore, in our example, you should add the line:

::

    elif [[ $machine == "samba" ]]; then
        echo -e "from main_launch import *\nx = Calculo('$here' ,'$work', '$database', '$xc', $encut)\nx.launch_samba($structure)" > launch_$structure.py

Next, under the ``LAUNCH`` commentary, you need to add another ``elif``:

::

    elif [[ $machine == "samba" ]]; then
        echo "Samba!"
        sbatch run_$structure.slurm > output_sbatch; awk '{ print $4 }' output_sbatch > $here/jobid_$structure

.. note::
    
    I'm considering that your machine have a single partition. If that's not the case, you can add a variable in the ``launch_samba`` function with the partition and call it using a bash variable in the vulcan.sh script with the name ``partition``. Follow the ``explor`` and ``occigen`` examples if you want to add yours.

This is the command to submit your calculation to the queue. In my case, I launch with the **sbatch** command, but you modify it according to your machine.
We're done with the permanent configuration! Let's move on to the last part (finally!).

.. _variables:

Variables on the vulcan.sh
--------------------------

The last thing you need to configure are the variables at the top of the file. I'll put a standand configuration and explain it.

::

    #-------Parameters------#
    PROJECT="VULCAN_PROJECT"
    namework="tutorial_vulcan"
    machine="samba"
    partition_explor="std" # On Explor --"std", "sky", "mysky"
    partition_occigen="HSW24" # On Occigen -- "BDW28", "HSW24"
    xc="pbe"

    vasp="True"
    qe="False"
    #------Calculation------#
    ncore_test="False"

    slab="False"
    bader="False"

    #########
    md="False"
    gamma="False"
    relax="True"
    dos="False"
    stm="False"
    ##
    cohp="False"
    nbands="500"
    ##

    ##
    adsorption="False"
    plane_of_separation="9.2" # In the z direction
    calculate_chgdiff="False"
    ##

    ##########
    #-------Parameters------#

The variables ``PROJECT`` and ``namework`` will define the directory's names in the $SCRATCH to go in and launch the calculation. We set ``vasp="True"`` (I'm working on adding Quantum Espresso to work with this code as well).
Next, in the "Calculations" section, there are several really important keys for the calculation. Each of them will activate one part of the ``run_step_relax`` function in the vulcan.py script, so you can choose what type of calculation you want to do. This code can do multiple types of series calculations for each structure, following the order that the variables appears. 
In this example, I'm activating only the ``relax`` flag, meaning that I want to do a simple relaxation in the Cu-FCC structure.

.. note::

    The ``run_step_relax`` function is faily optimized, because it will make backup files of the important files on each type of calculation. If you want to restart your calculation when it stops, it will verify what are the flags that finished well and continue exactly from where you stoped. This way you don't need to be afraid to re-launch all the calculations.

And that's it! Just watch it run.

.. _backup:

Automated backup
================

At the end of the calculation, there are several functions to save the important files as a backup. The way the script is structured now can be a bit heavy on the backup (I'm tired of losing files), but if you need it is possible to modify it easily on the ``run_step_relax`` function.

The most important backup that the script does is to save the current state of each structure to a SQLite3 database (``.db format``), this way you'll have all your results in a single place.
The second layer of the backup is done by copying several important files the each structure and saving it to a **CALC_STATES** folder in the same directory that you launched your calculations. In this example it will be in the "~/tutorial-vulcan" folder.
