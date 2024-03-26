ABFE
====

.. image:: https://img.shields.io/badge/License-MIT--0-teal
   :target: https://opensource.org/license/mit-0
   :alt: License

This package has been developed as a specific protocol for calculating ABFE, employing the double-system/single-box (DSSB) AbsoluteDG method for calculating absolute binding free energies through non-equilibrium switching. It involves two distinct simulations: the first focuses on the transition of a ligand from a bound state in the enzyme's active site to an unbound state within the simulation box. Conversely, the second simulation observes the ligand moving from an unbound state to a bound state within the active site. Throughout the calculations, the ligand is maintained in the active site using a system of restraints.

The package is based on tools implemented in `pmx <https://github.com/deGrootLab/pmx/tree/abfe_dev>`_ for the protocol of calculating absolute binding free energy and `sire <https://github.com/OpenBioSim/sire>`_ for calculating the correction introduced by restraints into the system.

In this protocol, all necessary tools are adapted to use a scheme of holding ligands using only pairwise distance restraints, rather than Boresh restraints, as was done in the basic approach.

For detailed technical information, you may refer to the ``ABFE_cli.py`` and ``ABFE_cli.sh`` scripts in ``ABFE/docker_build/abfe_pipeline`` directory for the main simulation launch protocol and calculation organization. In the same folder you can find ``AbsoluteDG_child.py`` for the pairwise distance restraints generation block, and ``StandardState.py`` for calculating the correction to the final free energy due to the error introduced by the presence of restraints in the system.

The package itself consists of files for assembling a Docker image ready for running simulations, as well as an example of execution.

Usage
=====

First, download the repository.

.. code-block:: bash

    git clone https://github.com/insilicomedicine/ABFE.git

Then, build the Docker image:

.. code-block:: bash

    cd ./ABFE/docker_build
    DOCKER_BUILDKIT=1 docker buildx build --progress=plain -t abfe:user .

After that, you can run the image:

.. code-block:: bash

     docker run --name abfe -it -v /path/to/your/data:/home/jovyan/data abfe:user

Inside a container, you may need to change the ownership of the ``/home/jovyan/data`` directory.

.. code-block:: bash

    cd /home/jovyan/data
    sudo chown -R jovyan:users . 

After finishing all your calculations, you will need to change the ownership back to work with the files outside of the container.

Example
=======

You can find example data for processing in 'example' folder. Just copy both sub folders (pdbs and scripts) to path/to/your/data

.. code-block:: bash

    #From root folder
    cd ./ABFE/example
    cp -r * path/to/your/data
    cd path/to/your/data

Then run the script in scripts folder

.. code-block:: bash

    cd ./scripts
    bash run_bench.sh

Data
=======

A list of 97 kinase complexes employed in the original research, along with their corresponding experimental energies, is available in  ``data/kinase_systems.csv``. 

Additionally, Jupyter Notebook with example code for calculating the stability of ligand positioning during the equilibration phase can be found in ``data/Restraint_stability.ipynb``. 
