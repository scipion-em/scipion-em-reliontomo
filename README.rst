=================
Reliontomo plugin
=================

This plugin provide wrappers around several programs of `RELION <https://www3.mrc-lmb.cam.ac.uk/relion/index.php/Main_Page>`_ software suite (version 3.0) for its use in Tomography.


Current development
-------------------

This plugin is currently in **BETA** mode.



Installation
------------

You will need to use `3.0 <https://scipion-em.github.io/docs/docs/scipion-modes/how-to-install.html>`_ version of Scipion to be able to run these protocols.



Protocols
-----------

* **Extract coordinates from pseudo-subtomograms** : Protocol to extract a set of 3D coordinates from a set of pseudo-subtomograms.
* **Import coordinates 3D from a star file** : Protocol to import a set of 3D coordinates from a star file.
* **Import subtomograms from a star file** : Protocol to import a set of subtomograms from a star file.
* **3D Classification of subtomograms** : 3D Classification of subtomograms.
* **Tomo CTF refine** : Tomo CTF refine
* **De novo 3D initial model** : Generate a de novo 3D initial model from the pseudo-subtomograms.
* **Apply operation to Relion particles** : Operate on the particles star file.
* **Make pseudo-subtomograms** : Make pseudo-subtomograms.
* **Sharpen a 3D reference maps** : Sharpen a 3D reference map and estimate the gold-standard FSC curves for subtomogram averaging.
* **Prepare data for Relion 4** : Prepare data for Relion 4
* **Reconstruct particle from tilt series** :  Reconstruct particle from the original tilt series images.
* **Auto-refinement of subtomograms** :  Auto-refinement of subtomograms.
* **Rec. particle averaging subtomograms** : This protocol reconstructs a volume using Relion. Reconstruct a volume from a given set of particles. The alignment parameters will be converted to a Relion star file and used as direction projections to reconstruct.
* **Tomo frame align** : Tomo frame align
* **Reconstruct tomograms from prepare data prot** : This protocol reconstructs a single tomogram using Relion. It is very useful  to check if the protocol "Prepare data" has been applied correctly (in terms of flip  options, for example).
* **Matching coordinates** : Protocol to generate a set of pseudosubtomograms taking into account the intersection of  a set of 3D coordinates and a set of pseudosubtomograms as input.

**Latest plugin versions**
==========================

If you want to check the latest version and release history go to `CHANGES <https://github.com/scipion-em-reliotomo/reliontomo/blob/master/CHANGES.txt>`_


**Installing the plugin**
=========================

In order to install the plugin follow these instructions:

.. code-block::

    scipion installp -p scipion-em-reliontomo


or through the **plugin manager** by launching Scipion and following **Configuration** >> **Plugins**


**To install in development mode**

Clone or download the plugin repository

.. code-block::

    git clone https://github.com/scipion-em/scipion-em-reliontomo.git

Install the plugin in developer mode.

.. code-block::

    scipion installp -p local/path/to/scipion-em-reliontomo --devel


To check the installation, simply run one of the tests. A complete list of tests can be displayed by executing ``scipion3 tests --grep reliontomo --show``

Supported versions
------------------

Relion 3.0, 4.0


References
----------

1. Scheres et al., JMB, 2012 
2. Scheres et al., JSB, 2012 
3. Kimanius et al., eLife, 2016 
4. Zivanov et al., eLife, 2018


Buildbot status
---------------

Status devel version:

.. image:: http://scipion-test.cnb.csic.es:9980/badges/tomo_devel.svg


Status production version:

.. image:: http://scipion-test.cnb.csic.es:9980/badges/tomo_prod.svg

