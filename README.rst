====================================
Scipion plugin for Relion Tomography
====================================

.. image:: https://img.shields.io/pypi/v/scipion-em-reliontomo.svg
        :target: https://pypi.python.org/pypi/scipion-em-reliontomo
        :alt: PyPI release

.. image:: https://img.shields.io/pypi/l/scipion-em-reliontomo.svg
        :target: https://pypi.python.org/pypi/scipion-em-reliontomo
        :alt: License

.. image:: https://img.shields.io/pypi/pyversions/scipion-em-reliontomo.svg
        :target: https://pypi.python.org/pypi/scipion-em-reliontomo
        :alt: Supported Python versions

.. image:: https://img.shields.io/pypi/dm/scipion-em-reliontomo
        :target: https://pypi.python.org/pypi/scipion-em-reliontomo
        :alt: Downloads

This plugin provide wrappers around several programs of `RELION <https://relion.readthedocs.io/en/release-5.0/index.html>`_ software suite for its use in Tomography.


Installation
------------

You will need to use `3.0 <https://scipion-em.github.io/docs/release-3.0.0/docs/scipion-modes/how-to-install.html>`_ version of Scipion to be able to run these protocols.


Protocols
-----------

VERSION 4.0:

* **Import coordinates 3D from a star file** : Protocol to import a set of 3D coordinates from a star file.
* **Import subtomograms from a star file** : Protocol to import a set of subtomograms from a star file.
* **3D Classification of subtomograms** : 3D Classification of subtomograms.
* **Tomo CTF refine** : Tomo CTF refine
* **De novo 3D initial model** : Generate a de novo 3D initial model from the pseudo-subtomograms.
* **Apply operation to Relion particles** : Operate on the particles star file.
* **Make pseudo-subtomograms** : Make pseudo-subtomograms.
* **Sharpen a 3D reference map** : Sharpen a 3D reference map and estimate the gold-standard FSC curves for subtomogram averaging.
* **Prepare data for Relion 4** : Prepare data for Relion 4
* **Reconstruct particle from tilt series** :  Reconstruct particle from the original tilt series images.
* **Auto-refinement of subtomograms** :  Auto-refinement of subtomograms.
* **Rec. particle averaging subtomograms** : This protocol reconstructs a volume using Relion. Reconstruct a volume from a given set of particles. The alignment parameters will be converted to a Relion star file and used as direction projections to reconstruct.
* **Tomo frame align** : Tomo frame align
* **Reconstruct tomograms from prepare data prot** : This protocol reconstructs a single tomogram using Relion. It is very useful  to check if the protocol "Prepare data" has been applied correctly (in terms of flip  options, for example).

VERSION 5.0:

* **Import subtomograms from a star file** : Protocol to import a set of subtomograms from a star file.
* **3D Classification of subtomograms** : 3D Classification of subtomograms.
* **Tomo CTF refine** : apply the refinement of the different parameters related to the CTF estimation.
* **De novo 3D initial model** : Generate a de novo 3D initial model from the pseudo-subtomograms.
* **Apply operation to Relion particles** : Operate on the particles star file: re-center and duplicates removal.
* **Sharpen a 3D reference map** : Sharpen a 3D reference map and estimate the gold-standard FSC curves for subtomogram averaging.
* **Extract subtomos** : extract the subtomograms from the tilt-series.
* **Reconstruct particle** : generate a reference map from a particle set by averaging the particles using their positions and orientations
* **Auto-refinement of subtomograms** :  Auto-refinement of subtomograms.
* **Bayesian polishing** : refine the projections that map 3D space onto the images of the tilt series.
* **Tomogram reconstruction** : reconstruct tomograms.
* **Motion correction** : motion correction of the angular frame stacks using Relion's own implementation.


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

Relion 4.0, 5.0


References
----------

1. Zivanov et al., eLife, 2022
2. Burt et al., bioRxiv, 2024


