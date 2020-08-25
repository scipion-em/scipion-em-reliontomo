=================
Reliontomo plugin
=================

This plugin provide wrappers around several programs of `RELION <https://www3.mrc-lmb.cam.ac.uk/relion/index.php/Main_Page>`_ software suite for its use in Tomography.

Installation
------------

You will need to use `3.0 <https://scipion-em.github.io/docs/docs/scipion-modes/how-to-install.html>`_ version of Scipion to be able to run these protocols. To install the plugin, you have two options:

a) Stable version

   .. code-block::

      scipion3 installp -p scipion-em-reliontomo

b) Developer's version

   * download repository

   .. code-block::

      git clone https://github.com/scipion-em/scipion-em-reliontomo.git

   * install

   .. code-block::

      scipion3 installp -p path_to_scipion-em-reliontomo --devel

To check the installation, simply run one of the tests. A complete list of tests can be displayed by executing ``scipion3 tests --grep reliontomo --show``

Supported versions
------------------

Relion 3.0

Protocols
---------

* ctf 3D estimation
* subtomogram 3D classification
* subtomogram 3D refinement
* subtomogram reconstruction

References
----------

1. Scheres et al., JMB, 2012 
2. Scheres et al., JSB, 2012 
3. Kimanius et al., eLife, 2016 
4. Zivanov et al., eLife, 2018
