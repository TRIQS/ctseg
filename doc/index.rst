.. _welcome:

The hybridization-expansion segment solver
******************************************

.. sidebar:: ctseg 3.1.0

   This is the homepage of ctseg v3.1.0.
   For changes see the :ref:`changelog page <changelog>`.
      
      .. image:: _static/logo_github.png
         :width: 75%
         :align: center
         :target: https://github.com/triqs/ctseg


The :ref:`TRIQS-based <triqslibs:welcome>` hybridization-expansion segment
solver allows to solve the problem of a **quantum impurity** embedded in a
conduction bath. It is restricted to **purely density-density interaction** 
vertices, but is very much optimized for this situation. Moreover, these interactions may be **dynamical** in the charge and longitudinal spin channel (even in a multiorbital context), as well as in the **transverse** spin channel (only for a single-orbital model). Finally, the code supports **off-diagonal** hybridization terms.

The code can be used both in C++ and Python to measure the **one-particle Green's function**, **density-density correlation functions**, **three- and four-point correlation functions** (and the corresponding **vertices**), various **improved estimators** corresponding to these correlation functions, as well as various histograms.

You can learn how to use it in the documentation:

    
.. toctree::
   :maxdepth: 2

   tour/contents
   reference/contents

To install the code, report issues or get more information, check out:

.. toctree::
   :maxdepth: 2
   :hidden:

   install
   issues
   about

.. note::

  For impurity models with static, non-density-density interactions and/or static hopping terms on the impurity, we refer the user to the TRIQS/cthyb (`<https://triqs.github.io/cthyb>`_) solver.
