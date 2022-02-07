.. module:: ctseg

.. _accessors:

Solver Accessors (Python)
===========================


Instances of the class Solver have a number of accessors, some of which are used as inputs to the calculation, some of which are used to retrieve the outputs of the calculation.

Input accessors
----------------
The input accessors are the following:

    +-------------+---------------------------------------------------------------+-----------------+
    | Member      | Description                                                   | Type            |
    +=============+===============================================================+=================+
    | G0          | Bare Green's function                                         | BlockGfImFreq   |
    +-------------+---------------------------------------------------------------+-----------------+
    | K_tau       | Retarded interaction kernel                                   | GfImTime        |
    +-------------+---------------------------------------------------------------+-----------------+
    | Kprime_tau  | Derivative of retarded interaction kernel                     | GfImTime        |
    +-------------+---------------------------------------------------------------+-----------------+
   
For calculations with static interactions, only ``G0`` is needed.

Output accessors
------------------

The output accessors are the following (they are read-only):

    +-------------+---------------------------------------------------------------+-----------------+
    | Member      | Description                                                   | Type            |
    +=============+===============================================================+=================+
    | nn_tau      | Density-density correlation function in imaginary time        | GfImTime        |
    +-------------+---------------------------------------------------------------+-----------------+
    | nn_omega    | Density-density correlation function on Matsubara frequencies | GfImFreq        |
    +-------------+---------------------------------------------------------------+-----------------+
    | G_tau       | Imaginary-time Green's function                               | BlockGfImTime   |
    +-------------+---------------------------------------------------------------+-----------------+
    | F_tau       | Improved estimator in imaginary time                          | BlockGfImTime   |
    +-------------+---------------------------------------------------------------+-----------------+
    | G_legendre  | Imaginary-time Green's function                               | BlockGfLegendre |
    +-------------+---------------------------------------------------------------+-----------------+
    | F_legendre  | Improved estimator in imaginary time                          | BlockGfLegendre |
    +-------------+---------------------------------------------------------------+-----------------+
    | G_omega     | Matsubara Green's function                                    | BlockGfImFreq   |
    +-------------+---------------------------------------------------------------+-----------------+
    | F_omega     | Improved estimator on Matsubara frequencies                   | BlockGfImFreq   |
    +-------------+---------------------------------------------------------------+-----------------+
    | Sigma_omega | Self-energy on Matsubara frequencies                          | BlockGfImFreq   |
    +-------------+---------------------------------------------------------------+-----------------+
    | nn          | Density-density static correlations                           | numpy.array     |
    +-------------+---------------------------------------------------------------+-----------------+
    | hist        | Histogram of perturbation order                               | numpy.array     |
    +-------------+---------------------------------------------------------------+-----------------+
    | state_hist  | Histogram of perturbation order resolved by state             | numpy.array     |
    +-------------+---------------------------------------------------------------+-----------------+
      
 The data contained by the output accessors is meaningful only if the corresponding measurement (see :doc:`options <options>`) has been turned on. Here is a list of the measures and the corresponding outputs:

    +-------------------+------------------------+-------------+
    | Measure name      | Measured observable    | Byproduct   |
    +===================+========================+=============+
    | measure_gt        | G_tau, F_tau           |             |
    +-------------------+------------------------+-------------+
    | measure_gl        | G_legendre, F_legendre |             |
    +-------------------+------------------------+-------------+
    | measure_gw        | G_omega, F_omega       | Sigma_omega |
    +-------------------+------------------------+-------------+
    | measure_nn        | nn                     |             |
    +-------------------+------------------------+-------------+
    | measure_nnt       | nn_tau                 |             |
    +-------------------+------------------------+-------------+
    | measure_nnw       | nn_omega               |             |
    +-------------------+------------------------+-------------+
    | measure_hist      | hist                   |             |
    +-------------------+------------------------+-------------+
    | measure_statehist | state_hist             |             |
    +-------------------+------------------------+-------------+

 Hence, if you e.g. switch on ``measure_gw``, you should not attempt to access ``G_tau`` : the corresponding container will be empty.

