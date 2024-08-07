.. _configfile:

Configuration File
==================

The configuration file provides important specifications about how the `resolvedvelocities` code will be run.  The following table summarizes the categories and parameters that should be in the configuration file.  Example values are provided where appropriate.

+----------+---------------------------------------+--------------------------------+---------------+--------------+---------------+
| Category | Parameter                             | Description                    | PFISR         | RISR-N       | RISR-C        |
+==========+=======================================+================================+===============+==============+===============+
| FILEIO   | :ref:`DATAFILE <datafile>`            | Input file name                | 20190328.006_lp_1min-fitcal.h5               |
|          +---------------------------------------+--------------------------------+---------------+--------------+---------------+
|          | :ref:`OUTPUT_NAME <outfilename>`      | Output file name               | 20190328.006_lp_1min-fitcal_vvels.h5         |
|          +---------------------------------------+--------------------------------+---------------+--------------+---------------+
|          | :ref:`OUTPUT_PATH <outfilepath>` \+   | Output file path               | /home/users/vvels/output/                    |
+----------+---------------------------------------+--------------------------------+---------------+--------------+---------------+
| OPTIONS  | :ref:`CHIRP <chirp>` \+               | Correction for frequency chirp | 0.0 \*                                       |
|          +---------------------------------------+--------------------------------+---------------+--------------+---------------+
|          | :ref:`NELIM <nelim>` \+               | Density limits                 | 2.0e9, 1.0e13 \*                             |
|          +---------------------------------------+--------------------------------+---------------+--------------+---------------+
|          | :ref:`CHI2LIM <chi2lim>` \+           | :math:`\chi^2` limits          | 0.1, 10. \*                                  |
|          +---------------------------------------+--------------------------------+---------------+--------------+---------------+
|          | :ref:`GOODFITCODE <goodfitcode>` \+   | Fit codes that can be used     | 1, 2, 3, 4 \*                                |
|          +---------------------------------------+--------------------------------+---------------+--------------+---------------+
|          | :ref:`COVAR <covar>`                  | A priori covariance matrix     | 9.0e6,9.0e6,2.5e3 \*                         |
|          +---------------------------------------+--------------------------------+---------------+--------------+---------------+
|          | :ref:`INTTIME <inttime>` \+           | Post-integration time          | 180                                          |
|          +---------------------------------------+--------------------------------+---------------+--------------+---------------+
|          | :ref:`USE_BEAMS <use_beams>` \+       | Beam codes to use              | 64016, 64157, 64964                          |
+----------+---------------------------------------+--------------------------------+---------------+--------------+---------------+
| VVELSLAT | :ref:`BINMLATDEF <binmlatdef>`        | MLAT steps for regular bins    | see :ref:`Velocity Bins <bins>`              |
|          +---------------------------------------+--------------------------------+---------------+--------------+---------------+
|          | :ref:`BINVERTDEF <binvertdef>`        | Vertices of each mlat bin      | see :ref:`Velocity Bins <bins>`              |
|          +---------------------------------------+--------------------------------+---------------+--------------+---------------+
|          | :ref:`MARPROT <marprot>` \+           | MARP rotation coordinates      | 0.0, 0.0, 0.0 | 74.7, -94.4, 26.0            |
|          +---------------------------------------+--------------------------------+---------------+--------------+---------------+
|          | :ref:`ALTLIM <altlim>` \+             | Altitude limits                | 150., 400. \*                                |
|          +---------------------------------------+--------------------------------+---------------+--------------+---------------+
|          | :ref:`OUTALTS <outalts>`              | Altitudes of geodetic output   | 200.,250.,300.,350.,400.,450.,500.           |
|          +---------------------------------------+--------------------------------+---------------+--------------+---------------+
|          | :ref:`UPB_BEAMCODE <upb_beamcode>` \+ | Beam code of "up-B" beam       | 64157 \*      |              | 65426 \*      |
|          +---------------------------------------+--------------------------------+---------------+--------------+---------------+
|          | :ref:`IONUP <ionup>` \+               | Ion upflow correction method   | UPB \*        |              | UPB \*        |
+----------+---------------------------------------+--------------------------------+---------------+--------------+---------------+
| VVELSALT | :ref:`ALTBIN <altbin>`                | Altitude bins                  | 80,150,4.5,9;150,300,20,20                   |
|          +---------------------------------------+--------------------------------+---------------+--------------+---------------+
|          | :ref:`MLATLIM <mlatlim>` \+           | Magnetic latitude limits       |               |              |               |
+----------+---------------------------------------+--------------------------------+---------------+--------------+---------------+
| PLOTTING | :ref:`PLOTSAVEDIR <plotsavedir>`      | Directory to save plots to     | /home/user/vvels/plots                       |
+----------+---------------------------------------+--------------------------------+---------------+--------------+---------------+
| PLOTTING | :ref:`PLOTPREFIX <plotprefix>`        | Prefix for summary plot names  | 20190328.006_lp_1min-fitcal                  |
+----------+---------------------------------------+--------------------------------+---------------+--------------+---------------+

\+ Optional parameter

\* Value is recommended for most normal use cases

Example configuration file for PFISR and RISR-N are provided in the root of the `source code directory <https://github.com/amisr/resolvedvelocities>`_.  These can be copied and adapted as necessary.


Detailed Descriptions
---------------------

.. _datafile:

DATAFILE: The name of the input processed AMISR data file the algorithm should be run on.  This can include a path if the file is not in the run directory.  This should be a hdf5 file output from the standard AMISR fitting routine.

.. _outfilename:

OUTFILENAME: The name of the output resolved velocities hdf5 data file.

.. _outfilepath:

OUTFILEPATH: (Optional) The path where the output resolved velocities hdf5 data file should be saved.  If not provided, the file will be saved in the run directory.

.. _chirp:

CHIRP: (Optional) The line-of-sight velocity correction that should be added to the measured line-of-sight velocity to correct for a frequency chirp in the radar transmit pulse.  For the AMISRs, this correction is only necessary for some RISR-N experiments.

.. _nelim:

NELIM: (Optional) A lower and upper limit of electron density for data that should be included in the velocity reconstruction.  In general, anomalously high or low density values are an indication that the ISR fitting procedure failed and the velocity value should not be trusted.  Data points with density values outside these limits will be filtered out.

.. _chi2lim:

CHI2LIM: (Optional) The lower and upper limits of :math:`\chi^2` values from the processed AMISR data file that should be included in the velocity reconstruction.  The :math:`\chi^2` parameter indicates goodness of fit of the original ISR fitting procedure.  Large values indicate excessively large errors while small values suggest "overfitting" and that the values should not be trusted even if errors are small.  Data points with :math:`\chi^2` values outside this range will be filtered out.

.. _goodfitcode:

GOODFITCODE: (Optional) The list of fit codes (assigned by the AMISR fitter) that indicate a successful fit.  Generally, fit codes :math:`\le 0` or :math:`\ge 5` mean the fit failed and that point should be excluded.

.. _covar:

COVAR: The a priori covariance matrix.  This is the expected variance of velocity in each of the three apex directions (e1,e2,e3) in the ionosphere (described as :math:`\Sigma_\nu` in Heinselman and Nicolls, 2008).

.. _inttime:

INTTIME: (Optional) Post-integration period (in seconds) that should be used to reconstruct vectors for a single time stamp in the output file.  Because the input processed data has discrete time stamps (usually on a cadence greater than one minute), the output file will typically not have exactly this resolution unless an integration time is selected that is an exact multiple of the cadence of the input data file.  Instead, the code will post-integrate the smallest number of input time stamps that cover a period greater than or equal to the specified integration time.  If this parameter is omitted, the native times of the input file are used (no post-integration).

.. _use_beams:

USE_BEAMS: (Optional) A list of specific beam codes that should be used for the velocity reconstruction.  If this field is omitted, all available beams from the experiment will be used.

.. _binmlatdef:

BINMLATDEF: The description of magnetic latitude spacing for regular mlat bins.  The format is "step, stride" where the start and stop points are determined dynamically by the available data and bins are centered around the radar site.  See the :ref:`Velocity Bins <bins>` page for more details.


.. _binvertdef:

BINVERTDEF: The list of vertices that define each bin.  These should be given in magnetic coordinates.  See the :ref:`Velocity Bins <bins>` page for more details about how to specify these.

.. _marprot:

MARPROT: (Optional) The rotation coordinates for the MARP coordinate system (geodetic latitude of the new Null Island, geodetic longitude of the new Null Island, and azimuth of the new North direction).  Typically, the site latitude and longitude and the radar boresight direction will be appropriate.  This option is generally only used for the Resolute Bay radars.

.. _altlim:

ALTLIM: (Optional) The altitude limits (in kilometers) of data that should be included in the velocity reconstruction.  The algorithm assumes that the line-of-sight velocity measured by the radar is a component of the :math:`\vec{E}\times\vec{B}` plasma drift velocity with minimal plasma motion along the field line.  This is roughly true for plasma motion in the main F-region, but ion-neutral collisions in the E-region rotate the plasma motion direction and ion upflow along the field lines becomes important at high altitudes, so it is important to limit the range of altitudes considered.

.. _outalts:

OUTALTS: The list of altitudes (in kilometers) at which geodetic components of the output velocity and electric field should be calculated at. The output files contain the Apex components, which are constant along magnetic field lines, but for convenience, geodetic components are also included for a discrete grid defined by the locations of the bin centers and this output altitude array.

.. _upb_beamcode:

UPB_BEAMCODE: (Optional) Beam code for the beam pointing directly up the magnetic field.  Some radar modes may not include an "up-B" beam.  The RISR-N field-of-view is tilted too far North for any beam to be oriented along the magnetic field, so config files written for RISR-N experiments should NEVER include this field.

.. _ionup:

IONUP: (Optional) The method used to calculate and correct for any ion upflow component of the velocity.

.. _altbin:

ALTBIN: The description of altitude bins.  Bins are regularly spaced between a start and stop altitude based on the prescribed step size and stride, but multiple groups of bins are allowed.  The format is "start,stop,step,stride;" with different groups separated by a semicolon.  In general, it makes sense for sorter steps in the E-region and long in the F-region.

.. _mlatlim:

MLATLIM: (Optional) The limits in magnetic latitude of which points should be considered.  This can help control which data points contribute to the reconstruction, specifically removing points far from the center of the FoV where the ionospheric velocity may not be the same.

.. _plotsavedir:

PLOTSAVEDIR: The path for a directory output summary plots should be saved in.

.. _plotprefix:

PLOTPREFIX: (Optional) The string to prefix standard summary plot filenames with. This can help identify the particular experiment file the summary plots are associated with.
