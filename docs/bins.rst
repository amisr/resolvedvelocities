.. _bins:

Velocity Bins
=============

The Heinselman and Nicolls algorithm works by defining several spatial regions where the :math:`\vec{E}\times\vec{B}` velocity is assumed to be uniform and each region of uniform flow is intersected by multiple radar beams.  The Apex components of the :math:`\vec{E}\times\vec{B}` velocity are constant along magnetic field lines, so the full 3D vector can be resolved from multiple independent line-of-sight velocity measurements with different look directions.  The performance of this algorithm is highly dependent on how these bins are selected.  Very large bins, particularly during periods of dynamic flows, may violate the assumption that the plasma drift velocity is constant within the bin and artificially smooth the results.  Bins that are too small may not have a sufficient diversity of look directions to unambiguously resolve the 3D velocity and result in large errors.  Overlapping bins can be specified to improve the spatial resolution, but the resulting :math:`\vec{E}\times\vec{B}` velocity estimates will have significant covariance.

Specifying Bins
---------------

There are two ways to specify bins int the :ref:`configuration file <configfile>`.

Regular magnetic latitude bins can be specified through the :ref:`BINMLATDEF <binmlatdef>` parameter. This field takes two parameters, `step, stride`. The `step` parameter refers to the distance between the start of each bin while the `stride` parameter is the width of each bin. Bins start at the magnetic latitude of the radar and spread both north and south to cover the entire range of data available for that paritular experiment. If `step < stride`, the bins will overlap, which is normal and helps achieve a higher spatial resolution with lower error on each bin.  Below is an example of a relatively standard setup with only one bin in magnetic longitude and several overlapping bins in magnetic latitude::

  BINMLATDEF = 0.25, 0.5

Arbitrary vertices for bins can be specified through the :ref:`BINVERTDEF <binvertdef>` parameter.  Vertices should be given in magnetic coordinates.  It is not necessary for bins to form regular rectangles and any number of vertices may be provided for each bin.  This gives maximum flexibility at the expense of an unwieldy field in the config file.  Use the syntax of a Python nested list, where each element of the list represents a single bin and each element of that list list is the magnetic latitude and longitude of a vertex::

  BINVERTDEF = [[[65.0,-95.],[65.0,-88.],[65.5,-88.],[65.5,-95.]],
   [[65.5,-95.],[65.5,-88.],[66.0,-88.],[66.0,-95.]],
   [[66.0,-95.],[66.0,-88.],[66.5,-88.],[66.5,-95.]],
   [[66.5,-95.],[66.5,-88.],[67.0,-88.],[67.0,-95.]],
   [[67.0,-95.],[67.0,-88.],[67.5,-88.],[67.5,-95.]]]

Auroral Zone
------------

In the auroral zone (PFISR), it is reasonable to assume most structures and velocity shears are oriented along lines of constant magnetic latitude.  Therefore, at PFISR it makes sense to select bins that are elongated in Apex latitude across the radar's field of view.  Typically, bins are selected that are 0.5 degrees wide in magnetic latitude and overlap by 0.25 degrees.

Polar Cap
---------

Choosing appropriate bins in the polar cap (RISR-N and RISR-C) can be much more challenging because it is rarely appropriate to assume flows align along lines of magnetic latitude in regions dominated by magnetosphere-driven convection on open field lines.  In fact, the proximity of Resolute Bay to the magnetic pole introduces additional complications purely due to the singularity in the coordinate system.  To resolve this, the resolved velocities algorithm is usually run in the Modified Apex Rotated Poles (MARP) coordinate system, which is an analogous magnetic coordinate system to the Modified Apex coordinate system used at PFISR but rotated so that Resolute Bay is located on the equator of the coordinate system.

When running `resolvedvelocities` on RISR-N or RISR-C data, make sure to specify rotation coordinates in the :ref:`MARPROT <marprot>` parameter in the :ref:`configuration file <configfile>` (setting MARPROT to the geodetic coordinates of Resolute Bay and the radar boresight generally works well).  Then, define bin vertices in MARP coordinate.  In this coordinate system, most bins should be very close to the equator.  When this option is used, bin coordinates are transformed back to standard Apex before being saved to the output datafile to simplify the data product for the end user.
