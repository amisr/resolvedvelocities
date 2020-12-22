.. _bins:

Velocity Bins
=============

The Heinselman and Nicolls algorithm works by defining several spatial region where the :math:`\vec{E}\times\vec{B}` velocity is assumed to be uniform but contain multiple AMISR beams.  This gives multiple independent line-of-sight velocity estimates from different look-angles, which lets the full 3D vector be resolved.  The performance of this algorithm is highly dependent on how these bins are selected.  Very large bins, particularly during periods of dynamic flows, may violate the assumption that the plasma drift velocity is constant over the bin and artificially smooth the results.  Small bins may not have a sufficient diversity of look directions to unambiguously resolve the 3D velocity and result in large errors.  Overlapping bins can be specified, but the resulting :math:`\vec{E}\times\vec{B}` velocity estimates will have significant covariance.

Specifying Bins
---------------

The vertices of velocity bins are specified through the :ref:`BINVERT <binvert>` parameter in the :ref:`configuration file <configfile>`.  Vertices should be given in magnetic coordinates.  Although any vertices can technically be selected, it usually makes sense to form regular rectangles in magnetic coordinates (latitude/longitude pairs), and makes the resulting velocities much easier to interpret.  Use the syntax of a Python nested list, where each element of the list represents a single bin and each element of that list list is the magnetic latitude and longitude of a vertex::

  BINVERT = [[[65.0,-95.],[65.0,-88.],[65.5,-88.],[65.5,-95.]],
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

When running `resolvedvelocities` on RISR-N or RISR-C data, make sure to specify rotation coordinates in the :ref:`MARPROT <marprot>` parameter in the :ref:`configuration file <configfile>` (setting MARPROT to the Apex coordinates of Resolute Bay generally work well).  Then, define bin vertices in MARP coordinate.  In this coordinate system, most bins should be very close to the equator.
