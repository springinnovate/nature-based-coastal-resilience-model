History
=======

Unreleased
----------

Added
~~~~~

- Added an AOI coastline classification utility that clips local coastline
  linework to a polygon AOI, segments it by Voronoi cells generated from nearby
  GCC points, assigns each line segment to its nearest GCC point, and writes a
  classified shoreline vector output.
- Documented classifier usage, required local inputs, output fields, and
  geospatial Python dependencies.
