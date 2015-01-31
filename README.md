# grid-physics
Exploring the possibilities for flexible movement on a grid.

<p align="center"><img src="grid_physics_screenshot.png" width="400"></p>

Problem:
--------
Previous models have typically provided for each atom to move individually but that leads to limited movement on a grid.

Solutions:
----------
  1. By searching for all subgraphs of a bonded set of atoms and allowing each subgraph to move as a whole we can achieve much better flexibility. This is implemented to investigate but too expensive, since there are many subgraphs for large molecules.
  2. By moving whole blocks of space around we can get some flexibility. We call this MPEG physics since it is similar to the motion compensation used in video compression. Molecules tend to clump together under this scheme.
  3. Another possibility is a combination of the two above ideas. We want to partition each molecule into sections that are likely to be able to move independently. We know that the molecules are planar graphs with edges only between neighboring lattice sites. Thus we can use MPEG-like blocks but only consider one molecule at a time to avoid the clumping.
