This folder contains all the .exelem and .exnode files that are outputted by OpenCMISS when running a steady-state Simplex problem.
Edited so that certain elements of your choosing are removed from the cuboid (at line 87).
Note the code isn't perfectly functioning. Can't remove so many elements in one region such that no element is attached to one particular node.
For example, if elements 75, 78, 34 and 38 were all attached to node 24 (this is a random example, not real values), and all 4 elements were removed, then an error would result.
