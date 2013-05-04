pimaimSupercellGenerator
========================

Generator of a supercell for PIMAIM. Creates the file restart.dat.

Compilation
-----------

Install SCONS on your computer. Then type:
$ scons
Everything will be done automatically. SCONS is not yet very well known but it's an improved configure+make.

Exemple of Input file named "in"
--------------------------------

3         Number of species in the liquid, NaF ZrF4 and YF3
0.5       Molecular fraction of molecule 1 (ZrF4)
0.5       Molecular fraction of molecule 2 (NaF)
0.0	  Molecular fraction of molecule 3 (YF3)
5
F
F
F
F
Zr
2
F
Na
4
F
F
F
Y
43.0
