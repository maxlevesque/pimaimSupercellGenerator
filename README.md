pimaimSupercellGenerator
========================

pimaimSupercellGenerator generates simulation supercells for the molecular dynamic package PIMAIM.
It generates the input file *restart.dat*.


Compilation
-----------

First, you need **SCONS** on your computer. SCONS is not yet very well known, but let's say it is an improved autotools (configure+make). SCONS can be found in the repositories of almost every linux distribution. Feedback about macOS is welcome.

Once SCONS is installed on your computer, start a new Terminal emulator, then:
```
$ scons
```
The files SConstruct and SConscript that can be found in the pimaimSupercellGenerator package are used by SCONS.


Exemple of Input file named "in"
--------------------------------
```
3     Number of species in the liquid, NaF ZrF4 and YF3
0.5   Molecular fraction of molecule 1 (ZrF4)
0.5   Molecular fraction of molecule 2 (NaF)
0.0   Molecular fraction of molecule 3 (YF3)
5     There are 5 atoms in ZrF4, whose symbols are
F
F
F
F
Zr
2     There are 2 atoms in NaF
F
Na
4     There are 4 atoms in YF3
F
F
F
Y
43.0  The generated cubic supercell will be 43 Bohr wide.
6     Number of molecule in each direction of the supercell.
```


Future optimizations
--------------------
[] The supercell can only be cubic. The number of molecules in the supercell can only be around i^3, with i a positive integer.
