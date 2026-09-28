# 2EVQ-derived beta-hairpin example

`angles.txt` contains 12 backbone \(\phi,\psi,\omega\) triplets derived from Model 1 of PDB entry 2EVQ.

The first-residue \(\phi\) angle is undefined for an isolated peptide and is represented by a reference value for coordinate construction. At the C terminus, the values required by the recursive construction are represented using the terminal geometry used for this example.

The generated PDB contains the backbone atoms N, CA, C, and O. The current program writes generic `ALA` residue labels, so the file represents the generated backbone geometry rather than the complete chemical identity of 2EVQ.
