# Bump mesh generation

1. `generate_mesh.m` : get the nodes and triangle points for the first unstructured mesh. This uses the distmesh package by Per-Olof Persson. Elements are clusters around the bump.
2. `gri_file.py` : generate .gri file. Containing nodes, triangles, and boundaries information.
3. `process_grid.py` : Process the .gri file and generates the (text files)matrices:
    - I2E: mapping from interior faces to elements. Contains 4 integers `elemL`, `faceL`, `elemR`, `faceR`. Containing the elements adjacent to the face and the local numbering of their corresponding face. The left (L) element is arbitrarily chosen as the one with the lower element index.
    - B2E: Mapping from boundary faces to elements and boundary group. Contains three integers `elem`, `face`, `bgroup`.
    - In: Normal vectors for interior faces.
    - Bn: Normal vector for boundary faces.
    - Area: Area of each element.
  The matrices are all 1-based. This file also does a sanity check:
    $$\sum\limits_{i=1}^3 l_i\vec{n}_{i,out}=0$$
4. `refine.py` : File to do uniform mesh refinement, where element is split into four subelements by means of edge bisection. This generates a new set of elements of nodes that can be processed with `gri_file.m` and `process_grid.py`. 
