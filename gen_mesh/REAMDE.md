# Bump mesh generation

1. generate_mesh.m : get the nodes and triangle points for the first unstructured mesh. This uses the distmesh package by Per-Olof Persson. Elements are clusters around the bump.
2. gri_file.py : generate .gri file. Containing nodes, triangles, and boundaries information.
3. process_grid.py : Process the .gri file and generates the (text files)matrices:
    - I2E: mapping from interior faces to elements.
    - B2E: Mapping from boundary faces to elements and boundary group.
    - In: Normal vectors for interior faces.
    - Bn: Normal vector for boundary faces.
    - Area: Area of each element.   
  This file also does a sanity check:
    $$\sum\limits_{i=1}^3 l_i\vec{n}_{i,out}=0$$
4. 
