# Discontinous-Galerkin solver for Euler equations

To compile the solver please run: g++ *.cpp -o main.exe -O2 -larmadillo
The armadillo package must be installes in order to compute some matrix operations (inverse of the mass matrix)

DG.h defines all the functions used by the solver 
flux_BC.cpp includes the flux(Fx,Fy,and Roe) and the boundary conditinos 
quad.cpp includes the basis functions, 1d/2d quad points, Jacobian calculations, Mass Matrix for cuved elements
readgrid.cpp is the function used to read the curved mesh. It just works for this particular case.
residuals.cpp includes 4 functions to compute the interior and edge residual contributions
solver.cpp put everything togheter and run the time marching scheme RK4
main.cpp call the solver and there you can set the mesh, p, tolerance, and CFL
The python files are for post-processing

Inside the folder Meshes you can find the meshes used and the code to curve them.

## Bump mesh generation

Files to generate the meshes in `gen_mesh` 

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


