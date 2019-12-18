# Grid-Tool
Utilities for handling block-structured grid with conformal interfaces.

## Mesh I/O
Read and write some widely-used mesh files(e.g. *.__fmt__, *.__msh__).  
It's a self-contained toolkit with operations that are easy to use.  
This utility is typically designed for a 3D CFD solver.  

## Glue
Given block connectivity information, it converts multi-block structured grid into unstructured format.  
Especially, it functions as __Plot3D + NMF -> Fluent Mesh__.  
This utility is typically designed for optimization.  