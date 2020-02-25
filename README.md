# Grid-Tools
Utilities for handling block-structured grid with conformal interfaces.

## Mesh I/O
Read and write some widely-used mesh files.  
Currently supported formats:
> * PLOT3D: *.__fmt__ or  *.__xyz__ 
> * FLUENT: *.__msh__

It aims to be a self-contained toolkit with operations that are easy to use.  
This utility is typically designed for a 3D CFD solver.  

## Block-Glue
Given block connectivity information, it converts multi-block structured grid into unstructured format.  
The block connectivity info is written in a `Neutral Map File`, where topology structure can be identified.  
The cartesian coordinates are stored in a `PLOT3D` file, whose format is classical and easy to understand. It should be noted that the "`IBLANK`" info within a PLOT3D grid will __NOT__ be used.  
In short, it functions as __PLOT3D + NMF -> FLUENT__.  
This utility is typically designed for optimization.  
