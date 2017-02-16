The .in files stored in here are read into the python files stored in IdealisedMeshes.
Note that separate python files had to be created for each mesh because of the format of the mesh.
Reading the node blocks corresponding to the inlet or outlet of each mesh often differed because they would start with different titles, or end on different lines.
For example, the line indicating the start of the nodes could either read 'CMBLOCK', 'MOUTH' OR 'CMBLOCK', 'INLET'. Same goes for the end of the node block, which could end in 'CMBLOCK' for the next section, or 'GOLIST\n'.
