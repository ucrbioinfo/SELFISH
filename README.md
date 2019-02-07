# SELFISH
SELFISH: Discovery of Differential Chromatin Interactions via a Self-Similarity Measure

<h3>Dependency </h3>
MATLAB 8.5 and newer.
<h3>Syntax</h3>
`orthoMCL.txt`
 `[X,Y,P] = SELFISH_DCI(contc1,contc2,norm1,norm2,THRESHOLD,RESOLUTION,INTERVAL)` returns the coordinates of differential chromatin interactions [X,Y] and their correponding corrected p-values in P for raw contact maps read from files contc1 and contc2 with normalization vectors read from files norm1 and norm2. 
 `[X,Y,P] = SELFISH_DCI(contc1,contc2,[],[],THRESHOLD,RESOLUTION,INTERVAL)` returns the coordinates of differential chromatin interactions [X,Y] and their correponding corrected p-values in P for contact maps read from files contc1 and contc2. 

<h3>Parameters</h3>
SELFISH_DCI doesn't take any parameters as input but it has two predefined parameters work best for 5kb resoultion.
1. `sigma0` first Gaussian filter sigma ($r_0 =  2*ceil(2*sigma0)+1$).
2. `s` the number of Gaussian radii used in the modeling.
<h3>Inputs</h3>
`contc1`           -   Hi-C contact map filename 1
`norm1`            -   normalization vector filename 1
`contc2`           -   Hi-C contact map filename 2
`norm2`            -   normalization vector filename 2
`THRESHOLD`        -   THRESHOLD at which DCIs are return (e.g. 10^-4)
`RESOLUTION`       -   Data resoultion in bp (e.g. 5000)
`INTERVAL`         -   The interval in bp for which DCIs are detected.
