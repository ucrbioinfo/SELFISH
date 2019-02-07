# SELFISH
SELFISH (Discovery of Differential Chromatin Interactions via a Self-Similarity Measure) is a tool for finding differential chromatin interactions between two Hi-C contact maps. It uses self-similarity to model interactions in a robust way. For more information read the full paper: <a href="https://www.biorxiv.org/content/10.1101/540708v1?rss=1" target="_blank">**Selfish: Discovery of Differential Chromatin Interactions via a Self-Similarity Measure**</a>. 

<img src="https://raw.githubusercontent.com/ucrbioinfo/SELFISH/master/ESvsNPC_final.png" width="500" height="500">

<h3>Dependency </h3>

- MATLAB 8.5 and newer.

<h3>Syntax</h3>

- ```[X,Y,P] = SELFISH_DCI(contc1,contc2,norm1,norm2,THRESHOLD,RESOLUTION,INTERVAL)``` returns the coordinates of differential chromatin interactions [X,Y] and their correponding corrected p-values in P for raw contact maps read from files contc1 and contc2 with normalization vectors read from files norm1 and norm2. 

- ```[X,Y,P] = SELFISH_DCI(contc1,contc2,[],[],THRESHOLD,RESOLUTION,INTERVAL)``` returns the coordinates of differential chromatin interactions [X,Y] and their correponding corrected p-values in P for contact maps read from files contc1 and contc2. 

<h3>Parameters</h3>
SELFISH_DCI doesn't take any parameters as input but it has two predefined parameters work best for 5kb resoultion.

- `sigma0 = 1.6` first Gaussian filter sigma (r0 = 2*ceil(2*sigma0)+1).

- `s = 10` the number of Gaussian radii used in the modeling.
<h3>Inputs</h3>

- `contc1`              Hi-C contact map filename 1

- `norm1`               normalization vector filename 1

- `contc2`              Hi-C contact map filename 2

- `norm2`               normalization vector filename 2

- `THRESHOLD`           p-value threshold used for calling DCIs (e.g. 10^-4)

- `RESOLUTION`          data resoultion in bp (e.g. 5000)

- `INTERVAL`            the interval in bp for which DCIs are detected.
