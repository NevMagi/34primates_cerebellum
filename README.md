# 34primates-cerebellum
Phylogenetic comparative analyses of primate cerebellar and ansiform area volumes.

In the current project we investigated the evolutionary dynamics of cerebellar and ansiform area (crura I-II) in a 34-primate dataset.

All custom code used in this project can be accessed in this repository. We used R version 4.1.0 (2021-05-18) for all our analyses. 
We also upload input files, that we obtained from the BrainBox API through Python Colab scripts. Many transformed dataframes can be found and may be used for own analyses.
Lastly, all outputs including (intermediate) figures, text files, and tables with full outcomes are uploaded. 
Final figures for publication were created in Inkscape.


First, place all input files in your home directory, and run script 0.cleanInput to get volumetric tables in the correct format. 
This step is necessary to proceed with all the other scripts, which may otherwise have input files missing.
Otherwise, all scripts should work after simply changing the working directory to your own.


Scripts are quite richly annotated, and should facilitate understanding of the meaning of each file in this folder.
For any questions about the code, or troubleshooting, you can contact the authors directly (see the publication at: doi...)


In this main folder you will find:
- Scripts (scr)
- Input: all input files (used and not used in the current project).
- Output: intermediate and final files for all analyses (foldered and number 0 through 7).
- CERES: volumetry reports for the 10 human brains in our study.
- Colab: contains the Python Colaboratoy script that was used to retrieve volumes from the BrainBox API.
- PublicationFigures: figures (including Inkscape files) for the publication main text and supplementary data.
- RevisionCode: A folder similar to this main folder, with input data, scripts, and output. These data all correspond to changes made based on Reviewer's comments.

![alt text](https://github.com/NevMagi/34primates-evo-cerebellum/blob/main/Hamadryas.png?raw=true)


