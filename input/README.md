This folder contains all files used in analysis.
These files are necessary to perform code in this repository. 
Put them in your working directory and make sure to specify the path at the beginning of each script.

Volumetric data come from Heuer et al. (2019); https://github.com/neuroanatomy/34primates/blob/master/data/derived/stats/stats.csv and by obtaining volumes from the current study from the BrainBox API through Python Colab scripts. 
The BrainBox API is only accesible from within the project. 
Contact the authors to be added to the BrainBox project.

Phylogenetic trees were obtained from the 10kTrees website (Arnold et al. (2010), https://onlinelibrary.wiley.com/doi/full/10.1002/evan.20251). 
Volumes for cerebellar and cerebral volume from the Stephan dataset (Stephan et al. (1981); https://pubmed.ncbi.nlm.nih.gov/7014398/) were manually copied into a csv file. 
Species look-up tables were obtained from Heuer et al. (2019); https://github.com/neuroanatomy/34primates/tree/master/src)

The folder structure is somewhat messy due to data wrangling activities, but we decided to keep all files in becasue they are rather small and may be useful for your own research purposes.
Only script 0.cleaningInput.R and 1.raw_stats.R should require input from this folder.
If you want to minimally reproduce our results, one should only keep the files that these scripts need.

We made a short-list for the necessary data with descriptions (please check manually if this list is complete):

Input data description (alphabetical):

- 10kTrees_34PrimateSpecies_adapted.csv: Slightly adapted (from Heuer et al. (2019)) look-up table for appending species information such as Wilson & Reeder, English, and Wikipedia names, as well as clade and ape membership to the correct specimen IDs.
- bodyweight-isleretal2008.csv: bodyweight place-holder values obtained from Isler et al. (2008), which originally come largely from Stephan et al. (1981). We decided against using these data, as described in the main text.
- catn.R: printing function for phylogenetic analyses.
- consensusTree_10kTrees_Primates_Version3_13.nex: consensus phylogenetic tree for the 13 species for which the current study obtained complete data (including ansiform area (Crura I-II) volumes).
- consensusTree_10kTrees_Primates_Version3_15.nex: consensus phylogenetic tree for the 15 species for which both the current study and the Stephan collection (Stephan et al. (1981)) reported cerebellar and cerebral volumes.
- consensusTree_10kTrees_Primates_Version3_34.nex: consensus phylogenetic tree for the 34 species included in the current study.
- dfCerebellum: cerebellar volumetric data obtained from the BrainBox API.
- dfCrusMahta: Crura I+II volumetric data, as segmented by Mahta Abbaspour, obtained from the BrainBox API.
- dfCrusNeville: Crura I+II volumetric data, as segmented by Neville Magielse, obtained from the BrainBox API.
- dfCrusVanessa: Crura I+II volumetric data, as segmented by Vanessa Steigauf, obtained from the BrainBox API.
- manual_automated_segmentations.csv: contains both manual segmentations and segmentations from CERES (reports from CERES are uploaded to a separate folder in this GitHub repository) added manually.
- manual_automated_segmentations_withMEAN.numbers: convenient file to read the mean values and compare them. Not used in analyses.
- stats2.csv: Neueroanatomical measurements obtained in Heuer et al. (2019).

As an important note, data from the Stephan collection have not been published open access.
Therefore, we are not able to share them here. 
You can still access our scripts, but will have to obtain the volumetric data from the Stephan et al. 1981 (DOI: 10.1159/000155963) yourselves.
Additionally, we did not want to give away the species in the database without explicit consent of the authors/publishers, so the phylogenetic trees have to also be obtained from 10kTrees before being able to use our 7.robustness script.

