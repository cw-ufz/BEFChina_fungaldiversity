Journal publication
"Experimental evidence of functional group-dependent effects of tree diversity on soil fungi in subtropical forests"

Christina Weißbecker1,2*, Tesfaye Wubet1,3*, Guillaume Lentendu4, Peter Kühn5, Thomas Scholten5, Helge Bruelheide6,3, François Buscot1,3 
1UFZ - Helmholtz Centre for Environmental Research, Department of Soil Ecology, Theodor-Lieser-Straße 4, 06120 Halle (Saale), Germany; 
2University of Leipzig, Institute of Biology, Johannis-Allee 3-5, 04103 Leipzig, Germany; 
3German Centre for Integrative Biodiversity Research (iDiv) Jena – Halle – Leipzig, Deutscher Platz 5e, 04103 Leipzig, Germany; 
4University of Kaiserslautern, Erwin-Schroedinger Straße, 67663 Kaiserslautern, Germany, 
5Eberhard Karls University Tübingen, Chair of Soil Science and Geomorphology, Rümelinstraße 19-23, 72070 Tübingen, Germany; 
6Martin Luther University Halle-Wittenberg, Institute of Biology / Geobotany and Botanical Garden, Am Kirchtor 1, 06108 Halle (Saale), Germany

Author for correspondence: Christina Weißbecker, Email: christina.weissbecker@ufz.de

Abstract
Deconvoluting relative contributions of specific biotic and abiotic drivers on soil fungal community compositions facilitates predictions about 
ecosystems’ functional responses to environmental changes, such as losses of plant diversity, but is hindered by the complex interactions involved.
Experimental assembly of tree species allows separation of the respective effects of plant community composition (biotic components) and soil 
conditions (abiotic components) enabling a much greater statistical power than can be achieved in observational studies. Therefore, we analyzed 
these contributions by assessing fungal communities via pyrotag sequencing of the internal transcribed spacer (ITS2) rDNA region in young subtropical 
forest plots included in a large experiment on effects of tree species richness. Spatial and soil-related factors were the main drivers of soil
 fungal alpha and beta-diversity which implies strong early-stage environmental filtering and dispersal limitation. Tree-related variables, such as
 tree community composition, significantly affected arbuscular mycorrhizal fungal and pathogen community structure, while differences in tree host 
species and host abundance affected ectomycorrhizal fungal community composition. Due to persisting legacy effects and the limited size of tree saplings, 
providing low carbon inputs to the ecosystem by rhizodeposits and leaf litter at this early stage of the experiment, we expect to find even stronger tree 
related effects on fungal community composition with proceeding forest development. 

In this README you find the related information where to find the associated data:
	1) Fungal ITS2 DNA amplicon 454 raw sequences
	2) Bioinformatic workflow description
	3) Metadata and processed sequence data
	4) R analyses scripts and in-between results
	5) Publication

1) Fungal ITS2 DNA amplicon 454 raw sequences

Raw sequences and some associated metadata was deposited at the European Nucleotide Archive under study accession number PRJEB12020.

2) Bioinformatic workflow description

configuration.BEF_vipA_454_ITS_1283866816.tsv			# Initial setting for in-house bioinformatic program pipeline amp_pipe version 1.5
BEF_vipA_454_ITS_1283866816.documentation.txt			# Overview of bioinformtatic analysis (also some previous settings for a former draft), summary stats are unfortunately incomplete
BEF_vipA_454_ITS_1283866816.read_counts.tsv			# Overview read counts per sample at each analysis step

3) Metadata and processed sequence data

This part of the data is found at Zenodo repository (all versions: DOI: 10.5281/zenodo.1215504) (https://zenodo.org/record/1215505#.W30mD6LsQis)
which is linked to this GitHub repository (all versions: DOI 10.5281/zenodo.1215486) (https://zenodo.org/record/1215487#.W30nlqLsQis)

BEF_vipA_454_ITS_1283866816.json_Nov2017.biom			# Final bioinformatic output
BEF_vipA_454_ITS_2011_curated_taxonomy_guild.csv		# Manually refined taxonomic and FunGuild annotation

Further metadata:
	BEF_vipA_2011_local_9_trees.csv				# tree species identity of sampling tree and surrounding eight tree neighbors
	BEF_vipA_2011_soil_typo_metadata.csv			# soil and topographic measurments investigated for the individual samples
	BEF_vipA_2011_sample_location.csv			# GPS coordinates for each sampling location
	BEF_vipA_2011_Tree_mycotype_composition.csv		# Mycorrhizal Type of each sampling tree individual
	BEF_vipA_2011_Sampling_tree.csv				# Sampling tree species identity

4) R analysis scripts

There are submitted a total of 13 numbered analysis scripts and 3 associated R functions.


5) Publication

Publication can be found here: 





