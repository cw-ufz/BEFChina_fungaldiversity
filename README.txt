Journal publication
"Experimental evidence of functional group-dependent effects of tree diversity on soil fungi in subtropical forests"

Christina Weißbecker1,2*, Tesfaye Wubet1,3*, Guillaume Lentendu4, Peter Kühn5, Thomas Scholten5, Helge Bruelheide6,3, François Buscot1,3 
1UFZ - Helmholtz Centre for Environmental Research, Department of Soil Ecology, Theodor-Lieser-Straße 4, 06120 Halle (Saale), Germany; 2University of Leipzig, Institute of Biology, Johannis-Allee 3-5, 04103 Leipzig, Germany; 3German Centre for Integrative Biodiversity Research (iDiv) Jena – Halle – Leipzig, Deutscher Platz 5e, 04103 Leipzig, Germany; 4University of Kaiserslautern, Erwin-Schroedinger Straße, 67663 Kaiserslautern, Germany, 5Eberhard Karls University Tübingen, Chair of Soil Science and Geomorphology, Rümelinstraße 19-23, 72070 Tübingen, Germany; 6Martin Luther University Halle-Wittenberg, Institute of Biology / Geobotany and Botanical Garden, Am Kirchtor 1, 06108 Halle (Saale), Germany

Author for correspondence: Christina Weißbecker, Email: christina.weissbecker@ufz.de

Abstract
•	Unraveling the independent relative contributions of biotic and abiotic drivers on soil fungal community composition enables predictions about functioning of ecosystems and how they react to environmental changes such as plant diversity loss. We analyzed these independent contributions in the framework of a large experiment on forest tree species richness, thus attaining a statistical power that is much larger than that of observational studies.
•	Soil fungal communities in young planted subtropical forests were assessed by pyrotag sequencing of the fungal internal transcribed spacer (ITS2) rDNA region. 
•	Spatial and soil related factors constituted the dominant drivers of soil fungal alpha and beta-diversity variation. Tree related variables such as tree beta-diversity significantly affected AM fungal and pathogen beta-diversity, while differences in tree host species and host abundance affected EcM fungal beta-diversity.
•	Effects of tree community changes and tree species loss on the soil fungal community might be undetectable considering the fungal community as a whole (as dominated by saprotrophs) but might critically impact the fungal functional subgroups of pathogens, EcM and AM fungi and related ecosystem functions.

In this README you find the related information where to find the associated data:
	1) Fungal ITS2 DNA amplicon 454 raw sequences
	2) Bioinformatic workflow description
	3) Metadata and processed sequence data
	4) R analyses scripts and inbetween results
	5) Publication

1) Fungal ITS2 DNA amplicon 454 raw sequences

Raw sequences and some associated metadata is deposited at the European Nucleotide Archive under study accession number PRJEB12020.

2) Bioinformatic workflow description

configuration.BEF_vipA_454_ITS_1283866816.tsv			# Initial setting for in-house bioinformatic program pipeline amp_pipe version 1.5
BEF_vipA_454_ITS_1283866816.documentation.txt			# Overview of analysis (also some previous settings for a former draft), summary stats are unfortunately incomplete
BEF_vipA_454_ITS_1283866816.read_counts.tsv				# Overview read counts per sample at each analysis step

3) Metadata and processed sequence data

This part of the data is found at Zenodo repository DOI: 
which is linked to this GitHub repository

BEF_vipA_454_ITS_1283866816.json_Nov2017.biom			# Final bioinformatic output
BEF_vipA_454_ITS_2011_curated_taxonomy_guild.csv		# Manually refined taxonomic and FunGuild annotation
Further metadata:
	BEF_vipA_2011_local_9_trees.csv						# tree species identity of sampling tree and surrounding eight tree neighbors
	BEF_vipA_2011_soil_typo_metadata.csv				# soil and topographic measurments investigated for the individual samples
	BEF_vipA_2011_sample_location.csv					# GPS coordinates for each sampling location
	BEF_vipA_2011_Tree_mycotype_composition.csv			# Mycorrhizal Type of each sampling tree individual
	BEF_vipA_2011_Sampling_tree.csv						# Sampling tree species identity

4) R analysis scripts

There are submitted a total of 13 numbered analysis scripts and 3 associated R functions.


5) Publication

Publication can be found here:





