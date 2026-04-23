=======================================
## SYGQ5 Research project 1 
# Reconstructing the historical species ranges of Madagascar's birds
=======================================

Using the ~1,100 undigitized specimens from the 1929-1931 Archbold-Vernay Expedition to Madagascar,
this project looks to complate the following aims:

1. To further develop and test an AI-assisted workflow to extract location information from the NHM-AVE historical museum
   specimens. 
2. To illustrate how extracted historic location data from the AVE collection can be used to generate historical species
   ranges. 
3. To explore how sampling effort from historical expeditions might be quantified and what impact this might have on
   resulting modelled species distributions.

   Method:
   1. Data collection
   2. Specimen label imaging
   3. LLM model testing and fuzzy testing
   4. API-LLM data extraction
   5. Geo-coding specimen localities
   6. Post-processing geo-coded specimens
   7. Species distribution modelling
  
   Listed below are the code file names in this repository, Research Project 1, with the method step in which they are
   utilised (as per the project report).

   Specimen label cropping and mergimg.ipynp - Step 3.
   API-LLM data extraction.ipynb - Step 3. and 4.
   Fuzzy-matching test.ipynb - Step 3
   SDM predictor layers.R - Step 7
   Final_maxent_workflow.R - Step 7
   One-way_ANOVA_test_code.R - Step 7
