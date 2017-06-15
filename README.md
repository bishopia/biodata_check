# biodata_check
biodata_check() is a function that checks user taxonomy against the BioData Taxonomic System

Some colleagues have recently expressed an interest in being able to quickly compare lists of 
taxonomic names that they use in the their own work with USGS BioData Taxonomic System. I 
agree that this would be a useful tool and that it would help promote taxonomic consistency 
in North America. I hope you find this useful and would appreciate feedback if one is inclined 
to do so.

This repository concerns the "biodata_check()" function. In addition to the code itself, I have 
provided an R Markdown script with notes on proper data formatting, required R packages, and 
several examples of its use.

So far, this function can be used in two ways. First, a user can input a list of names 
and the function will return that list and the Biodata Taxon Name thats correspond to it. Second, a 
user can input a species abundance dataset (site by species matrix) and the function will return 
an updated matrix with updated nomenclature based on corresponding Biodata Taxon Names. Newly formed 
duplicates columns will be merged together.


-ian
