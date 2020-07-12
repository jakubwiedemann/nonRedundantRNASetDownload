#README

##Description

The purpose of this tool is to get the current non-redundant set of RNAs or identify which structures has been added/removed since last update.

The workflow consist of:

1. Downloading non-redundant set from http://rna.bgsu.edu/rna3dhub/ and saving it in RNA_SETS folder.
2. If RNA_SETS folder contains only 1 file that means that there is no initial set and it needs to be created:
2.1. The latest file is parsed and a new file 'init_set.txt' is being created. This file contain all of the structures from the latest release.

