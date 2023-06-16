# Materials Cloud 3-dimensional structure IDs

The current format of `mc3d.yaml` is a dictionary with for each database (crystallographic open database ([COD](http://www.crystallography.net/cod/)), inorganic crystal structures database ([ICSD](https://www.fiz-karlsruhe.de/en/produkte-und-dienstleistungen/inorganic-crystal-structure-database-icsd)) and materials platform for data science ([MPDS](https://mpds.io/#start))) a list of {‘id’: …, ‘version’: …} dictionaries. 
This is the clearest and most succinct way to identify the whole dataset without providing the actual structures, which we can’t be shared for the MPDS and ICSD due to licensing constraints.
The file contains 79854 inorganic structures. 