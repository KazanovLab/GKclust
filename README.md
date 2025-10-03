<img src="docs/genaclust_logo.jpg" alt="logo" title="GENAclust logo" height="400" align="right" />

# GENAclust

GENAclust is a computational tool designed to detect clusters of mutations in complete genomes. 

## Installation

To install GENAclust from GitHub navigate to a directory where you want to install the program and clone the repository:
```
git clone https://github.com/KazanovLab/GENAclust
```

To compile GENAclust and install it system-wide navigate into the cloned directory and run:
```
cd GENAclust
make
sudo make install
```

## Quick start

GENAclust can generate output in several formats simultaneously.
* The short output (-s option) includes only mutations that are part of the identified clusters.
* The full output (-f option) contains all mutations from the input VCF files.
* The list of cluster (-m option) contains only list of clusters for particular VCF file.

Other options include:
* The input directory containing VCF file (-i) or the list of full paths to VCF files (-l) in a separate text file.
* The reference genome location (-g).
* The output directory (-o).
* The p-value for determining the between-mutation distance threshold (-t, default: 0.01).

Examples:

Short format output:
```
genaclust -i /inputdir/ -o /outputdir/ -g /refgenome/hg19.fa
```

Full format output:
```
genaclust -f /inputdir/ -o /outputdir/ -g /refgenome/hg19.fa
```


## Reporting Bugs and Feature Requests
Please use the [GitHub issue tracker](https://github.com/KazanovLab/SBSClust/issues) to report bugs or suggest features.

## Citing
Ponomarev GV, Kurtoğlu B, Kazanov MD. "GENAclust: A Bioinformatics Tool for Detecting Spatially Clustered Mutations in Complete Genomes". *to be submitted*

## Funding
This study was supported by the Scientific and Technological Research Council of Turkey (TUBITAK) under Grant Number 123E476. The authors thank TUBITAK for their support. 

