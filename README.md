<img src="docs/sbsclust_logo.jpg" alt="logo" title="Secondary structure" height="400" align="right" />

# SBSClust

SBSClust is a computational tool designed to detect clusters of mutations in complete genomes. 

## Installation

To install SBSClust from GitHub navigate to a directory where you want to install the program and clone the repository:
```
git clone https://github.com/KazanovLab/SBSClust
```

To compile SBSClust and install it system-wide navigate into the cloned directory and run:
```
cd SBSClust
make
sudo make install
```

## Quick start

SBSClust can generate output in short, long, or both formats simultaneously.
* The short output (-s option) includes only mutations that are part of the identified clusters.
* The full output (-f option) contains all mutations from the input VCF files.

Other options include:
* The input directory containing VCF file (specified after the -s/-f options)
* The reference genome location (-g).
* The output directory (-o).
* The p-value for determining the between-mutation distance threshold (-t, default: 0.01).

