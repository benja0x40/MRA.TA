MRA.TA - Multi-Resolution Analysis of Tiling Arrays
================================================================================

MRA.TA is an R package for multi-resolution representation and segmentation of
genomic profiles from tiling arrays such as ChIP-on-chip and Chromosome
Conformation Capture on Chip (4C) experimental data.

### Main features ###

#### Domaingram representation ####

#### Multi-resolution segmentation ####

### Installation ###

#### Prerequisites ####

  - R environment version 3.x
  - R packages: `devtools`, `stringr`, `getopt`, `plotrix`
  - [Bioconductor](http://www.bioconductor.org/) packages: `Biostrings`, `GenomicRanges`
  
#### Installing MRA.TA ####

```R
library("devtools")
install_github("benja0x40/MRA.TA")
```

### References ###

<a name="1"></a>1. Leblanc B., Comet I., Bantignies F., and Cavalli G., *Chromosome Conformation Capture on Chip (4C): data processing.* Book chapter to appear in *Polycomb Group Proteins.* Lanzuolo C., Bodega B. editors, Methods in Molecular Biology (2016).  
[publisher](https://www.springer.com/gp/book/9781493963782)

<a name="2"></a>2. de Wit E., Braunschweig U., Greil F., Bussemaker H. J. & van Steensel B. *Global chromatin domain organization of the Drosophila genome.* PLoS Genetics 4, e1000045 (2008).  
[publisher](http://dx.doi.org/10.1371/journal.pgen.1000045) | [pubmed](https://www.ncbi.nlm.nih.gov/pubmed/18369463)
