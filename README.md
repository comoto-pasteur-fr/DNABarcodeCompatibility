

About
=================

DNABarcodeCompatibility: R-package to find least-redundant sets of compatible barcodes for multiplex experiments performed on Illumina sequencing platforms.



Installation 
================

* First install [R](https://www.r-project.org/) if not yet installed.

* Within a R console, type in the following commands:
    
```
# Install the devtools package
install.packages("devtools")

# Enable both CRAN and Bioconductor repositories
setRepositories(ind=1:2)

# Install DNABarcodeCompatibility
devtools::install_github("comoto-pasteur-fr/DNABarcodeCompatibility")
```
* Note:

DNABarcodeCompatibility uses the rJava package and therefore requires [java (JDK)](http://www.oracle.com/technetwork/java/javase/downloads/jdk8-downloads-2133151.html) to be installed on your system.

Documentation
================

[Introduction](https://comoto-pasteur-fr.github.io/DNABarcodeCompatibility/)

[API documentation](https://comoto-pasteur-fr.github.io/DNABarcodeCompatibility/DNABarcodeCompatibility-manual.pdf)

Support
=========

Please use the [github ticket system](https://github.com/comoto-pasteur-fr/DNABarcodeCompatibility/issues) to report issues or suggestions. 
We also welcome pull requests.



Reference
==========

Céline Trébeau, Jacques Boutet de Monvel, Virginie Wong Jun Tai, Raphaël Etournay. (2018, May 31). comoto-pasteur-fr/DNABarcodeCompatibility: First stable release (Version v1.0.0). Zenodo. [http://doi.org/10.5281/zenodo.1256863](http://doi.org/10.5281/zenodo.1256863)



