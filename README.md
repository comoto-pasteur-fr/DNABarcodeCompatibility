

About
=================

DNABarcodeCompatibility: R-package to find least-redundant sets of compatible barcodes for multiplex experiments.



Installation 
================

* Requirements
    + Install [R](https://www.r-project.org/) if not yet installed (R >= 3.2 is required).


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

In addition, DNABarcodeCompatibility can be used from the DNABarcodeCompatibility [graphical user interface](https://github.com/comoto-pasteur-fr/DNABarcodeCompatibility_GUI). In such a case, additional dependencies must be installed: [Java (JDK 8)](http://www.oracle.com/technetwork/java/javase/downloads/jdk8-downloads-2133151.html) and `rJava` R-package (`install.packages("rJava")`). 


Documentation
================

[Introduction](https://comoto-pasteur-fr.github.io/DNABarcodeCompatibility/)

[API documentation](https://comoto-pasteur-fr.github.io/DNABarcodeCompatibility/DNABarcodeCompatibility-manual.pdf)

Related tool: see the [graphical user interface (GUI)](https://github.com/comoto-pasteur-fr/DNABarcodeCompatibility_GUI) of DNABarcodeCompatibility.

Support
=========

Please use the [github ticket system](https://github.com/comoto-pasteur-fr/DNABarcodeCompatibility/issues) to report issues or suggestions. 
We also welcome pull requests.



Reference
==========

Céline Trébeau, Jacques Boutet de Monvel, Fabienne Wong Jun Tai, Raphaël Etournay. (2018, May 31). comoto-pasteur-fr/DNABarcodeCompatibility: First complete release (Version v0.0.0.9000). Zenodo. [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1256863.svg)](https://doi.org/10.5281/zenodo.1256863)



