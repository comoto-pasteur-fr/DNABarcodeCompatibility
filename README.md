

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


Documentation
================

[Introduction](https://comoto-pasteur-fr.github.io/DNABarcodeCompatibility/)

[API documentation](https://comoto-pasteur-fr.github.io/DNABarcodeCompatibility/DNABarcodeCompatibility-manual.pdf)

* Related tool:

In addition, DNABarcodeCompatibility can be used from a graphical user interface. We now provide two alternative interfaces:

1) The [Java based graphical user interface](https://github.com/comoto-pasteur-fr/DNABarcodeCompatibility_GUI). In such a case, additional dependencies must be installed: [Java (JDK 8 - 64 bits)](http://www.oracle.com/technetwork/java/javase/downloads/jdk8-downloads-2133151.html) and `rJava` R-package (`install.packages("rJava")`).

2) The [Shiny web-based graphical user interface](https://github.com/comoto-pasteur-fr/DNABarcodeCompatibility_Shiny). This interface can be run locally within a web browser or deployed on a web server.



Support
=========

Please use the [github ticket system](https://github.com/comoto-pasteur-fr/DNABarcodeCompatibility/issues) to report issues or suggestions. 
We also welcome pull requests.



Reference
==========

Céline Trébeau, Jacques Boutet de Monvel, Fabienne Wong Jun Tai, Raphaël Etournay. (2018, May 31). comoto-pasteur-fr/DNABarcodeCompatibility: First complete release (Version v0.0.0.9000). Zenodo. [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1256863.svg)](https://doi.org/10.5281/zenodo.1256863)



