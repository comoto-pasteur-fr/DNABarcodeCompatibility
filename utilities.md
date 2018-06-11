Utilities
===========

### Create and use gh-pages branch to store html file for display through git io 

Note that html files are displayed as raw code on the standard git repository.
Once the new gh-pages branch is created, go on GitHub > Settings and declare the gh-pages branch for the git io.

1. First create an orphan branch (nothing in common with master) locally:

```
git checkout --orphan gh-pages
git rm --cached -r .
git rm -rf .
```

2. Push changes to github:

```
git push origin gh-pages
```


3. Switch back to master branch locally:

```
git checkout master
```

4. Create a separate local repository folder for the gh-pages branch for clarity:

```
cd $HOME; mkdir DNABarcodeCompatibility_static; cd DNABarcodeCompatibility_static
```

```
git clone -b gh-pages --single-branch https://github.com/comoto-pasteur-fr/DNABarcodeCompatibility.git .
```

5. Solve merging conflicts:

```
# Select master branch in which to merge another "doc" branch
git checkout master

# Merger the doc branch into the master branch
git merge doc

# Solve conflicts manually using a GUI
git mergetool -t diffmerge
```

### R installation on Linux: dependencies

```
sudo apt-get install -y libcurl4-openssl-dev libssl-dev libssh2-1-dev openjdk-11-jdk
```

### Handling R library paths

* Find where local R libraries are installed (R console)

```
.libPaths()
```

or in a unix shell: 

```
echo $R_LIBS_USER
```

* Append a new library path (R console)

```
newPath="absolute/path/to/library"
.libPaths(newPath)
```

### Add the Illumina dataset to the package

* Store the datasets in the data/ folder
```
IlluminaIndexesRaw <- read.table("/Users/retourna/Index_Illumina.txt")
devtools::use_data(IlluminaIndexesRaw)
```

```
export_dataset_to_file = function(dataset) {
  if (class(dataset)=="data.frame") {
    write.table(dataset,
                textfile <- tempfile(), row.names = FALSE, col.names = FALSE, quote=FALSE)
    print("Dataset successfully exported into a file.")
    return(textfile)
  } else print("The input dataset isn't a data.frame: NOT exported into file")
}

IlluminaIndexes <- file_loading_and_checking(export_dataset_to_file(DNABarcodeCompatibility::IlluminaIndexesRaw))
devtools::use_data(IlluminaIndexes, overwrite = TRUE)
```



* Document the dataset in R/data.R
```
#' Barcode dataset from Illumina.
#'
#' 48 barcodes from Illumina TruSeq Small RNA kits
#'
#' @format A data frame with 48 rows and 2 variables:
#' \describe{
#'   \item{V1}{barcode identifier}
#'   \item{V2}{DNA sequence}
#'   ...
#' }
#'
"IlluminaIndexes"
```

