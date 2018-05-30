Utilities
===========

Create and use gh-pages branch to store html file for display through git io.

1. First create an orphan branch (nothing in common with master) locally:

```
git checkout --orphan gh-pages
git rm --cached -r .
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