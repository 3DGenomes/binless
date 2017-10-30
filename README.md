# binless
## Resolution-independent normalization of Hi-C data

### Installation

* To get started, download this repository to a fresh folder on your computer
* Go to that folder on your machine, and hit `R CMD INSTALL --preclean binless` in a shell.
* If it complains that some packages are not installed, you must install them in R using install.packages.

We plan to put this in Bioconductor to make the installation easier in the future

### Example dataset

* Try out the file `fast_binless.R`, executing each line step by step in rstudio. In the `example/` folder, we provide files to play with, taken from publically available data (Rao *et al.*, 2014).

* Alternatively you can use your own data. Start with something not too large, for example 10Mb.

* **This is a fast and approximate version of the full algorithm, so you will not get statistically significant output**, and it might not look as *smooth* as the full-blown algorithm. We will release it soon, sit tight.
