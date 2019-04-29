# binless

## Resolution-independent normalization of Hi-C data

Refer to the [CHANGELOG](CHANGELOG.md) for the latest changes

### Installation

#### The easy way
Let `devtools` do the installation for you. Type the following in an R shell

```
install.packages("devtools")
devtools::install_github("3DGenomes/binless",subdir="binless")
```

Installation should take about 10 minutes.

#### Manual installation
You can also install it manually as follows:

* Download the latest release [here](../../releases/latest)
* Unpack it in a clean folder on your machine, and hit `R CMD INSTALL --preclean binless`
  in a shell.
* If it complains that some packages are not installed, you must install them in
  R using `install.packages`

#### Requirements
The minimum required R version is `3.2.0`
Binless uses the following packages: `data.table`, `Hmisc`, `foreach`,
`doParallel`, `MASS`, `matrixStats`, `ggplot2`, `dplyr`, `Matrix`, `quadprog`,
`scales`, `utils`

Also, you will need the GNU Scientific Library (GSL) available on your computer.
If you use a mac and homebrew, `brew install gsl` should be enough.

Binless has been developed and tested on a MacBook Pro (2015) and on a CentOS 7 
linux workstation with 128Gb RAM and 32 cores. Resource usage can go from modest 
(fast binless will run on a laptop for loci <10Mb) to huge (fast binless on
human chromosome 1 at 5kb base resolution requires about 500Gb of RAM).

### How does it work?

In the `example/` folder, we provide plots and files to perform a normalization,
taken from publicly available data (Rao *et al.*, 2014). Alternatively you can
use your own data.  Start with something not too large, for example 1Mb. If you
want a quick and dirty overview, skip to the *Fast binless* section. Otherwise,
read on.

### Preprocessing

Try out the file `preprocessing.R`, executing each line step by step (e.g.
using rstudio). *Binless* takes mapped read pairs as input, in
[TADbit](https://3dgenomes.github.io/TADbit/index.html) `.tsv` format (see below
for file format descriptions). *Binless* normalization is performed on multiple
datasets simultaneously. All `.tsv` files are pre-processed and put into a
single `CSnorm` object that will hold all relevant information you need. At this
point, you will also be able to visualize your data using our **base-resolution
(arrows) representation**.

### Optimized binless

See `optimized_binless.R`. The resolution-independent normalization is performed
on the `CSnorm` object you built at the previous step. Once normalized, datasets
can be combined, and signal and difference detection can be performed.  **This
is the full-blown version of the algorithm, with statistically
significant output**. Note that this is a beta version, so check for updates
regularly.

### Fast binless

See `fast_binless.R`. Here, we implemented a fast approximation with fixed
fusion penalty and approximate decay and bias terms. You can either use a `CSnorm` object
produced at the preprocessing stage, or directly provide the binned raw matrix.
**This is a fast and approximate version of the full algorithm, so you will not
get statistically significant output**. But you can try out a whole chromosome ;)

### Chromosome-wide normalization

Binless is meant for locus-specific normalization. If however you need
binless on a whole chromosome, the script `chromosome_binless.R` can be used to
combine optimized and fast binless in order to get a good first guess. See the
paper for the procedure used to estimate the required parameters.

### Base-resolution (arrow) plots

If you just want to see how your data looks like with the arrow plots and are
not interested in binless normalization, check out the short tutorial
`arrow_plot.R`.

### Questions? Problems?

Check out the [FAQ](https://github.com/3DGenomes/binless/issues?q=label%3AFAQ)
in the GitHub issues tracker. If you don't find anything, just file an issue by
clicking on the "New issue" button. You can use it to report bugs, ask questions
or suggest new features. We are looking forward to your feedback!

### Citation

If you use this software, please acknowledge the following paper

Spill YG, Castillo D, Vidal E, Marti-Renom MA, "Binless normalization of Hi-C data
provides significant interaction and difference detection independent of resolution",
Nature Communications, 10:1938 (2019) [doi:10.1038/s41467-019-09907-2](https://rdcu.be/byHAs) 

### File format description

[TADbit](https://3dgenomes.github.io/TADbit/index.html) `.tsv`: tab-separated
text file containing paired-end mapped reads for a single experiment. All lines
starting with # will be discarded. The order of the reads is unimportant. It has
the following columns
1. `id`: ID of the mapped pair (alphanumeric without spaces)
1. `chr1`: chromosome of read 1 (alphanumeric without spaces)
1. `pos1`: position of the 5' end of read 1 (integer)
1. `strand1`: whether the read maps on the forward (1) or reverse (0)
     direction
1. `length1`: how many bases were mapped (integer)
1. `re.up1`: upstream restriction site closest to pos1 (integer <= pos1)
1. `re.dn1`: downstream restriction site closest to pos1 (integer > pos1)
1. `chr2`
1. `pos2`
1. `strand2`
1. `length2`
1. `re.up2`
1. `re.dn2`

Binned raw matrix (used for fast binless): tab, comma or space-separated text file
containing multiple datasets. The first line is a header that must start with
`"name" "bin1" "pos1" "bin2" "pos2" "distance" "observed" "nobs"`. Optionally, more columns
can be added but make sure their column names are different.
1. `name`: The name of the dataset. Put in "" if you use spaces
1. `bin1`: the label for the bin on the x axis. For example `"[begin1,end1)"`
1. `pos1`: The position in bases of the center of the bin on the x axis
1. `bin2`: The label for the bin on the y axis
1. `pos2`: The position on the y axis
1. `distance`: The average distance (in bases) between the two bins
1. `observed`: How many paired-end reads were mapped within (bin1,bin2)
1. `nobs`: The number of observables, i.e. four times the number of cut site intersections in that bin 

Note that `name` will be converted to R factors, so you an also
provide them as integers starting at 1 (i.e. use 1 for the first dataset, 2 for the second etc.).
Also, **you must have pos2 >= pos1, and the data must be sorted by name, pos1 and pos2, in that order**.


