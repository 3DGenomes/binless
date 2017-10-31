# binless
## Resolution-independent normalization of Hi-C data

### Installation

* Download this repository to a fresh folder on your computer
* Go to that folder on your machine, and hit `R CMD INSTALL --preclean binless`
  in a shell.
* If it complains that some packages are not installed, you must install them in
  R using `install.packages`

We plan to put this in Bioconductor to make the installation easier in the
future.

### How does it work?

In the `example/` folder, we provide files to play with, taken from publicly
available data (Rao *et al.*, 2014). Alternatively you can use your own data.
Start with something not too large, for example 10Mb. If you want to try it out
right away, skip to the *Fast binless* section. Otherwise, read on.


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
is the full-blown version of the algorithm, so you will get statistically
significant output**. Note that this is a beta version, so check for updates
frequently.

### Fast binless

See `fast_binless.R`. Here, we implemented a fast approximation with fixed
fusion penalty and an approximate decay. You can either use a `CSnorm` object
produced at the preprocessing stage, or directly provide the binned raw matrix.
**This is a fast and approximate version of the full algorithm, so you will not
get statistically significant output**, and it might not look as *smooth* as the
full-blown algorithm. But you can try out a whole chromosome ;)

### Questions? Problems?

Please use the GitHub issues tracker to report anything, rather than
sending me an e-mail.

### Citation

If you use this software, please acknowledge the following paper

Spill YG, Castillo D, Marti-Renom MA, "Binless normalization of Hi-C data
provides significant interaction and difference detection independently of
resolution", *to be submitted*

### File format description

[TADbit](https://3dgenomes.github.io/TADbit/index.html) `.tsv`: tab-separated
text file containing paired-end mapped reads for a single experiment. All lines
starting with # will be discarded. It has the following columns
    . `id`: ID of the mapped pair (alphanumeric without spaces)
    . `chr1`: chromosome of read 1 (alphanumeric without spaces)
    . `pos1`: position of the 5' end of read 1 (integer)
    . `strand1`: whether the read maps on the forward (1) or reverse (0)
    direction
    . `length1`: how many bases were mapped (integer)
    . `re.up1`: upstream restriction site closest to pos1 (integer <= pos1)
    . `re.dn1`: downstream restriction site closest to pos1 (integer > pos1)
    . `chr2`
    . `pos2`
    . `strand2`
    . `length2`
    . `re.up2`
    . `re.dn2`


Binned raw matrix (used for fast binless): tab or space-separated text file
containing multiple datasets. The first line is a header that must start with
`"name" "bin1" "bin2" "observed"`. Optionally, more columns can be added.
    . `name`: The name of the dataset. Put in "" if you use spaces
    . `bin1`: the label for the bin on the x axis. Put in "" if you use spaces,
    for example "[begin1,end1)"
    . `bin2`: The label for the bin on the y axis
    .`observed`: How many paired-end reads were mapped within (bin1,bin2)

