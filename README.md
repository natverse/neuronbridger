<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->
[![natverse](https://img.shields.io/badge/natverse-Part%20of%20the%20natverse-a241b6)](https://natverse.github.io)
[![Docs](https://img.shields.io/badge/docs-100%25-brightgreen.svg)](http://natverse.github.io/neuronbridger/reference/)
[![Travis build
status](https://travis-ci.com/natverse/neuronbridger.svg?branch=main)](https://travis-ci.com/natverse/neuronbridger)
<img src="man/figures/logo.svg" align="right" height="139" /> [![Codecov
test
coverage](https://codecov.io/gh/natverse/neuronbridger/branch/main/graph/badge.svg)](https://codecov.io/gh/natverse/neuronbridger?branch=main)
[![R build
status](https://github.com/natverse/neuronbridger/workflows/R-CMD-check/badge.svg)](https://github.com/natverse/neuronbridger/actions)
<!-- [![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://www.tidyverse.org/lifecycle/#stable) -->
<!-- [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3843544.svg)](https://doi.org/10.5281/zenodo.3843544) -->
<!-- badges: end -->

neuronbridger
=============

The goal of `neuronbridger` is to provide R client utilities for
interacting with the [NeuronBridge](https://neuronbridge.janelia.org/)
neuron matching service. This can help a user design sparse genetic
driver lines for *D. melanogaster* neurons (i.e.
[split-GAL4](https://www.janelia.org/lab/rubin-lab/our-research/gal4-driver-lines/split-gal4-lines)
lines) and discover what connectome-derived neuronal cell types are
targetted by which extant genetic driver lines. This is becoming
increasingly important as wet-lab biologists attempt to make sense of
connectomic data [(Bates et
al. 2019)](https://doi.org/10.1016/j.conb.2018.12.012).

The [NeuronBridge website](https://neuronbridge.janelia.org/) provides
and interface for loading and examining neuron matching results for *D.
melanogaster* brains. Specifically, its aim is to enable one to search
light (LM) and electron microscopy (EM) data sets of the *D.
melanogaster* nervous system provided by the
[FlyLight](https://www.janelia.org/project-team/flylight) and
[FlyEM](https://www.janelia.org/project-team/flyem) projects at [Janelia
Research Campus](https://www.janelia.org/).

You can find similar neurons based on shape regardless of data set. The
[NeuronBridge](https://neuronbridge.janelia.org/) tool was built by
[FlyLight](https://www.janelia.org/project-team/flylight) and scientific
computing at [Janelia Research Campus](https://www.janelia.org/), see
acknowledgements.

Using this R package in concert with the
[natverse](https://github.com/natverse/natverse) ecosystem is
recommended. For example, the package
[vfbr](https://github.com/jefferis/vfbr) grants R side client utilities
for
[virtualflybrain.org](https://v2.virtualflybrain.org/org.geppetto.frontend/geppetto),
allowing users to pull registered, fly brain confocal stacks for GAL4
lines. The package
[hemibrainr](https://github.com/flyconnectome/hemibrainr) contains more
specific tools for annotations and analysis of the [hemibrain
connectome](https://www.janelia.org/project-team/flyem/hemibrain).

Data sets
---------

The major EM dataset at the time when this package was built was the
[hemibrain
connectome](https://www.janelia.org/project-team/flyem/hemibrain).
Connectome data can be seen using the [neurPrint
website](https://neuprint.janelia.org/help/videos?dataset=hemibrain) and
accessed programmatically in R using
[neuprintr](https://github.com/natverse/neuprintr).

Users of [NeuronBridge](https://neuronbridge.janelia.org/) can discover
which [Genaration 1 GAL4
lines](https://gen1mcfo.janelia.org/cgi-bin/gen1mcfo.cgi) and
[split-GAL4 lines](https://splitgal4.janelia.org/cgi-bin/splitgal4.cgi)
contain which hemibrain connectome neurons. GAL4 lines are genetic
driver lines that can be used to express transgenes (e.g. flourescent
proteins) in relatively small numbers of neurons in live flies. Linking
them to results from the connectome enables researchers to
experimentally manipulate components of the fly brain’s ‘wiring
diagram’.

How have these searches been performed?
---------------------------------------

[NeuronBridge](https://neuronbridge.janelia.org/) uses results from a
‘colour depth mask search’. The methodology has recently been published
[(Otsuna et
al. 2018)](https://www.biorxiv.org/content/10.1101/318006v1). This
method represents 3d voxel space in a 2d image by encoding depth as
colour, and allows for a fast pixel-based comparison across specimens.
Both LM image volumes and EM reconstructions can be represented in this
space, leading to efficient LM-&gt;EM and EM-&gt;LM searching.

Images are represented as colour ‘maximum intensity projection images’
(MIPs). A MIP is there not a ‘stack’ of images, but a single pane.

Here is an example of a MIP for a hemibrain neuron
([542634818](https://neuprint.janelia.org/view?bodyid=542634818), the
DM1 olfactory uPN):

![mip\_em\_example](https://raw.githubusercontent.com/natverse/neuronbridger/main/inst/images/mip_em_example.png)

And for a GAL4 lines that seems to contains that neuron
([R84D10](https://v2.virtualflybrain.org/org.geppetto.frontend/geppetto?id=VFBexp_FBtp0063448),
data from stochastic labelling (MCFO) of line [(Meissener et
al. 2020)](https://doi.org/10.1101/2020.05.29.080473)):

![mip\_gmr\_example](https://raw.githubusercontent.com/natverse/neuronbridger/main/inst/images/mip_gmr_example.png)

[NeuronBridge](https://neuronbridge.janelia.org/) loads precomputed
matches between FlyLight Split GAL4 and MCFO vs Hemibrain 1.1. See its
[release
notes](https://neuronbridge.janelia.org/releasenotes/NEURONBRIDGE) and
[data version release
notes](https://neuronbridge.janelia.org/releasenotes/DATA).

The data is also available on
[quilt](https://open.quiltdata.com/b/janelia-flylight-color-depth/tree).

[NBLAST](https://github.com/natverse/nat.nblast) can also be used to
match EM and LM neurons. However, it works best with ‘skeletonised’
neurons. Using ‘colour depth mask search’ circumvents this
‘skeletonisation’ step for both volumetric EM data and LM data. This
step is still non-trivial. We have used NBLAST to match skeletons
between two EM dataset, see: [neuron matching with
hemibrainr](https://flyconnectome.github.io/hemibrainr/articles/match_making.html).

Why would I need to look at these results in R?
-----------------------------------------------

The [NeuronBridge website](https://neuronbridge.janelia.org/) is great
for examining results. However, if we want to consider many neurons are
once it is easiest to do so in an environment like R. A common goal
these days is to use tools such as NeuronBridge to design ‘split-GAL4’
lines to sparsely target neurons of experimental interest to an
investigator, ideally an individual cell type. A cell type may contain
as little as 1-5 neurons on average [(Bates et
al. 2019)](https://pubmed.ncbi.nlm.nih.gov/30703584/).

With `neuronbridger` one thing we can try to do is search for lines that
are a hit for many neurons we want (e.g. all members of a cell type from
the hemibrain connectome). Another, is to try to avoid lines that seem
to contain neurons we do not want (e.g. similar looking hemibrain
neurons, or common line contaminant you do not want, for some
experimental reason). This could be critical for some experiments.
Ultimately clean lines enable one to see neurons better during, for
example, calcium imaging, or allow experimenters to only manipulate
those cells when performing, for example, optogenetic behavioural
experiments.

Tutorial
========

Let us to grips with `neuronbridger`.

Installation
------------

``` r
# install
if (!require("remotes")) install.packages("remotes")
remotes::install_github("natverse/neuronbridger")

# use 
library(neuronbridger)
```

Example
-------

Now we can have a look at what is available. Let us be interested in the
hemibrain neuron 542634818, the DM1 uPN. This neuron is an olfactory
projection neurons from the antennal lobe glomerulus ‘DM1’ and sends
signals into deeper regions of the fly brain in response to odours in
the fly’s environment, such as apple cider vinegar and other food-like
smells.

### Get MIP information for a connectome neuron

``` r
# What neurons can be searched using NeuronBridge?
all.nb.ids = neuronbridge_ids()
length(neuronbridge_ids)
## So there is over 36k neurons/lines that can be searched

# I am interested in the hemibrain neuron, 542634818. 
## Do we have it?
id = "542634818"
id%in%all.nb.ids
## Yes!

# What do we have on it?
nb.info = neuronbridge_info(id)
View(nb.info)
## So here you see one entry because only one MIP file matches the given ID
## In some cases, esp. for GAL4 line, there may be multiple entries because multiple images for stochastic labelling have been taken.
## In this data frame, you find interesting fields such as the dataset the ID if from (libraryName), the locations at which NeuronBridge saves
## MIP files related to this neuron, the gender of the fly brain from which this data item came.
## Importantly there is an internal NeuronBridge ID for this data item (nb.id, 2820604089967050763).
## Every MIP has its own nb.id because multiple MIPs could relate to the same GAL4 line.
## see ?neuronbridge_info

# Let us now see the related MIP file
dm1.mip = neuronbridge_mip(id)
## This gets every MIP file associated with id

# Plot the MIP image in an rgl viewer
plot_mip(dm1.mip)
```

### Searching for ‘hits’ among genetic driver lines

Now that we know that our hemibrain connectome neuron, the DM1 uPN
(542634818), is in the NeuronBridge database, and we have seen its
colour MIP, we may want to know what its matches are. Let us go for it:

``` r
# In order to look for hits, we need the nb.id not the id
## Why? Well in this case it is moot, but for some data items, i..e GAL4 liens rather than hemibrain neurons, there are multiple MIP files
## Each one may be from a different MCFO experiment, i.e. a different random subset of neurons contained in the full line.
## The id would specify all of these files. The nb.id specifies just the one with which you wish to search.
nb.id = nb.info$nb.id
length(nb.id)
# We just have one

# Find information on hits
nb.hits = neuronbridge_hits(nb.id=nb.id)
## This data frame is already ranked so that the top hits are at the top.

# We can see the top 10 hits
View(nb.hits[1:10,])
## Looks like the Gen1 GMR GAL4 lines R84D10 and R26B04 target the DM1 uPN.
## We could design a split-GAL4 lien using hemidrivers from these lines.
## Though they are likely to also give us similar looking PNs that are not
## the DM1 uPN, in addition to the DM1 uPN ....

# We can take a quick look at the score distribution
hist(nb.hits$normalizedScore)
## As expected, most lines are not going to contain our neuron!

# Let us scan through our top 10 hits and see what they look like
scan_mip(mips = nb.hits, no.hits = 10)

# Good, some nice looking possibilities there.
## We have saved a lot of MIP files in our local directory.
## Which is good if we wish to keep acccessing them, but we can clear with:
neuronbridge_remove_mips()

# This said, best practice would be to set the option 'neuronbridger' to a location
## on your computer, where MIP files cna be downloaded and view at your leisure. I.e.
# options(neuronbridger="PATH/TO/MY/MIPS")
```

### Filtering for the best hits

However, there are a lot of neurons that look like the DM1 uPNs. In
fact, we can see a load of them like this:

``` r
if (!require("hemibrainr")) remotes::install_github("flyconnectome/hemibrainr")

# Load package to examine hemibrain data
library(neuprintr)

# Get all olfactory PNs
all.pns = neuprintr::neuprint_read_neurons(hemibrainr::pn.ids)
## There's a lot of neurons to grab, it might take a while

# Now plot!
nat::nopen3d()
hemibrainr::hemibrain_view()
plot3d(all.pns[id], col = "black", lwd = 5, soma = 1000) # in black, the DM1 uPN
plot3d(all.pns, soma = 500, lwd = 0.5, alpha = 0.25)

# See the problem?
## These are all different neurons, and many different neuron types. There is only oneDM1 uPN. 
## However, there is a lot that looks similar that might also be labelled by the same lines.
```

![em\_dm1\_upn](https://raw.githubusercontent.com/natverse/neuronbridger/main/inst/images/em_dm1_upn.png)

So how can we make sure we only get the PN we want, and not these other
PNs?

``` r
# So this is the one we want
search = "542634818"

# And we do not want these other PNs
avoid = setdiff(hemibrainr::pn.ids,search)

# Search! ASnd get info on those neurons we do not like as well.
nb.a = neuronbridge_avoid(search = search, avoid = avoid)
## This may take a long time as 'avoid' has many IDs
### Luckily, we have a loading bar!

# Examine the results
View(nb.a)
## A lot of our top hits have neurons we would rather not have huh
```

### Which connectome neurons are inside each genetic driver line?

We may also want to predict which neurons are in a given genetic driver
line. Let us consider the case of MB543B, a generation 1 GMR line that
is pretty ‘sparse’ and seems to label neurons from the mushroom body
[(Aso et al., 2014)](https://elifesciences.org/articles/04577). There
are 5 ‘MCFO’ images for this neuron. This means that the FlyLight team
has stochastically labelled a subset of neurons in this line five times
[(Meissener et al. 2020)](https://doi.org/10.1101/2020.05.29.080473).
Each time they took a confocal stack, which was later converted into a
colour MIP. We can examine the ‘MIP search matches’ for each of these
images to try to determine the ‘contents’ of this line i.e. which
neurons from the hemibrain connectome it might target.

These are the individual MIPs for different MCFO experiments for the
GAL4 line MB543B:

<img width="800" alt="MB543B#1" src="https://raw.githubusercontent.com/natverse/neuronbridger/main/inst/images/MB543B_1.png">  
<img width="800" alt="MB543B#2" src="https://raw.githubusercontent.com/natverse/neuronbridger/main/inst/images/MB543B_2.png">
<img width="800" alt="MB543B#3" src="https://raw.githubusercontent.com/natverse/neuronbridger/main/inst/images/MB543B_3.png">
<img width="800" alt="MB543B#4" src="https://raw.githubusercontent.com/natverse/neuronbridger/main/inst/images/MB543B_4.png">
<img width="800" alt="MB543B#5" src="https://raw.githubusercontent.com/natverse/neuronbridger/main/inst/images/MB543B_5.png">

Let us see which connectome neurons have high matching scores:

``` r
# Interesting line that labels mushroom body neurons
## But which ones?
line = "MB543B"

# Let's find out
contents = neuronbridge_line_contents(line = line, threshold = 23000)
## Note! Choosing the right threshold for the nomarlsied MIP-comparison score
### Is critical. This may take some titring for lines you are really interested in
#### Though this value is a good first pass.

# So this is what is likely in this line!
## Note we have meta-data from neuprint! So we can see neuron types and names!!
View(contents)

# Now let's check what these neurons look like
if (!require("neuprintr")) remotes::install_github("natverse/neuprintr")

# Load package
library(neuprintr)
## Note you need to 'login' to neuPrint through R
### Look at the package README and/or examine: ?neuprint_login

# Let's see the EM neurons all together
cont.neurons = neuprintr::neuprint_read_neurons(unique(contents$publishedName))

# Plot in 3d!
hemibrainr::hemibrain_view()
plot3d(cont.neurons, lwd  = 2, soma = 500)
plot3d(hemibrainr::hemibrain.surf, col = "grey", alpha = 0.1)

# And compare with the LM data:
open3d()
mips = neuronbridge_mip(line)
scan_mip(mips,type="images", sleep = 5)
```

![em\_hits](https://raw.githubusercontent.com/natverse/neuronbridger/main/inst/images/em_hits.png)

We can use this to try to work out how to design a ‘split’ line. The
split-GAL4 method allows an experimenter to ‘intersect’ the expression
of two GAL4 lines, and in so doing produce a ‘split-GAL4’ line that only
targets those neurons in both patterns. Note that some MIPs in the
NeuronBridge database already represent split-GAL4 line.

So returning to our DM1 uPN issue, let’s see what we might hope to get
from combining the GAL4 lines that appear to be our bets hit:

``` r
# Based on the results above
line1 = "R84D10" # this line seems to be the top hit for the EM DM1 uPN
line2 = "VT033006" # this is not the secodn highest, but has fewer of the neurons we do not like if we inspect nb.a above.

# What neurons might we expect?
potential.split = neuronbridge_predict_split(line1 = line1, line2 = line2,
                                             threshold1 = 23000, threshold2 = 23000)
message("Potentially ~", nrow(potential.split), " neurons in this predicted split")
View(potential.split)
```

### Conclusions

Bear in mind that this function and some others are very sensitive to
the value you give the threshold argument!! You will have to empirically
determine a good value for this - it might differ given the kind of
neurons you are looking at for example.

Of course - do not go planning any experiments based on the results of
this package without looking at the 3D stacks for GAL4 lines first! The
results from these functions are only as good as the LM-EM matching
scores from the colour MIP searches performed by FlyLight. In addition,
the more stochastic labelling (MCFO) that has been performed on a line,
the more certain we can be of its contents - but many lines need a lot
more MCFO.

Acknowledgments
===============

Data and tools
--------------

[NeuronBridge](https://neuronbridge.janelia.org/) was developed by
scientific advisors (Geoffrey Meissner, Wyatt Korff, Gudrun Ihrke), data
scientists (Hideo Otsuna) and software developers (Jody Clements,
Cristian Goina, Rob Svirskas, Konrad Rokicki, Antje Kazimiers) at
[Janelia Research Campus](https://www.janelia.org/).

Please examine the [NeuronBridge usage
terms](https://neuronbridge.janelia.org/usage). If you use results from
this package, you should cite this package, [Clements et
al. 2020](https://janelia.figshare.com/articles/NeuronBridge_Codebase/12159378/1)
for the NeuronBridge codebase, [Ostsuna et
al. 2018](https://www.biorxiv.org/content/10.1101/318006v1) for the
methodology, and the [relevant
publications](https://neuronbridge.janelia.org/usage) for any data sets
used. See below.

This package was created by [Alexander Shakeel
Bates](https://scholar.google.com/citations?user=BOVTiXIAAAAJ&hl=en)
while in the group of [Rachel
Wilson](https://neuro.hms.harvard.edu/faculty-staff/rachel-wilson). Note
the package version mirrors the [version of NeuronBridge
data]((https://neuronbridge.janelia.org/releasenotes/DATA)) it was made
to work with. You can cite this package as:

``` r
citation(package = "neuronbridger")
```

**Bates AS** (2020). *neuronbridger: R client utilities for interacting
with the neuronbridge matching service.* **R package** version 2.1.1.
<a href="https://github.com/natverse/neuronbridger" class="uri">https://github.com/natverse/neuronbridger</a>

Citations
---------

-   **The hemibrain connectome (hemibrain:v1.1)**: Scheffer, L.K., Xu,
    C.S., Januszewski, M., Lu, Z., Takemura, S.-Y., Hayworth, K.J.,
    Huang, G.B., Shinomiya, K., Maitlin-Shepard, J., Berg, S., et al.
    (2020). A connectome and analysis of the adult Drosophila central
    brain. Elife 9. [doi:
    https://doi.org/10.1101/2020.05.29.080473](https://doi.org/10.1101/2020.05.29.080473)

-   **Gen 1 GAL4 line MCFO**: Meissner, G.W., Dorman, Z., Nern, A.,
    Forster, K., Gibney, T., Jeter, J., Johnson, L., He, Y., Lee, K.,
    Melton, B., et al. (2020). An image resource of subdivided
    Drosophila GAL4-driver expression patterns for neuron-level
    searches. bioRxiv. [doi:
    https://doi.org/10.1101/2020.05.29.080473](https://doi.org/10.1101/2020.05.29.080473)

-   **Gen 1 VT GAL4 lines and hemidrivers**: Tirian, L., and Dickson, B.
    (2017). The VT GAL4, LexA, and split-GAL4 driver line collections
    for targeted expression in the Drosophila nervous system. bioRxiv.
    [doi:
    https://doi.org/10.1101/198648](https://doi.org/10.1101/198648)

-   **Gen 1 GMR GAL4 lines**: Jenett, A., Rubin, G.M., Ngo, T.-T.B.,
    Shepherd, D., Murphy, C., Dionne, H., Pfeiffer, B.D., Cavallaro, A.,
    Hall, D., Jeter, J., et al. (2012). A GAL4-Driver Line Resource for
    Drosophila Neurobiology. Cell Rep. 2, 991–1001. [doi:
    https://doi.org/10.1016/j.celrep.2012.09.011](https://doi.org/10.1016/j.celrep.2012.09.011)

-   **GMR hemidrivers**: Dionne, H., Hibbard, K.L., Cavallaro, A., Kao,
    J.-C., and Rubin, G.M. (2018). Genetic Reagents for Making
    Split-GAL4 Lines in Drosophila. Genetics 209, 31–35. [doi
    :https://doi.org/10.1534/genetics.118.300682](https://doi.org/10.1534/genetics.118.300682)

-   **JRC2018F brain and VNS templates**: Bogovic, J.A., Otsuna, H.,
    Heinrich, L., Ito, M., Jeter, J., Meissner, G.W., Nern, A.,
    Colonell, J., Malkesman, O., Ito, K., et al. (2018). An unbiased
    template of the Drosophila brain and ventral nerve cord. bioRxiv.
    [doi:
    https://doi.org/10.1101/376384](https://doi.org/10.1101/376384)

-   **Colour MIP search tool**: Otsuna, H., Ito, M., and Kawase, T.
    (2018). Color depth MIP mask search: a new tool to expedite
    Split-GAL4 creation. bioRxiv. [doi:
    https://doi.org/10.1101/318006](https://doi.org/10.1101/318006)

-   **NeuronBridge codebase**: Clements, J., Goina, C., Kazimiers A.,
    Otsuna, H., Svirskas, R., Rokicki K. (2020) NeuronBridge Codebase.
    [software](https://janelia.figshare.com/articles/NeuronBridge_Codebase/12159378/1)

-   **Janelia split-GAL4 lines**:
    [various](https://neuronbridge.janelia.org/about)

-   **What is a fly brain cell type**: Bates, A.S., Janssens, J.,
    Jefferis, G.S., and Aerts, S. (2019). Neuronal cell types in the
    fly: single-cell anatomy meets single-cell genomics. Curr. Opin.
    Neurobiol. 56, 125–134. [doi:
    https://doi.org/10.1016/j.conb.2018.12.012](https://doi.org/10.1016/j.conb.2018.12.012)
