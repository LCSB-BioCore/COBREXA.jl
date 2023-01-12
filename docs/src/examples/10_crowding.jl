# # FBA with ???

# Let's starting with loading the models and packages.

!isfile("e_coli_core.xml") &&
    download("http://bigg.ucsd.edu/static/models/e_coli_core.xml", "e_coli_core.xml")

using COBREXA, Tulip

model = load_model("e_coli_core.xml")

import Random
Random.seed!(1) # for repeatability of random numbers below

# TODO add crowding example via smoment
