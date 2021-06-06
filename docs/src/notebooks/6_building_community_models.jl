# # Building and analysing a small community model

#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/notebooks/@__NAME__.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/notebooks/@__NAME__.ipynb)

# Here we will use `COBREXA` to build and analyze a small community model  
# consisting of three *E. coli* mutants using the `CoreModel`. We will use an
# objective function that enforces equal growth rates.

# ## Load the base model

!isfile("e_coli_core.json") &&
    download("http://bigg.ucsd.edu/static/models/e_coli_core.json", "e_coli_core.json");

using COBREXA
using Tulip

model = load_model(CoreModel, "e_coli_core.json")
