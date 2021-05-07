# # Building and analysing a community model

#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/notebooks/@__NAME__.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/notebooks/@__NAME__.ipynb)

# In this tutorial we will use `COBREXA` to build and analyze a community model  
# consisting of multiple variants of *E. coli* knockouts using the `CoreModel`.
# Here each knockout will only be able to metabolize one sugar. 

# If it is not already present, load a large scale *E. coli* model.

!isfile("iML1515.json") &&
    download("http://bigg.ucsd.edu/static/models/iML1515.json", "iML1515.json")
#
using COBREXA

m = load_model(StandardModel, "iML1515.json")

knockouts = ["EX_glc__D_e",
            "EX_"]
