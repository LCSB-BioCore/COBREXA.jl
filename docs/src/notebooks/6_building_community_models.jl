# # Building and analysing a community model

# In this tutorial we will use `COBREXA` to build and analyze a community model  
# consisting of multiple variants of *E. coli* knockouts using the `CoreModel`.

# If it is not already present, load a large scale *E. coli* model.

!isfile("iML1515.json") &&
    download("http://bigg.ucsd.edu/static/models/iML1515.json", "iML515.json")
#
using COBREXA

m = load_model(CoreModel, "iML1515.json")
