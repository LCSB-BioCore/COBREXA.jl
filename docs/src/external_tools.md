# External Tools
Numerous external databases exist that can add functionality to constraint based models. 
Currently, `CobraTools.jl` tries to make it easier to incorporate two such databases in analyses.
The first is Brenda, a database that maps enzymes to kinetic data.
The second is Equilibrator, a database that allows one to calculate thermodynamic information, e.g. ``\Delta G``, from reactions.
`CobraTools.jl` includes lightweight interfaces to these sources.

## Brenda Interface
To use the Brenda interface functions, you will need to download the database as a txt file [available here](https://www.brenda-enzymes.org/download_brenda_without_registration.php) (~250 MB).
Once the database has been downloaded, it can be parsed by `parse_brenda`.
```@docs
parse_brenda
```
This function returns an array of `BrendaEntry` structs, which are composed of `EnzymeParams` for each field extracted from the Brenda database.
Currently, only the ID (=EC number), TN (=turn over number), KM (=Michaelis-Menten constant, ``K_M``), KI (=Inhibition term for Michaelis-Menten kintics), and KKM (=ratio of TN/KM) numbers for each enzyme class (ID, or EC number) are extracted. All the structs have pretty printing enabled.
```@docs
CobraTools.BrendaEntry
CobraTools.EnzymeParams
```
```@setup brenda
brenda_loc = model_location = joinpath("..","..", "test", "data", "small_brenda.txt")
```
```@example brenda
using CobraTools

brenda_data = parse_brenda(brenda_loc)
```
```@example brenda
brenda_data[1]
```
```@example brenda
brenda_data[1].TN
```
## Equilibrator Interface
The Equilibrator interface requires that the Equilibrator-API has been installed and can be accessed through Julia's PyCall package. Refer to the [Equilibrator-API website](https://gitlab.com/equilibrator/equilibrator-api) for installation instructions. Within Julia, if you can call `pyimport("equilibrator_api")` successfully, then you will be able to use the functions exposed here. To actually use the functions insert `using PyCall` in your main level script (before or after `using CobraTools`).
```@docs
map_gibbs_rxns
map_gibbs_external
map_gibbs_internal
```
