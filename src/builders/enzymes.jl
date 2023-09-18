
# the design space is:
# - total enzyme consumption y/n? (can be done very easily manually given the individual enzyme masses are present)
# - isozymes y/n (separates smoment from gecko)
# - mass groups y/n (this is basically summing up enzymes/isozymes by a label)
# - allow overproduction of enzymes (i.e., have extra variables for enzymes/isozymes to allow some slack)
# - anyone is probably able to do a partial sum over the above things themselves, we should make sure they squash well

#TODO: this does not need the variables allocated, it's just a fancy product
enzyme_mass_sum(; forward_fluxes, reverse_fluxes, enzymes, reaction_enzyme_association) = missing

isozyme_mass_mapping(; forward_fluxes, reverse_fluxes, isozymes, ...) = missing
isozyme_mass_group_mapping(; forward_fluxes, reverse_fluxes, isozymes, ...) = missing
