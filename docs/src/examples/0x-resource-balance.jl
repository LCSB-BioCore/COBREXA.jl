
# Copyright (c) 2021-2024, University of Luxembourg                         #src
# Copyright (c) 2021-2024, Heinrich-Heine University Duesseldorf            #src
#                                                                           #src
# Licensed under the Apache License, Version 2.0 (the "License");           #src
# you may not use this file except in compliance with the License.          #src
# You may obtain a copy of the License at                                   #src
#                                                                           #src
#     http://www.apache.org/licenses/LICENSE-2.0                            #src
#                                                                           #src
# Unless required by applicable law or agreed to in writing, software       #src
# distributed under the License is distributed on an "AS IS" BASIS,         #src
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.  #src
# See the License for the specific language governing permissions and       #src
# limitations under the License.                                            #src

# # Resource balance analysis

using COBREXA
using SBtabFBCModels

import Downloads: download

!isfile("bsubtilis_rba_model.tsv") &&
    download("https://rba.inrae.fr/models/Bacillus-subtilis-168-WT/data/rba-sbtab_model.tsv", "bsubtilis_rba_model.tsv")

!isfile("bsubtilis_rba_parameters.tsv") &&
    download("https://rba.inrae.fr/models/Bacillus-subtilis-168-WT/data/rba-sbtab_parameters.tsv", "bsubtilis_rba_parameters.tsv")

md, tabs = SBtabFBCModels.load_tables("bsubtilis_rba_model.tsv");
p_md, p_tabs = SBtabFBCModels.load_tables("bsubtilis_rba_parameters.tsv");

using DataFrames
using ConstraintTrees
using JSON

todict(s) = JSON.parse(replace(s, "'" => "\""))

summary = DataFrame(tabs[2]) # summary
mets = DataFrame(tabs[3])
rxns = DataFrame(tabs[4])
enzs = DataFrame(tabs[5])
prots = DataFrame(tabs[6])
macro_mols = DataFrame(tabs[7])
compartments = DataFrame(tabs[8])
cell_modules = DataFrame(tabs[9]) # empty
processes = DataFrame(tabs[10]) # machines
cell_targets = DataFrame(tabs[11])
met_cons = DataFrame(tabs[12])
dens_cons = DataFrame(tabs[13])
machine_cap_con = DataFrame(tabs[14])
enz_cap_con = DataFrame(tabs[15])

medium_comp = DataFrame(p_tabs[1]) # mmol/l
enz_cap_f = DataFrame(p_tabs[2]) # 1/h
enz_cap_b = DataFrame(p_tabs[3]) # 1/h
machine_cap = DataFrame(p_tabs[4]) # 1/h
comp_cap = DataFrame(p_tabs[5]) # mmol/gDW
target_vals = DataFrame(p_tabs[6]) # no unit

# uniprot
using BioSequences, FASTX

bsub_proteome = Dict()
FASTAReader(open("Bsub_proteome.fasta")) do reader
    for record in reader
        id = split(identifier(record), "|")[2]
        seq = sequence(record)
        c = split(seq,"")
        uc = String.(unique(c)) 
        bsub_proteome[id] = Dict(u => count(==(u), c) for u in uc)
    end
end

uni_id = Dict(k1 => get(todict(k2), "UniprotID", nothing) for (k1, k2) in zip(prots.ID, prots.Annotation) if !ismissing(k2) && startswith(k2, "{"))

# metabolite mass balances

SBtabFBCModels.parse_reaction_formula(rxns[1,:].ReactionFormula)

processes
todict(processes[1,:].ProcessComponents)
todict(processes[1,:].MachineSubunits)

processes[1,:].InitiationCofactors
