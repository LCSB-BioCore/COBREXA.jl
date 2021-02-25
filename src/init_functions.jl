
function __init__()
   @require PyCall="438e738f-606a-5dbb-bf0a-cddfbfd45ab0" include("equilibrator_tools.jl") # only add the equilibrator functions if PyCall is brought into scope
end
