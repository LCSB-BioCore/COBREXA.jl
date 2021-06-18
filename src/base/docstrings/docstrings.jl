macro _change_bound!(model_type, index_type, example)
    return :(
        """
        change_bound!(
            model::$model_type,
            rxn::$index_type;
            lower_bound = -_constants.default_reaction_bound,
            upper_bound = _constants.default_reaction_bound,
        )
    
    Change the bounds of a reaction with `rxn` in `model` in-place. Note
    that if the bound argument is not supplied then a default (see
    `_constants.default_reaction_bound`) is used. 
    
    See also: [`change_bound`](@ref), [`change_bounds!`](@ref), [`change_bounds!`](@ref) 
    
    # Example
    ```
    change_bound!(model, $example; lower_bound=-10, ub=10)
    change_bound!(model, $example; lower_bound=-10.2) # upper_bound is set to _constants.default_reaction_bound
    change_bound!(model, $example; upper_bound=10)
    ```
    """
    )
end
