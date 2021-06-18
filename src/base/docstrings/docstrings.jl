macro _change_bound_s_bang(model_type, index_name, index_type, example, isplural, isinplace)
    if isplural
        if isinplace
            add_bang = "!"
            descr_sentence = "    Change the bounds of all reactions with `$(index_name)` in `model` in-place."
        else
            add_bang = ""
            descr_sentence = "    Return a shallow copy of the `model` where the bounds of all reactions with
            `$(index_name)` have been changed."
        end
        add_s = "s"
        bounds1 = "lower_bounds=[-10, -20], upper_bounds=[10, 22]"
        bounds2 = "lower_bounds=[-10.2, -14.3]"
        bounds3 = "upper_bounds=[10.2, 23]"
        lbs = "fill(-_constants.default_reaction_bound, length($(index_name)))"
        ubs = "fill(_constants.default_reaction_bound, length($(index_name)))"
    else
        if isinplace
            add_bang = "!"
            descr_sentence = "    Change the bounds of a reaction `$(index_name)` in `model` in-place."
        else
            add_bang = ""
            descr_sentence = "    Return a shallow copy of the `model` where the bounds of reaction with
            `$(index_name)` has been changed."
        end
        add_s = ""
        bounds1 = "lower_bound=-10, upper_bound=10"
        bounds2 = "lower_bound=-10.2"
        bounds3 = "upper_bound=10"
        lbs = "-_constants.default_reaction_bound"
        ubs = "_constants.default_reaction_bound"
    end
    return :(
        """
        change_bound$($add_s)$($add_bang)(
            model::$($model_type),
            $($index_name)::$($index_type);
            lower_bound$($add_s) = $($lbs)
            upper_bound$($add_s) = $($ubs),
        )
    
    $($descr_sentence) Note
    that if the bound argument is not supplied then a default (see
    `_constants.default_reaction_bound`) is used. 
    
    See also: [`change_bound`](@ref), [`change_bounds!`](@ref), [`change_bound!`](@ref), [`change_bounds!`](@ref)
    
    # Example
    ```
    change_bound$($add_s)$($add_bang)(model, $($example); $($bounds1))
    change_bound$($add_s)$($add_bang)(model, $($example); $($bounds2)) # upper_bound$($add_s) set to _constants.default_reaction_bound
    change_bound$($add_s)$($add_bang)(model, $($example); $($bounds3))
    ```
    """
    )
end
