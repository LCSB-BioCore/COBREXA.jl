"""
$(TYPEDSIGNATURES)

A helper function to quickly create a sparse matrix from a dictionary that
describes it. Reverse of [`make_mapping_dict`](@ref).
"""
function make_mapping_mtx(
    vars::Vector{String},
    semantics::Vector{String},
    var_sem_val::Dict{String,Dict{String,Float64}},
)::Types.SparseMat
    # TODO: move this to general utils or so
    rowidx = Dict(vars .=> 1:length(vars))
    colidx = Dict(semantics .=> 1:length(semantics))
    n = sum(length.(values(var_sem_val)))
    R = Vector{Int}(undef, n)
    C = Vector{Int}(undef, n)
    V = Vector{Float64}(undef, n)
    i = 1
    for (cid, col_val) in var_sem_val
        for (rid, val) in col_val
            R[i] = rowidx[rid]
            C[i] = colidx[cid]
            V[i] = val
            i += 1
        end
    end
    sparse(R, C, V, length(vars), length(semantics))
end

"""
$(TYPEDSIGNATURES)

A helper function to quickly create a sparse matrix from a dictionary that
describes it. Reverse of [`make_mapping_mtx`](@ref).
"""
function make_mapping_dict(
    vars::Vector{String},
    semantics::Vector{String},
    mtx::Types.SparseMat,
)::Dict{String,Dict{String,Float64}}
    # TODO: move this to general utils
    Dict(
        sid => Dict(vars[vidx] => val for (vidx, val) in zip(findnz(mtx[:, sidx])...)) for
        (sidx, sid) in enumerate(semantics)
    )
end

"""
$(TYPEDEF)

A structure of function that stores implementation of a given variable
semantic. Use [`semantics`](@ref) for lookup. If you want to create and
register new semantics, use [`@make_variable_semantics`](@ref).
"""
Base.@kwdef struct Semantics
    # TODO: move this to types?
    """
    Returns a vector identifiers that describe views in the given semantics.
    """
    ids::Function

    """
    Returns the size of the vector of the identifiers, in a possibly more
    efficient way than measuring the length of the ID vector.
    """
    count::Function

    """
    Returns a mapping of semantic values to variables IDs in the model.
    """
    mapping::Function

    """
    Same as `mapping` but returns a matrix (with variables in rows and the
    semantic values in columns), which is possibly more efficient or handy in
    specific cases.
    """
    mapping_matrix::Function

    """
    Returns either `nothing` if the semantics does not have specific
    constraints associated, or a vector of floating-point values that represent
    equality constraints on semantic values, or a tuple of 2 vectors that
    represent lower and upper bounds on semantic values.
    """
    bounds::Function
end

const variable_semantics = Dict{Symbol,Semantics}()

"""
$(TYPEDSIGNATURES)

Get all available semantics in a symbol-indexed dictionary.
"""
function get_semantics()::Dict{Symbol,Semantics}
    variable_semantics
end

"""
$(TYPEDSIGNATURES)

Get a tuple of functions that work with the given semantics, or `nothing` if
the semantics doesn't exist.
"""
function get_semantics(semantics::Symbol)::Types.Maybe{Semantics}
    get(variable_semantics, semantics, nothing)
end

"""
$(TYPEDSIGNATURES)

Like [`get_semantics`](@ref) but throws a `DomainError` if the semantics is not
available.
"""
function semantics(semantics::Symbol)::Semantics
    res = get_semantics(semantics)
    isnothing(res) && throw(DomainError(semantics, "unknown semantics"))
    res
end

"""
$(TYPEDSIGNATURES)

Inject a new functionality for variable semantics defined by `sym` into
`themodule` (which should ideally be COBREXA.Accessors).

`name` is a human readable description of the semantic object. `example` is a
string that closer describes the semantics, which is inserted into the main
semantic-accessing function.
"""
function make_variable_semantics(
    themodule::Module,
    source,
    sym::Symbol,
    name::String,
    example::String,
)
    haskey(themodule.Internal.variable_semantics, sym) && return

    ids = Symbol(sym, :_ids)
    count = Symbol(sym, :_count)
    mapping = Symbol(sym, :_variables)
    mapping_mtx = Symbol(sym, :_variables_matrix)
    bounds = Symbol(sym, :_bounds)

    idsfn = Expr(
        :macrocall,
        Symbol("@doc"),
        source,
        Expr(
            :string,
            :TYPEDSIGNATURES,
            """

List the semantics for model variables that can be interpreted as $name.

$example

Use [`$count`](@ref) to quickly determine the amount of $name in the
model. See the documentation of [`$mapping`](@ref) for closer
definition of the correspondence of $name and model variables.
""",
        ),
        :(function $ids(a::AbstractMetabolicModel)::Vector{String}
            String[]
        end),
    )

    countfn = Expr(
        :macrocall,
        Symbol("@doc"),
        source,
        Expr(
            :string,
            :TYPEDSIGNATURES,
            """

Count of $name that the model describes, should be equal to the length of
vector returned by [`$ids`].
""",
        ),
        :(function $count(a::AbstractMetabolicModel)::Int
            length($ids(a))
        end),
    )

    mappingfn = Expr(
        :macrocall,
        Symbol("@doc"),
        source,
        Expr(
            :string,
            :TYPEDSIGNATURES,
            """

Bipartite mapping of $name described by the model to the actual
variables in the model. Returns a dictionary of $name assigned to the
variable IDs and their linear coefficients. See the documentation of
[`$ids`](@ref) for semantics.

To improve the performance, you may want to use [`$mapping_mtx`](@ref).
""",
        ),
        :(function $mapping(a::AbstractMetabolicModel)::Dict{String,Dict{String,Float64}}
                Dict()
            end),
    )

    mtxfn = Expr(
        :macrocall,
        Symbol("@doc"),
        source,
        Expr(
            :string,
            :TYPEDSIGNATURES,
            """

Bipartite mapping of $name described by the model to the actual
variables in the model, described as a sparse matrix mapping with rows
corresponding to model variables and columns corresponding to $name.

By default, this is derived from [`$mapping`](@ref) in all models. For
safety reasons, this is never automatically inherited by wrappers.
""",
        ),
        :(function $mapping_mtx(a::AbstractMetabolicModel)::SparseMat
            make_mapping_mtx(variable_ids(a), $ids(a), $mapping(a))
        end),
    )

    boundsfn = Expr(
        :macrocall,
        Symbol("@doc"),
        source,
        Expr(
            :string,
            :TYPEDSIGNATURES,
            """

Bounds for $name described by the model. Either returns `nothing` if there are
no bounds, or a vector of floats with equality bounds, or a tuple of 2 vectors
with lower and upper bounds.
""",
        ),
        :(function $bounds(a::AbstractMetabolicModel)
            $(sym == :metabolite ? :(spzeros($count(a))) : nothing)
        end),
    )

    Base.eval.(Ref(themodule), [idsfn, countfn, mappingfn, mtxfn, boundsfn])

    # extend the AbstractModelWrapper defaults
    Base.eval(themodule, :(function $ids(w::AbstractModelWrapper)::Vector{String}
        $ids(unwrap_model(w))
    end))

    Base.eval(themodule, :(function $count(w::AbstractModelWrapper)::Int
        $count(unwrap_model(w))
    end))

    Base.eval(
        themodule,
        :(function $mapping(w::AbstractModelWrapper)::Dict{String,Dict{String,Float64}}
            $mapping(unwrap_model(w))
        end),
    )

    Base.eval(themodule, :(function $bounds(w::AbstractModelWrapper)
        $bounds(unwrap_model(w))
    end))

    # TODO here we would normally also overload the matrix function, but that
    # one will break once anyone touches variables of the models (which is quite
    # common). We should have a macro like @model_does_not_modify_variable_set
    # that adds the overloads. Or perhaps AbstractModelWrapperWithSameVariables?
    #
    # The same probably goes for other semantics;
    # AbstractModelWrapperThatOnlyTouchesSemantics(...) ? (Which has an
    # alternative in forcing people to overload all semantic functions in all
    # cases of adding semantics, which might actually be the right way.)

    themodule.Internal.variable_semantics[sym] =
        Semantics(Base.eval.(Ref(themodule), (ids, count, mapping, mapping_mtx, bounds))...)
end

"""
$(TYPEDSIGNATURES)

Convenience macro for running [`make_variable_semantics`](@ref).
"""
macro make_variable_semantics(sym, name, example)
    src = __source__
    quote
        $make_variable_semantics($Accessors, $src, $sym, $name, $example)
    end
end

"""
$(TYPEDSIGNATURES)

Convenience helper -- many models carry no other variable semantics than the
reactions; this macro declares precisely the same about the model type.
"""
macro all_variables_are_reactions(mt)
    m = esc(mt)
    quote
        $Accessors.reaction_ids(model::$m) = $Accessors.variable_ids(model)
        $Accessors.reaction_count(model::$m) = $Accessors.variable_count(model)
        $Accessors.reaction_variables(model::$m) =
            Dict(var => Dict(var => 1.0) for var in $Accessors.variable_ids(model))
        $Accessors.reaction_variables_matrix(model::$m) =
            $SparseArrays.spdiagm(fill(1.0, $Accessors.variable_count(model)))
    end
end
