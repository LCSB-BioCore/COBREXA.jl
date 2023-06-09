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
    Dict(
        sid => Dict(vars[vidx] => val for (vidx, val) in zip(findnz(mtx[:, sidx])...)) for
        (sidx, sid) in enumerate(semantics)
    )
end

const Semantic = Tuple{Function,Function,Function,Function}

const variable_semantics = Dict{Symbol,Semantic}()

"""
$(TYPEDSIGNATURES)

Get a tuple of functions that work with the given semantics, or `nothing` if
the semantics doesn't exist.
"""
function get_semantics(semantics::Symbol)::Types.Maybe{Semantic}
    get(variable_semantics, semantics, nothing)
end

"""
$(TYPEDSIGNATURES)

Like [`get_semantics`](@ref) but throws a `DomainError` if the semantics is not
available.
"""
function semantics(semantics::Symbol)::Semantic
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

    plural = Symbol(sym, :s)
    count = Symbol(:n_, plural)
    mapping = Symbol(sym, :_variables)
    mapping_mtx = Symbol(sym, :_variables_matrix)

    pluralfn = Expr(
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
        :(function $plural(a::AbstractMetabolicModel)::Vector{String}
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
vector returned by [`$plural`].
""",
        ),
        :(function $count(a::AbstractMetabolicModel)::Int
            0
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
[`$plural`](@ref) for semantics.

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
            make_mapping_mtx(variables(a), $plural(a), $mapping(a))
        end),
    )

    Base.eval.(Ref(themodule), [pluralfn, countfn, mappingfn, mtxfn])

    Base.eval(themodule, :(function $plural(w::AbstractModelWrapper)::Vector{String}
        $plural(unwrap_model(w))
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

    # TODO here we would normally also overload the matrix function, but that
    # one will break once anyone touches variables of the models (which is
    # common). We should have a macro like @model_does_not_modify_variable_set
    # that adds the overloads. Or perhaps AbstractModelWrapperWithSameVariables?
    #
    # The same probably goes for other semantics;
    # AbstractModelWrapperThatOnlyTouchesSemantics(...) ? (Which has an
    # alternative in forcing people to overload all semantic functions in all
    # cases of adding semantics, which might actually be the right way.)

    themodule.Internal.variable_semantics[sym] =
        Base.eval.(Ref(themodule), (plural, count, mapping, mapping_mtx))
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
        $Accessors.reactions(model::$m) = $Accessors.variables(model)
        $Accessors.n_reactions(model::$m) = $Accessors.n_variables(model)
        $Accessors.reaction_variables(model::$m) =
            Dict(var => Dict(var => 1.0) for var in $Accessors.variables(model))
        $Accessors.reaction_variables_matrix(model::$m) =
            $SparseArrays.spdiagm(fill(1.0, $Accessors.n_variables(model)))
    end
end

"""
$(TYPEDSIGNATURES)

Convenience helper -- for many models boundary variables are exchange reactions.
"""
macro all_boundary_variables_are_exchanges(mt)
    m = esc(mt)
    quote
        $Accessors.exchanges(model::$m) = [
            vid for vid in $Accessors.variables(model) if
            length(keys($Accessors.reaction_stoichiometry(model, vid))) == 1
        ]
        $Accessors.n_exchanges(model::$m) = length($Accessors.exchanges(model))
        $Accessors.exchange_variables(model::$m) = Dict(
            vid => Dict(vid => 1.0) for vid in $Accessors.variables(model) if
            length(keys($Accessors.reaction_stoichiometry(model, vid))) == 1
        )
    end
end
