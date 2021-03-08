"""
EnzymeParams

This struct has `val`, `substrate`, `ph`, and `temperature` fields, extrated from the database.
"""
struct EnzymeParams
    val :: Union{Nothing, Float64}
    substrate :: Union{Nothing, String}
    ph :: Union{Nothing, Float64}
    temperature :: Union{Nothing, Float64}
end

function Base.show(io::IO, m::EnzymeParams)
    println(io, "Value: ", m.val)
    println(io, "Substrate: ", m.substrate)
    println(io, "pH: ", m.ph)
    println(io, "Temperature: ", m.temperature)
end

"""
BrendaEntry

This struct has ID, TN, KI, KM, KKM fields, corresponding to the fields in the database.
"""
struct BrendaEntry
    ID :: String # EC number
    TN :: Array{EnzymeParams, 1} # turnover number [[TN, substrate, temp, ph], ...]
    KI :: Array{EnzymeParams, 1} # turnover number [[KI, substrate, temp, ph], ...]
    KM :: Array{EnzymeParams, 1} # turnover number [[KM, substrate, temp, ph], ...]
    KKM :: Array{EnzymeParams, 1} # turnover number [[KKM, substrate, temp, ph], ...]

    # Keywords below are not implemented (yet)...
    # AC	# activating compound
    # AP	# application
    # CL	# cloned
    # CR	# crystallization
    # EN	# engineering
    # EXP	# expression
    # GI	# general information on enzyme
    # GS	# general stability
    # IN	# inhibitors
    # LO	# Localization
    # ME	# Metals/Ions
    # MW	# molecular weight
    # NSP	# natural substrates/products	reversibilty information in {...}
    # OS	# oxygen stability
    # OSS	# organic solvent stability
    # PHO	# pH-optimum
    # PHR	# pH-range
    # PHS	# pH stability
    # PI	# isoelectric point
    # PM	# posttranslation modification
    # PR	# protein
    # PU	# purification
    # RE	# reaction catalyzed
    # RF	# references
    # REN	# renatured
    # RN	# accepted name (IUPAC) 
    # RT	# reaction type
    # SA	# specific activity
    # SY	# synonyms
    # SP	# substrates/products	reversibilty information in {...}
    # SS	# storage stability
    # ST	# source/tissue
    # SU	# subunits
    # SY	# systematic name 
    # TO	# temperature optimum
    # TR	# temperature range
    # TS  # temperature	stability    
end

function Base.show(io::IO, m::BrendaEntry)
    println(io, "EC: ", m.ID)
end

"""
    parse_brenda(brenda_txt_file_path)

Parse the Brenda txt file line by line. 
Download it from https://www.brenda-enzymes.org/download_brenda_without_registration.php.
Returns an array with Brenda 
"""
function parse_brenda(brenda_loc)
    brenda_data = Array{BrendaEntry, 1}()
    open(brenda_loc) do io
        ID = ""
        TN = Array{EnzymeParams, 1}()
        KI = Array{EnzymeParams, 1}()
        KM = Array{EnzymeParams, 1}()
        KKM = Array{EnzymeParams, 1}()
        firstline = true
        for ln in eachline(io)
            # println(ID)
            if startswith(ln, "ID\t")
                if firstline || ID == "0.0.0.0"
                    firstline = false
                else
                    push!(brenda_data, BrendaEntry(ID, TN, KI, KM, KKM))
                end
                ID = split(split(ln, "\t")[2], " ")[1]
                TN = Array{EnzymeParams, 1}()
                KI = Array{EnzymeParams, 1}()
                KM = Array{EnzymeParams, 1}()
                KKM = Array{EnzymeParams, 1}()
            end

            if startswith(ln, "TN\t")
                tn_data = get_enzyme_data(ln)
                if !isnothing(tn_data.val)
                    push!(TN, tn_data)
                end
            end
            if startswith(ln, "KI\t")
                ki_data = get_enzyme_data(ln)
                if !isnothing(ki_data.val)
                    push!(KI, ki_data)
                end
            end
            if startswith(ln, "KM\t")
                km_data = get_enzyme_data(ln)
                if !isnothing(km_data.val)
                    push!(KM, km_data)
                end
            end
            if startswith(ln, "KKM\t")
                kkm_data = get_enzyme_data(ln)
                if !isnothing(kkm_data.val)
                    push!(KKM, kkm_data)
                end
            end
        end
    end

    return brenda_data
end

"""
[value, substrate, ph, temperature] = get_enzyme_data(ln::String)

Return enzyme data from line in Brenda txt file.
"""
function get_enzyme_data(ln::String)
    prts = split(ln, "\t")[2]
    
    val_ = length(split(prts," ")) < 2 ? nothing : split(prts," ")[2]
    val = isnothing(val_) || contains(val_, "-") ? nothing : parse(Float64, val_)
    
    substrate_ = match(r"\{{1}\d*\D*\d*\}{1}", prts)
    substrate = isnothing(substrate_) ? nothing : substrate_.match[2:end-1]

    ph_ = match(r"(pH\s){1}\d{1,2}\.\d{1}", prts)
    ph = isnothing(ph_) ? nothing : parse(Float64, split(ph_.match, " ")[2])

    temp_ = match(r"\d+\.?\d{0,2}°C", prts)
    temp = isnothing(temp_) ? nothing : parse(Float64, split(temp_.match,"°")[1])

    return EnzymeParams(val, substrate, ph, temp)
end