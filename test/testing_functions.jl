"""
Test if model is the same after it was read in, saved, and then re-read.
"""
function read_write_read_test(model, format)
    tmpfile = joinpath("data", "temp." * format)

    write_model(model, tmpfile)
    tmpmodel = read_model(tmpfile, StandardModel)

    try
        rm(tmpfile)
    catch
        @warn "The file $tmpfile could not be removed"
    end
    model_comparison_test(model, tmpmodel)
end
