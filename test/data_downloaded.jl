
function check_data_file_hash(path, expected_checksum)
    actual_checksum = bytes2hex(sha256(open(path)))
    if actual_checksum != expected_checksum
        @error "The downloaded data file `$path' seems to be different from the expected one. Tests will likely fail." actual_checksum expected_checksum
    end
end

function download_data_file(url, path, hash)
    if isfile(path)
        check_data_file_hash(path, hash)
        @info "using cached `$path'"
        return path
    end

    Downloads.download(url, path)
    check_data_file_hash(path, hash)
    return path
end


isdir("downloaded") || mkdir("downloaded")
df(s) = joinpath("downloaded", s)

model_paths = Dict{String,String}(
    "iJO1366.json" => download_data_file(
        "http://bigg.ucsd.edu/static/models/iJO1366.json",
        df("iJO1366.json"),
        "9376a93f62ad430719f23e612154dd94c67e0d7c9545ed9d17a4d0c347672313",
    ),
    "iJO1366.mat" => download_data_file(
        "http://bigg.ucsd.edu/static/models/iJO1366.mat",
        df("iJO1366.mat"),
        "b5cfe21b6369a00e45d600b783f89521f5cc953e25ee52c5f1d0a3f83743be30",
    ),
    "iJO1366.xml" => download_data_file(
        "http://bigg.ucsd.edu/static/models/iJO1366.xml",
        df("iJO1366.xml"),
        "d6d9ec61ef6f155db5bb2f49549119dc13b96f6098b403ef82ea4240b27232eb",
    ),
    "ecoli_core_model.xml" => download_data_file(
        "http://systemsbiology.ucsd.edu/sites/systemsbiology.ucsd.edu/files/Attachments/Images/downloads/Ecoli_core/ecoli_core_model.xml",
        df("ecoli_core_model.xml"),
        "78692f8509fb36534f4f9b6ade23b23552044f3ecd8b48d84d484636922ae907",
    ),
    "e_coli_core.json" => download_data_file(
        "http://bigg.ucsd.edu/static/models/e_coli_core.json",
        df("e_coli_core.json"),
        "7bedec10576cfe935b19218dc881f3fb14f890a1871448fc19a9b4ee15b448d8",
    ),
    "e_coli_core.xml" => download_data_file(
        "http://bigg.ucsd.edu/static/models/e_coli_core.xml",
        df("e_coli_core.xml"),
        "b4db506aeed0e434c1f5f1fdd35feda0dfe5d82badcfda0e9d1342335ab31116",
    ),
    "e_coli_core.mat" => download_data_file(
        "http://bigg.ucsd.edu/static/models/e_coli_core.mat",
        df("e_coli_core.mat"),
        "478e6fa047ede1b248975d7565208ac9363a44dd64aad1900b63127394f4175b",
    ),
    "iJR904.mat" => download_data_file(
        "http://bigg.ucsd.edu/static/models/iJR904.mat",
        df("iJR904.mat"),
        "d17be86293d4caafc32b829da4e2d0d76eb45e1bb837e0138327043a83e20c6e",
    ),
)
