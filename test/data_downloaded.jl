
# TODO this should be downloaded by documentation scripts

isdir("downloaded") || mkdir("downloaded")
df(s) = joinpath("downloaded", s)
const dl = A.download_data_file

model_paths = Dict{String,String}(
    "iJO1366.json" => dl(
        "http://bigg.ucsd.edu/static/models/iJO1366.json",
        df("iJO1366.json"),
        "9376a93f62ad430719f23e612154dd94c67e0d7c9545ed9d17a4d0c347672313",
    ),
    "iJO1366.mat" => dl(
        "http://bigg.ucsd.edu/static/models/iJO1366.mat",
        df("iJO1366.mat"),
        "b5cfe21b6369a00e45d600b783f89521f5cc953e25ee52c5f1d0a3f83743be30",
    ),
    "iJO1366.xml" => dl(
        "http://bigg.ucsd.edu/static/models/iJO1366.xml",
        df("iJO1366.xml"),
        "d6d9ec61ef6f155db5bb2f49549119dc13b96f6098b403ef82ea4240b27232eb",
    ),
    "ecoli_core_model.xml" => dl(
        "http://systemsbiology.ucsd.edu/sites/systemsbiology.ucsd.edu/files/Attachments/Images/downloads/Ecoli_core/ecoli_core_model.xml",
        df("ecoli_core_model.xml"),
        "78692f8509fb36534f4f9b6ade23b23552044f3ecd8b48d84d484636922ae907",
    ),
    "e_coli_core.json" => dl(
        "http://bigg.ucsd.edu/static/models/e_coli_core.json",
        df("e_coli_core.json"),
        "7bedec10576cfe935b19218dc881f3fb14f890a1871448fc19a9b4ee15b448d8",
    ),
    "e_coli_core.xml" => dl(
        "http://bigg.ucsd.edu/static/models/e_coli_core.xml",
        df("e_coli_core.xml"),
        "b4db506aeed0e434c1f5f1fdd35feda0dfe5d82badcfda0e9d1342335ab31116",
    ),
    "e_coli_core.mat" => dl(
        "http://bigg.ucsd.edu/static/models/e_coli_core.mat",
        df("e_coli_core.mat"),
        "478e6fa047ede1b248975d7565208ac9363a44dd64aad1900b63127394f4175b",
    ),
    "iJR904.mat" => dl(
        "http://bigg.ucsd.edu/static/models/iJR904.mat",
        df("iJR904.mat"),
        "d17be86293d4caafc32b829da4e2d0d76eb45e1bb837e0138327043a83e20c6e",
    ),
    "Recon3D.json" => dl(
        "http://bigg.ucsd.edu/static/models/Recon3D.json",
        df("Recon3D.json"),
        "aba925f17547a42f9fdb4c1f685d89364cbf4979bbe7862e9f793af7169b26d5",
    ),
    "yeast-GEM.mat" => dl(
        "https://github.com/SysBioChalmers/yeast-GEM/raw/v8.6.2/model/yeast-GEM.mat",
        df("yeast-GEM.mat"),
        "c2587e258501737e0141cd47e0f854a60a47faee2d4c6ad582a00e437676b181",
    ),
    "yeast-GEM.xml" => dl(
        "https://github.com/SysBioChalmers/yeast-GEM/raw/v8.6.2/model/yeast-GEM.xml",
        df("yeast-GEM.xml"),
        "c728b09d849b744ec7640cbf15776d40fb2d9cbd0b76a840a8661b626c1bd4be",
    ),
)
