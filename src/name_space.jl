function name_space_mapping()
    metanetx_url = "https://www.metanetx.org/cgi-bin/mnxget/mnxref/reac_xref.tsv"
    nspfile = download(metanetx_url, joinpath("..", ))
end