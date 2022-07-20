
# Faster model distribution

When working with a large model, it may happen that the time required to send
the model to the worker nodes takes a significant portion of the total
computation time.

You can use Julia serialization to prevent this. Because the shared filesystem
is a common feature of most HPC installations around, you can very easily
utilize it to broadcast a serialized version of the model to all worker nodes.
In COBREXA, that functionality is wrapped in a [`Serialized`](@ref) model type, which provides a tiny abstraction around this functionality:

- you call [`serialize_model`](@ref) before you start your analysis to place
  the model to the shared storage (and, possibly, free the RAM required to hold
  the model)
- the freshly created [`Serialized`](@ref) model type is tiny -- it only stores
  the file name where the model data can be found
- all analysis functions automatically call [`precache!`](@ref) on the model to
  get the actual model data loaded before the model is used, which
  transparently causes the loading from the shared media, and thus the fast
  distribution

The use is pretty straightforward. After you have your model loaded, you simply
convert it to the small serialized form:

```julia
model = load_model(...)
# ... prepare the model ...

cachefile = tempname(".") 
sm = serialize_model(model, cachefile)
```

Then, the analysis functions is run as it would be with a "normal" model:
```julia
screen(sm, ...)
```

After the analysis is done, it is useful to remove the temporary file:
```julia
rm(cachefile)
```

!!! warn "Caveats of working with temporary files"
    Always ensure that the temporary filenames are unique -- if there are two
    jobs running in parallel and both use the same filename to cache and
    distribute a model, both are likely going to crash, and almost surely
    produce wrong results.  At the same time, avoid creating the temporary
    files in the default location of `/tmp/*`, as the temporary folder is
    local, thus seldom shared between HPC nodes.
