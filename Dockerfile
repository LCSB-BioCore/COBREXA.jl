FROM julia

RUN julia -e 'import Pkg; Pkg.add("COBREXA"); Pkg.resolve(); Pkg.status(); Pkg.instantiate(); Pkg.precompile()'

CMD ["julia"]
