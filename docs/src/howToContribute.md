# How to contribute to/develop COBREXA

If you want to contribute to the `COBREXA` package, please fork the present
repository by following [these instructions](https://docs.github.com/en/github/getting-started-with-github/fork-a-repo).

## Step 1: Retrieve a local version of COBREXA

There are two ways that you can retrieve a local copy of the package: one is to
manually clone the forked repository, and the second one is to use the
integrated Julia package manager.

### Option 1: Manually clone your fork

!!! warning "Warning" 
    Please make sure you have _forked_ the repository, as described above.

You can do this as follows from the command line:

```bash
$ git clone git@github.com:yourUsername/COBREXA.jl.git COBREXA.jl
$ cd COBREXA.jl
$ git checkout -b yourNewBranch origin/develop
```

where `yourUsername` is your Github username and `yourNewBranch` is the name of a new branch.

Then, in order to develop the package, you can install your cloned version as
follows (make sure you are in the `COBREXA.jl` directory):

```julia
(v1.1) pkg> add .
```

(press `]` to get into the packaging environment)

This adds the `COBREXA.jl` package and all its dependencies. You can verify
that the installation worked by typing:

```julia
(v1.1) pkg> status
```

If everything went smoothly, this should print something similar to:

```julia
(v1.1) pkg> status
    Status `~/.julia/environments/v1.1/Project.toml`
    [babc4406] COBREXA v0.1.0 #yourNewBranch (.)
```

Now, you can readily start using the `COBREXA` module:

```julia
julia> using COBREXA
```

### Option 2: Use the Julia package manager

When you are used to using the Julia package manager for developing or
contributing to packages, you can type:

```julia
(v1.1) pkg> dev COBREXA
```

This will install the `COBREXA` package locally and check it out for
development. You can check the location of the package with:

```julia
(v1.1) pkg> status
    Status `~/.julia/environments/v1.4/Project.toml`
  [a03a9c34] COBREXA v0.0.5 [`~/.julia/dev/COBREXA`]
```

The default location of the package is `~/.julia/dev/COBREXA`.

You can then set your remote by executing these commands in a regular shell:

```bash
$ cd ~/.julia/dev/COBREXA
$ git remote rename origin upstream # renames the origin as upstream
$ git remote add origin git@github.com:yourUsername/COBREXA.jl.git
$ git fetch origin
```

where `yourUsername` is your Github username.

!!! warning "Warning" 
    Make sure that your fork is located at `github.com/yourUsername/COBREXA.jl`.

Then, checkout a branch `yourNewBranch`:

```bash
$ cd ~/.julia/dev/COBREXA
$ git checkout -b yourNewBranch origin/develop
```

Then, you can readily use the `COBREXA` package:

```julia
julia> using COBREXA
```

After making changes, precompile the package:

```julia
(v1.1) pkg> precompile
```

## Step 2: Activate COBREXA

!!! warning "Warning" 
    Please note that you cannot use the dependencies of COBREXA directly,
    unless they are installed separately or the environment has been activated:

```julia
(v1.1) pkg> activate .
(COBREXA) pkg> instantiate
```

Now, the environment is activated (you can see it with the prompt change
`(COBREXA) pkg>`). Now, you can use the dependency. For instance:

```julia
julia> using JuMP
```

!!! warning "Warning"
    If you do not  `activate` the environment before using any of the
    dependencies, you will see a red error messages prompting you to install the
    dependency explicitly.

## Step 3: Contribute!

!!! tip "Tip"
    Adding `[skip ci]` in your commit message will skip CI from automatically
    testing that commit. 
