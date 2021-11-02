# Contributing to COBREXA.jl

:+1::tada: Thanks for taking the time to contribute to
[COBREXA.jl](https://github.com/LCSB-BioCore/COBREXA.jl)! :tada::+1:

## How to report a bug or suggest an enhancement

Please use the [GitHub issue
tracker](https://github.com/LCSB-BioCore/COBREXA.jl/issues) to report any
problems with the software, and discuss any potential questions about COBREXA
use.

Before creating bug reports, please check [the open
issues](https://github.com/LCSB-BioCore/COBREXA.jl/issues), you might find
out that the issue is already reported and known.

General guidelines for reporting issues:

- If creating a bug report, include a complete description of how we can
  reproduce the bug, including e.g. links to datasets and any external scripts
  used. Ideally, try to create a code snippet that causes the problem on a
  fresh installation of COBREXA.jl (often called the "minimal crashing
  example")
- If possible, use the supplied issue templates and fill in all fields.
- If your issue is already described in an issue that is "closed", do not
  reopen it. Instead, open a new issue and include a link to the original
  issue. (The fact that the original issue might have been mistakenly closed
  may be an issue on its own.)
- Enhancement proposals should refer a viable way for implementing the
  enhancement. If there are multiple possibilities for implementation, we will
  welcome a discussion about which one is optimal for COBREXA.jl.

## How to test a development version of the package?

### Step 1: Load COBREXA.jl from the source from the git repository

There are two ways that you can retrieve a local copy of the development repo:
you can either clone the repository manually, or use Julia package manager to
get a development version for you.

#### Option 1: Using Julia package manager

When you are used to using the Julia package manager for developing or
contributing to packages, you can type:

```julia
(v1.6) pkg> dev COBREXA
```

This will install the `COBREXA` package locally and check it out for
development. You can check the location of the package with:

```julia
(v1.6) pkg> status
    Status `~/.julia/environments/v1.4/Project.toml`
  [a03a9c34] COBREXA v0.0.5 [`~/.julia/dev/COBREXA`]
```

The default location of the package is `~/.julia/dev/COBREXA`.

#### Option 2: Cloning with `git` manually

You can use `git` to get the sources as follows:

```bash
$ git clone git@github.com:LCSB-BioCore/COBREXA.jl.git
```

When the cloning process finishes, you shold see the package cloned in a new
directory `COBREXA.jl`. To install this version to your Julia, change to the
directory first, and start Julia:

```bash
$ cd COBREXA.jl
$ julia
```

With Julia, you can install the development version of the package from the
directory as follows:

```julia
(v1.6) pkg> add .
```

(press `]` to get into the packaging environment)

This adds the `COBREXA.jl` package and all its dependencies. You can verify
that the installation worked by typing:

```julia
(v1.6) pkg> status
```

If you are planning to develop the package, it is often easier to install the
package in development mode, with `dev` command:

```julia
(v1.6) pkg> dev .
```

That causes the package to always load with whatever code changes that you
added to the source directory.

#### Finally: load COBREXA.jl

With both of above options, you should get COBREXA.jl installed, which means
that the following command should, without errors, load the package and make
COBREXA.jl functions available for testing:

```julia
julia> using COBREXA
```

You may now freely modify the code and test the result.

Remember that if you want to work in the environment of the package, you need
to *activate* it. That causes, among other, that the additional dependencies
specified with packaging `add` command will be written automaticaly to
`Project.toml` file of your local COBREXA.jl clone, not to your global
environment. Activation is simple: when in the directory of the package, just
type the command into the packaging shell:

```julia
(v1.6) pkg> activate
```

### Step 2: Publish your changes

You are expected to make a fork of the main COBREXA.jl repository, and open a
pull request from that one to the `develop` branch of the main repository.
For creating the fork, just hit the "Fork" button on GitHub.

After that, change the directory to your repository and adjust the remotes:

```bash
$ cd ~/.julia/dev/COBREXA             # or any other directory, as needed
$ git remote rename origin upstream   # renames the origin (the main COBREXA.jl repo) to upstream
$ git remote add origin git@github.com:yourUsername/COBREXA.jl.git  # adds the link to your clone as new origin
$ git fetch origin                    # fetches the refs from your repo
```

In the above code, change `yourUsername` is your GitHub username.

When the renaming is done, start a new branch at `upstream/master`. In the code
snippet, substitute `yn` for your initials (Your Name here) and give the new
feature a better name than `somefeature`:
```bash
$ git checkout -b yn-somefeature origin/master
```

Commit any changes and features that you like to the new branch. When the
commits look complete to you, push the branch to your repository fork:

```bash
$ git push -u origin yn-somefeature
```

This makes your changes visible in your repository. After that, you can
navigate to [GitHub's pull request
page](https://github.com/LCSB-BioCore/COBREXA.jl/pulls), where you should
immediately see a big green button that helps you to create a pull request for
this branch. Read the section below for precise details and guidelines on
submitting the pull requests.

## How to submit a pull request (PR) with your modification/enhancement?

1. **Make a fork of the repository**, commit the modifications in a **separate
   branch** and push the branch to your fork.
2. Make a pull request where you describe the motivation and expected outcome
   for the users. Specifically, consider any possible incompatibilities, and the
   necessity to increment the version number after your changes are applied.
   Set the target branch to `develop`.
3. After submitting the pull request, verify that all status checks (tests,
   documentation) are passing. Make sure any new contribution is properly
   documented and tested (you may want to check with coverage tools, using
   `test --coverage` from the Julia packaging shell)

After you submitted a pull request, a label might be assigned that allows us
to track and manage issues and pull requests.

### Code culture and style recommendations

Follow basic rules for software maintainability and extensibility:
- Do not reimplement functionality that is available in other packages, unless
  the reimplementation is either trivial and short, or there is a grave need to
  do so because the other implementations are deficient in some manner.
- Try to keep the function names and interfaces consistent with ecosystem
  standards and the other functions in the package. Consistency reduces the
  amount of surprise on the user side, thus lowers the need to reach for
  documentation, and in turn makes the software much easier and faster to use.
- Code less. Shorter code is almost always better unless demonstrated
  otherwise, e.g. with a benchmark. Avoid repetitive boilerplate (there should
  be ways to generate it, if needed).
- Keep the functionality "open" and composable. In particular, avoid all
  unnecessarily opaque and leaky abstractions (common in object-oriented
  programming).
- Avoid producing lots of "informative" text side-output by default, unless
  that is what the user asked for.
- Adhere to the code formatting rules defined by
  [JuliaFormatter](https://github.com/domluna/JuliaFormatter.jl). We usually
  have a bot running that checks all PRs and reports whether the code is
  properly formatted.

Follow the common rules for making easily mergable and reviewable PRs:
- Create one PR for each logical "feature" you want to merge. If your change is
  more complex and contains multiple "stages", open multiple PRs.
- Keep the test coverage reasonably high.
- If you commit many small, partial changes in a PR, you may help us save
  energy by prefixing your commit names with `[skip ci]`, which deactivates the
  CI trigger on that commit. With each skipped CI, you save a few watt-hours of
  energy. Testing just the "final" commit of the pull-request branch is
  sufficient.

## For developers: What is the expected branch management/workflow?

The workflow is based on [GitLab
flow](https://docs.gitlab.com/ee/topics/gitlab_flow.html), i.e., a `develop`
branch with `feature` branches being merged into the `develop` branch, all
periodically merged to `master` branch. Depending on your access rights, you
may open the `feature` branch in this repository, on in your fork.

The guidelines can be summarized as such:

- when making a contribution, create one new branch and open one new PR for
  each new independent feature or bugfix
- do not push to another branch unless it is your own
- try to get a review before merging unless the change is trivial and
  non-impacting
- consider prefixing your branch names with your initials, so that one can
  easily see who owns which branch (e.g. `ad-somefeature` would be committed by
  Arthur Dent)
