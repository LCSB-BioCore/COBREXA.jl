# Contributing to COBREXA.jl

:+1::tada: Thanks for taking the time to contribute to
[COBREXA.jl](https://github.com/LCSB-BioCore/COBREXA.jl)! :tada::+1:

## How to test a development version of the package?

You can use the "devel" feature of Julia packaging system. One possible way is here:

1. Clone the repository as usual, `cd` into the directory
2. Start `julia`, type `]` to switch to the packaging interface and type
   `develop --local .`
3. The package will be always recompiled and loaded from your local repository.

Alternatively, you can checkout the whole repository directly using Julia, by
typing:
```
develop https://github.com/LCSB-BioCore/COBREXA.jl
```

That will create the repo for you in `~/.julia/dev/`, which you can use for
development and pushing/pulling.


## How to report a bug or suggest an enhancement

Please use the [GitHub issue
tracker](https://github.com/LCSB-BioCore/COBREXA.jl/issues) to report any
problems with the software, and discuss any potential questions about COBREXA
use.

Before creating bug reports, please check [the open
issues](https://github.com/LCSB-BioCore/COBREXA.jl/issues) as you might find
out that you don't need to create one. When you are creating a bug report,
please include as many details as possible. Fill out the required template, the
information it asks for helps us resolve issues faster.

> If you find a Closed issue that seems like it is the same thing that
  you're experiencing, open a new issue and include a link to the original issue
  in the body of your new one.

If reporting issues, please do not forget to include a description of
conditions under which the issue occurs; preferably a code snippet that can be
run separately (sometimes termed "minimal crashing example").

## How to submit a pull request (PR) with your modification/enhancement?

1. Make a fork of the repository, commit the modifications in a separate branch.
2. Make a PR, follow instructions in the PR template (if available). Describe
   the motivation and expected outcome for the users. Specifically, consider
   any possible incompatibilities, and the necessity to increment the version
   number after your changes are applied.
3. Submit your pull request
4. Verify that all status checks (tests, documentation) are passing. Make sure
   any new contribution is properly documented and tested (you may want to
   check with coverage tools, using `test --coverage` from the Julia packaging
   shell)

After you submitted a pull request, a label might be assigned that allows us
to track and manage issues and pull requests.

## What is the workflow?

The workflow is based on [GitLab
flow](https://docs.gitlab.com/ee/topics/gitlab_flow.html), i.e., a `master`
branch with `feature` branches being merged into the `master` branch. Depending
on your access rights, you may open the `feature` branch in this repository, on
in your fork.

The guidelines can be summarized as such:

- when making a contribution, create one new branch and open one new PR for
  each new independent feature or bugfix
- do not push to another branch unless it is your own
- try to get a review before merging unless the change is trivial and
  non-impacting
- consider prefixing your branch names your initials (e.g. `ad-somefeature` by
  Arthur Dent)
