---
name: Crash / error / performance problem report
about: Something that should work but does not, works very slowly, or produces wrong results.
---

## Minimal code example to reproduce the problem

> Include a reproducible sample of code that causes the problem. Others should
> be able to reproduce your problem by just copypasting your code into Julia
> shell. Use `download` to get any necessary datasets from the internet.
> Preferably, keep the sample code as minimalistic as possible.

```
using COBREXA
download("https://badmodels.org/model1.xml", "model.xml")
load_model("model.xml")
# ...
```

## Expected result

> Describe what you expected to happen if the above code worked correctly, and
> why you think that behavior is correct.
>
> If possible, copypaste the expected result directly to the issue description.

## Actual behavior

> What was the unexpected thing that happened instead?
>
> If possible, copypaste the error message or the wrong result directly into
> the issue.
>
> In case you are reporting a performance problem, include some timing
> information, using e.g. `@time` or `@benchmark`.

## Optional: Suggestions for fixing

> Feel free to suggest whatever change in the package that might remove the
> described problematic behavior. If you already have a fix, you may want to
> open a pull request instead of the issue.
>
> Erase the section if the cause of the problem is unclear.

## Optional: Environment

> If relevant to the problem, it is beneficial to provide all information about
> your environment, including the version of Julia. If possible, paste your
> package status here (type `]status` into Julia shell).
