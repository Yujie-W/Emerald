# Emerald

**E**arth **M**odeling with **E**nhanced **R**epresentations of **A**reo-Aqua and **L**and **D**ynamics

## Model components
- `EmeraldLand` Land Surface Model based on CliMA Land (Remastered)
- `EmeraldOcean` Ocean Model (To be implemented)

## Installation
`Emerald` is not a registered package, to use `Emerald`, one need to install it manually through Julia REPL:
```julia
using Pkg;
Pkg.add(PackageSpec(url="https://github.com/Yujie-W/Emerald.git"));
```

To use more testing features, you may install a specific branch, say `wyujie`
```julia
using Pkg;
Pkg.add(PackageSpec(url="https://github.com/Yujie-W/Emerald.git", rev="wyujie"));
```

## Tutorials
For the details about how to use Emerald, you may look into the `test/tutorial` folder, and subscribe [`Bilibili`](https://space.bilibili.com/22576885) and [`Bilibili Live`](http://live.bilibili.com/14243417).
