# Emerald

## Model components
- `Emerald` Land Surface Model based on CliMA Land (Remastered)
- `Sapphire` Ocean Model based on CliMA Ocean (To be implemented)
- `Diamond` Atmosphere Model based on CliMA Atmos (To be implemented)

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
For the details about how to use Emerald, you may look into the `test/tutorial` folder, and subscribe YouTube channel [`Plant Form & Function`](https://www.youtube.com/channel/UCTxiti8ntVAiU4sJ5N2kTnQ).

## Change logs
| Version Tag | Linear Change History                     |
|:------------|:------------------------------------------|
| Note        | Historty in ResearchProjects              |
| A1          | Version used for initial run              |
| A2          | Use minimum beta instead of mean beta     |
|             | Use moving average optimal T from CLM5    |
| A3          | Add vertical Vcmax profile                |
| A4          | Extend PAR definition to 400-750 nm       |
| A5          | Fix SIF yield scaling bug                 |
|             | Add ePAR scenario                         |
|             | Add prescibed CI option                   |
| Note        | History in ClimaLand-0.1                  |
| A6          | Fix a bug in Medlyn model (coeff 1.6)     |
| A7          | Add quadratic colimitation to C3 J (0.7)  |
| Note        | History in ClimaLand-0.2                  |
| B1          | Use Land v0.2                             |
| B2          | Add temperature dependency for K_D (VJP)  |
| B3          | Add limits to dgdt                        |
|             | Fix a bug in critical flow (not used)     |
|             | Fix a bug in dif rad SIF (not in v0.1)    |
| Note        | Moving to Emerald from CliMA Land         |
| B4*         | Add SZA < 89 degree limitation            |
|             | Swtich to global run per time step        |
|             | Add LAI = 0 case for SPAC                 |
|             | Add state field to SPAC                   |
|             | Pass water and enerby budgets             |
|             | Scale beta based on roots with + flow     |
|             | Add root disconnection feature (LAI => 0) |
|             | Use exp functions in radiation extinction |
|             | Use δLAI per canopy layer (still uniform) |
|             | Add beta function type LIDF               |
|             | Fix a bug in Verhoef LIDF (not used)      |
|             | Add soil trace gas diffusion (e.g. H₂O)   |
| 2023-06-14  | Run soil only RT as well in global model  |
|||

`*` for current development branch
