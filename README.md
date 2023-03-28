# Emerald

## Model components
- `Emerald` Land Surface Model based on CliMA Land (Remastered)
- `Sapphire` Ocean Model based on CliMA Ocean (To be implemented)
- `Diamond` Atmosphere Model based on CliMA Atmos (To be implemented)

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
| 2023-Mar-28 | Add root disconnection feature (LAI => 0) |
|||

`*` for current development branch
