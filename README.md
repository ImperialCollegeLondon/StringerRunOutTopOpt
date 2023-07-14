[![DOI](https://zenodo.org/badge/584917542.svg)](https://zenodo.org/badge/latestdoi/584917542)
[![](https://img.shields.io/badge/DOI-10.1016/j.compstruct.2023.117314-green)](https://doi.org/10.1016/j.compstruct.2023.117314)

# StringerRunOutTopOpt

:loudspeaker: Multi-scale topology optimisation method based on a continuum shell formulation with floating node method (FNM) and level-set method (LSM) capabilities for the design of stringer run-out regions in the presence of kissing bond damage.

## Dependencies

- Abaqus 2021
- Intel Fortran 19.1

## Usage


```
abaqus interactive job=_ input=sro_wing_fnm_disbond.inp global=wing_global_abq21.odb double
```

## Manuscript

The manuscript is available at [https://doi.org/10.1016/j.compstruct.2023.117314](https://doi.org/10.1016/j.compstruct.2023.117314).


## Authors  

Rui O. S. S. da Costa ([r.costa18@imperial.ac.uk](mailto:r.costa18@imperial.ac.uk))  
Silvestre T. Pinho ([silvestre.pinho@imperial.ac.uk](mailto:silvestre.pinho@imperial.ac.uk))

Department of Aeronautics  
Imperial College London  
South Kensington Campus  
SW7 2AZ, London  
United Kingdom

## License

This project is licensed under the BSD-3-Clause License - see the [LICENSE](LICENSE) file for details.

To use the code, please cite the repository according to [CITATION.cff](CITATION.cff) and the [manuscript](https://doi.org/10.1016/j.compstruct.2023.117314).
