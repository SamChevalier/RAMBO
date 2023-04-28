# RAMBO
Repository for the paper [Scalable Bilevel Optimization for Generating Maximally Representative OPF Datasets](https://arxiv.org/abs/2304.10912), submitted the ISGT Europe 2023.

This repository collects maximaly representative OPF data using ```JuMP.jl``` and a custom extension of ```PowerModels.jl``` in ```Julia-1.8```.

Before running the code, make sure to activate the virtual environment from ```Project.toml```, e.g., by running

```
julia> ]
pkg> activate .
```

Instruction for running the bilevel solver are given in ```generate_OPF_points.jl```. Instructuins for generating randomly sampled benchmark datasets are given in ```generate_OPF_random_points.jl```.
