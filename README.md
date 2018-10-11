# Velocity Gradient Technique: a simple Demo

## What is Velocity Gradient Technique?



## Onilne Demo

The demo is hosted on our [laztech Github page](https://github.com/GindaChen/laztech-page/blob/master/VGT-demo/demo-VGT.ipynb) for Jupyter Notebook version. We also provide the Julia script for terminal user to run on the terminal environment.

## Run on Local Machine

*The simulation data used on the demo can be found here: *

### Prerequisite

1. **Julia**. VGT currently support [Julia 0.6](https://julialang.org/downloads/oldreleases.html). 
2. (Optional) **Jupyter Notebook**. After you download Julia, you could install Jupyter Notebook easily.
3. **Julia Packages**: `HDF5, PyCall, PyPlot`. This little script might be helpful:

```julia
# Package Install
# In your Julia terminal, type:
Pkg.add("HDF5");
Pkg.add("PyCall");
Pkg.add("PyPlot");

# You might also want to install Jupyter Notebook for the 
Pkg.add("IJulia")
```

To open Jupyter Notebook in your computer, simply type in your terminal:

```shell
jupyter-notebook
```

and navigate to the folder where you installed.



### Quick Start

1. Download the [VGT-demo](https://github.com/GindaChen/laztech-page/VGT-demo) folder to your local machine.
2. Open `demo-VGT.jl` if you're using Julia terminal or `demo-VGT.ipynb` if you're using Jupyter Notebook.
3. If you're using terminal, make sure you have the data cube in the same folder and change the corresponding data name in the code:

```julia
# Line 6~8 in demo-VGT.jl 
### Choose your cube here ###
dataname= "run-0-ppv.h5";   # <----- Put your cube name here
#---------------------------#
```

