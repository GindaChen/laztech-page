# Velocity Gradient Technique: a simple Demo

*If you're looking a demo page, come [here](https://github.com/GindaChen/laztech-page/blob/master/VGT-demo/demo-VGT.ipynb) . The cube we're using in the demo is a bit too large for file transfer... We'll still post it on site.*

### 

## Run on Local Machine

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
3. 