"""
#### written by Aaron Watters 3/15/2018
# dc_network
Dendritic call regulatory network explorer

[![Binder](https://mybinder.org/badge.svg)](https://mybinder.org/v2/gh/flatironinstitute/dc_network/master)

This repository uses Jupyter notebooks and encrypted data.

You must provide a decryption key when running the `notebooks//Decrypt data.ipynb`
notebook  to
decrypt the data before viewing the data using the
`notebooks//Visualize network.ipynb` to view the network.

The repository is designed to use 
[`repo2docker`](https://repo2docker.readthedocs.io/en/latest/)
to build a `docker` container which includes
the code for running the analysis and all needed dependencies.
The `repo2docker` tool requires the `docker` infrastructure
and Python 3 to run.  See the 
[installation instructions](https://repo2docker.readthedocs.io/en/latest/install.html).

To build the docker image from the command line

```bash
jupyter-repo2docker --no-run \
    --image-name dc_network \
    https://github.com/flatironinstitute/FREYA
```

To run the image after it has been built use `docker run`:

```bash
docker run -it --rm -p 8888:8888 dc_network:latest
```

"""