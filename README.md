# alphadia-validate
A set of tools to validate DIA data.


## Run locally
1. create a conda environment, install AlphaDIA (cf. https://github.com/MannLabs/alphadia?tab=readme-ov-file#pip-installation)
```bash
pip install "alphadia[stable]==1.9.3-dev2""
```
2. (non-Windows) install `mono`
3. Install jupyter
4. open `notebooks/showcase.ipynb` and `notebooks/shared_fragment_histogram.ipynb`
Please ignore all the other notebooks, they were work in progress.


## Run in Docker
This is a convenient solution to try this out without any modication to your system (except for Docker).


1. Install [Docker](https://docs.docker.com/engine/install/ubuntu/).
2. Create the following directories:
```bash
mkdir -p ~/alphadia-validate/data/data
mkdir -p ~/alphadia-validate/data/output
```
2. Download the raw data (`dia_data` object) as a `pkl` file [here](https://datashare.biochem.mpg.de/s/pckjZUEBChOvA9v)
to `~/data/alphadia-validate/data/output`.
This is currently required because we didn't get mono to run in the container.

3. Clone this repository, `cd` into it and build the container
```bash
docker build -f Dockerfile --progress=plain -t alphadia-validiate .
```

4. Run the container
```bash
BASE_FOLDER=~/alphadia-validate
docker run -p 8887:8888 -v $BASE_FOLDER:/app/base/ -it alphadia-validiate
```

5. Open the Jupyter notebook in your browser: http://localhost:8887/notebooks/notebooks/showcase.ipynb
