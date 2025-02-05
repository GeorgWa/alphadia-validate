#FROM jupyter/minimal-notebook:latest

# https://github.com/michaelosthege/pythonnet-docker
FROM --platform=linux/amd64 mosthege/pythonnet:python3.10.10-mono6.12-pythonnet3.0.1

# Prevents Python from writing pyc files.
ENV PYTHONDONTWRITEBYTECODE=1

# Keeps Python from buffering stdout and stderr to avoid situations where
# the application crashes without emitting any logs due to buffering.
ENV PYTHONUNBUFFERED=1


WORKDIR /app

RUN git clone https://github.com/MannLabs/alphadia.git \
    && cd alphadia \
    && git checkout main

RUN --mount=type=cache,target=/root/.cache/pip \
    cd alphadia && pip install ".[stable]"

#RUN #cd alphadia && pip install -r requirements/requirements.txt
# TODO fix the alphadia version


RUN pip install jupyter

# Make port 8888 available to the world outside this container
EXPOSE 8888

COPY notebooks/showcase /app/notebooks/showcase
COPY notebooks/showcase.ipynb /app/notebooks
# add all other directories here
COPY notebooks/xic /app/notebooks/xic
COPY notebooks/magnus_utils /app/notebooks/magnus_utils
COPY notebooks/mirror_plotting.py /app/notebooks



CMD ["jupyter", "notebook", "--ip=0.0.0.0", "--port=8888", "--no-browser", "--allow-root", \
"--NotebookApp.token=''", "--NotebookApp.password=''"]

#docker build -f Dockerfile  --progress=plain --build-arg="ALPHABASE_REF=latest" -t validiate .

# run bash:
# DATA_FOLDER=.
# docker run -v $DATA_FOLDER:/app/data/ -it validiate bash
