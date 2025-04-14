# adapted from https://github.com/binder-examples/minimal-dockerfile
FROM python:3.11
# install the notebook package
RUN pip install --no-cache --upgrade pip && \
    pip install --no-cache notebook jupyterlab

RUN printf '[install]\ncompile = no\n[global]\nno-cache-dir = True' >> /etc/pip.conf

# create user with a home directory
ARG NB_USER
ARG NB_UID
ENV USER ${NB_USER}
ENV HOME /home/${NB_USER}

RUN adduser --disabled-password \
    --gecos "Default user" \
    --uid ${NB_UID} \
    ${NB_USER}
WORKDIR ${HOME}
USER ${USER}

# from here it's the same as the Dockerfile in the root folder

COPY misc/clean_reqs.sh .

RUN git clone https://github.com/MannLabs/alphadia.git \
    && cd alphadia \
    && git checkout slimmer_utils

RUN git clone https://github.com/MannLabs/alphabase.git \
    && cd alphabase \
    && git checkout v1.6.2

# reduce the dependencies of alphadia and alphabase for a smaller docker image
RUN cd alphabase && ../clean_reqs.sh requirements/requirements.txt psutil scikit-learn pydivsufsort pyarrow
RUN cd alphabase && pip install "."

RUN cd alphadia && ../clean_reqs.sh requirements/_requirements.freeze.txt alphabase alphatims peptdeep rocket_fft tokenizers huggingface-hub torchmetrics transformers directlfq pythonnet zstandard torch nvidia triton
RUN cd alphadia && ../clean_reqs.sh requirements/requirements.txt         alphabase alphatims peptdeep rocket_fft tokenizers huggingface-hub torchmetrics transformers directlfq pythonnet zstandard torch nvidia triton


RUN cd alphadia && pip install ".[stable]"
RUN pip install matplotlib progressbar altair seaborn
RUN pip install jupyter


# copy over the code
COPY notebooks/src /app/notebooks/src
COPY notebooks/showcase.ipynb /app/notebooks
COPY misc /app/misc

ENV BASE_FOLDER=/app/base

RUN rm -r alphadia alphabase

## change owner to be able to write to that folder
#COPY --chown=${USER}:${USER} example_nbs example_nbs

# build & debug:
# docker build --progress=plain --build-arg="NB_USER=user" --build-arg="NB_UID=1000" -f .binder/Dockerfile -t alphaquant-binder .
# docker run -it alphaquant-binder bash
