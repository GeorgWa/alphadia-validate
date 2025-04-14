FROM python:3.11

# Prevents Python from writing pyc files.
ENV PYTHONDONTWRITEBYTECODE=1

# Keeps Python from buffering stdout and stderr to avoid situations where
# the application crashes without emitting any logs due to buffering.
ENV PYTHONUNBUFFERED=1

WORKDIR /app

# prepare pip
RUN pip install --upgrade pip
RUN printf '[install]\ncompile = no\n[global]\nno-cache-dir = True' >> /etc/pip.conf

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
RUN apt-get clean && \
  rm -rf /var/lib/apt/lists/*

# Make port 8888 available to the world outside this container
EXPOSE 8888

CMD ["jupyter", "notebook", "--ip=0.0.0.0", "--port=8888", "--no-browser", "--allow-root", \
"--NotebookApp.token=''", "--NotebookApp.password=''"]
#
##docker build -f Dockerfile  --progress=plain  --build-arg="ALPHABASE_REF=latest" -t validiate .
#
## run bash:
## DATA_FOLDER=.
## docker run -v $DATA_FOLDER:/app/data/ -it validiate bash
