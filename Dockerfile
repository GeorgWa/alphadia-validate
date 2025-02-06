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


# we are installing alpharaw and alphadia from special hackathon branches
RUN git clone https://github.com/MannLabs/alpharaw.git \
    && cd alpharaw \
    && git checkout tmp_hackathon

RUN cd alpharaw && pip install .

RUN git clone https://github.com/MannLabs/alphadia.git \
    && cd alphadia \
    && git checkout tmp_hackathon

RUN git clone https://github.com/MannLabs/alphabase.git \
    && cd alphabase \
    && git checkout v1.4.2

# reduce the dependencies of alphadia and alphabase
COPY misc/requirements_alphabase.txt alphabase/requirements.txt
COPY misc/requirements_alphadia.txt alphadia/requirements/requirements.txt
COPY misc/requirements_alphadia.txt alphadia/requirements/requirements_loose.txt

RUN cd alphadia && pip install ".[stable]"

RUN pip install matplotlib progressbar altair seaborn
RUN pip install jupyter

# Make port 8888 available to the world outside this container
EXPOSE 8888

#RUN git clone https://github.com/GeorgWa/alphadia-validate.git adv

# copy over the code
COPY notebooks/src /app/notebooks/src
COPY notebooks/showcase.ipynb /app/notebooks
COPY misc /app/misc

ENV BASE_FOLDER=/app/base

RUN rm -r alphadia alphabase alpharaw

CMD ["jupyter", "notebook", "--ip=0.0.0.0", "--port=8888", "--no-browser", "--allow-root", \
"--NotebookApp.token=''", "--NotebookApp.password=''"]

#docker build -f Dockerfile  --progress=plain --build-arg="ALPHABASE_REF=latest" -t validiate .

# run bash:
# DATA_FOLDER=.
# docker run -v $DATA_FOLDER:/app/data/ -it validiate bash
