FROM 812206152185.dkr.ecr.us-west-2.amazonaws.com/latch-base:02ab-main

# Its easy to build binaries from source that you can later reference as
# subprocesses within your workflow.
RUN apt-get update -y &&\
    apt-get install -y wget parallel ncbi-blast+

ENV CONDA_DIR /opt/conda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
	/bin/bash ~/miniconda.sh -b -p /opt/conda
ENV PATH=$CONDA_DIR/bin:$PATH
RUN conda install -c conda-forge -c bioconda mmseqs2

ENV VIRTUAL_ENV=/opt/venv
RUN python3 -m venv $VIRTUAL_ENV
ENV PATH="$VIRTUAL_ENV/bin:$PATH"

# Install dependencies:
COPY /data/requirements.txt .
RUN /opt/venv/bin/python3 -m pip install --upgrade pip && pip install -r requirements.txt

# You can use local data to construct your workflow image.  Here we copy a
# pre-indexed reference to a path that our workflow can reference.

# STOP HERE:
# The following lines are needed to ensure your build environement works
# correctly with latch.
WORKDIR /root
COPY data /root/data
COPY wf /root/wf
COPY wf/phamlite.py /root/wf/phamlite.py
ARG tag
ENV FLYTE_INTERNAL_IMAGE $tag
ENV LATCH_AUTHENTICATION_ENDPOINT https://nucleus.latch.bio
RUN  sed -i 's/latch/wf/g' flytekit.config
RUN python3 -m pip install --upgrade latch
