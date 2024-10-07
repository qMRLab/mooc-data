FROM jupyter/base-notebook:8ccdfc1da8d5

USER root

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        build-essential=12.4ubuntu1 \
        emacs \
        git \
        inkscape \
        jed \
        libsm6 \
        libxext-dev \
        libxrender1 \
        lmodern \
        netcat \
        unzip \
        nano \
        curl \
        wget \
        gfortran \
        cmake \
        libjpeg-dev \
        bsdtar \
        rsync \
        imagemagick \
        gnuplot-x11 \
        libopenblas-base \
        octave \
        liboctave-dev \
        octave-info \
        octave-parallel \
        octave-struct \
        octave-io \
        octave-statistics \
        octave-optim \
        octave-image \
        python3-dev \
        ttf-dejavu && \
    apt-get clean && \
    apt-get autoremove && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN cd $HOME/work;\
    pip install scipy \
                plotly \
                dash \
                dash_core_components \
                dash_html_components \
                dash_dangerously_set_inner_html \
                dash-renderer \
                sh \
                flask \
                oct2py \
                ipywidgets \
                Markdown \
                nbconvert \
                jupyterlab \
                jupytext\
                pandas \
                numpy \
                datascience \
                folium \
                random2 \
                matplotlib \
                sklearn \
                nilearn ; \
    git clone https://github.com/jvelazquez-reyes/qMT_tutorial-ISMRM2022.git;           \
    cd qMT_tutorial-ISMRM2022;\
    git clone https://github.com/neuropoly/qMRLab.git; \
    cd qMRLab; \
    git checkout d15a553f9d93457c3ed59861380852c54458c2b4; \
    cd ..; \
    chmod -R 777 $HOME/work/qMT_tutorial-ISMRM2022; \
    octave --eval "cd qMRLab; \
                      startup; \
                      pkg list;"

WORKDIR $HOME/work/qMT_tutorial-ISMRM2022

USER $NB_UID
