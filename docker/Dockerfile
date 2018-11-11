FROM ubuntu:16.04

LABEL name="coolbox"
LABEL base.image="ubuntu:16.04"
LABEL version="1"
LABEL software="CoolBox"
LABEL software.version="latest"
LABEL website="https://github.com/Nanguage/CoolBox"
LABEL documentation="https://github.com/Nanguage/CoolBox/wiki"
LABEL license="https://github.com/Nanguage/CoolBox/blob/master/LICENSE"
LABEL tags="Bioinformatics,Genomics,Hi-C,Visualization"

MAINTAINER nanguage@yahoo.com

ENV MINICONDA https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh

ENV LD_ALL C.UTF-8
ENV LANG C.UTF-8

# add apt source and install requirement packages
RUN mv /etc/apt/sources.list /etc/apt/sources.list.bkp && \
    bash -c 'echo -e "deb mirror://mirrors.ubuntu.com/mirrors.txt xenial main restricted universe multiverse\n\
deb mirror://mirrors.ubuntu.com/mirrors.txt xenial-updates main restricted universe multiverse\n\
deb mirror://mirrors.ubuntu.com/mirrors.txt xenial-backports main restricted universe multiverse\n\
deb mirror://mirrors.ubuntu.com/mirrors.txt xenial-security main restricted universe multiverse\n\n" > /etc/apt/sources.list' && \
    cat /etc/apt/sources.list.bkp >> /etc/apt/sources.list && \
    cat /etc/apt/sources.list
RUN apt-get clean all && \
    apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y  \
        autotools-dev   \
        automake        \
        cmake           \
        curl            \
        grep            \
        sed             \
        dpkg            \
        fuse            \
        git             \
        wget            \
        zip             \
        openjdk-8-jre   \
        build-essential \
        pkg-config      \
        python          \
	python-dev      \
        python-pip      \
        bzip2           \
        ca-certificates \
        libglib2.0-0    \
        libxext6        \
        libsm6          \
        libxrender1     \
        git             \
        mercurial       \
        subversion      \
        zlib1g-dev &&   \
        apt-get clean && \
        apt-get purge && \
        rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# install miniconda
RUN echo 'export PATH=/opt/conda/bin:$PATH' > /etc/profile.d/conda.sh && \
    wget --quiet ${MINICONDA} -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh

ENV PATH=/opt/conda/bin:$PATH

# add channels
RUN conda config --add channels r
RUN conda config --add channels bioconda
RUN conda upgrade conda

# install dependency
RUN conda install --yes pybigwig
RUN pip install --upgrade pip
RUN pip install jupyter ipywidgets
RUN jupyter nbextension enable --py widgetsnbextension

# install DLO-HiC-Tools
RUN pip install coolbox


CMD ["/bin/bash"]

WORKDIR /data
