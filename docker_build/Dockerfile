# Define the build argument for GPU support
ARG BASE_IMAGE="nvidia/cuda:12.0.1-devel-ubuntu22.04"
#ARG GPU_SUPPORT=false

# Set the base image depending on GPU_SUPPORT
ARG ROOT_CONTAINER
ARG CUDA_VERSION=12.0.1
ARG UBUNTU_VERSION=22.04
FROM ${BASE_IMAGE}

ARG NB_USER="jovyan"
ARG NB_UID="1000"
ARG NB_GID="100"

# Fix: https://github.com/hadolint/hadolint/wiki/DL4006
# Fix: https://github.com/koalaman/shellcheck/wiki/SC3014
SHELL ["/bin/bash", "-o", "pipefail", "-c"]

USER root

# Install all OS dependencies for notebook server that starts but lacks all
# features (e.g., download as all possible file formats)
ENV DEBIAN_FRONTEND noninteractive
RUN apt-get update --yes && \
    # - apt-get upgrade is run to patch known vulnerabilities in apt-get packages as
    #   the ubuntu base image is rebuilt too seldom sometimes (less than once a month)
    apt-get upgrade --yes && \
    apt-get install --yes --no-install-recommends \
    # - bzip2 is necessary to extract the micromamba executable.
    bzip2 \
    ca-certificates \
    locales \
    sudo \
    # - tini is installed as a helpful container entrypoint that reaps zombie
    #   processes and such of the actual executable we want to start, see
    #   https://github.com/krallin/tini#why-tini for details.
    tini \
    wget && \
    apt-get clean && rm -rf /var/lib/apt/lists/* && \
    echo "en_US.UTF-8 UTF-8" > /etc/locale.gen && \
    locale-gen

# Configure environment
ENV CONDA_DIR=/opt/conda \
    SHELL=/bin/bash \
    NB_USER="${NB_USER}" \
    NB_UID=${NB_UID} \
    NB_GID=${NB_GID} \
    LC_ALL=en_US.UTF-8 \
    LANG=en_US.UTF-8 \
    LANGUAGE=en_US.UTF-8
ENV PATH="${CONDA_DIR}/bin:${PATH}" \
    HOME="/home/${NB_USER}"

# Copy a script that we will use to correct permissions after running certain commands
COPY fix-permissions /usr/local/bin/fix-permissions
RUN chmod a+rx /usr/local/bin/fix-permissions

# Enable prompt color in the skeleton .bashrc before creating the default NB_USER
# hadolint ignore=SC2016
RUN sed -i 's/^#force_color_prompt=yes/force_color_prompt=yes/' /etc/skel/.bashrc && \
   # Add call to conda init script see https://stackoverflow.com/a/58081608/4413446
   echo 'eval "$(command conda shell.bash hook 2> /dev/null)"' >> /etc/skel/.bashrc

# Create NB_USER with name jovyan user with UID=1000 and in the 'users' group
# and make sure these dirs are writable by the `users` group.
RUN echo "auth requisite pam_deny.so" >> /etc/pam.d/su && \
    sed -i.bak -e 's/^%admin/#%admin/' /etc/sudoers && \
    sed -i.bak -e 's/^%sudo/#%sudo/' /etc/sudoers && \
    useradd -l -m -s /bin/bash -N -u "${NB_UID}" "${NB_USER}" && \
    echo "${NB_USER} ALL=(ALL) NOPASSWD: ALL" > /etc/sudoers.d/${NB_USER} && \
    chmod 0440 /etc/sudoers.d/${NB_USER} && \
    mkdir -p "${CONDA_DIR}" && \
    chown "${NB_USER}:${NB_GID}" "${CONDA_DIR}" && \
    chmod g+w /etc/passwd && \
    fix-permissions "${HOME}" && \
    fix-permissions "${CONDA_DIR}"

ENV DEBIAN_FRONTEND=noninteractive
ENV GROMACS_VERSION=2023
ENV PLUMED_KERNEL=/usr/local/lib/libplumedKernel.so
# Update and install required dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    cmake \
    curl \
    ca-certificates \
    libfftw3-dev \
    libgomp1 \
    liblapack3 \
    openmpi-common \
    libopenmpi-dev \
    tmux \
    htop \
    nano \
    less \
    htop \
    openmpi-bin \
    git \
    vim \
    && rm -rf /var/lib/apt/lists/*

# Download and install PLUMED
RUN cd / \
    && curl -sSLO https://github.com/plumed/plumed2/releases/download/v2.9.0/plumed-2.9.0.tgz \
    && tar -xzf plumed-2.9.0.tgz \
    && cd /plumed-2.9.0 \
    && ./configure --prefix=/usr/local \
    && make -j$(nproc) \
    && make install \
    && cd /

# Download and install GROMACS with MPI support

RUN cd / \
    && curl -sSLO https://ftp.gromacs.org/gromacs/gromacs-${GROMACS_VERSION}.tar.gz \
    && tar -xvzf gromacs-${GROMACS_VERSION}.tar.gz \
    && cd /gromacs-${GROMACS_VERSION} \
    && . /plumed-2.9.0/sourceme.sh \
    && echo "4" | plumed patch -p \
    && cd / \
    && mkdir /gromacs-${GROMACS_VERSION}-build \
    && cd /gromacs-${GROMACS_VERSION}-build \
    && cmake /gromacs-${GROMACS_VERSION} \
        -DGMX_MPI=ON \
        -DGMX_GPU=CUDA \
	-DCUDAToolkit_ROOT=/usr/local/cuda-12.0 \
        -DCUDA_TOOLKIT_ROOT_DIR=/usr/local/cuda-12.0 \	
        -DGMX_DOUBLE=OFF \
        -DGMX_EXTERNAL_BLAS=OFF \
        -DGMX_EXTERNAL_LAPACK=OFF \
        -DCMAKE_INSTALL_PREFIX=/usr/local \
    && make -j$(nproc) \
    && make install \
    && mkdir /gromacs-${GROMACS_VERSION}-nompibuild \
    && cd /gromacs-${GROMACS_VERSION}-nompibuild \
    && cmake /gromacs-${GROMACS_VERSION} \
        -DGMX_MPI=OFF \
        -DGMX_GPU=CUDA \
        -DCUDAToolkit_ROOT=/usr/local/cuda-12.0 \
        -DCUDA_TOOLKIT_ROOT_DIR=/usr/local/cuda-12.0 \	
        -DGMX_DOUBLE=OFF \
        -DGMX_EXTERNAL_BLAS=OFF \
        -DGMX_EXTERNAL_LAPACK=OFF \
        -DCMAKE_INSTALL_PREFIX=/usr/local \
    && make -j$(nproc) \
    && make install \
    && ldconfig \
    && cd / \
    && rm -rf /gromacs-${GROMACS_VERSION}* \
    && rm -rf /gromacs-${GROMACS_VERSION}-build \
    && rm -rf /gromacs-${GROMACS_VERSION}-nompibuild; 

USER ${NB_UID}

# Pin python version here, or set it to "default"
ARG PYTHON_VERSION=3.10

# Setup work directory for backward-compatibility
RUN mkdir "/home/${NB_USER}/work" && \
    fix-permissions "/home/${NB_USER}"


# Install Python, Mamba
# Cleanup temporary files and remove Micromamba
# Correct permissions
# Do all this in a single RUN command to avoid duplicating all of the
# files across image layers when the permissions change
COPY --chown="${NB_UID}:${NB_GID}" initial-condarc "${CONDA_DIR}/.condarc"
WORKDIR /tmp
RUN set -x && \
    arch=$(uname -m) && \
    if [ "${arch}" = "x86_64" ]; then \
        # Should be simpler, see <https://github.com/mamba-org/mamba/issues/1437>
        arch="64"; \
    fi && \
    wget -qO /tmp/micromamba.tar.bz2 \
        "https://micromamba.snakepit.net/api/micromamba/linux-${arch}/latest" && \
    tar -xvjf /tmp/micromamba.tar.bz2 --strip-components=1 bin/micromamba && \
    rm /tmp/micromamba.tar.bz2 && \
    PYTHON_SPECIFIER="python=${PYTHON_VERSION}" && \
    if [[ "${PYTHON_VERSION}" == "default" ]]; then PYTHON_SPECIFIER="python"; fi && \
    # Install the packages
    ./micromamba install \
        --root-prefix="${CONDA_DIR}" \
        --prefix="${CONDA_DIR}" \
        --yes \
        "${PYTHON_SPECIFIER}" \
        'mamba' && \
    rm micromamba && \
    # Pin major.minor version of python
    mamba list python | grep '^python ' | tr -s ' ' | cut -d ' ' -f 1,2 >> "${CONDA_DIR}/conda-meta/pinned" && \
    mamba clean --all -f -y && \
    fix-permissions "${CONDA_DIR}" && \
    fix-permissions "/home/${NB_USER}"

# Configure container startup
ENTRYPOINT ["tini", "-g", "--"]
CMD ["start.sh"]

# Copy local files as late as possible to avoid cache busting
COPY start.sh /usr/local/bin/

# Switch back to jovyan to avoid accidental container runs as root
USER ${NB_UID}

WORKDIR "${HOME}/work"

# Download and install custom pmx fork
RUN mamba install MDAnalysis mdtraj acpype rdkit openbabel tqdm -y
RUN cd /home/jovyan \
    && git clone --branch abfe_dev https://github.com/deGrootLab/pmx.git \
    && cd pmx \
    && git checkout 3b135ad 

COPY --chown=jovyan:100 ./abfe_pipeline /home/${NB_USER}/work/
COPY --chown=jovyan:100 ./AbsoluteDG.py ./DoubleBox.py ./Restraints.py ./gmx.py /home/${NB_USER}/pmx/src/pmx

RUN cd /home/jovyan \
    && cd pmx \
    && python setup.py develop

# CMD ["-   -input_file", "none", "--ligand_name", "LIG", "--processes", "32", "--lig", "ligand", "--apoCase", "apo", "--replicas", "1"]
