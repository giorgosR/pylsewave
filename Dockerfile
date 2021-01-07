## ------------------------------------------------------------------ ##
##  Dockerfile for creating env with jupyter notebooks for pylsewave  ##
## ------------------------------------------------------------------ ##
FROM ubuntu:18.04

ENV PATH="/root/miniconda3/bin:${PATH}"
ARG PATH="/root/miniconda3/bin:${PATH}"

RUN apt-get update

RUN apt-get install -y wget && rm -rf /var/lib/apt/lists/*

RUN wget \
    https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && mkdir /root/.conda \
    && bash Miniconda3-latest-Linux-x86_64.sh -b \
    && rm -f Miniconda3-latest-Linux-x86_64.sh 

RUN conda --version

#RUN mkdir /usr/giorgosr && cd /junbs

WORKDIR /usr/giorgosr/junbs

COPY ./JupyterNbs/*.ipynb ./

RUN mkdir ./results

COPY ./.conda/pylsewave.yml ./pylsewave.yml

ARG conda_env=pylsewave

RUN conda env create -n ${conda_env} -f pylsewave.yml

# RUN echo ". /root/miniconda3/etc/profile.d/conda.sh" >> ~/.bashrc && \ 
#     echo "conda activate $(head -1 ${environment} | cut -d' ' -f2)" >> ~/.bashrc

# Pull the environment name out of the environment.yml
RUN echo "source activate $(head -1 ./pylsewave.yml | cut -d' ' -f2)" >> ~/.bashrc
ENV PATH /root/miniconda3/envs/$conda_env/bin:$PATH
ENV CONDA_DEFAULT_ENV $conda_env

# # ENV TINI_VERSION v0.6.0
# # ADD https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini /usr/bin/tini
# # RUN chmod +x /usr/bin/tini
# # ENTRYPOINT ["/usr/bin/tini", "--"]

EXPOSE 8888

ENTRYPOINT ["jupyter", "notebook", "--port=8888", "--no-browser", "--ip=0.0.0.0", "--allow-root"]
