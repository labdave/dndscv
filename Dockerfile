FROM bioconductor/release_base2

MAINTAINER Slav <smk70@duke.edu>

# Helps clean up Docker images
RUN rm -rf /var/lib/apt/lists/*

# Retrieve scripts and support files from github
RUN git clone https://github.com/labdave/dndscv.git

# add repository to SYSPATH
ENV PATH /dndscv:$PATH

# change the permission of the repo
RUN chmod 777 -R /dndscv


# invalidates cache every 24 hours
ADD http://master.bioconductor.org/todays-date /dndscv/

RUN R -f /dndscv/install_dndscv.R

# END HERE





