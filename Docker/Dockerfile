FROM ubuntu:16.04

MAINTAINER PhenoMeNal-H2020 Project ( phenomenal-h2020-users@googlegroups.com )


# Install R and Bruker2h5
RUN apt-get update && apt-get install -y --no-install-recommends r-base r-base-dev \
                              libcurl4-openssl-dev libssl-dev libhdf5-dev git && \
    echo 'options("repos"="http://cran.rstudio.com", download.file.method = "libcurl")' >> /etc/R/Rprofile.site && \
    R -e "install.packages(c('h5','devtools', 'getopt','optparse'))" && \
    R -e "library(devtools); install_github('jianlianggao/h5splitter4BATMAN/'); remove.packages('devtools')" && \
    apt-get purge -y r-base-dev git libcurl4-openssl-dev libssl-dev && \
    apt-get -y clean && apt-get -y autoremove && rm -rf /var/lib/{cache,log}/ /tmp/* /var/tmp/*

ADD runBruker2h5.R /usr/local/bin
RUN chmod 0755 /usr/local/bin/runBruker2h5.R

ADD runTest1.sh /usr/local/bin
RUN chmod 0755 /usr/local/bin/runTest1.sh


# Define entry point, useful for generale use
ENTRYPOINT ["runTest1.sh"]

