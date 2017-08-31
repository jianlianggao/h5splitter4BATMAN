FROM ubuntu:16.04

MAINTAINER PhenoMeNal-H2020 Project ( phenomenal-h2020-users@googlegroups.com )


# Install R and Bruker2h5
RUN apt-get update && apt-get install -y --no-install-recommends r-base r-base-dev \
                              libcurl4-openssl-dev libssl-dev libhdf5-dev git && \
    echo 'options("repos"="http://cran.rstudio.com", download.file.method = "libcurl")' >> /etc/R/Rprofile.site && \
    R -e "install.packages(c('h5','devtools','getopt','optparse'))" && \
    R -e "library(devtools); install_github('jianlianggao/h5splitter4BATMAN/h5split4batman'); remove.packages('devtools')" && \
    apt-get purge -y r-base-dev git libcurl4-openssl-dev libssl-dev && \
    apt-get -y clean && apt-get -y autoremove && rm -rf /var/lib/{cache,log}/ /tmp/* /var/tmp/*


ADD runh5split4batman.R /usr/local/bin
RUN chmod 0755 /usr/local/bin/runh5split4batman.R

ADD runh5splitter.sh /usr/local/bin
RUN chmod 0755 /usr/local/bin/runh5splitter.sh


# Define entry point, useful for generale use
ENTRYPOINT ["runh5splitter.sh"]

