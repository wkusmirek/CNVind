FROM ubuntu:xenial
MAINTAINER Wiktor Kuśmirek <kusmirekwiktor@gmail.com>

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update  && \
    apt-get install -y software-properties-common  && \
    apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9  && \
    add-apt-repository 'deb [arch=amd64,i386] https://cran.rstudio.com/bin/linux/ubuntu xenial/'  && \
    apt-get install -y apt-transport-https

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y r-base libssl-dev libssh2-1-dev libxml2-dev libcurl4-openssl-dev libpq-dev wget

RUN Rscript -e "install.packages(c('roxygen2','devtools','DBI','getopt'))"

RUN Rscript -e "source('https://bioconductor.org/biocLite.R');biocLite('CODEX')"

RUN wget https://cran.r-project.org/src/contrib/Archive/optparse/optparse_1.4.4.tar.gz && \
    Rscript -e "install.packages('optparse_1.4.4.tar.gz', repos=NULL, type='source')"

RUN Rscript -e "install.packages('CODEXCOV', repos = 'http://zsibio.ii.pw.edu.pl/nexus/repository/r-all', dependencies=FALSE)"

ADD run_CODEXCOV_from_Yhat.R /opt
