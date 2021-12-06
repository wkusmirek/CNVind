FROM biodatageeks/cnv-opt-codexcov
MAINTAINER Wiktor Kuśmirek <kusmirekwiktor@gmail.com>

RUN Rscript -e "install.packages('CODEXCOV', repos = 'http://zsibio.ii.pw.edu.pl/nexus/repository/r-all')"

ADD run_CODEXCOV_from_Yhat.R /opt
