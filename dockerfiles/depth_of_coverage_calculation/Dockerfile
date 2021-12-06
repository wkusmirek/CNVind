FROM biodatageeks/cnv-opt-codexcov
MAINTAINER Wiktor Kuśmirek <kusmirekwiktor@gmail.com>

RUN Rscript -e "install.packages('CODEXCOV', repos = 'http://zsibio.ii.pw.edu.pl/nexus/repository/r-all')"

ADD count_coverage_for_single_sample_by_CODEX.R /opt
