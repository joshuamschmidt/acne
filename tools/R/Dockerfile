FROM rocker/r-ver:4.2.1
# image for data.table
# rocker/r-ver is based on debian

MAINTAINER Joshua Schmidt joshmschmidt1@gmail.com

RUN install2.r --error --skipinstalled --ncpus -1 \
  data.table \
  R.utils \
  && rm -rf /tmp/downloaded_packages\
  && strip /usr/local/lib/R/site-library/*/libs/*.so

WORKDIR /app
