FROM python:3.11.4

MAINTAINER Joshua Schmidt joshmschmidt1@gmail.com

WORKDIR /usr/src/app

COPY requirements.txt ./

RUN pip install --no-cache-dir -r requirements.txt

COPY . .

WORKDIR /data

CMD ["/bin/bash"]

