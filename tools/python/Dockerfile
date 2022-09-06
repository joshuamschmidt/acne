FROM python:3.10.6

MAINTAINER Joshua Schmidt joshmschmidt1@gmail.com

WORKDIR /usr/src/app

COPY requirements.txt ./

RUN pip install --no-cache-dir -r requirements.txt

COPY . .

WORKDIR /data

CMD ["/bin/bash"]

COPY ../../bin/makeBAF.py /usr/local/bin

COPY ../../bin/splitGS.py /usr/local/bin

RUN chmod +x /usr/local/bin/makeBAF.py
RUN chmod +x /usr/local/bin/splitGS.py