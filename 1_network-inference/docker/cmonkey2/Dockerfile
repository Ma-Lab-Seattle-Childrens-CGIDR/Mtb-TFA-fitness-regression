FROM weiju/cmonkey2:latest

RUN apt-get update
RUN apt-get install -y datamash sqlite3

RUN mkdir -p /cmonkey/mount
RUN chmod 777 /cmonkey/ /cmonkey/mount
VOLUME /cmonkey/mount

WORKDIR /cmonkey

ADD cmonkey_execute.sh cmonkey_execute.sh
ADD cmonkey_process.sh cmonkey_process.sh
ADD cmonkey_process.py cmonkey_process.py

ENTRYPOINT ["./cmonkey_execute.sh"]
