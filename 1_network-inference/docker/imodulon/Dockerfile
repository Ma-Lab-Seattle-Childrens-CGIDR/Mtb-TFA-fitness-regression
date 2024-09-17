FROM avsastry/modulome-workflow:v1.0

USER root

RUN mkdir /root/mount
VOLUME /root/mount

WORKDIR /home/jovyan
ADD imodulon_execute2.py imodulon_execute2.py

ENTRYPOINT ["python3", "imodulon_execute2.py"]
