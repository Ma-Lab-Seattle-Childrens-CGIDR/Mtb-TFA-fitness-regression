FROM avsastry/modulome-workflow:v1.0

USER root

RUN apt-get update
RUN apt-get -y install datamash

USER $NB_USER

RUN mkdir /home/$NB_USER/mount
VOLUME /home/$NB_USER/mount

ADD imodulon_execute.sh imodulon_execute.sh

ENTRYPOINT ["./imodulon_execute.sh"]
