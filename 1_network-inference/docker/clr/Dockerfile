FROM rhub/r-minimal:4.4.0

RUN apk add --no-cache libgomp

RUN installr -d optparse parmigene tidyr larsgr/BSplineMI

WORKDIR /root

ADD clr_execute.R clr_execute.R

ENTRYPOINT ["Rscript", "--vanilla", "clr_execute.R"]

CMD ["--help"]

# set number of threads to use at runtime like `--env 8`
