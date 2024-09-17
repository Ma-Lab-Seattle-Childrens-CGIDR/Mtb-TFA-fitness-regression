FROM r-minimal:devel

RUN installr -d doParallel doRNG GENIE3 optparse

WORKDIR /root

RUN mkdir /root/mount
VOLUME /root/mount

ADD genie3_execute.R genie3_execute.R

ENTRYPOINT ["Rscript", "--vanilla", "genie3_execute.R"]
CMD ["--help"]

# Run with `docker run -it --volume ./mount:/root/mount genie3 mount/expression.tsv mount/regulators.tsv mount/output.txt`
# `mount` should be a directory, and should contain the relevant input files,
# which must also be provided to the script.
# Script output will also be written to the `mount` directory, at the location specified.
# Run `docker run -it genie3 --help` to see additional command line options
