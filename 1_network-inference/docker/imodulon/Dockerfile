FROM python:3.12.2-alpine3.19 as build

RUN apk add --no-cache gcc g++ libc-dev

RUN python -m pip install 'numpy >= 1.19' \
                          'pandas >= 2.2, < 3' \
                          'scikit-learn == 1.4.0' \
                          'scipy >= 1.6.0'

# ---------------------------------------- #

FROM python:3.12.2-alpine3.19

COPY --from=build /usr/local/lib/python3.12/site-packages /usr/local/lib/python3.12/site-packages
COPY --from=build /usr/include/c++ /usr/include/c++
COPY --from=build /usr/lib /usr/lib

WORKDIR /root

ADD imodulon_execute3.py imodulon_execute.py
ADD imodulon_process.py imodulon_process.py

RUN mkdir /root/mount
VOLUME /root/mount

ENTRYPOINT ["python", "imodulon_execute.py"]
