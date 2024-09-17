FROM python:3.12.2-alpine3.19 as build

RUN apk add --no-cache gcc g++ libc-dev

RUN python -m pip install 'joblib >= 1.2.0' \
                          'numpy >= 1.19' \
                          'pandas >= 2.2, < 3' \
                          'scikit-learn == 1.4.0' \
                          'scipy >= 1.6.0' \
                          'threadpoolctl >= 2.0.0'

# ---------------------------------------- #

FROM python:3.12.2-alpine3.19

COPY --from=build /usr/local/lib/python3.12/site-packages /usr/local/lib/python3.12/site-packages
COPY --from=build /usr/include/c++ /usr/include/c++
COPY --from=build /usr/lib /usr/lib

WORKDIR /root

ADD linear_model.py linear_model.py

RUN mkdir /root/mount
VOLUME /root/mount

ENTRYPOINT ["python", "linear_model.py"]
