FROM ipython/scipystack

MAINTAINER Edward J. Kim <jkim575@illinois.edu>

VOLUME /notebooks
WORKDIR /notebooks

EXPOSE 8889

# You can mount your own SSL certs as necessary here
ENV PEM_FILE /key.pem
# $PASSWORD will get `unset` within notebook.sh, turned into an IPython style hash
ENV PASSWORD Dont make this your default
ENV USE_HTTP 0

ADD notebook.sh /
RUN chmod u+x /notebook.sh

RUN apt-get install -y -q \
    python-pyfits \
    python-mpi4py \
    python3-pyfits \
    python-mpi4py \
    libcfitsio3 \
    python-tables \
    python3-tables

RUN curl -O http://mirrors.kernel.org/ubuntu/pool/universe/h/healpix-cxx/libhealpix-cxx0_3.11.2-7_amd64.deb \
    && dpkg -i libhealpix-cxx0_3.11.2-7_amd64.deb

RUN curl -O http://mirrors.kernel.org/ubuntu/pool/universe/h/healpy/python-healpy_1.7.4-2_amd64.deb \
    && dpkg -i python-healpy_1.7.4-2_amd64.deb

RUN curl -O https://pypi.python.org/packages/source/M/MLZ/MLZ-1.1.tar.gz \
    && tar xvzf MLZ-1.1.tar.gz \
    && cd MLZ-1.1 \
    && python setup.py install

CMD ["/notebook.sh"]
