Bootstrap: docker
From: ubuntu

%files
    environment.yml /root/
    scripts /root/

%appenv preprocessing
    export PATH=$PATH:/var/cellranger-5.0.1/bin

%appinstall preprocessing
    apt update && apt install wget bash -y
    cd /var/
    wget -O cellranger-5.0.1.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-5.0.1.tar.gz?Expires=1614525742&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cHM6Ly9jZi4xMHhnZW5vbWljcy5jb20vcmVsZWFzZXMvY2VsbC1leHAvY2VsbHJhbmdlci01LjAuMS50YXIuZ3oiLCJDb25kaXRpb24iOnsiRGF0ZUxlc3NUaGFuIjp7IkFXUzpFcG9jaFRpbWUiOjE2MTQ1MjU3NDJ9fX1dfQ__&Signature=kv50KOHd9EG5ybOFPLgGcyESYWcv8u0j8uN22NVjTPz3~lKlxmIKid2KCsjb8t2Q5Vw7t7NxrFXk55zGWUMKHOm9b68Nn7tbREPYLksFW5e6Cg9AE3rS4iF2U9P94a940NyEJRp4yU1p9BLNnc2ljUfLTAusgZdsYmc5O2luyKq-lEsefg9VyLbmV7jAqD5~TBEtsB3O8sl4t4xiqaWghB8IYkFgN~txyO~2Ji6tD7vBAJ7FicQVM77DAEt218Z4LMklvXFD5pSUd7Z~t2fI5ral6VXmnHQyErrspSrSVYUJutJB0vZmqPLlBllcpcia80p-0tt50B5WKnCMaBVz2w__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA"
    wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz
    tar xvf cellranger-5.0.1.tar.gz
    tar xvf refdata-gex-GRCh38-2020-A.tar.gz
    rm -f *.gz
    echo "export PATH=$PATH:/var/cellranger-5.0.1/bin" >> ~/.bashrc

%appenv downstream
    export PATH=/root/miniconda/bin:$PATH

%appinstall downstream
    cd /root/
    export PATH=/root/miniconda/bin:$PATH
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
        /bin/bash Miniconda3-latest-Linux-x86_64.sh -b -p /root/miniconda && \
        rm -f Miniconda3-latest-Linux-x86_64.sh && \
        conda env update --file environment.yml --name base

%apprun preprocessing
    bash scripts/run_cellranger.sh $@

%apprun downstream
    python scripts/firstline_downstream_analysis.py $@