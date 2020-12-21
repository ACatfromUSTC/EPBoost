FROM python:3.6.8

RUN apt update && \ 
    apt install -y gcc python3-dev python3-pip libxml2-dev libxslt1-dev zlib1g-dev g++ git cmake build-essential vim bedtools
RUN pip install pip -U && \ 
    pip config set global.index-url https://pypi.tuna.tsinghua.edu.cn/simple && \
    pip install numpy pandas scikit-learn scikit-image matplotlib tqdm catboost