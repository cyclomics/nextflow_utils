FROM nfcore/base:2.1

WORKDIR /app
COPY environment.yml /tmp/environment.yml
RUN conda env update -n base -f /tmp/environment.yml

# Add code form the cycas git submodule 
COPY Cycas/cycas/ /cycas
