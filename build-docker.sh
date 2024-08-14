DOCKER_TAG="0.1.1"
docker build . -t cyclomics/cyclomicsseq-dev:${DOCKER_TAG}
docker tag cyclomics/cyclomicsseq-dev:${DOCKER_TAG} cyclomics/cyclomicsseq-dev:latest
docker push cyclomics/cyclomicsseq-dev:${DOCKER_TAG}
docker push cyclomics/cyclomicsseq-dev:latest
