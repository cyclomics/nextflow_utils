docker build . -t cyclomics/cyclomicsseq-dev:0.1.1
docker tag cyclomics/cyclomicsseq-dev:0.1.1 cyclomics/cyclomicsseq-dev:latest
docker push cyclomics/cyclomicsseq-dev:0.1.1
docker push cyclomics/cyclomicsseq-dev:latest
