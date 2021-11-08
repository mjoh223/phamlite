#!/usr/bin/bash
wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.12.0+-x64-linux.tar.gz
tar xvzf ncbi-blast-2.12.0+-x64-linux.tar.gz
cp ncbi-blast-2.12.0+/bin/* /app/.apt/usr/bin
