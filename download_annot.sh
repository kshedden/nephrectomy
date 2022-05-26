#!/bin/bash

BASE=/home/kshedden/data/Markus_Bitzer

rm -rf $(BASE)/tmp
mkdir $(BASE)/tmp

dropbox_uploader.sh download "CSCAR-Nephrology Collaboration/Kerby .xml (2022)"
mv "Kerby .xml (2022)" tmp/k1

dropbox_uploader.sh download "CSCAR-Nephrology Collaboration/XML_updated"
mv "XML_updated" tmp/k3

rm -rf tmp/final
mkdir tmp/final

cp tmp/k3/*xml tmp/final
cp -u tmp/k1/*xml tmp/final

cd tmp/final
tar -cvf annotations.tar *xml
gzip -f annotations.tar
mv annotations.tar.gz $(BASE)/Annotations
