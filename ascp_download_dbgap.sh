#!/usr/bin/bash

cd /fh/fast/stanford_j/Xiaoyu/QTL/data/TBD

/app/bin/ascp -QTr -l 300M -k 1 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh -W A7D295310CA25469A9EB7025F8675A34BA6AFD9A8EFA1427383709308585BD6871F4AFECC1D7F4F69A25A45A41416B229D dbtest@gap-upload.ncbi.nlm.nih.gov:data/instant/jamesdai/60547 .

decrypt the files:
1.Configure the tookit: /fh/fast/dai_j/CancerGenomics/Tools/sratoolkit.2.8.2-1-ubuntu64/bin/vdb-config -i
In the toolkit configuration console (opened in Step 3), import the repository key (dbGaP download page), and create a project file workspace. see an example in /fh/fast/stanford_j/Xiaoyu/QTL/data/TBD/60547/vdb_config.png
2. In that workspace you can decrypt the files by using vdb-dycrypt from the toolkit: 
cd /fh/fast/stanford_j/Xiaoyu/QTL/data/TBD/60547/data/dbGaP-16806
/fh/fast/dai_j/CancerGenomics/Tools/sratoolkit.2.8.2-1-ubuntu64/bin/vdb-decrypt ../../../60547

