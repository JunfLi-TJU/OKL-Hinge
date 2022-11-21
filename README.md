# OKL-Hinge
Source codes of algorithms and datasets for our paper "Improved Kernel Alignment Regret Bound for Online Kernel Learning", accepted in AAAI 2023.

We implement all algorithms with R on a Windows machine with 2.8 GHz Core(TM) i7-1165G7 CPU. execute each experiment 10 times with random permutation of all datasets and average all of the results.

To run the code, you must set the paths following the code, or set new paths.
The default path of codes is "D:/experiment/Conference Paper/ECML/ECML2022". 
The path of datasets is "D:/experiment/AAAI2023/dataset". 
The store path is "D:/experiment/AAAI2023/Result/". 

The baseline algorithms include: OKS, RF-OKS, Raker and LKMBooks. Our algorithms include: OKS++, IOKS, RF-OKS++ and RF-IOKS.

The datasets are downloaded from: https://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/ and http://archive.ics.uci.edu/ml/datasets.php

binary classification datasets: magic04 (Num:19020, Fea: 10), phishing (Num: 11055, Fea: 68), a9a (Num: 32561, Fea: 123), SUSY (Num: 20000, Fea: 18).

regression datasets bank (Num:8192, Fea: 32), elevators (Num:16599, Fea: 18), ailerons (Num:13750, Fea: 40), Hardware (Num:28179, Fea: 96)
