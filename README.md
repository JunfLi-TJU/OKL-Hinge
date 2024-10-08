# OKL-Hinge
Source codes of algorithms and datasets for our paper "Improved Kernel Alignment Regret Bound for Online Kernel Learning", accepted in AAAI 2023.

We implement all algorithms with R on a Windows machine with 2.8 GHz Core(TM) i7-1165G7 CPU. execute each experiment 10 times with random permutation of all datasets and average all of the results.

To run the code, you must set the paths following the code, or set new paths.
The default path of codes is "D:/experiment/AAAI2023/code". 
The path of datasets is "D:/experiment/AAAI2023/dataset". 
The store path is "D:/experiment/AAAI2023/Result/". 

The baseline algorithms include: FOGD, NOGD, SkeGD and BAO2KS. Our algorithm is POMDR.

We are grateful to Xiao Zhang [1] for providing the code of SkeGD.
To run SkeGD (or OKS-SkeGD), the datasets should be located in the "c_data" folder.
The main file is named "experiment_OKS.R".

The datasets are downloaded from: https://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/ and http://archive.ics.uci.edu/ml/datasets.php

binary classification datasets: 
w8a (Num:49749, Fea:300), magic04 (Num:19020, Fea:10), mushrooms (Num:8124, Fea:112), a9a (Num:48842, Fea:123), 
SUSY (Num:50000, Fea:18), ijcnn1 (Num:141691, Fea:22), minist12 (Num:12700, Fea:780)

Reference

[1] Xiao Zhang and Shizhong Liao. Incremental Randomized Sketching for Online Kernel Learning. ICML,2019.
