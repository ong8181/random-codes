# Echo State Network (ESN)
Echo Stata Network (ESN) is an approach to Recurrent Neural Network (RNN) training. Here ESN is implemented with Python according to a tutrial described by Herbert Jaeger (http://minds.jacobs-university.de/uploads/papers/ESNTutorialRev.pdf).

Two argorithms are implemented: One for time series prediction and the other for MNIST image classification.

For time series prediction, time series data to be predicted is chaotic model time series analyzed in Deyle et al. (2016) _Proceeding of the Royal Society B_ (https://doi.org/10.1098/rspb.2015.2258) and Chang et al. (2017) _Ecological Research_ (https://doi.org/10.1007/s11284-017-1469-9). Data can be downloaded by clicking the following link: [Supplementary data in Chang et al. (2017)](https://esj-journals.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1007%2Fs11284-017-1469-9&file=ere0785-sup-0006.csv)). 

For MNIST image classification, the algorithm is similar to that presented in Schaetti et al. (2016) https://hal.archives-ouvertes.fr/hal-02131170. MNIST image data can be downloaded by excuting ```mnist = datasets.fetch_openml('mnist_784', version=1,)```.

