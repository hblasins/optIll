# Designing illuminant spectral power distribution for surface classification

There are many application areas where users have full control over the illumination, and have the ability to shape its spectral power distribution. In those cases spectrally optimal light can accentuate features in images and hence improve image classification tasks.

<p align="center"> 
<img width="500px" src="https://github.com/hblasins/optIll/blob/master/Figures/Overview.png">
</p>

In this project we propose two approaches that estimate the optimal spectral power distribution of the illuminant. The unsupervised approach uses nonnegative sparse Principal Component Analysis to derive optimal, and physically realizable illuminant spectra. The supervised approach incorporates linear image formation model directly into classification algorithm and uses alternating minimization to simultaneously search for classifier decision boundaries and optimal lights.

If you use these tools please cite
```
@inproceedings{Blasinski_2017_CVPR,
    author = {Blasinski, Henryk and Farrell, Joyce and Wandell, Brian},
    title = {Designing Illuminant Spectral Power Distributions for Surface Classification},
    booktitle = {The IEEE Conference on Computer Vision and Pattern Recognition (CVPR)},
    month = {July},
    year = {2017}
}
```

The [manuscript](http://openaccess.thecvf.com/content_cvpr_2017/papers/Blasinski_Designing_Illuminant_Spectral_CVPR_2017_paper.pdf and the [supplement](http://openaccess.thecvf.com/content_cvpr_2017/supplemental/Blasinski_Designing_Illuminant_Spectral_2017_CVPR_supplemental.pdf) ara available on the [CVPR2017 website](http://openaccess.thecvf.com/CVPR2017.py).

## Dependencies

To succesffuly run the scripts in this project the following dependencies need to be 
installed.

* MATLAB - we tested the code with MATLAB 2015 and 2016. The illuminant selection algorithms do not require additional toolboxes, but some evaluation scripts do, specifically: Image Processing and Statistics and Machine Learning.
* [ISET](http://imageval.com) - Image Systems Engineering Toolbox for camera sensor simulations. A light version of ISET is provided with this code repository.
* [CVX](http://cvxr.com/) - a Matlab toolbox for convex optimization. As an alternative MATLAB Optimization Toolbox can also be used. **As of 06/2017 CVX is only supported on MATLAB 2016b and older.**

## Data

Sample data used in our experimens can be downloaded from the [Stanford Digital Repository](https://purl.stanford.edu/rq453qp3526). The `Images.zip` file contains all images captured using the experimental setup. The `Results.zip` contains pixel classification results, this data should be reproducible by running the appropriate `classify***.m` scripts. Note that given the extensive cross-validation these results will take long to generate, which is why we are also releasing them. Once you download the files, unzip them directly into the root folder of this repository.

## Getting started

Once you start MATLAB please run the `install.m` script from the `./Code` directory. This script adds all the relevant project sub-directories to MATLAB path.

To have a better sense of what this project is about and how optimal illuminants can be used run the `s_syntheticExample.m` script. It generates a toy problem where a target feature is invisible under broadband illuminant and shows how to use the algorithms we describe to derive the optimal illuminant. The feature becomes more visible when optimal light is used. 

To see the algorithms in action you can have a look at `s_simulationExample.m` or `s_simulationNumChannels.` scripts. These use the ISET toolbox to simulate camera captures and compare conventional RGB cameras and broadband illuminants to using monochrome cameras and optimal lights. 

The supervised and unsupervised illuminant selection algorithms are in `./Code/Algorithms` folder which also contains test scripts `t_***.m` that illustrate how these functions should be used and where possible (`sparsePCA.m`) compare the iterative ADMM implementation with that of `cvx`.


