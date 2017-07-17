# Designing illuminant spectral power distribution for surface classification

There are many application areas where users have full control over the illumination, and have the ability to shape its spectral power distribution. In those cases spectrally optimal light can accentuate features in images and hence improve image classification tasks.

In this project we propose two approaches that estimate the optimal spectral power distribution of the illuminant. The unsupervised approach uses nonnegative sparse Principal Component Analysis to derive optimal, and physically realizable illuminant spectra. The supervised approach incorporates linear image formation model directly into classification algorithm and uses alternating minimization to simultaneously search for classifier decision boundaries and optimal lights.

If you use these tools please cite
```
@inproceedings{blasinski2017optIll,
    title={Designing illuminant spectral power distribution for surface classification},
    author={Blasinski, Henryk and Farrell, Joyce and Wandell, Brian},
    booktitle={IEEE Computer Vision and Pattern Recognition, CVPR},
    year={2017},
    organization={IEEE}
}
```

## Dependencies

To succesffuly run the scripts in this project the following dependencies need to be 
installed.

* MATLAB - we tested the code with MATLAB 2015 and 2016. The illuminant selection algorithms do not require additional toolboxes, but some evaluation scripts do, specifically: Image Processing and Statistics and Machine Learning.
* [ISET](http://imageval.com) - Image Systems Engineering Toolbox for camera sensor simulations. A light version of ISET is provided with this code repository.
* [CVX](http://cvxr.com/) - a Matlab toolbox for convex optimization. As an alternative MATLAB Optimization Toolbox can also be used. **As of 06/2017 CVX is only supported on MATLAB 2016b and older.**

## Getting started



