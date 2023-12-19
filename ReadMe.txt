[![Video Preview](https://github.com/zkamranian/Video-Object-Segmentation/blob/master/video-object-segmentation.avi)](https://github.com/zkamranian/Video-Object-Segmentation/blob/master/video-object-segmentation.avi)

Joint Motion Boundary Detection and CNN-based Feature Visualization for Video Object Segmentation 

Co-segmentation via Visualization

Version 1
---------------------------------------------------------------------------------------------------------------------------------
The source code is prepared for both image co-segmentation and video object segmentation [I,II]. 
The code can be run on the related images or videos which their objects have been already trained in VGG-16. However, it can be applied to any image groups or videos, provided that a new CNN is trained on a dataset that includes its common object class.
Using the code:

1- Please download, install and compile Matconvnet from (http://www.vlfeat.org/matconvnet/). To speed up code execution, please compile the software in GPU support version.

2- Please download 'imagenet-vgg-verydeep-16.mat' [1] from (http://www.vlfeat.org/matconvnet/pretrained/), and load the model.

[1] Simonyan K, Zisserman A (2014) Very deep convolutional networks for  large-scale image recognition. arXiv preprint arXiv:14091556.
 
3- For oversegmentation, gPb-UCM code is used [2,3]. Please download and install mcg-2.0 code from (https://www2.eecs.berkeley.edu/Research/Projects/CS/vision/grouping/resources.html)

[2]P. Arbelaez, J. Pont-Tuset, J. T. Barron, F. Marques, J. Malik, Multiscale combinatorial grouping, in: Computer Vision and Pattern Recognition (CVPR), Conference on, IEEE, 2014, pp. 328-335.
[3] P. Arbelaez, M. Maire, C. Fowlkes, J. Malik, Contour detection and hierarchical image segmentation, Pattern Analysis and Machine Intelligence, IEEE transactions on 33 (2011) 898-916.

4- Please download and compile Matlab wrappers from (http://calvin.inf.ed.ac.uk/software/fast-video-segmentation/)
to find motion boundary detection [4].

[4] Papazoglou A, Ferrari V (2013) Fast object segmentation in unconstrained video. In: Computer Vision (ICCV) International Conference, IEEE, pp 1777–1784. 

5- run main_demo.m


Please note that the code is the first version. 
In this version, 
- some functions are simpler than original ones.
- Inspired by [5], we use spline regression to learn the local energy. 
Thus, some results may are different from what is reported on the paper [II]. 

 
[5] X. Dong, J. Shen, L. Shao, M.-H. Yang, Interactive cosegmentation using global and local energy optimization, Image Processing, IEEE Transactions on 24 (2015) 3966{3977.

For further information, please do not hesitate to contact zahra_kamranian@yahoo.com

If you use this software for academic research, please cite the following publication:

[I] Kamranian Z, Tombari F, Nilchi ARN, Monadjemi A, Navab N (2018) Co-segmentation via visualization. Journal of Visual Communication and Image Representation 55:201–214.
[II] Kamranian Z, Nilchi ARN, Sadeghian H, Tombari F, Navab N (2019)Joint Motion Boundary Detection and CNN-based Feature Visualization for Video Object Segmentation. Neural Computing with Applications.
