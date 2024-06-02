MATLAB 2021a Scripts for the Analysis of Microtubule Dynamics Visualized with Dark-Field Microscopy

These scripts allow for:

•	Tracking the position of the MT tip over time.

•	Quantifying the protofilament configuration at the microtubule tip over time.

For optimal performance of the algorithm, the dark-field intensity needs to be calibrated. The calibration procedure is described at the end of these instructions.


HOW TO USE THE SCRIPTS

I. Add the Scripts to MATLAB Path

Add the scripts in the subfolders 'B4D_v3p2' and 'steps_bumps_toolkit' to the MATLAB Path.

II. PreProcess.m Script

The PreProcess.m script is used to crop the TIFF file with microscopy data, adjust the brightness and contrast, and denoise the data. Use the following sequence of steps:
1.	Edit lines 5 and 6 of PreProcess.m to specify the folder and name of the TIFF file with dark-field microscopy data.
2.	Run the PreProcess.m file.
3.	The middle frame of the TIFF file will be displayed.

•	To change contrast: press the left mouse button and move the cursor horizontally while keeping the button pressed.

•	To change brightness: press the left mouse button and move the cursor vertically while keeping the button pressed.

4.	Move and resize the Region Of Interest (ROI) to contain the microtubule of interest throughout the whole movie. Ensure that the tip of the microtubule is separated from any edge of the ROI by at least 5 pixels. Note: the larger the size of the ROI, the longer it takes the script to process the data, but the higher the quality of the tracking.
5.	Press Enter.
6.	In the command line, specify the number of the first and last frames to analyze. Wait for the message 'video has been prepared for MT tracking' to appear.

III. DF_MT_tracker.m Script

The DF_MT_tracker.m script is used to track the position of the microtubule tip and quantify the protofilament configuration at the tip. Use the following sequence of steps:

FOR CALIBRATION:
1.	Edit lines 2-15 of the DF_MT_tracker.m file to set tracking parameters. Typically, only the following parameters need to be modified:

•	calibration_on = 1 (should be set to one)

•	show = 1 (display tracking results at each frame)

•	pix (the size of the pixel in meters)

2.	Run the DF_MT_tracker.m file. The first frame of the pre-selected video fragment is displayed after (left half) and before (right half) denoising.
3.	On the denoised image, left-click on the point where the microtubule seed starts and then right-click on the point where the microtubule seed ends (this point defines the dynamic end to be tracked). This launches the calibration procedure. The program reports the calibrated 'tube' parameter and the corresponding dependence of the PF number vs. time.
4.	Record the 'tube' parameter if satisfied.

FOR TRACKING:
1.	Edit lines 2-15 of the DF_MT_tracker.m file to set tracking parameters. Typically, only the following parameters need to be modified:

•	calibration_on = 0 (should be set to zero)

•	show = 1 (display tracking results at each frame)

•	pix (the size of the pixel in meters)

•	tube (use the value determined during the calibration above)

2.	Run the DF_MT_tracker.m file. The first frame of the pre-selected video fragment is displayed after (left half) and before (right half) denoising.
3.	On the denoised image, left-click on the point where the microtubule starts and then right-click on the point where the microtubule ends (this point defines the dynamic end to be tracked). This launches the tracking procedure, which can take several minutes.
4.	Two plots are displayed:

•	The position of the microtubule end vs. frame number (blue) overlaid on the tip extension vs. frame number (orange).

•	Mean protofilament number in each of the zones at the microtubule tip vs. frame number. The size of the zones is defined by the parameter rtip (line 12), and the number of zones is defined by the parameter level_dist (line 11).

OUTPUT FILES:

• ".... MT_length.csv"-	The position of the microtubule end vs. frame number

•	".... Tip_extension.csv"- The tip extension vs. frame number

•	".... PF_number.csv" - Mean protofilament number in each of the zones at the microtubule tip vs. frame number

REFERENCES:

[1] M. Maggioni, V. Katkovnik, K. Egiazarian, A. Foi, "A Nonlocal 
    Transform-Domain Filter for Volumetric Data Denoising and 
    Reconstruction", IEEE Trans. Image Process., vol. 22, no. 1,
	pp. 119-133, Jan. 2013. doi:10.1109/TIP.2012.2210725

[2] M. Maggioni, A. Foi, "Nonlocal Transform-Domain Denoising of 
    Volumetric Data With Groupwise Adaptive Variance Estimation", 
    Proc. SPIE Electronic Imaging 2012, San Francisco, CA, USA, Jan. 2012

[3] R. Vincent, "Brainweb:  Simulated  brain  database", online at
    http://mouldy.bic.mni.mcgill.ca/brainweb/, 2006.

[4] Z. Wang, A. Bovik, H. Sheikh, E. Simoncelli, "Image quality 
    assessment: from error visibility to structural similarity",
    IEEE Trans. Image Process., vol. 13, no. 4, pp. 600-612, April 2004.

[5] J. V. Manjon, P. Coupe, A. Buades, D. L. Collins, M. Robles, 
    "New methods for MRI denoising based on sparseness and self-similarity", 
    Medical Image Analysis, vol. 16, no. 1, pp. 18-27, January 2012

[6] M. A. Little and N. S. Jones, "Sparse Bayesian step-filtering for high-throughput 
analysis of molecular machine dynamics," 2010 IEEE International Conference on Acoustics, 
Speech and Signal Processing, Dallas, TX, USA, 2010, pp. 4162-4165, 
doi: 10.1109/ICASSP.2010.5495722.

[7] M.A. Little, B.C. Steel, F. Bai, Y. Sowa, T. Bilyard, D.M. Mueller, R.M. Berry,
N.S. Jones (2010) "Steps and Bumps: Precision Extraction of Discrete States of
Molecular Machines using Physically-based, High-throughput Time Series Analysis",
arXiv:1004.1234v1 [q-bio.QM]

[8] B. Kalafut, K. Visscher (2008) "An objective, model-independent method for
detection of non-uniform steps in noisy signals", Comp. Phys. Comm.,
179(2008):716-723.

[9] S.H. Chung, R.A. Kennedy (1991), "Forward-backward non-linear filtering
technique for extracting small biological signals from noise", J. Neurosci. Methods. 
40(1):71-86.



