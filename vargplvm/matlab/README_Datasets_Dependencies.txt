# 1) VARGPLVM0p13 is the software toolbox developed for the method described in the paper.
     This software depends on various other toolboxes, which can be found in:
      http://staffwww.dcs.shef.ac.uk/people/N.Lawrence/software.html
      Specifically, the additional toolboxes needed are the following:
	DATASETS0p1371
	GP0p136
	KERN0p225
	MLTOOLS
	MOCAP0p136
	NDLUTIL0p162
	netlab
	NOISE0p141
	OPTIMI0p132
	PRIOR0p22

# 2) All the demos used to produce the paper's results can be found in the file: 
	'demosDynamics.m'. Alternatively, the submitted .tex file also contains the same code
	and can be run as a regular .m file using MATweave:
	http://staffwww.dcs.shef.ac.uk/people/N.Lawrence/matweave.html

	
# 3) The datasets are not included because they are very large in size. 
	Below you can find instructions on how you can obtain them and process
	them in the same way we did. Note, however, that the videos have to be
	processed using WINDOWS (unless you do not wish to process them with the
	provided script), because aviread does not support video compression under
	LINUX.

	
####### Mocap dataset #####
Can be downloaded from:
http://mocap.cs.cmu.edu/
The subject / motion used and the preproprecssing is described in the paper.



########    Obtaining and processing the VIDEO datasets #############
a) 'dog' dataset:
http://www.fitfurlife.com/
OR
http://www.youtube.com/watch?v=K4Dm2HcfUKQ

We only kept the 61 frames showing the first dog on the treadmill
(around 0.02 - 0.04). This cropping can be done in MATLAB or using
a free program such as http://www.virtualdub.org/download.html, as
we did.

The .avi can be converted to a .mat file using our script:
[Y, height, width] = preprocessVideo(<your_avi_file>, 1, 1, 'avi', []);
save 'DATA_DogTr.mat' 'Y' 'height' 'width'
and place the .mat file into the DEPENDENCIES folder.

b) For the 'ocean' demo you can work as in a), after downloading the
video from:
http://www.youtube.com/watch?v=12QU41ZfqHc
We downloaded the HD version of that video. If you wish to do the same, 
you can work as in a) but maybe you can first split the .avi file into
pieces as we did, because our machine could not load the whole dataset at once.

c) similarly the 'missa' dataset can be download from:
http://www.cipr.rpi.edu/resource/sequences/sequence01.html
Ideally, you can transform that into avi and use the script above to
obtain a .mat file. If you are unable to do that, you can still
use the script "preprocessVideo" but you have to give the correct dimensions
for the frames and set the 'form' argument accordingly.

For compatibility with the rest of the scripts, save the .mat file with the name
'miss-americaHD.mat' and place it into the DEPENDENCIES folder.


