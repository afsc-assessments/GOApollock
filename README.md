# GOApollock
Files for the Gulf of Alaska pollock assessment. Initialized from the accepted model in 2020 (pk20_8.tpl).

This repository contains three main components. First is the R package containing functions to interact with the assessment model (reading, writing, processing outputs, making figures, etc.). Second is the folder "markdown" which contains the RMD files needed to build the SAFE. The SAFE is not reproducible from this repo because it is missing some data inputs and outputs which are not currently shared publicly (but available upon request). Finally, the model source code for the ADMB and TMB models are in the 'source' folder. 

Releases are done at the end of each year when the assessment and SAFE is finalized. Thus one can reproduce an older assessment by installing a previous version of this package and then running the scripts and such in the appropriate folder. 
The 'TMB' folder contains some bridging materials and will eventually be deleted.
