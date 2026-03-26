# neutralwinds
Estimate the profile of neutral winds in the lower ionosphere (E region) using the method outlined in [Heinselman and Nicholls, 2008](https://amisr.com/publications/media/pub-pdfs/2008_Heinselman_10.10292007RS003805.pdf).  This is related to the standard AMISR [resolvedvelocities](https://github.com/amisr/resolvedvelocities) technique, but does not depend on that code.

## Usage Instructions
This code is presently written in python 2, which was depricated in 2020.  It is HIGHLY recomended that you set up an independent virtual environment or use your favorite flavor of containerization as it requires installing old versions of packages that will not be compatable with newer python code.

1. Clone the repository locally.
```
git clone https://github.com/amisr/neutralwinds.git
cd neutralwinds
```

2. Version pinned requisite packages are available in `requirements.txt`.  These can be installed with pip.
```
pip install -r requirements.txt
```

3. Enter the `ProcessingCode` directory.
```
cd ProcessingCode
```

4. Modify the config file `default.ini` so that `RootFileDirectory` points to where you would like to save the output on the local machine.  You can change any other parameters in the config file as well based on how you would like the code to run.  Feel free to copy this config file and save it elsewhere.

5. Run the script `ProcessEregionNeutralWinds.py` with command line arguments of the input alternating code file, the input longpulse file, and the config file.
```
python ProcessEregionNeutralWinds.py /path/to/ac_file.h5 /path/to/lp_file.h5 /path/to/config_file.ini
```

## Acknowledgments
This software was originally provided by S. Kaeppler.
