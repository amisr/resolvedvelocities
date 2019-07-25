Instructions for the synthetic data module

The synthetic data module is designed to help evaluate how well the resolved velocities algorithm is reconstructing an existing, known field.  An arbitrary field in the ionosphere can be constructed, an then a synthetic data set is created by calculating what line-of-site velocity an AMISR with a particular beam pattern would measure.  The synthetic data set is then run through the standard resolved velocities algorithm and the output is compared with the known input ionosphere.  This helps evaluate:

- How well the resolved velocities algorithm can reconstruct features such as flow shears and vorticies
- If certain beam patterns or configurations are better for reconstructing the true velocity field.

To run:
1. Specify the "true" field by modifying the [FIELD] section of synth_config_example.ini which contains both an array of geodetic coordinates and the geodetic East, North, and Up components of the plasma drift velocity at those coordinates.  The script create_synth_field.py can be modified and run to help generate these arrays for certain simple field configurations.

2. Specify the radar location and beam configuration in the [RADAR] section of synth_config_example.ini.

3. Modify vvels_config_example.ini to specify what options you would like the resolved velocities algorithm to use.  This is just a standard vvels config file and all the normal options apply.

4. Run run_sythetic_data.py with the two config files as command line input to generate the synthetic data set, resolve vectors from this synthetic data set, and compare the resolved vectors with the "true" field values.

	$ python run_synthetic_data.py synth_config_example.ini vvels_config_example.ini
