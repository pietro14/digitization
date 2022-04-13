MC_data_gen.py
===============
This simple script applies both the smearing due to diffusion and the electronic noise background to a MC sample track (GEANT4 output).
To run, the script need to be in the same location where ConfigFile_new.txt is.

The output file is a `.root` file containing all the TH2F histograms generated.

USAGE
-----
The usage of the script is the following:

`python MC_data_new.py ConfigFile_new.txt [options]`

With the following options available:

```Javascript
 Identifier               Option                         Default value
 
     -I           `<path_to_inputfolder>`         `<path_to_current_directory>+src/`
     -O           `<path_to_outputfolder>`        `<path_to_current_directory>+out/`
     
```
Given the input folder, the script will run over all the .root files it can find in that location.
The output file contains also a TDirectoryFile used as a storage for values imported from the `ConfigFile_new.txt`. They are stored in single-binned histogram, so you can easily access it using
```Javascript
param_dir->cd()
'histogram_name'->GetBinContent(1)
```
Also the type of particle, and its initial energy, for each event are stored in a subfolder. You can access it using

```Javascript
event_info->cd()
'histogram_name'->GetBinContent(1)
```

EXAMPLE
--------
Here an example is provided.

+ First of all, download the repository with `git clone git@github.com:CYGNUS-RD/digitization.git`
+ You want to specify the folder in which your GEANT4 simulation output is. If you don't have any MC output file, you can download one [here](https://drive.google.com/open?id=1hut-cRycXGwYfO5eJLUXaKKzAwQU_i0p)
+ Run the script with the following command line: `python MC_data_new.py ConfigFile_new.txt -I <path_to_input_folder>`

You will find the output in the default `out/` folder.

You can draw the image opening the output in an interactive ROOT session. To make the image similar to the experimental data, we advice to use the following commands

```Javascript
gStyle->SetPalette(kGreyScale)
gStyle->SetOptStat(0)
```
and to set properly the z-axis scale once the TH2F has been written with `COLZ` option.

Run in batch
-------------
To run in batch using PBS queue system you can use the script `submit_digi_batch.py`

Example command:

```
python scripts/submit_digi_batch.py `pwd` --inputdir /nfs/cygno/CYGNO-MC-data/pbs_outputs/CYGNO_60_40_ER_6_keV/ --outdir /nfs/cygno2/users/$USER/digitization-out/ --tag LIMEsaturation_test --conf ConfigFile_new_saturation.txt
```
If you want easily submit multiple jobs, you can use a similar script to `scripts/run_batch.sh` 

You can check the status of the jobs you submitted with:

```
qstat | grep $USER
```

And you can cancel a job with:

```
qdel <job_number>.xmaster
```

Suggestions for debugging and contributing
------------
If you have made minor changes to the code, and the physical model has not changed, the output should be the same (except for statistical fluctuations). 
Once you have set the same seed for random distributions, you can use the script compare_digitizations.py to easily compare the output of two simulations. For instance: 

```Javascript
python3 compare_digitizations.py output1.root output2.root
```

Note, the two root file should have the same number of events.

Wiki page and documentation
------------
You can find more information about the digitization and CYGNO Montecarlo on the [Wiki-page](https://github.com/CYGNUS-RD/WIKI-documentation/wiki/Digitization) 


Work in progress
------------
+ Add an option in `ConfigFile.txt` to choose between different detectors and geometries, in order to simulate other setups without manually changing the parameters
+ Parallelize background generation to make the script run faster
+ Speed up the simulation with saturation effect, maybe with a parametrization (done)
+ Find a way to apply saturation effect to non-spot tracks (maybe already done)
