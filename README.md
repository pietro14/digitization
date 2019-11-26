MC_data_gen.py
===============
This simple script applies both the smearing due to diffusion and the electronic noise background to a MC sample track (GEANT4 output).
To run, the script need to be in the same location where ConfigFile.txt is. This file contains all the parameters that can be manually set as the user prefers.
The output file is a `.root` file containing all the TH2F histograms generated.

USAGE
-----
The usage of the script is the following:

`python img_generator.py ConfigFile.txt [options]`

With the following options available:

```Javascript
 Identifier               Option                         Default value
 
     -I           `<path_to_inputfolder>`         `<path_to_current_directory>+src/`
     -O           `<path_to_outputfolder>`        `<path_to_current_directory>+out/`
     
```
Given the input folder, the script will run over all the .root files it can find in that location.
The output file contains also a TDirectoryFile used as a storage for values imported from the `ConfigFile.txt`. They are stored in single-binned histogram, so you can easily access it using
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
+ Run the script with the following command line: `python MC_data_gen.py ConfigFile.txt -I <path_to_input_folder>`

You will find the output in the default `out/` folder.

You can draw the image opening the output in an interactive ROOT session. To make the image similar to the experimental data, we advice to use the following commands

```Javascript
gStyle->SetPalette(kGreyScale)
gStyle->SetOptStat(0)
```
and to set properly the z-axis scale once the TH2F has been written with `COLZ` option.

Work in progress
------------
+ Add an option in `ConfigFile.txt` to choose between different detectors and geometries, in order to simulate other setups without manually changing the parameters
+ Parallelize background generation to make the script run faster

