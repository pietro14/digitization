#!/bin/bash

# This script is for submitting to the queue several jobs at once, for instance when doing scans over one o more parameters (HV, z_gem, etc...)
# Run the script from the main directory as: `source scripts/run_batch.sh`
# Currently it works for ER root files with name CYGNO_60_40_ER_<energy>_keV, you can adapt it for NR and other input files

# The script creates a new config file for each run and saves it in ./config. So, creates config dir, if it doesn't exist  
mkdir -p ./config

# set input folders 
rootfile_inputdir="/nfs/cygno/CYGNO-MC-data/pbs_outputs/"    # for 6 keV ER (Fe55)
#rootfile_inputdir="/nfs/cygno2/CYGNO-MC-data/pbs_outputs/"  # for other energies (Cu, Rb, Mo, Ag, Ba, Tb)

z_gem=250 
saturation=True
bckg=True
noiserun=4432
events=100
energy=6
tag_config="\'Data\'"   

GEM1_HV=440
GEM2_HV=440
GEM3_HV=440

diff_const_sigma0T="0.1225"
diff_coeff_T="0.01232"
diff_const_sigma0L="0.0676"
diff_coeff_L="0.00978"

camera_aperture=0.95
ion_pot=0.0462

A=1
beta="1.0e-5" 
absorption_l=1000

# example: scan over GEM1_HV
for GEM1_HV in 440 350 386 406 420 431 440 
do
	#example: scan over z_gem
	for z_gem in 250 
	do

		#set current input dir
		inputdir=$(echo "$rootfile_inputdir"/CYGNO_60_40_ER_"$energy"_keV/|tr "." "p")
		#set current output idr 
		outdir="/nfs/cygno/users/pmeloni/CYGNO/feb22_lug22/digitization-out"
		tag=$(echo LIMEsaturation_abs"$absorption_l"_beta"$beta"_sT0350_"$z_gem"mm_HV"$GEM1_HV"_pedrun"$noiserun"_"$energy"keV_"$A"A|tr "." "p" | sed 's/e-5//g')

		conf="config/config_"$tag".txt"

		# creating new config file in ./config dir
		cp ConfigFile_new.txt "$conf" 

		# setting parameters in the current config file
		sed -i 's/'\''z_gem'\'' *\t*: .*/'\'z_gem\'\ :\ 255.+$z_gem\.\,'/'  "$conf"
		sed -i 's/'\''beta'\'' *\t*: .*/'\'beta\'\ :\ $beta\,'/'  "$conf"
		sed -i 's/'\''ion_pot'\'' *\t*: .*/'\'ion_pot\'\ :\ $ion_pot\,'/'  "$conf"
		sed -i 's/'\''camera_aperture'\'' *\t*: .*/'\'camera_aperture\'\ :\ $camera_aperture\,'/'  "$conf"
		sed -i 's/'\''absorption_l'\'' *\t*: .*/'\'absorption_l\'\ :\ $absorption_l\.\,'/'  "$conf"
		sed -i 's/'\''GEM1_HV'\'' *\t*: .*/'\'GEM1_HV\'\ :\ $GEM1_HV\.\,'/' "$conf"
		sed -i 's/'\''GEM2_HV'\'' *\t*: .*/'\'GEM2_HV\'\ :\ $GEM2_HV\.\,'/' "$conf"
		sed -i 's/'\''GEM3_HV'\'' *\t*: .*/'\'GEM3_HV\'\ :\ $GEM3_HV\.\,'/' "$conf"
		sed -i 's/'\''saturation'\''\t* *: .*,/'\'saturation\'\ :\ $saturation\,'/' "$conf"
		sed -i 's/'\''bckg'\''\t* *: .*,/'\'bckg\'\ :\ $bckg\,'/' "$conf"
		sed -i 's/'\''noiserun'\''\t* *: .*,/'\'noiserun\'\ :\ $noiserun\,'/' "$conf"
		sed -i 's/'\''events'\''.*/'\'events\'\ :\ $events\,'/' "$conf"
		sed -i 's/'\''tag'\''.*/'\'tag\'\ :\ $tag_config\,'/' "$conf"
		sed -i 's/'\''diff_const_sigma0T'\''.*/'\'diff_const_sigma0T\'\ :\ $diff_const_sigma0T\,'/' "$conf"
		sed -i 's/'\''diff_coeff_T'\''.*/'\'diff_coeff_T\'\ :\ $diff_coeff_T\,'/' "$conf"
		sed -i 's/'\''diff_const_sigma0L'\''.*/'\'diff_const_sigma0L\'\ :\ $diff_const_sigma0L\,'/' "$conf"
		sed -i 's/'\''diff_coeff_L'\''.*/'\'diff_coeff_L\'\ :\ $diff_coeff_L\,'/' "$conf"

		# if you want to submit the jobs 
		python3 scripts/submit_digi_batch.py `pwd` --inputdir $inputdir --outdir $outdir --tag $tag --conf $conf

		# if you want to run without submitting to the queue
		#python MC_data_new.py $conf -I $inputdir -O "$outdir/out/$tag"

		# if you want just to see the config files
		#cat $conf

	done
done


