from ROOT import *
import multiprocessing as mp
import os
import math
import optparse
import time
import sys
import numpy as np

## FUNCTIONS DEFINITION

def smearing(z_hit, y_hit, energyDep_hit, options):
	Z=list(); Y=list()
	Z*=0; Y*=0
	Z.append(np.random.normal(loc=z_hit, scale=options.diff_param, size=int(energyDep_hit*opt.Conversion_Factor)))
	Y.append(np.random.normal(loc=y_hit, scale=options.diff_param, size=int(energyDep_hit*opt.Conversion_Factor)))
	return Z, Y


def AddBckg(histo, options):
	for s in range(0, options.z_pix):
		for t in range(0, options.y_pix):
			r=np.random.normal(loc=options.noise_mean, scale=options.noise_sigma)
                	histo.SetBinContent(s,t,(histo.GetBinContent(s,t)+r))
	return None 


def AddTrack(Tree, hist_list, smear_func):
	for i in range(0, Tree.numhits):
		S=smear_func(Tree.z_hits[i], Tree.y_hits[i], Tree.energyDep_hits[i], opt)
		for j in range(0, len(S[0])):
			for t in range(0, len(S[0][j])):
				hist_list[entry].Fill(S[0][j][t],S[1][j][t])
				
	return None

def SaveValues(par, out):

	out.cd()
	fold=TDirectoryFile('fold', '')
	fold.cd()

	for k,v in par.items():
		h=TH1F(k, '', 1, 0, 1)
		h.SetBinContent(1, v)
		h.Write()

	#fold=None

	return None


## MAIN EXECUTION

if __name__ == "__main__":

	#CALL: python img_generator.py ConfigFile.txt -I <<INPUT_FOLDER>> -F MC_img_runs.txt (optional: -S 9xxxx)

	parser = optparse.OptionParser("usage: %prog [options] arg1 arg2")

	parser.add_option("-I", "--inputfolder", dest="infolder", default=os.getcwd()+"src", help="specify the folder containing input files")
	parser.add_option("-O", "--outputfolder", dest="outfolder", defaul=os.getcwd()+"out", help="specify the output destination folder")
	parser.add_option("-F", "--run_number_file", dest="runs", default="MC_runs.txt", help="specify the .txt file in which the last number of run used is")
	parser.add_option("-S", "--set_run_number", dest="def_Nrun", default=None, help="manually sets the run number")	

	(opt, args) = parser.parse_args()

	config = open(args[0], "r") 		#GET CONFIG FILE
	params = eval(config.read()) 		#READ CONFIG FILE

	for k,v in params.items():
		setattr(opt, k, v)

#### CODE EXECUTION ####

	for infile in os.listdir(opt.infolder):
		
		if infile.endswith('.root'):		
			if opt.def_Nrun == None:

	                	run_number = open(opt.runs, "r")        ## GET INCREMENTAL NUMBER OF FILE
        	        	Nrun=run_number.read()                  ## READ LAST RUN NUMBER

                		run_number = open(opt.runs, "w")
                		run_number.write(str(int(Nrun)+1))      ## WRITE ON FILE NEW RUN NUMBER
                		run_number.close()                      ## CLOSE THE FILE

        		else:

                		Nrun=opt.def_Nrun


			rootfile=TFile.Open(opt.infolder+infile)
			tree=rootfile.Get('nTuple')
			outfile=TFile('out/histograms_Run'+Nrun+'.root', 'RECREATE')
               		SaveValues(params, outfile) ## SAVE PARAMETER OF THE RUN

			final_imgs=list(); nobckg_imgs=list()

			for num in range(0, tree.GetEntries()-9998):
			
				final_imgs.append(TH2F()) ## HISTOGRAMS LIST SETUP
			
			for entry in range(0, tree.GetEntries()-9998):
		
				print('number of tracks processed: %d'%(entry+1))
				tree.GetEntry(entry)

				final_imgs[entry]=TH2F('pic_run'+Nrun+'_ev'+str(entry), '', opt.z_pix, -opt.z_dim*0.5, opt.z_dim*0.5, opt.y_pix, -opt.y_dim*0.5, opt.y_dim*0.5) #smeared track with background

				AddTrack(tree, final_imgs, smearing)

				if opt.bckg:
					AddBckg(final_imgs[entry], opt)
				
		## WRITE EACH TRACK ON THE ROOT FILE CORRESPONDING TO THE ENERGY

				final_imgs[entry].Write()	

			outfile.Close()


