from ROOT import *
import multiprocessing as mp
import os
import math
import optparse
import time
import sys
import numpy as np
import root_numpy as rn

## FUNCTIONS DEFINITION

def smearing(z_hit, y_hit, energyDep_hit, options):
	Z=list(); Y=list()
	Z*=0; Y*=0
	Z.append(np.random.normal(loc=z_hit, scale=options.diff_param, size=int(energyDep_hit*opt.Conversion_Factor)))
	Y.append(np.random.normal(loc=y_hit, scale=options.diff_param, size=int(energyDep_hit*opt.Conversion_Factor)))
	return Z, Y


def AddBckg(options):
	
	bckg_array=np.zeros((options.z_pix,options.y_pix))
	for s in range(0, options.z_pix):
		for t in range(0, options.y_pix):
			bckg_array[s][t]=np.random.normal(loc=options.noise_mean, scale=options.noise_sigma)
               
	return bckg_array 


def AddTrack(Tree, hist_list, smear_func):
	for i in range(0, Tree.numhits):
		S=smear_func(Tree.z_hits[i], Tree.y_hits[i], Tree.energyDep_hits[i], opt)
		for j in range(0, len(S[0])):
			for t in range(0, len(S[0][j])):
				hist_list[entry].Fill(S[0][j][t],S[1][j][t])
	signal_array=rn.hist2array(hist_list[entry])				
	
	return signal_array

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

	parser.add_option("-I", "--inputfolder", dest="infolder", default=os.getcwd()+"/src", help="specify the folder containing input files")
	parser.add_option("-O", "--outputfolder", dest="outfolder", default=os.getcwd()+"/out", help="specify the output destination folder")
	
        (opt, args) = parser.parse_args()

	config = open(args[0], "r") 		#GET CONFIG FILE
	params = eval(config.read()) 		#READ CONFIG FILE

	for k,v in params.items():
		setattr(opt, k, v)

#### CODE EXECUTION ####

        run_count=90000
	t0=time.time()
        if not os.path.exists(opt.outfolder):
            os.makedirs(opt.outfolder)
            
        for infile in os.listdir(opt.infolder):
		
		if infile.endswith('.root'):		

			rootfile=TFile.Open(opt.infolder+infile)
			tree=rootfile.Get('nTuple')
			outfile=TFile(opt.outfolder+'/histograms_Run'+str(run_count)+'.root', 'RECREATE')
               		SaveValues(params, outfile) ## SAVE PARAMETER OF THE RUN

			final_imgs=list(); nobckg_imgs=list()
			
			for entry in range(0, tree.GetEntries()):

				tree.GetEntry(entry)

				final_imgs.append(TH2F('pic_run'+str(run_count)+'_ev'+str(entry), '', opt.z_pix, -opt.z_dim*0.5, opt.z_dim*0.5, opt.y_pix, -opt.y_dim*0.5, opt.y_dim*0.5)) #smeared track with background

				signal=AddTrack(tree, final_imgs, smearing)

				if opt.bckg:
					background=AddBckg(opt)
				
				total=signal+background
				final_imgs[entry-1]=rn.array2hist(total, final_imgs[entry-1])

		## WRITE EACH TRACK ON THE ROOT FILE CORRESPONDING TO THE ENERGY

				print('%d images generated'%(entry+1))
				final_imgs[entry].Write()

                        run_count+=1
			outfile.Close()
	t1=time.time()
	print(t1-t0)

