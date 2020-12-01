from ROOT import *
import multiprocessing as mp
import os
import math
import optparse
import time
import sys
import numpy as np
import root_numpy as rn
from scipy.stats import expon
from scipy.stats import poisson

#sys.path.append("../reconstruction")
import swiftlib as sw


## FUNCTIONS DEFINITION


def smearing(z_hit, y_hit,x_hit, energyDep_hit, options):
    Z=list(); Y=list()
    Z*=0; Y*=0
    ## in Geant4 x is the drift axis
    #Z.append(np.random.normal(loc=z_hit, scale=np.sqrt(options.diff_const_sigma0+options.diff_coeff_B*(np.abs(x_hit-options.x_gem))), size=int(energyDep_hit*opt.Conversion_Factor)))
    #Y.append(np.random.normal(loc=y_hit, scale=np.sqrt(options.diff_const_sigma0+options.diff_coeff_B*(np.abs(x_hit-options.x_gem))), size=int(energyDep_hit*opt.Conversion_Factor)))
    nph=Nphotons(energyDep_hit, options)
    nph2=int(energyDep_hit*opt.Conversion_Factor) #not used, keeping here for back compatibility
    Z.append(np.random.normal(loc=(z_hit+0.5*options.z_dim)*options.z_pix/options.z_dim, scale=np.sqrt(options.diff_const_sigma0+options.diff_coeff_B*(np.abs(x_hit-options.x_gem))/10.)*options.z_pix/options.z_dim, size=int(nph))) #use nph2 here not to use gain fluctuations
    Y.append(np.random.normal(loc=(y_hit+0.5*options.y_dim)*options.y_pix/options.y_dim, scale=np.sqrt(options.diff_const_sigma0+options.diff_coeff_B*(np.abs(x_hit-options.x_gem))/10.)*options.y_pix/options.y_dim, size=int(nph))) #use nph2 here not to use gain fluctuations

    #print("distance from gem = "+str(np.abs(x_hit-options.x_gem))+" mm")
    return Z, Y

def Nphotons(energyDep_hit, options):

        # compute the number of ionization (primary) electrons) with a poisson distribution
        n_ioniz_el_mean=round(energyDep_hit/options.ion_pot)   # mean number of ionization electrons
        primary=poisson(n_ioniz_el_mean)           # poisson distribution for primary electrons
        n_ioniz_el=primary.rvs()                   # number of primary electrons
        # compute the number of secondary electrons (multiplication) considering fluctuations in the first GEM foil only
        n_el_oneGEM=0                              # number of secondary electrons
        gain2=expon(loc=0,scale=options.GEM_gain)  # exponential distribution for the GAIN in the first GEM foil
        #print ('--ioniz el= %d'%(n_ioniz_el))
        for k in range(0,n_ioniz_el):
            nsec = gain2.rvs()                     # number of secondary electrons in the first GEM multiplication for each ionization electron
            #nsec = options.GEM_gain               # use this line and comment out the previous one to eliminate gain fluctuations
            n_el_oneGEM += nsec
            #print ('--   loop on ioniz el, k= %d - nsec= %d - nel_onegem= %d'%(k,nsec,n_el_oneGEM))

        n_tot_el=n_el_oneGEM*pow(options.GEM_gain,2)  # total number of secondary electrons considering the gain in the 2nd and 3rd GEM foils
        nmean_tot_ph=n_tot_el*options.photons_per_el       # mean total number of photons
        photons=poisson(nmean_tot_ph)                    # poisson distribution for photons
        n_tot_ph=photons.rvs()                   # number of total photons

        # compute the number of photons hitting the sensor
        demag=options.y_dim/options.sensor_size  # optical de-magnification (in the config there are y_dim and z_dim. Which to use? In principle they should be the same --> a single number in the config)
        a=options.camera_aperture               # camera aperture
        omega=1./math.pow((4*(demag+1)*a),2)   # solid angle ratio
        n_photons=n_tot_ph*omega               # number of photons on the sensor
        #print ('  --> Energy= %f - Ion el= %d - tot photons= %d - number of photons= %d'%(energyDep_hit,n_ioniz_el,n_tot_ph,n_photons))
        return n_photons

def AddBckg(options, i):
    bckg_array=np.zeros((options.z_pix,options.y_pix))
    if options.bckg:
        #for s in range(0, options.z_pix):
        #    for t in range(0, options.y_pix):
        #        bckg_array[s][t]=np.random.normal(loc=options.noise_mean, scale=options.noise_sigma)
        if sw.checkfiletmp(int(options.noiserun)):
            options.tmpname = "/tmp/histograms_Run%05d.root" % int(options.noiserun)
        else:
            print ('Downloading file: ' + sw.swift_root_file(options.tag, int(options.noiserun)))
            options.tmpname = sw.swift_download_root_file(sw.swift_root_file(options.tag, int(options.noiserun)),int(options.noiserun))
        tmpfile = TFile.Open(options.tmpname)
        tmphist = tmpfile.Get("pic_run%05d_ev%d"% (int(options.noiserun),i))
        bckg_array = rn.hist2array(tmphist) 
    #print(bckg_array)
    
    return bckg_array


def AddTrack(Tree, histo, smear_func):
    for i in range(0, Tree.numhits):
        S=smear_func(Tree.z_hits[i], Tree.y_hits[i], Tree.x_hits[i], Tree.energyDep_hits[i], opt)
        #print(S)
        for j in range(0, len(S[0])):
            for t in range(0, len(S[0][j])):
                                #hist_list[entry].Fill(S[0][j][t],S[1][j][t])
                histo.Fill(S[0][j][t],S[1][j][t])
                #print(str(S[0][j][t])+"  "+str(S[1][j][t]))

#    signal_array=rn.hist2array(hist_list[entry])
    signal_array=rn.hist2array(histo)
    return signal_array

def AddTrackGen(numhits, xx_hits, yy_hits, zz_hits, edep_hits, histo, smear_func):
    for i in range(0, numhits):
        S=smear_func(xx_hits[i], yy_hits[i], zz_hits[i], edep_hits[i], opt)
        #print(S)
        for j in range(0, len(S[0])):
            for t in range(0, len(S[0][j])):
                histo.Fill(S[0][j][t],S[1][j][t])

    #signal_array=rn.hist2array(hist_list[lastentry])                
    signal=rn.hist2array(histo)
    return signal

def SaveValues(par, out):

    out.cd()
    out.mkdir('param_dir')
    #gDirectory.pwd()
    gDirectory.ls()
    out.cd('param_dir')
    
    for k,v in par.items():
        if (k!='tag'):
            h=TH1F(k, '', 1, 0, 1)
            h.SetBinContent(1, v)
            h.Write()
    out.cd()

    return None

def SaveEventInfo(info_dict, folder, out):
    
    out.cd()
    #gDirectory.pwd()
    out.cd('event_info')

    for k,v in info_dict.items():
        h=TH1F(k, '', 1,0,1)
        h.SetBinContent(1,v)
        h.Write()
    out.cd()
    info_dict.clear()

    return None

######################################### MAIN EXECUTION ###########################################

if __name__ == "__main__":

    #CALL: python img_generator.py ConfigFile.txt -I <<INPUT_FOLDER>>

    parser = optparse.OptionParser("usage: %prog [options] arg1 arg2")

    parser.add_option("-I", "--inputfolder", dest="infolder", default=os.getcwd()+"/src", help="specify the folder containing input files")
    parser.add_option("-O", "--outputfolder", dest="outfolder", default=os.getcwd()+"/out", help="specify the output destination folder")
    
    (opt, args) = parser.parse_args()

    config = open(args[0], "r")         #GET CONFIG FILE
    params = eval(config.read())         #READ CONFIG FILE

    for k,v in params.items():
        setattr(opt, k, v)

#### CODE EXECUTION ####
    run_count=1
    t0=time.time()
        
    if not os.path.exists(opt.outfolder): #CREATING OUTPUT FOLDER
        os.makedirs(opt.outfolder)
            
    for infile in os.listdir(opt.infolder): #READING INPUT FOLDER
            
        if opt.rootfiles==True:
                # code to be used with input root files from Geant
            if infile.endswith('.root'):    #KEEPING .ROOT FILES ONLY
                rootfile=TFile.Open(opt.infolder+infile)
                tree=rootfile.Get('nTuple')            #GETTING NTUPLES
            
                infilename=infile[:-5]    
                outfile=TFile(opt.outfolder+'/'+infilename+'_Run'+str(run_count)+'.root', 'RECREATE') #OUTPUT NAME
                #fname=opt.outfolder+'/'+infilename+'_Run'+str(run_count)+'.root'
                #fname=outfile.
                #print('opening file %s'%(fname))
                outfile.mkdir('event_info')
                SaveValues(params, outfile) ## SAVE PARAMETERS OF THE RUN
                    
                final_imgs=list();
                
                if opt.events==-1:
                    totev=tree.GetEntries()
                else:
                    
                    if opt.events<=tree.GetEntries():
                        totev=opt.events
                    else:
                        totev=tree.GetEntries()
                    
                for entry in range(0, totev): #RUNNING ON ENTRIES
                    tree.GetEntry(entry)
                    event_dict={'partID_'+str(entry): tree.pdgID_hits[0], 'E_init_'+str(entry): tree.ekin_particle[0]*1000}
                    #print(str(tree.pdgID_hits[0])+" "+str(tree.ekin_particle[0]*1000))
                    SaveEventInfo(event_dict, 'event_info', outfile)
                    #final_imgs.append(TH2F('pic_run'+str(run_count)+'_ev'+str(entry), '', opt.z_pix, -opt.z_dim*0.5, opt.z_dim*0.5, opt.y_pix, -opt.y_dim*0.5, opt.y_dim*0.5)) #smeared track with background
                    #final_image=TH2I('pic_run'+str(run_count)+'_ev'+str(entry), '', opt.z_pix, -opt.z_dim*0.5, opt.z_dim*0.5, opt.y_pix, -opt.y_dim*0.5, opt.y_dim*0.5) #smeared track with background
                    final_image=TH2I('pic_run'+str(run_count)+'_ev'+str(entry), '', opt.z_pix, 0, opt.z_pix-1, opt.y_pix, 0, opt.y_pix-1) #smeared track with background

                    #signal=AddTrack(tree, final_imgs, smearing)
                    signal=AddTrack(tree, final_image, smearing)
                    background=AddBckg(opt,entry+1)
                    total=signal+background

                    #final_imgs[entry]=rn.array2hist(total, final_imgs[entry])
                    final_image=rn.array2hist(total, final_image)
                    print('%d images generated'%(entry+1))
                    #final_imgs[entry].Write()
                    outfile.cd()
                    final_image.Write()            

                print('COMPLETED RUN %d'%(run_count))
                run_count+=1
                #outfile.Close()

        if opt.rootfiles==False:    
                # code to be used with input txt files from SRIM
                if infile.endswith('.txt'):    #KEEPING part.txt FILES ONLY NB: one single file with a specific name for the moment. there are other .txt files in the folder...this has to be fixed...

                    textfile=open(opt.infolder+infile, "r")
            
                    infilename=infile[:-4]    
                    outfile=TFile(opt.outfolder+'/'+infilename+'_Run'+str(run_count)+'.root', 'RECREATE') #OUTPUT NAME
                    #fname=opt.outfolder+'/'+infilename+'_Run'+str(run_count)+'.root'
                    #fname=outfile.
                    #print('opening file %s'%(fname))
                    outfile.mkdir('event_info')
                    SaveValues(params, outfile) ## SAVE PARAMETERS OF THE RUN

                    final_imgs=list();

                    zvec=list(); yvec=list(); xvec=list(); evec=list()
                    zvec*=0; yvec*=0; xvec*=0; evec*=0
                    lastentry=0
                    nhits=0
                    z_dist=TH1F('z_dist','z_dist; z [mm]; Entries/0.1mm',1000,-50,50)
                    y_dist=TH1F('y_dist','y_dist; y [mm]; Entries/0.1mm',1000,-50,50)
                    x_dist=TH1F('x_dist','x_dist; x [mm]; Entries/0.1mm',1000,-50,50)
                    de_dist=TH1F('de_dist','de_dist; dE_hit [keV]; Entries/0.1keV',1100,-10,100)
                    content = textfile.readlines() #READ IN TXT FILE
        
                    for nlines,line in enumerate(content):           #LOOP OVER LINES
                        myvars=line.split()
                        entry=int(myvars[0])-1
                        if entry!=lastentry or nlines==len(content)-1 : #NEW EVENT FOUND (OR LAST LINE IN THE FILE) - STORE INFORMATIONS ON PREVIOUS ONE
                            #event_dict={'partID_'+str(entry): tree.pdgID_hits[0], 'E_init_'+str(entry): tree.ekin_particle[0]*1000}
                            #SaveEventInfo(event_dict, 'event_info', outfile)
                            #final_imgs.append(TH2F('pic_run'+str(run_count)+'_ev'+str(lastentry), '', opt.z_pix, -opt.z_dim*0.5, opt.z_dim*0.5, opt.y_pix, -opt.y_dim*0.5, opt.y_dim*0.5)) #smeared track with background
                            #final_image=TH2I('pic_run'+str(run_count)+'_ev'+str(lastentry), '', opt.z_pix, -opt.z_dim*0.5, opt.z_dim*0.5, opt.y_pix, -opt.y_dim*0.5, opt.y_dim*0.5)#smeared track with background
                            final_image=TH2I('pic_run'+str(run_count)+'_ev'+str(lastentry), '', opt.z_pix, 0., opt.z_pix, opt.y_pix, 0., opt.y_pix)#smeared track with background
                           
                            signal=AddTrackGen(nhits, xvec, yvec, zvec, evec, final_image, smearing)
                            background=AddBckg(opt,entry+1)
                            total=signal+background
                            #final_imgs[lastentry]=rn.array2hist(total, final_imgs[lastentry])
                            final_image=rn.array2hist(total, final_image)
                            print('%d images generated with new code'%(lastentry))
                            #final_imgs[lastentry].Write()
                            outfile.cd()
                            final_image.Write()
                            zvec*=0; yvec*=0; xvec*=0; evec*=0 #RESET LISTS
                            lastentry=entry #END OF THE EVENT

                            if lastentry>=opt.events and opt.events!=-1:
                                break
                        nhits=int(myvars[1])
                        zvec.append(float(myvars[4]))
                        z_dist.Fill(float(myvars[4]))
                        yvec.append(float(myvars[3]))
                        y_dist.Fill(float(myvars[3]))
                        xvec.append(float(myvars[2]))
                        x_dist.Fill(float(myvars[2]))
                        evec.append(float(myvars[5]))
                        de_dist.Fill(float(myvars[5]))
                                
            
                    print('COMPLETED RUN %d'%(run_count))
                    run_count+=1
                    z_dist.Write()
                    y_dist.Write()
                    x_dist.Write()
                    de_dist.Write()
                    #outfile.Close()


    t1=time.time()
    if opt.donotremove == False:
        sw.swift_rm_root_file(opt.tmpname)
    print('\n')
    print('Generation took %d seconds'%(t1-t0))

