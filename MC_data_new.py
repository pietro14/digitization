import ROOT as rt
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

def NelGEM2(energyDep,z_hit,options):
    n_ioniz_el=energyDep/options.ion_pot
    drift_l = np.abs(z_hit-options.z_gem)
    n_ioniz_el_mean = np.abs(n_ioniz_el*np.exp(-drift_l/options.absorption_l)) 
    primary=poisson(n_ioniz_el_mean)           # poisson distribution for primary electrons
    n_ioniz_el=primary.rvs()                   # number of primary electrons
    n_el_oneGEM=0                              # number of secondary electrons
    gain1=expon(loc=0,scale=GEM1_gain)  # exponential distribution for the GAIN in the first GEM foil
    for k in range(0,n_ioniz_el):
        nsec = gain1.rvs()*extraction_eff_GEM1      # number of secondary electrons in the first GEM multiplication for each ionization electron
        n_el_oneGEM += nsec

    # total number of secondary electrons considering the gain in the 2nd GEM foil
    n_tot_el=n_el_oneGEM*GEM2_gain*extraction_eff_GEM2


    return n_tot_el

def cloud_smearing3D(x_hit,y_hit,z_hit,energyDep_hit,options):
    X=list(); Y=list(); Z=list() 
    X*=0; Y*=0; Z*=0
    nel=NelGEM2(energyDep_hit,z_hit,options)

    ## arrays of positions of produced electrons after GEM2
    X=(np.random.normal(loc=(x_hit), scale=np.sqrt(options.diff_const_sigma0T+options.diff_coeff_T*(np.abs(z_hit-options.z_gem))/10.), size=int(nel)))
    Y=(np.random.normal(loc=(y_hit), scale=np.sqrt(options.diff_const_sigma0T+options.diff_coeff_T*(np.abs(z_hit-options.z_gem))/10.), size=int(nel)))
    Z=(np.random.normal(loc=(z_hit-z_ini), scale=np.sqrt(options.diff_const_sigma0L+options.diff_coeff_L*(np.abs(z_hit-options.z_gem))/10.), size=int(nel)))   
    #print("distance from the GEM : %f cm"%((np.abs(z_hit-opt.z_gem))/10.))
    return X, Y, Z

def ph_smearing2D(x_hit,y_hit,z_hit,energyDep_hit,options):
    X=list(); Y=list()
    X*=0; Y*=0
    ## electrons in GEM2
    nel = NelGEM2(energyDep_hit,z_hit,options)
    ## photons in GEM3
    nph = nel * GEM3_gain *omega * options.photons_per_el * options.counts_per_photon
    ## arrays of positions of produced photons
    X=(np.random.normal(loc=(x_hit), scale=np.sqrt(options.diff_const_sigma0T+options.diff_coeff_T*(np.abs(z_hit-options.z_gem))/10.), size=int(nph)))
    Y=(np.random.normal(loc=(y_hit), scale=np.sqrt(options.diff_const_sigma0T+options.diff_coeff_T*(np.abs(z_hit-options.z_gem))/10.), size=int(nph)))
    return X, Y

# Nph_saturation() is not needed anymore: now we use Nph_saturation_vectorized() 
#
#def Nph_saturation(histo_cloud,options,xmin_vox,xmax_vox,ymin_vox,ymax_vox,zmin_vox,zmax_vox):
#    Nph_array = np.zeros((histo_cloud.GetNbinsX(),histo_cloud.GetNbinsY()))
#    Nph_tot = 0
#    for i in range(xmin_vox, xmax_vox):
#        for j in range(ymin_vox,ymax_vox):
#            hout = 0
#            for k in range(zmin_vox,zmax_vox):
#                hin = histo_cloud.GetBinContent(i,j,k)
#                nel_in = hin
#                hout += (nel_in * options.A * GEM3_gain)/(1 + options.beta * GEM3_gain  * nel_in) 
#                
#
#            nmean_ph= hout * omega * options.photons_per_el * options.counts_per_photon     # mean total number of photons
#            photons=poisson(nmean_ph)                    # poisson distribution for photons
#            n_ph=photons.rvs()  
#            Nph_array[i-1][j-1] = n_ph
#            Nph_tot += Nph_array[i-1][j-1]
#            #if hout>0:
#            #    print("number final electrons per voxel: %f"%hout)
#            
#    return Nph_tot, Nph_array


    
def Nph_saturation_vectorized(histo_cloud,options):
    Nph_array = np.zeros((histo_cloud.shape[0],histo_cloud.shape[1]))
    Nph_tot = 0

    nel_in=histo_cloud
    hin=(nel_in  * options.A * GEM3_gain)/(1 + options.beta * GEM3_gain  * nel_in) 
    hout=np.sum(hin,axis=(2))
    nmean_ph= hout * omega * options.photons_per_el * options.counts_per_photon     # mean total number of photons
    poisson_distr = lambda x: poisson(x).rvs()
    n_ph=poisson_distr(nmean_ph) 
    Nph_array=n_ph
    Nph_tot=np.sum(n_ph)
            
    return Nph_tot, Nph_array


def AddBckg(options, i):
    bckg_array=np.zeros((options.x_pix,options.y_pix))
    if options.bckg:
        if sw.checkfiletmp(int(options.noiserun)):
            options.tmpname = "/tmp/histograms_Run%05d.root" % int(options.noiserun)
            #FIXME
            #options.tmpname = "/nfs/cygno/users/dimperig/CYGNO/CYGNO-tmp/histograms_Run%05d.root" % int(options.noiserun)
        else:
            print ('Downloading file: ' + sw.swift_root_file(options.tag, int(options.noiserun)))
            options.tmpname = sw.swift_download_root_file(sw.swift_root_file(options.tag, int(options.noiserun)),int(options.noiserun))
        tmpfile =rt.TFile.Open(options.tmpname)
        tmphist = tmpfile.Get("pic_run%05d_ev%d"% (int(options.noiserun),i))
        bckg_array = rn.hist2array(tmphist) 
    #print(bckg_array)
    
    return bckg_array



def SaveValues(par, out):

    out.cd()
    out.mkdir('param_dir')
    #gDirectory.pwd()
    rt.gDirectory.ls()
    out.cd('param_dir')
    
    for k,v in par.items():
        if (k!='tag'):
            h=rt.TH1F(k, '', 1, 0, 1)
            h.SetBinContent(1, v)
            h.Write()
    out.cd()

    return None

def SaveEventInfo(info_dict, folder, out):  # THIS FUNCTION IS CURRENTLY NOT USED, MAYBE IT SHOULD BE REMOVED 
    
    out.cd()
    #gDirectory.pwd()
    out.cd('event_info')

    for k,v in info_dict.items():
        h=rt.TH1F(k, '', 1,0,1)
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
    params = eval(config.read())        #READ CONFIG FILE

    for k,v in params.items():
        setattr(opt, k, v)

    ## fit from Fernando Amaro's single GEM gain measurement
    GEM1_gain = 0.0347*np.exp((0.0209)*opt.GEM1_HV)
    GEM2_gain = 0.0347*np.exp((0.0209)*opt.GEM2_HV)
    GEM3_gain = 0.0347*np.exp((0.0209)*opt.GEM3_HV)
    print("GEM1_gain = %d"%GEM1_gain)
    print("GEM2_gain = %d"%GEM2_gain)
    print("GEM3_gain = %d"%GEM3_gain)
    
    ## dividing Fernando's to Francesco&Karolina's single GEM gain measurement
    extraction_eff_GEM1 = 0.87319885*np.exp(-0.0020000000*opt.GEM1_HV)
    extraction_eff_GEM2 = 0.87319885*np.exp(-0.0020000000*opt.GEM2_HV)
    extraction_eff_GEM3 = 0.87319885*np.exp(-0.0020000000*opt.GEM3_HV)
    print("extraction eff GEM1 = %f" % extraction_eff_GEM1 )
    print("extraction eff GEM2 = %f" % extraction_eff_GEM2 )
    print("extraction eff GEM3 = %f" % extraction_eff_GEM3 )

    demag=opt.y_dim/opt.sensor_size
    a=opt.camera_aperture
    omega=1./math.pow((4*(demag+1)*a),2)   # solid angle ratio
    #print(omega)

#### CODE EXECUTION ####
    run_count=1
    t0=time.time()
   
    # Set fixed seed for random distributions for debugging purposes
    np.random.seed(seed=0)


    eventnumber = np.array([-999], dtype="int")
    particle_type = np.array([-999], dtype="int")
    energy_ini = np.array([-999], dtype="float32")
    theta_ini = np.array([-999], dtype="float32")
    phi_ini = np.array([-999], dtype="float32")


    if not os.path.exists(opt.outfolder): #CREATING OUTPUT FOLDER
        os.makedirs(opt.outfolder)
            
    for infile in os.listdir(opt.infolder): #READING INPUT FOLDER
            
            if infile.endswith('.root'):    #KEEPING .ROOT FILES ONLY
                
                #FIXME
                z_ini = 255.
                zbins = int(opt.zcloud/opt.z_vox_dim)
                rootfile=rt.TFile.Open(opt.infolder+"/"+infile)
                tree=rootfile.Get('nTuple')            #GETTING NTUPLES
            
                infilename=infile[:-5]    
                outfile=rt.TFile('%s/histograms_Run%05d.root' % (opt.outfolder,run_count), 'RECREATE') #OUTPUT NAME (only run number)
                outfile.mkdir('event_info')
                SaveValues(params, outfile) ## SAVE PARAMETERS OF THE RUN
                outtree = rt.TTree("info_tree", "info_tree")
                outtree.Branch("eventnumber", eventnumber, "eventnumber/I")
                outtree.Branch("particle_type", particle_type, "particle_type/I")
                outtree.Branch("energy_ini", energy_ini, "energy_ini/F")
                outtree.Branch("theta_ini", theta_ini, "theta_ini/F")
                outtree.Branch("phi_ini", phi_ini, "phi_ini/F")

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
                    eventnumber[0] = tree.eventnumber
                    #FIXME
                    if (opt.NR==True): 
                        energy_ini[0] = tree.ekin_particle
                        particle_type[0] = tree.particle_type
                    else: 
                        energy_ini[0] = tree.ekin_particle[0]*1000
                        particle_type[0] = 0
                    phi_ini[0] = -999.
                    theta_ini[0] = -999.
                    phi_ini[0] = np.arctan2( (tree.y_hits[1]-tree.y_hits[0]),(tree.z_hits[1]-tree.z_hits[0]) )
                    theta_ini[0] = np.arccos( (tree.x_hits[1]-tree.x_hits[0]) / np.sqrt( np.power((tree.x_hits[1]-tree.x_hits[0]),2) + np.power((tree.y_hits[1]-tree.y_hits[0]),2) + np.power((tree.z_hits[1]-tree.z_hits[0]),2)) )
                    outtree.Fill()

                    ## with saturation
                    if (opt.saturation):

                        S3D_x=np.array([])
                        S3D_y=np.array([])
                        S3D_z=np.array([])

                        for ihit in range(0,tree.numhits):
                            #print("Processing hit %d of %d"%(ihit,tree.numhits))

                            ## here swapping X with Z beacuse in geant the drift axis is X
                            if (opt.NR == True):
                                S3D = cloud_smearing3D(tree.x_hits[ihit],tree.y_hits[ihit],tree.z_hits[ihit],tree.energyDep_hits[ihit],opt)
                            else:
                                S3D = cloud_smearing3D(tree.z_hits[ihit],tree.y_hits[ihit],tree.x_hits[ihit],tree.energyDep_hits[ihit],opt)

                            S3D_x=np.append(S3D_x, S3D[0])
                            S3D_y=np.append(S3D_y, S3D[1])
                            S3D_z=np.append(S3D_z, S3D[2])

                        # compute max and min for x,y,z. Those values define the smallest volume that contains all electrons
                        # +/-2 is needed for electrons close to the hedge
                        xmax=+2+int(round(max( (0.5*opt.x_dim+S3D_x)*opt.x_pix/opt.x_dim))) 
                        xmin=-2+int(round(min( (0.5*opt.x_dim+S3D_x)*opt.x_pix/opt.x_dim)))
                        ymax=+2+int(round(max( (0.5*opt.y_dim+S3D_y)*opt.y_pix/opt.y_dim) ))
                        ymin=-2+int(round(min( (0.5*opt.y_dim+S3D_y)*opt.y_pix/opt.y_dim))) 
                        zmax=+2+int(round(max( (0.5*zbins*opt.z_vox_dim+S3D_z)/opt.z_vox_dim))) 
                        zmin=-2+int(round(min( (0.5*zbins*opt.z_vox_dim+S3D_z)/opt.z_vox_dim)))

                        # numpy histo is faster than ROOT histo
                        histo_cloud_entries=np.array(
                                [(0.5*opt.x_dim+S3D_x)*opt.x_pix/opt.x_dim ,
                                 (0.5*opt.y_dim+S3D_y)*opt.y_pix/opt.y_dim , 
                                 (0.5*zbins*opt.z_vox_dim+S3D_z)/opt.z_vox_dim]).transpose()

                        histo_cloud, hedge = np.histogramdd(
                                histo_cloud_entries,
                                bins=(xmax-xmin,ymax-ymin,zmax-zmin),
                                range=([xmin,xmax],[ymin,ymax],[zmin,zmax]),
                                normed=None, weights=None, density=None)

                        # apply saturation vectorized function
                        result_GEM3 = Nph_saturation_vectorized(histo_cloud,opt)   
                        array2d_Nph = result_GEM3[1]
                        tot_ph_G3 = np.sum(array2d_Nph)
                        
                        # add empty block of pixels around the histo_cloud and get the dimensions of the final image
                        array2d_Nph=np.pad(array2d_Nph, ((xmin, opt.x_pix-xmax), (ymin, opt.y_pix-ymax)), 'constant', constant_values=0)
                           
                        #print("tot num of sensor counts after GEM3 including saturation: %d"%(tot_ph_G3))
                        #print("tot num of sensor counts after GEM3 without saturation: %d"%(opt.A*tot_el_G2*GEM3_gain* omega * opt.photons_per_el * opt.counts_per_photon))
                        #print("Gain GEM3 = %f   Gain GEM3 saturated = %f"%(GEM3_gain, tot_ph_G3/(opt.A * tot_el_G2*omega * opt.photons_per_el * opt.counts_per_photon) ))   

                    ## no saturation
                    else:
                        signal=rt.TH2I('sig_pic_run'+str(run_count)+'_ev'+str(entry), '', opt.x_pix, 0, opt.x_pix-1, opt.y_pix, 0, opt.y_pix-1) 
                        tot_ph_G3=0
                        for ihit in range(0,tree.numhits):
                            ## here swapping X with Z beacuse in geant the drift axis is X
                            if (opt.NR==True):
                                S2D = ph_smearing2D(tree.x_hits[ihit],tree.y_hits[ihit],tree.z_hits[ihit],tree.energyDep_hits[ihit],opt)
                            else:
                                S2D = ph_smearing2D(tree.z_hits[ihit],tree.y_hits[ihit],tree.x_hits[ihit],tree.energyDep_hits[ihit],opt)
                            
                            for t in range(0, len(S2D[0])):
                                tot_ph_G3+=1

                                signal.Fill((0.5*opt.x_dim+S2D[0][t])*opt.x_pix/opt.x_dim, (0.5*opt.y_dim+S2D[1][t])*opt.y_pix/opt.y_dim ) 
                        array2d_Nph=rn.hist2array(signal)
                        array2d_Nph = array2d_Nph 
                        #print("tot num of sensor counts after GEM3 without saturation: %d"%(tot_ph_G3))


                    background=AddBckg(opt,entry)
                    total=array2d_Nph+background

                    final_image=rt.TH2I('pic_run'+str(run_count)+'_ev'+str(entry), '', opt.x_pix, 0, opt.x_pix-1, opt.y_pix, 0, opt.y_pix-1) #smeared track with background
                    final_image=rn.array2hist(total, final_image)
                    outfile.cd()
                    final_image.Write()            

                outfile.cd('event_info') 
                outtree.Write()
                print('COMPLETED RUN %d'%(run_count))
                run_count+=1
                #outfile.Close()

    t1=time.time()
    if opt.donotremove == False:
        sw.swift_rm_root_file(opt.tmpname)
    print('\n')
    print('Generation took %d seconds'%(t1-t0))
