import ROOT as rt
import multiprocessing as mp
import os
import math
import optparse
import time
import sys
import numpy as np
import root_numpy as rn
from root_numpy import hist2array, array2hist
import random
from scipy.stats import expon
from scipy.stats import poisson
import math
#sys.path.append("../reconstruction")
import swiftlib as sw
import ctypes
#import quaternion
from scipy.spatial.transform import Rotation as R

import requests

from rebin import * 

#from memory_profiler import profile

## FUNCTIONS DEFINITION

def NelGM1_vectorized(N_ioniz_el):
    n_el_oneGEM=N_ioniz_el*0
    if isinstance(N_ioniz_el, int):   # in case there is onyl one hit (very low energy)
        for j in range(0,int(round(N_ioniz_el))):
                nsec = expon(loc=0,scale=GEM1_gain).rvs()*extraction_eff_GEM1
                n_el_oneGEM += nsec
                n_el_oneGEM=N_ioniz_el*0
    else:
        for i, n in enumerate(N_ioniz_el):
            for j in range(0,int(round(N_ioniz_el[i]))):
                nsec = expon(loc=0,scale=GEM1_gain).rvs()*extraction_eff_GEM1
                n_el_oneGEM[i] += nsec

    return n_el_oneGEM


def NelGEM2_vectorized(energyDep,z_hit,options):
    n_ioniz_el=energyDep/options.ion_pot
    drift_l = np.abs(z_hit-options.z_gem)
    n_ioniz_el_mean = np.abs(n_ioniz_el*np.exp(-drift_l/options.absorption_l))
    poisson_distr = lambda x: poisson(x).rvs()
    n_ioniz_el = poisson_distr(n_ioniz_el_mean)

    # total number of secondary electrons considering the gain in the 2nd GEM foil
    n_tot_el=NelGM1_vectorized(n_ioniz_el)*GEM2_gain*extraction_eff_GEM2

    return np.round(n_tot_el)

#@profile
def cloud_smearing3D_vectorized(x_hit,y_hit,z_hit,energyDep_hit, options):

    nel=NelGEM2_vectorized(energyDep_hit,z_hit, options)

    X=np.array([])
    Y=np.array([])
    Z=np.array([])

    sigma_x = np.sqrt(options.diff_const_sigma0T+options.diff_coeff_T*(np.abs(z_hit-options.z_gem))/10.)
    sigma_y = np.sqrt(options.diff_const_sigma0T+options.diff_coeff_T*(np.abs(z_hit-options.z_gem))/10.)
    sigma_z = np.sqrt(options.diff_const_sigma0L+options.diff_coeff_L*(np.abs(z_hit-options.z_gem))/10.)

    if isinstance(nel, float):
        X = np.concatenate([np.random.normal(loc=(x_hit), scale=sigma_x, size=int(nel)) ])
        Y = np.concatenate([np.random.normal(loc=(y_hit), scale=sigma_y, size=int(nel)) ])
        Z = np.concatenate([np.random.normal(loc=(z_hit-z_ini), scale=sigma_z, size=int(nel)) ])
    else:
        X = np.concatenate([np.random.normal(loc=(x_hit_i), scale=sigma_i, size=int(nel_i)) for i, (x_hit_i, sigma_i, nel_i) in enumerate(zip(x_hit, sigma_x, nel))])
        Y = np.concatenate([np.random.normal(loc=(y_hit_i), scale=sigma_i, size=int(nel_i)) for i, (y_hit_i, sigma_i, nel_i) in enumerate(zip(y_hit, sigma_y, nel))])
        Z = np.concatenate([np.random.normal(loc=(z_hit_i-z_ini), scale=sigma_i, size=int(nel_i)) for i, (z_hit_i, sigma_i, nel_i) in enumerate(zip(z_hit, sigma_z, nel))])

    return X, Y, Z

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
def Nph_saturation(histo_cloud,options,xmin_vox,xmax_vox,ymin_vox,ymax_vox,zmin_vox,zmax_vox):
    Nph_array = np.zeros((histo_cloud.GetNbinsX(),histo_cloud.GetNbinsY()))
    Nph_tot = 0
    for i in range(xmin_vox, xmax_vox):
        for j in range(ymin_vox,ymax_vox):
            hout = 0
            for k in range(zmin_vox,zmax_vox):
                hin = histo_cloud.GetBinContent(i,j,k)
                nel_in = hin
                hout += (nel_in * options.A * GEM3_gain)/(1 + options.beta * GEM3_gain  * nel_in) 
                

            nmean_ph= hout * omega * options.photons_per_el * options.counts_per_photon     # mean total number of photons
            photons=poisson(nmean_ph)                    # poisson distribution for photons
            n_ph=photons.rvs()  
            Nph_array[i-1][j-1] = n_ph
            Nph_tot += Nph_array[i-1][j-1]
            #if hout>0:
            #    print("number final electrons per voxel: %f"%hout)
            
    return Nph_tot, Nph_array


#@profile
def Nph_saturation_vectorized(histo_cloud,options):
    Nph_array = np.zeros((histo_cloud.shape[0],histo_cloud.shape[1]),dtype="int16")
    Nph_tot = 0

    nel_in=histo_cloud
    if nel_in.nbytes/1024./1024./1024.>4. :
        hout=0
        nslices=nel_in.shape[2]
        for zhit in range(nslices):
            hout+=(nel_in[:,:,zhit]  * options.A * GEM3_gain)/(1 + options.beta * GEM3_gain  * nel_in[:,:,zhit]) 
    else:
        hin=(nel_in  * options.A * GEM3_gain)/(1 + options.beta * GEM3_gain  * nel_in) 
        hout=np.sum(hin,axis=(2))


    #nel_in=np.sum(histo_cloud,axis=(2))
    #del histo_cloud
    #hout=(nel_in  * options.A * GEM3_gain)/(1 + options.beta * GEM3_gain  * nel_in) 
    
    nmean_ph= hout * omega * options.photons_per_el * options.counts_per_photon     # mean total number of photons
    poisson_distr = lambda x: poisson(x).rvs()
    n_ph=poisson_distr(nmean_ph) 
    Nph_array=n_ph
    Nph_tot=np.sum(Nph_array)
    return Nph_tot, Nph_array


def AddBckg(options, i):
    bckg_array=np.zeros((options.x_pix,options.y_pix))
    if options.bckg:
        if os.path.exists("%s/histograms_Run%05d.root" % (options.bckg_path,int(options.noiserun))):
            options.tmpname = "%s/histograms_Run%05d.root" % (options.bckg_path,int(options.noiserun))
        elif sw.checkfiletmp(int(options.noiserun)):
            #FIXME
            options.tmpname = "/tmp/histograms_Run%05d.root" % int(options.noiserun)
        else:
            print ('Downloading file: ' + sw.swift_root_file(options.tag, int(options.noiserun)))
            options.tmpname = sw.swift_download_root_file(sw.swift_root_file(options.tag, int(options.noiserun)),int(options.noiserun))
        tmpfile =rt.TFile.Open(options.tmpname)
        n_pics = len([k.GetName() for k in tmpfile.GetListOfKeys() if 'pic' in k.GetName()])
        tmphist = tmpfile.Get("pic_run%05d_ev%d"% (int(options.noiserun),i%n_pics))
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
        if (k!='tag' and k!='bckg_path' and k!='sample_file' and k!='presigned_URLs' and k!='SRIM_file'):
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

def round_up_to_even(f):
    return math.ceil(f / 2.) * 2

def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)

def angle_between(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2'::

            >>> angle_between((1, 0, 0), (0, 1, 0))
            1.5707963267948966
            >>> angle_between((1, 0, 0), (1, 0, 0))
            0.0
            >>> angle_between((1, 0, 0), (-1, 0, 0))
            3.141592653589793
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

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
   
    # check if 'fixed seed' for random distributions (for debugging purposes)
    if (opt.fixed_seed==True):
        np.random.seed(seed=0)


    eventnumber = np.array([-999], dtype="int")
    particle_type = np.array([-999], dtype="int")
    energy_ini = np.array([-999], dtype="float32")
    theta_ini = np.array([-999], dtype="float32")
    phi_ini = np.array([-999], dtype="float32")
            ### saving track length##
    #param_tree = rt.TTree("param_tree","param_tree") #creating a tree

    track_length_3D=np.empty((1),dtype="float32")
    x_vertex=np.empty((1),dtype="float32")
    y_vertex=np.empty((1),dtype="float32")
    z_vertex=np.empty((1),dtype="float32")
    x_vertex_end=np.empty((1),dtype="float32")
    y_vertex_end=np.empty((1),dtype="float32")
    z_vertex_end=np.empty((1),dtype="float32")
    energy = np.empty((1),dtype="float32")
    px = np.empty((1),dtype="float32")
    py = np.empty((1),dtype="float32")
    pz = np.empty((1),dtype="float32")
    proj_track_2D= np.empty((1),dtype="float32")
    theta = np.empty((1),dtype = "float32")
    phi = np.empty((1),dtype = "float32")
    nhits_og = np.empty((1),dtype = "int32")
    xhits_og = np.empty((10000),dtype = "float32")
    yhits_og = np.empty((10000),dtype = "float32")
    zhits_og = np.empty((10000),dtype = "float32")
    EDepHit_og = np.empty((10000),dtype = "float32")

    
    if not os.path.exists(opt.outfolder): #CREATING OUTPUT FOLDER
        os.makedirs(opt.outfolder)

    if opt.presigned_URLs!='':
        presignedURL_list = open(opt.presigned_URLs)
        dict_presignedURL = eval(presignedURL_list.read())

        
    for infile in os.listdir(opt.infolder): #path to geant4 simulation files
        if infile.endswith('.root'):
            
            samplefile = rt.TFile.Open(opt.infolder+'/'+infile)
            sampletree = samplefile.Get('nTuple')
               
            #outfile=rt.TFile('%s/histograms_Run%05d.root' % (opt.outfolder,run_count), 'RECREATE') #OUTPUT NAME (only run number)
            outfilename = '{}/digi_{}'.format(opt.outfolder,infile) #TO CHANGE WITH FILENAME
            #outfilename = '%s/histograms_Run%01d%04d.root' % (opt.outfolder,int(tree.particle_type),int(tree.ekin_particle))
            outfile=rt.TFile(outfilename, 'RECREATE') #OUTPUT NAME (only run number)
            #outfile.mkdir('event_info')
            SaveValues(params, outfile) ## SAVE PARAMETERS OF THE RUN
            outtree = rt.TTree("info_tree", "info_tree")
            outtree.Branch("eventnumber", eventnumber, "eventnumber/I")
            outtree.Branch("particle_type", particle_type, "particle_type/I")
            outtree.Branch("energy_ini", energy_ini, "energy_ini/F")
            outtree.Branch("theta_ini", theta_ini, "theta_ini/F")
            outtree.Branch("phi_ini", phi_ini, "phi_ini/F")
            #outtree.Branch('track_length_3D',track_length_3D,"track_length_3D/F")
            outtree.Branch('proj_track_2D',proj_track_2D,"proj_track_2D/F")
            outtree.Branch('x_vertex',x_vertex,"x_vertex/F")
            outtree.Branch('y_vertex',y_vertex,"y_vertex/F")
            outtree.Branch('z_vertex',z_vertex,"z_vertex/F")
            outtree.Branch('x_vertex_end',x_vertex_end,"x_vertex_end/F")
            outtree.Branch('y_vertex_end',y_vertex_end,"y_vertex_end/F")
            outtree.Branch('z_vertex_end',z_vertex_end,"z_vertex_end/F")
            #outtree.Branch('energy',energy,"energy/F")
            #outtree.Branch('px',px,"px/F")
            #outtree.Branch('py',py,"py/F")
            #outtree.Branch('pz',pz,"pz/F")
            #outtree.Branch('theta',theta,"theta/F")
            #outtree.Branch('phi',phi,"phi/F") 
            outtree.Branch('nhits_og',nhits_og,"nhits_og/I")
            #outtree.Branch('xhits_og',xhits_og,"xhits_og[nhits_og]/F")
            #outtree.Branch('yhits_og',yhits_og,"yhits_og[nhits_og]/F")
            #outtree.Branch('zhits_og',zhits_og,"zhits_og[nhits_og]/F")
            #outtree.Branch('EDepHit_og',EDepHit_og,"EDepHit_og[nhits_og]/F")  
            
            track_counter = 0
            
            final_imgs=list();
            
            if opt.NR==True: 
                rootfile=rt.TFile.Open(opt.SRIM_file) #SRIM file for NR
                tree=rootfile.Get('nTuple')            #GETTING NTUPLES from SRIM file
            else:
                tree = sampletree.CloneTree()
                
            
            for sample_entry in range(sampletree.GetEntries()):
                
                entry=0
                if opt.NR==True:
                    if track_counter>=opt.events: break
                    sampletree.GetEntry(sample_entry) 
                    requestE = sampletree.energyDep_NR
                    if sampletree.pdgID_hits[0]<1000 or requestE>=200. or len(sampletree.x_hits)<=1: continue
                    tree_entry=0   
                    nofoundflag = 0
                    #find entry in tree that match energy of sample_entry in sampletree
                    entry=0
                    for entry in range (tree.GetEntries()):
                        
                        tree.GetEntry(entry)
                        candidateE = tree.ekin_particle
                        requestparticle = str(sampletree.pdgID_hits[0]-1e9)[0]
                        candidateparticle = str(tree.particle_type)
                        
                        if abs(requestE-candidateE)<0.001*requestE and requestparticle==candidateparticle:
                            print('found a track in SRIM file')
                            tree_entry = entry
                            break
                        else: 
                            if entry<tree.GetEntries()-1: continue #look for another suitable SRIM event
                            else: nofoundflag = 1
                            
                    if nofoundflag or tree_entry==0: continue  
                    else: track_counter+=1
                    
                    print('entry {} tree_entry {}'.format(entry,tree_entry))
                    z_ini = 0.
                    tree.GetEntry(tree_entry)
                    print('Digitizing a NR with Z={} and energy {}, from GEANT4 the energy is {}'.format(tree.particle_type,candidateE, requestE))
                
                else: 
                    if track_counter>=opt.events: break
                    tree.GetEntry(sample_entry)
                    entry=sample_entry
                    z_ini = 0. #default: 255., ambe: 69
                    #if tree.energyDep_NR>0.: continue
                    #if tree.eventnumber!=3080: continue
                    if tree.pdgID_hits[0]>=1e9: continue
                    #if tree.energyDep<0. or tree.energyDep>50.: continue
                    if len(tree.x_hits)<=1: continue
                    track_counter+=1
    
                rotated_points = np.empty((tree.numhits,3),dtype="float32")
   
                # traslation of tracks (random or based on a sampled distribution)
                if opt.NR==True:
                    x_hits_tr = tree.x_hits
                    y_hits_tr = tree.y_hits
                    z_hits_tr = tree.z_hits
                else:
                    ##default:
                    #x_hits_tr = np.array(tree.z_hits)
                    #y_hits_tr = np.array(tree.y_hits)+20.
                    #z_hits_tr = np.array(tree.x_hits)+255.
                    ##AmBe (the center of the detector is translated):
                    x_hits_tr = np.array(tree.z_hits,dtype="float16")-100.
                    y_hits_tr = np.array(tree.y_hits,dtype="float16")+168.
                    z_hits_tr = np.array(tree.x_hits,dtype="float16") #+69., 174?
                if opt.randZ_range:
                    if opt.NR==True :
                        randx = random.randrange(-170.,170.)
                        randy = random.randrange(-170.,170.)
                        randz = random.randrange(0.,opt.randZ_range)
                        print('random z: {}'.format(randz))
                        for ihit in range(0,tree.numhits):
                            x_hits_tr[ihit]+=randx
                            y_hits_tr[ihit]+=randy
                            z_hits_tr[ihit]+=randz
                    else:
                        rand = (random.random()-0.5)*(opt.randZ_range)
                        for ihit in range(0,tree.numhits):
                            z_hits_tr[ihit]+=rand
                            
                #eventnumber[0] = tree.eventnumber
                eventnumber[0] = track_counter
                
                print(f'track number {track_counter} of energy {tree.energyDep} with vertex x={x_hits_tr[0]} y={y_hits_tr[0]} z={z_hits_tr[0]}')
                    
                #angle=angle between original direction (1,0,0) and desired direction
                v1 = (1,0,0)
                v2 = (sampletree.z_hits[1],sampletree.y_hits[1],sampletree.x_hits[1])
                angle = angle_between(v1,v2)
                axis = np.cross(v1,v2)
                #axis=cross product between the two vectors
                #matrix of rotation given angle and axis
                M=R.from_rotvec(angle*unit_vector(axis)).as_matrix()
                #M=R.random().as_matrix()
                if opt.NR==True:
                    x_hits_tr = M[0][0]*tree.x_hits+M[0][1]*tree.y_hits+M[0][2]*tree.z_hits+sampletree.z_hits[0]
                    y_hits_tr = M[1][0]*tree.x_hits+M[1][1]*tree.y_hits+M[1][2]*tree.z_hits+sampletree.y_hits[0]
                    z_hits_tr = M[2][0]*tree.x_hits+M[2][1]*tree.y_hits+M[2][2]*tree.z_hits+sampletree.x_hits[0]+255. #GEM coordinate
                    print('new coord {} {} {}'.format(sampletree.z_hits[0],sampletree.y_hits[0],sampletree.x_hits[0]+255.))

    
    
                #FIXME
                if opt.NR==True: 
                    proj_track_2D[0]=np.sum(np.sqrt(np.power(np.ediff1d(np.array(x_hits_tr)),2)+np.power(np.ediff1d(np.array(y_hits_tr)),2)))
                    energy_ini[0] = tree.ekin_particle
                    particle_type[0] = tree.particle_type
                    #phi_ini[0] = tree.phi_ini
                    #theta_ini[0] = tree.theta_ini
                    phi_ini[0] = np.arctan2( (y_hits_tr[1]-y_hits_tr[0]),(x_hits_tr[1]-x_hits_tr[0]) )
                    theta_ini[0] = np.arccos( (z_hits_tr[1]-z_hits_tr[0]) / np.sqrt( np.power((z_hits_tr[1]-z_hits_tr[0]),2) + np.power((y_hits_tr[1]-y_hits_tr[0]),2) + np.power((x_hits_tr[1]-x_hits_tr[0]),2)) )
                    x_vertex[0]= np.array(x_hits_tr)[0]
                    y_vertex[0]= np.array(y_hits_tr)[0]
                    z_vertex[0]= np.array(z_hits_tr)[0]
                    x_vertex_end[0]= np.array(x_hits_tr)[-1]
                    y_vertex_end[0]= np.array(y_hits_tr)[-1]
                    z_vertex_end[0]= np.array(z_hits_tr)[-1]
                    xhits_og = np.array(x_hits_tr)
                    yhits_og = np.array(y_hits_tr)
                    zhits_og = np.array(z_hits_tr)
                    EDepHit_og = np.array(tree.energyDep_hits)
                    nhits_og[0]=tree.numhits
                else: 
                    proj_track_2D[0]=np.sum(np.sqrt(np.power(np.ediff1d(np.array(x_hits_tr)),2)+np.power(np.ediff1d(np.array(y_hits_tr)),2)))
                    #energy_ini[0] = tree.ekin_particle[0]*1000 #if ER generated as primary particles
                    energy_ini[0] = tree.energyDep
                    particle_type[0] = tree.pdgID_hits[0]
                    phi_ini[0] = -999.
                    theta_ini[0] = -999.
                    phi_ini[0] = np.arctan2( (tree.y_hits[1]-tree.y_hits[0]),(tree.z_hits[1]-tree.z_hits[0]) )
                    theta_ini[0] = np.arccos( (tree.x_hits[1]-tree.x_hits[0]) / np.sqrt( np.power((tree.x_hits[1]-tree.x_hits[0]),2) + np.power((tree.y_hits[1]-tree.y_hits[0]),2) + np.power((tree.z_hits[1]-tree.z_hits[0]),2)) )
                    track_length_3D[0]=np.sum(np.array(tree.tracklen_hits))
                    xhits_og = np.array(x_hits_tr)
                    yhits_og = np.array(tree.y_hits)
                    zhits_og = np.array(z_hits_tr)
                    EDepHit_og = np.array(tree.energyDep_hits)
                    px[0]= np.array(tree.px_particle)[0]
                    py[0]= np.array(tree.py_particle)[0]
                    pz[0]= np.array(tree.pz_particle)[0]
                    x_vertex[0]= np.array(x_hits_tr)[0]
                    #y_vertex[0]= (np.array(tree.y_vertex_hits)[0]+0.5*opt.y_dim)*opt.y_pix/opt.y_dim
                    #z_vertex[0]= (np.array(tree.z_vertex_hits)[0]+0.5*opt.x_dim)*opt.x_pix/opt.x_dim
                    y_vertex[0]= np.array(y_hits_tr)[0]
                    z_vertex[0]= np.array(z_hits_tr)[0]
                    x_vertex_end[0]= np.array(x_hits_tr)[-1]
                    y_vertex_end[0]= np.array(y_hits_tr)[-1]
                    z_vertex_end[0]= np.array(z_hits_tr)[-1]
                    energy[0]=energy_ini[0]
                    theta[0]=theta_ini[0]
                    phi[0]=phi_ini[0]
                    nhits_og[0]=tree.numhits
                
                outtree.Fill()
                
    
                ## with saturation
                if (opt.saturation):
    
                    # non-vectorized smearing
                    #S3D_x=np.array([])
                    #S3D_y=np.array([])
                    #S3D_z=np.array([])
                   # for ihit in range(0,tree.numhits):
                    #    print("Processing hit %d of %d"%(ihit,tree.numhits))
                    #    ## here swapping X with Z beacuse in geant the drift axis is X
                    #   if (opt.NR == True):
                    #        S3D = cloud_smearing3D(x_hits_tr[ihit],y_hits_tr[ihit],z_hits_tr[ihit],tree.energyDep_hits[ihit],opt)
                    #    else:
                    #        S3D = cloud_smearing3D(x_hits_tr[ihit],y_hits_tr[ihit],z_hits_tr[ihit],tree.energyDep_hits[ihit],opt)
                    #    S3D_x=np.append(S3D_x, S3D[0])
                    #    S3D_y=np.append(S3D_y, S3D[1])
                    #    S3D_z=np.append(S3D_z, S3D[2])
                    #print('exited hits loop')
    
                    # vectorized smearing
                    # if ER file need to swapp X with Z beacuse in geant the drift axis is X
                    if (opt.NR == True):
                        S3D_x, S3D_y, S3D_z = cloud_smearing3D_vectorized(np.array(x_hits_tr),np.array(y_hits_tr),np.array(z_hits_tr),np.array(tree.energyDep_hits),opt)
                    else:
                        S3D_x, S3D_y, S3D_z = cloud_smearing3D_vectorized(x_hits_tr,y_hits_tr,z_hits_tr,np.array(tree.energyDep_hits,dtype="float16"),opt)
                
                    # if there are no electrons on GEM3, just use empty image 
                    if S3D_x.size == 0: 
                        array2d_Nph = np.zeros((opt.x_pix,opt.y_pix))
                    
                    #elif len(S3D_x)>500: #loop su più pezzi di traccia, e somma alla fine gli istogrammi
                    # if there are electrons on GEM3, apply saturation effect 
                    else:

                        xmax=max(S3D_x)
                        xmin=min(S3D_x)
                        ymax=max(S3D_y)
                        ymin=min(S3D_y)
                        zmax=max(S3D_z)
                        zmin=min(S3D_z)
                        # arrays with positions of all electrons in cloud

                        # numpy histo is faster than ROOT histo
                        histo_cloud_entries=np.array(
                                [S3D_x,
                                 S3D_y , 
                                 S3D_z],dtype="float16").transpose()
                        print(f'histo_cloud_entries = {histo_cloud_entries.nbytes/1024./1024./1024.} GB')
    
                        xbin_dim=opt.x_dim/opt.x_pix #opt.x_vox_dim 
                        ybin_dim=opt.y_dim/opt.y_pix #opt.y_vox_dim
                        zbin_dim=opt.z_vox_dim
    
    
                        x_n_bin=round_up_to_even((xmax-xmin)/xbin_dim)-1
                        y_n_bin=round_up_to_even((ymax-ymin)/ybin_dim)-1
                        z_n_bin=round_up_to_even((zmax-zmin)/zbin_dim)-1
    
     

                        histo_cloud, edge = np.histogramdd(
                                histo_cloud_entries,
                                bins=(x_n_bin, y_n_bin, z_n_bin),
                                range=([xmin,xmax],[ymin,ymax],[zmin,zmax]),
                                normed=None, weights=None, density=None)
    
                        #if histo_cloud.nbytes/1024/1024/1024>6. : continue
                        # apply saturation vectorized function

                        result_GEM3 = Nph_saturation_vectorized(histo_cloud,opt)  

                #####OPPURE 
                       # histo_cloud=rt.TH3I('test','test',,xmin,xmax,y_n_bin,ymin,ymax,z_n_bin,zmin,zmax)
                       # print(f'histo: ({x_n_bin},{y_n_bin},{z_n_bin}), array: {histo_cloud_entries.shape}')
                       # histo_cloud=rn.array2hist(histo_cloud_entries,histo_cloud)
                       # result_GEM3=Nph_saturation(histo_cloud,opt,xmin,xmax,ymin,ymax,zmin,zmax)
                        
                        
                        array2d_Nph = result_GEM3[1]

                        tot_ph_G3 = np.sum(array2d_Nph)
    
                        #x_n_bin2=round_up_to_even((xmax-xmin)/xbin_dim)
                        #y_n_bin2=round_up_to_even((ymax-ymin)/ybin_dim)
    
                        #xmax_bin2=round((xmax)/xbin_dim)+int(opt.x_pix/2.)
                        #xmin_bin2=round((xmin)/xbin_dim)+int(opt.x_pix/2.)
                        #ymax_bin2=round((ymax)/ybin_dim)+int(opt.y_pix/2.)
                        #ymin_bin2=round((ymin)/ybin_dim)+int(opt.y_pix/2.)
                        xmax_bin2=round((xmax)/xbin_dim)+int(opt.x_pix/2.)
                        xmin_bin2=round((xmin)/xbin_dim)+int(opt.x_pix/2.)
                        ymax_bin2=round((ymax)/ybin_dim)+int(opt.y_pix/2.)
                        ymin_bin2=round((ymin)/ybin_dim)+int(opt.y_pix/2.)
                        #print(ymin_bin2)
                        #print(ymax_bin2)
                        #print(xmin_bin2)
                        #print(xmax_bin2)
    
                        #x_n_bin2=xmax_bin2-xmin_bin2
                        #y_n_bin2=ymax_bin2-ymin_bin2

                        #print(x_n_bin2)
                        #print(y_n_bin2)
    
                        #xedges2 = np.linspace(xmin, xmax, num=x_n_bin2)
                        #yedges2 = np.linspace(ymin, ymax, num=y_n_bin2)
                        
                        array2d_Nph=np.around(array2d_Nph)

                        # for rebinning we are using the code in this repo: https://github.com/jhykes/rebin
                        # not sure if we have to add an acknowledgement in the README, or do something else to respect the copyright/license 
                        #array2d_Nph = rebin2d(edge[0], edge[1], array2d_Nph, xedges2,  yedges2, interp_kind=3)  
                        #array2d_Nph=np.around(array2d_Nph)
    
                        #array2d_Nph=np.pad(array2d_Nph, ( ( round((opt.x_pix-x_n_bin2)/2.), round((opt.x_pix-x_n_bin2)/2.+1)),  ( round((opt.y_pix-y_n_bin2)/2.), round((opt.y_pix-y_n_bin2)/2.+1))), 'constant', constant_values=0) 
                        #print('add {} bins to the left, {} bins to the right, {} bins on the bottom and {} on the top. Total bins x={} total bins y={}'.format(xmin_bin2+1,opt.x_pix-xmax_bin2,ymin_bin2+1,opt.y_pix-ymax_bin2,xmin_bin2+1+opt.x_pix-xmax_bin2+x_n_bin2,ymin_bin2+1+opt.y_pix-ymax_bin2+y_n_bin2))
                        array2d_Nph=np.pad(array2d_Nph, ( ( xmin_bin2 , opt.x_pix-x_n_bin-xmin_bin2 ),  ( ymin_bin2 , opt.y_pix-y_n_bin-ymin_bin2 ) ), 'constant', constant_values=0)
                       
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
                            S2D = ph_smearing2D(x_hits_tr[ihit],y_hits_tr[ihit],z_hits_tr[ihit],tree.energyDep_hits[ihit],opt)
                        else:
                            S2D = ph_smearing2D(x_hits_tr[ihit],y_hits_tr[ihit],z_hits_tr[ihit],tree.energyDep_hits[ihit],opt)
                        
                        for t in range(0, len(S2D[0])):
                            tot_ph_G3+=1
    
                            signal.Fill((0.5*opt.x_dim+S2D[0][t])*opt.x_pix/opt.x_dim, (0.5*opt.y_dim+S2D[1][t])*opt.y_pix/opt.y_dim ) 
                    array2d_Nph=rn.hist2array(signal)
                    array2d_Nph = array2d_Nph 
                    #print("tot num of sensor counts after GEM3 without saturation: %d"%(tot_ph_G3))
                

                background=AddBckg(opt,entry)
                total=array2d_Nph+background
    
                final_image=rt.TH2I('pic_run'+str(run_count)+'_ev'+str(track_counter), '', opt.x_pix, 0, opt.x_pix-1, opt.y_pix, 0, opt.y_pix-1) #smeared track with background
              
                final_image=rn.array2hist(total, final_image)

    
                outfile.cd()
                final_image.Write()            
                #param_tree.Fill()
                #final_image.Clear()

        #param_tree.Write()
            #outfile.cd('event_info') 
            outtree.Write()
            print('COMPLETED RUN %d'%(run_count))
            run_count+=1
            outfile.Close()
        
            #Write output to cloud storage
            if opt.presigned_URLs!='':
                with open(outfilename, 'rb') as f:
                    files = {'file': (outfilename, f)}
                    url_out = dict_presignedURL[infile]
                    http_response = requests.post(url_out['url'], data=url_out['fields'], files=files)

    t1=time.time()
    if opt.donotremove == False:
        sw.swift_rm_root_file(opt.tmpname)
    print('\n')
    print('Generation took %d seconds'%(t1-t0))
