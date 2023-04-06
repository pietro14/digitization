import ROOT as rt
import multiprocessing as mp
import os
import math
import optparse
import time
import sys
import numpy as np
import root_numpy as rn
import random
from scipy.stats import expon
from scipy.stats import poisson
import math

import diplib as dip  # to compute max and min faster 
import importlib
from scipy.spatial.transform import Rotation as R


#sys.path.append("../reconstruction")
import swiftlib as sw

#from rebin import * 


## FUNCTIONS DEFINITION

def NelGM1_vectorized(N_ioniz_el):
    n_el_oneGEM=N_ioniz_el*0
    if isinstance(N_ioniz_el, int):   # in case there is onyl one hit (very low energy)
        for j in range(0,int(round(N_ioniz_el))):
                nsec = expon(loc=0,scale=GEM1_gain).rvs()*extraction_eff_GEM1
                n_el_oneGEM += nsec
                n_el_oneGEM=N_ioniz_el
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
    n_tot_el=n_el_oneGEM*GEM2_gain*extraction_eff_GEM2*options.A

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
    #nmean_ph= hout * omega * options.photons_per_el * options.counts_per_photon     # mean total number of photons
    #poisson_distr = lambda x: poisson(x).rvs()
    #n_ph=poisson_distr(nmean_ph) 
    #Nph_array=n_ph
    #Nph_tot=np.sum(n_ph)
            
    return hout #Nph_tot, Nph_array
    
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
        if (k!='tag' and k!='bckg_path' and k!='Vignetting' and k!='Vig_Map' and k!='NR_list' and k!='translation'):
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

def TrackVignetting (arrTr,xpix,ypix,VignMap):

    for i in range(0,ypix):
        for j in range(0,xpix):
            if(arrTr[i][j]!=0):

                arrTr[i][j]=round(arrTr[i][j]*VignMap.GetBinContent(VignMap.GetXaxis().FindBin(j),VignMap.GetYaxis().FindBin(i)),0)

    return arrTr

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
    
    SimCentre=[0,0,0]
    if "translation" in params.keys(): SimCentre=opt.translation

    SRIM_events=[]
    if opt.NR==True and opt.NR_list!='':
        NR_list=importlib.import_module(opt.NR_list)
        SRIM_events=NR_list.ionlist

    if not os.path.exists(opt.outfolder): #CREATING OUTPUT FOLDER
        os.makedirs(opt.outfolder)
            
    for infile in os.listdir(opt.infolder): #READING INPUT FOLDER
            
            if infile.endswith('.root'):    #KEEPING .ROOT FILES ONLY
                
                #FIXME
                z_ini = 0.
                rootfile=rt.TFile.Open(opt.infolder+"/"+infile)
                tree=rootfile.Get('nTuple')            #GETTING NTUPLES
            
                infilename=infile[:-5]    
                #outfile=rt.TFile('%s/histograms_Run%05d.root' % (opt.outfolder,run_count), 'RECREATE') #OUTPUT NAME (only run number)
                outfilename = '{}/digi_{}'.format(opt.outfolder,infile)
                outfile=rt.TFile(outfilename, 'RECREATE')
                outfile.mkdir('event_info')
                SaveValues(params, outfile) ## SAVE PARAMETERS OF THE RUN

                track_length_3D=np.empty((1),dtype="float32")
                x_vertex=np.empty((1),dtype="float32")
                y_vertex=np.empty((1),dtype="float32")
                z_vertex=np.empty((1),dtype="float32")
                x_vertex_end=np.empty((1),dtype="float32")
                y_vertex_end=np.empty((1),dtype="float32")
                z_vertex_end=np.empty((1),dtype="float32")
                gem_distance=np.empty((1),dtype="float32")
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

                outtree = rt.TTree("info_tree", "info_tree")
                outtree.Branch("eventnumber", eventnumber, "eventnumber/I")
                outtree.Branch("particle_type", particle_type, "particle_type/I")
                outtree.Branch("energy_ini", energy_ini, "energy_ini/F")
                outtree.Branch("theta_ini", theta_ini, "theta_ini/F")
                outtree.Branch("phi_ini", phi_ini, "phi_ini/F")
                outtree.Branch('proj_track_2D',proj_track_2D,"proj_track_2D/F")
                outtree.Branch('x_vertex',x_vertex,"x_vertex/F")
                outtree.Branch('y_vertex',y_vertex,"y_vertex/F")
                outtree.Branch('z_vertex',z_vertex,"z_vertex/F")
                outtree.Branch('gem_distance',gem_distance,'gem_distance/F')
                outtree.Branch('x_vertex_end',x_vertex_end,"x_vertex_end/F")
                outtree.Branch('y_vertex_end',y_vertex_end,"y_vertex_end/F")
                outtree.Branch('z_vertex_end',z_vertex_end,"z_vertex_end/F")
                outtree.Branch('nhits_og',nhits_og,"nhits_og/I")

                ### saving track length##
                '''
                param_tree = rt.TTree("param_tree","param_tree") #creating a tree

                param_tree.Branch('track_length_3D',track_length_3D,"track_length_3D/F")
                param_tree.Branch('proj_track_2D',proj_track_2D,"proj_track_2D/F")
                param_tree.Branch('x_vertex',x_vertex,"x_vertex/F")
                param_tree.Branch('y_vertex',y_vertex,"y_vertex/F")
                param_tree.Branch('z_vertex',z_vertex,"z_vertex/F")
                param_tree.Branch('x_vertex_end',x_vertex_end,"x_vertex_end/F")
                param_tree.Branch('y_vertex_end',y_vertex_end,"y_vertex_end/F")
                param_tree.Branch('z_vertex_end',z_vertex_end,"z_vertex_end/F")
                param_tree.Branch('energy',energy,"energy/F")
                param_tree.Branch('px',px,"px/F")
                param_tree.Branch('py',py,"py/F")
                param_tree.Branch('pz',pz,"pz/F")
                param_tree.Branch('theta',theta,"theta/F")
                param_tree.Branch('phi',phi,"phi/F") 
                param_tree.Branch('nhits_og',nhits_og,"nhits_og/I")
                param_tree.Branch('xhits_og',xhits_og,"xhits_og[nhits_og]/F")
                param_tree.Branch('yhits_og',yhits_og,"yhits_og[nhits_og]/F")
                param_tree.Branch('zhits_og',zhits_og,"zhits_og[nhits_og]/F")
                param_tree.Branch('EDepHit_og',EDepHit_og,"EDepHit_og[nhits_og]/F")  
                '''                

                final_imgs=list();
                
                if opt.events==-1:
                    totev=tree.GetEntries()
                else:
                    
                    if opt.events<=tree.GetEntries():
                        totev=opt.events
                    else:
                        totev=tree.GetEntries()
                    
                VignMap=rt.TH2D()
                if(opt.Vignetting):
                    VignFile=rt.TFile.Open('VignettingMap/'+opt.Vig_Map,'read')
                    VignMap=VignFile.Get("normmap_lime")
                    VignMap.Smooth
                    print(f'Will apply vignetting from {opt.Vig_Map} map')
    
                for entry in range(0, totev): #RUNNING ON ENTRIES

                    tree.GetEntry(entry)
                    print("Entry %d of %d." % (entry, totev))#, end="\r")

                    if opt.NR==True:
                        if tree.ekin_particle>300: continue
                        print(f' Event energy {tree.ekin_particle} keV')
                        v1 = (1,0,0)
                        v2=tuple(map(lambda i, j: i - j, (SRIM_events[entry][3],SRIM_events[entry][5],SRIM_events[entry][7]), (SRIM_events[entry][2],SRIM_events[entry][4],SRIM_events[entry][6])))
                        angle = angle_between(v1,v2)
                        axis = np.cross(v1,v2)
                        #axis=cross product between the two vectors
                        #rotation matrix given angle and axis
                        M=R.from_rotvec(angle*unit_vector(axis)).as_matrix()
                    
                        x_hits_tr = M[0][0]*tree.x_hits+M[0][1]*tree.y_hits+M[0][2]*tree.z_hits+SRIM_events[entry][2]+SimCentre[0]
                        y_hits_tr = M[1][0]*tree.x_hits+M[1][1]*tree.y_hits+M[1][2]*tree.z_hits+SRIM_events[entry][4]+SimCentre[1]
                        z_hits_tr = M[2][0]*tree.x_hits+M[2][1]*tree.y_hits+M[2][2]*tree.z_hits+SRIM_events[entry][6]+SimCentre[2]
                        
                    else:
                        ##default:
                        #x_hits_tr = np.array(tree.z_hits)
                        #y_hits_tr = np.array(tree.y_hits)+20.
                        #z_hits_tr = np.array(tree.x_hits)+255.
                        ##AmBe (the center of the detector is translated):
                        print(f'Event energy {tree.energyDep} keV')
                        x_hits_tr = np.array(tree.z_hits)-100.
                        y_hits_tr = np.array(tree.y_hits)+160.
                        z_hits_tr = np.array(tree.x_hits)

                    # add random Z to tracks
                    #x_hits_tr = tree.x_hits
                    if opt.randZ_range:
                        rand = (random.random()-0.5)*(opt.randZ_range)
                        for ihit in range(0,tree.numhits):
                            x_hits_tr[ihit]+=rand

                    print(f'event position {x_hits_tr[0]} {y_hits_tr[0]} {z_hits_tr[0]}')
                    eventnumber[0] = tree.eventnumber
                    #FIXME
                    if (opt.NR==True): 
                        proj_track_2D[0]=np.sum(np.sqrt(np.power(np.ediff1d(x_hits_tr),2)+np.power(np.ediff1d(y_hits_tr),2)))
                        energy_ini[0] = tree.ekin_particle
                        particle_type[0] = tree.particle_type
                        phi_ini[0] = np.arctan2( (y_hits_tr[1]-y_hits_tr[0]),(x_hits_tr[1]-x_hits_tr[0]) )
                        theta_ini[0] = np.arccos( (z_hits_tr[1]-z_hits_tr[0]) / np.sqrt( np.power((z_hits_tr[1]-z_hits_tr[0]),2) + np.power((y_hits_tr[1]-y_hits_tr[0]),2) + np.power((x_hits_tr[1]-x_hits_tr[0]),2)) )
                    else: 
                        proj_track_2D[0]=np.sum(np.sqrt(np.power(np.ediff1d(np.array(tree.z_hits)),2)+np.power(np.ediff1d(np.array(tree.y_hits)),2)))
                        energy_ini[0] = tree.energyDep
                        particle_type[0] = tree.pdgID_hits[0]
                        if len(x_hits_tr)>1:
                            phi_ini[0] = np.arctan2( (y_hits_tr[1]-y_hits_tr[0]),(x_hits_tr[1]-x_hits_tr[0]) )
                            theta_ini[0] = np.arccos( (z_hits_tr[1]-z_hits_tr[0]) / np.sqrt( np.power((z_hits_tr[1]-z_hits_tr[0]),2) + np.power((y_hits_tr[1]-y_hits_tr[0]),2) + np.power((x_hits_tr[1]-x_hits_tr[0]),2)) )
                        else: 
                            phi_ini[0] = -999
                            theta_ini[0] = -999
                    #phi_ini[0] = np.arctan2( (tree.y_hits[1]-tree.y_hits[0]),(tree.z_hits[1]-tree.z_hits[0]) )
                    #theta_ini[0] = np.arccos( (tree.x_hits[1]-tree.x_hits[0]) / np.sqrt( np.power((tree.x_hits[1]-tree.x_hits[0]),2) + np.power((tree.y_hits[1]-tree.y_hits[0]),2) + np.power((tree.z_hits[1]-tree.z_hits[0]),2)) )
                    #track_length_3D[0]=np.sum(np.array(tree.tracklen_hits))
                    xhits_og = np.array(x_hits_tr)
                    yhits_og = np.array(y_hits_tr)
                    zhits_og = np.array(z_hits_tr)
                    #EDepHit_og = np.array(tree.energyDep_hits)
                    #px[0]= np.array(tree.px_particle)[0]
                    #py[0]= np.array(tree.py_particle)[0]
                    #pz[0]= np.array(tree.pz_particle)[0]
                    x_vertex[0]= np.array(x_hits_tr)[0]
                    y_vertex[0]= np.array(y_hits_tr)[0]
                    z_vertex[0]= np.array(z_hits_tr)[0]
                    x_vertex_end[0]= np.array(x_hits_tr)[-1]
                    y_vertex_end[0]= np.array(y_hits_tr)[-1]
                    z_vertex_end[0]= np.array(z_hits_tr)[-1]
                    gem_distance[0]= np.abs(z_hits_tr[0]-opt.z_gem) #distance from GEMs of interaction point (in mm)
                    #y_vertex[0]= (np.array(tree.y_vertex_hits)[0]+0.5*opt.y_dim)*opt.y_pix/opt.y_dim
                    #z_vertex[0]= (np.array(tree.z_vertex_hits)[0]+0.5*opt.x_dim)*opt.x_pix/opt.x_dim
                    #x_vertex_end[0]= np.array(x_hits_tr)[-1]
                    #y_vertex_end[0]= (np.array(tree.y_vertex_hits)[-1]+0.5*opt.y_dim)*opt.y_pix/opt.y_dim
                    #z_vertex_end[0]= (np.array(tree.z_vertex_hits)[-1]+0.5*opt.x_dim)*opt.x_pix/opt.x_dim
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
                        #for ihit in range(0,tree.numhits):
                        #    #print("Processing hit %d of %d"%(ihit,tree.numhits))
                        #    ## here swapping X with Z beacuse in geant the drift axis is X
                        #    if (opt.NR == True):
                        #        S3D = cloud_smearing3D(x_hits_tr[ihit],tree.y_hits[ihit],tree.z_hits[ihit],tree.energyDep_hits[ihit],opt)
                        #    else:
                        #        S3D = cloud_smearing3D(tree.z_hits[ihit],tree.y_hits[ihit],x_hits_tr[ihit],tree.energyDep_hits[ihit],opt)
                        #    S3D_x=np.append(S3D_x, S3D[0])
                        #    S3D_y=np.append(S3D_y, S3D[1])
                        #    S3D_z=np.append(S3D_z, S3D[2])

                        # vectorized smearing
                        # if ER file need to swapp X with Z beacuse in geant the drift axis is X
                        #if (opt.NR == True):
                        #    S3D_x, S3D_y, S3D_z = cloud_smearing3D_vectorized(np.array(x_hits_tr),np.array(tree.y_hits),np.array(tree.z_hits),np.array(tree.energyDep_hits),opt)
                        #else:
                        S3D_x, S3D_y, S3D_z = cloud_smearing3D_vectorized(x_hits_tr,y_hits_tr,z_hits_tr,np.array(tree.energyDep_hits),opt)


                        # if there are no electrons on GEM3, just use empty image 
                        if S3D_x.size == 0: 
                            array2d_Nph = np.zeros((opt.x_pix,opt.y_pix))

                        # if there are electrons on GEM3, apply saturation effect 
                        else:

                            xmin, xmax = dip.MaximumAndMinimum(S3D_x)   
                            ymin, ymax = dip.MaximumAndMinimum(S3D_y)
                            #print(xmin)
                            #print(xmax)
                            #print(ymin)
                            #print(ymax)
                            #print(opt.x_dim)


                            zmin, zmax = dip.MaximumAndMinimum(S3D_z)
                            deltaX=np.absolute(xmax-xmin)
                            deltaY=np.absolute(ymax-ymin)

                            

                            histo_cloud_entries=np.array(
                                    [S3D_x,
                                     S3D_y , 
                                     S3D_z]).transpose()

                            histo_cloud_entries=histo_cloud_entries[histo_cloud_entries[:, 2].argsort()]


                            # FIXME: create a function for the saturation loop
                            max_3Dhisto_volume=1*1e8      # (volume in number of voxels) that's around 0.5*1.6 GB of RAM
                            deltaZ=max(2*opt.z_vox_dim , opt.z_vox_dim * max_3Dhisto_volume/(deltaX/opt.x_vox_dim)/(deltaY/opt.y_vox_dim))
                            split_vals=np.arange(zmin, zmax, deltaZ)
                            split_at = histo_cloud_entries[:, 2].searchsorted(split_vals)
                            histo_cloud_entries_list = np.split(histo_cloud_entries, split_at)

                            i=0

                            xbin_dim=opt.x_vox_dim #opt.x_dim/opt.x_pix
                            ybin_dim=opt.y_vox_dim #opt.y_dim/opt.y_pix

                            x_n_bin=round_up_to_even((xmax-xmin)/xbin_dim)
                            y_n_bin=round_up_to_even((ymax-ymin)/ybin_dim)

                            hout=np.zeros(shape=(x_n_bin-1, y_n_bin-1))
                            for histo_entries in histo_cloud_entries_list:
                                if np.size(histo_entries)==0:
                                            continue
                                zmin, zmax = dip.MaximumAndMinimum(histo_entries.transpose()[2])



                                zbin_dim=opt.z_vox_dim

                                z_n_bin=max(2,round_up_to_even((zmax-zmin)/zbin_dim))  

                                histo_cloud, edge = np.histogramdd(
                                        histo_entries,
                                        bins=(x_n_bin-1, y_n_bin-1, z_n_bin-1),
                                        range=([xmin,xmax],[ymin,ymax],[zmin,zmax]),
                                        normed=None, weights=None, density=None)

                                
                                hout=hout+Nph_saturation_vectorized(histo_cloud,opt)   

                                i=i+1

                            hout = hout * omega * opt.photons_per_el * opt.counts_per_photon     # mean total number of photons
                            poisson_distr = lambda x: poisson(x).rvs()
                            n_ph=poisson_distr(hout) 
                            array2d_Nph = n_ph
                            #print("Shape sat img:", n_ph.shape)

                            
                            # FIXME Write a function padding()
                            # Define a translation vector
                            x_center_cloud=int(np.round(((xmax+xmin)/2)/opt.x_vox_dim))
                            y_center_cloud=int(np.round(((ymax+ymin)/2)/opt.y_vox_dim))
                            #print("x_center_cloud",x_center_cloud)
                            #print("y_center_cloud",y_center_cloud)
                            translation = np.array([x_center_cloud, y_center_cloud])
                            # Calculate the center position of the original array in the padded array
                            center = np.array([int(opt.x_pix/2), int(opt.y_pix/2)]) + translation
                            #print("Center:", center)
                            # Create the padded array
                            padded_array = np.zeros((opt.x_pix, opt.y_pix))
                            x_start = max(0, center[0] - array2d_Nph.shape[0]//2)
                            y_start = max(0, center[1] - array2d_Nph.shape[1]//2)
                            x_end = min(opt.x_pix, x_start + array2d_Nph.shape[0])
                            y_end = min(opt.y_pix, y_start + array2d_Nph.shape[1])
                            #print("PADDING [%d:%d,%d:%d]" %(x_start, x_end, y_start, y_end))
                            #print(" ")
                            padded_array[x_start:x_end, y_start:y_end] = array2d_Nph
                            array2d_Nph=padded_array
                            


                    ## no saturation
                    else:
                        signal=rt.TH2I('sig_pic_run'+str(run_count)+'_ev'+str(entry), '', opt.x_pix, 0, opt.x_pix-1, opt.y_pix, 0, opt.y_pix-1) 
                        tot_ph_G3=0
                        for ihit in range(0,tree.numhits):
                            ## here swapping X with Z beacuse in geant the drift axis is X
                            if (opt.NR==True):
                                S2D = ph_smearing2D(x_hits_tr[ihit],tree.y_hits[ihit],tree.z_hits[ihit],tree.energyDep_hits[ihit],opt)
                            else:
                                S2D = ph_smearing2D(tree.z_hits[ihit],tree.y_hits[ihit],x_hits_tr[ihit],tree.energyDep_hits[ihit],opt)
                            
                            for t in range(0, len(S2D[0])):
                                tot_ph_G3+=1

                                signal.Fill((0.5*opt.x_dim+S2D[0][t])*opt.x_pix/opt.x_dim, (0.5*opt.y_dim+S2D[1][t])*opt.y_pix/opt.y_dim ) 
                        array2d_Nph=rn.hist2array(signal)
                        array2d_Nph = array2d_Nph 
                        #print("tot num of sensor counts after GEM3 without saturation: %d"%(tot_ph_G3))


                    background=AddBckg(opt,entry)
                    if(opt.Vignetting):
                        array2d_Nph=TrackVignetting(array2d_Nph,opt.y_pix,opt.x_pix,VignMap)
                    total=array2d_Nph+background

                    final_image=rt.TH2I('pic_run'+str(run_count)+'_ev'+str(entry), '', opt.x_pix, 0, opt.x_pix-1, opt.y_pix, 0, opt.y_pix-1) #smeared track with background
                    final_image=rn.array2hist(total, final_image)


                    outfile.cd()
                    final_image.Write()            
                    #param_tree.Fill()

                #param_tree.Write()
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
