import os
import math
import time
import random
import optparse

import numpy as np
import root_numpy as rn
import ROOT as rt
from scipy.stats import expon, poisson
import diplib as dip  # to compiute max and min faster

# sys.path.append("../reconstruction")
import swiftlib as sw
from rebin import rebin2d

import csv

## FUNCTIONS DEFINITION

def generate_exposure_time_map(exposure_time, options):
    map = np.zeros((options.x_pix, options.y_pix), dtype=int)
    T = np.random.uniform(0, exposure_time + 180)
    if 0 < T < 180:
        row_threshold = int(T * options.y_pix / 180)
        map[:row_threshold, :] = 1
    elif 180 < T < exposure_time:
        map[:, :] = 1
    elif exposure_time < T < exposure_time + 180:
        row_threshold = int((T - (exposure_time)) * options.y_pix / 180)
        map[row_threshold:, :] = 1
    return map


def NelGM1_vectorized(N_ioniz_el):
    n_el_oneGEM = N_ioniz_el * 0
    if isinstance(N_ioniz_el, int):  # in case there is onyl one hit (very low energy)
        for j in range(0, int(round(N_ioniz_el))):
            nsec = expon(loc=0, scale=GEM1_gain).rvs() * extraction_eff_GEM1
            n_el_oneGEM += nsec
            n_el_oneGEM = N_ioniz_el
    else:
        for i, n in enumerate(N_ioniz_el):
            for j in range(0, int(round(N_ioniz_el[i]))):
                nsec = expon(loc=0, scale=GEM1_gain).rvs() * extraction_eff_GEM1
                n_el_oneGEM[i] += nsec

    return n_el_oneGEM


def NelGEM2_vectorized(energyDep, z_hit, options):
    n_ioniz_el = energyDep / options.ion_pot
    drift_l = np.abs(z_hit - options.z_gem)
    n_ioniz_el_mean = np.abs(n_ioniz_el * np.exp(-drift_l / options.absorption_l))
    poisson_distr = lambda x: poisson(x).rvs()
    n_ioniz_el = poisson_distr(n_ioniz_el_mean)

    # total number of secondary electrons considering the gain in the 2nd GEM foil
    n_tot_el = NelGM1_vectorized(n_ioniz_el) * GEM2_gain * extraction_eff_GEM2

    return np.round(n_tot_el)


def _compute_sigma(diff_const, diff_coeff, dz):
    return np.sqrt(diff_const + diff_coeff * dz / 10.0)


def _smear_vectorized(axis_hit, axis_sigma, nel):
    if isinstance(nel, float):
        return np.concatenate(
            [np.random.normal(loc=(axis_hit), scale=axis_sigma, size=int(nel))]
        )
    return np.concatenate(
        [
            np.random.normal(loc=(x_hit_i), scale=sigma_i, size=int(nel_i))
            for x_hit_i, sigma_i, nel_i in zip(axis_hit, axis_sigma, nel)
        ]
    )


def _smear(axis_hit, axis_sigma, nph):
    return np.random.normal(
        loc=(axis_hit),
        scale=axis_sigma,
        size=int(nph),
    )


def cloud_smearing3D_vectorized(x_hit, y_hit, z_hit, energyDep_hit, options):

    nel = NelGEM2_vectorized(energyDep_hit, z_hit, options)

    dz = np.abs(z_hit - options.z_gem)
    sigma_x = _compute_sigma(options.diff_const_sigma0T, options.diff_coeff_T, dz)
    sigma_y = _compute_sigma(options.diff_const_sigma0T, options.diff_coeff_T, dz)
    sigma_z = _compute_sigma(options.diff_const_sigma0L, options.diff_coeff_L, dz)

    X = _smear_vectorized(x_hit, sigma_x, nel)
    Y = _smear_vectorized(y_hit, sigma_y, nel)
    Z = _smear_vectorized(z_hit, sigma_z, nel)

    return X, Y, Z


def NelGEM2(energyDep, z_hit, options):
    n_ioniz_el = energyDep / options.ion_pot
    drift_l = np.abs(z_hit - options.z_gem)
    n_ioniz_el_mean = np.abs(n_ioniz_el * np.exp(-drift_l / options.absorption_l))
    primary = poisson(n_ioniz_el_mean)  # poisson distribution for primary electrons
    n_ioniz_el = primary.rvs()  # number of primary electrons
    n_el_oneGEM = 0  # number of secondary electrons
    gain1 = expon(
        loc=0, scale=GEM1_gain
    )  # exponential distribution for the GAIN in the first GEM foil
    for k in range(0, n_ioniz_el):
        nsec = (
            gain1.rvs() * extraction_eff_GEM1
        )  # number of secondary electrons in the first GEM multiplication for each ionization electron
        n_el_oneGEM += nsec

    # total number of secondary electrons considering the gain in the 2nd GEM foil
    n_tot_el = n_el_oneGEM * GEM2_gain * extraction_eff_GEM2

    return n_tot_el


def ph_smearing2D(x_hit, y_hit, z_hit, energyDep_hit, options):

    # Electrons in GEM2
    nel = NelGEM2(energyDep_hit, z_hit, options)

    # Photons in GEM3 (the factor A is added to be able to compare saturated and non-saturated results)
    nph = nel * options.A * GEM3_gain * omega * options.photons_per_el * options.counts_per_photon

    # Support values computed
    dz = np.abs(z_hit - options.z_gem)
    sigma = _compute_sigma(options.diff_const_sigma0T, options.diff_coeff_T, dz)

    # arrays of positions of produced photons
    X = _smear(x_hit, sigma, nph)
    Y = _smear(y_hit, sigma, nph)
    return X, Y

def Nph_saturation_vectorized(histo_cloud, options):
    Nph_array = np.zeros((histo_cloud.shape[0], histo_cloud.shape[1]))
    Nph_tot = 0
    nel_in = histo_cloud
    hin = (nel_in * options.A * GEM3_gain) / (1 + options.beta * GEM3_gain * nel_in)
    hout = np.sum(hin, axis=(2))
    return hout


def AddBckg(options, i):
    bckg_array = np.zeros((options.x_pix, options.y_pix))
    if options.bckg:
        if os.path.exists(
            "%s/histograms_Run%05d.root" % (options.bckg_path, int(options.noiserun))
        ):
            options.tmpname = "%s/histograms_Run%05d.root" % (
                options.bckg_path,
                int(options.noiserun),
            )
        elif sw.checkfiletmp(int(options.noiserun)):
            # FIXME
            options.tmpname = "/tmp/histograms_Run%05d.root" % int(options.noiserun)
        else:
            print(
                "Downloading file: "
                + sw.swift_root_file(options.tag, int(options.noiserun))
            )
            options.tmpname = sw.swift_download_root_file(
                sw.swift_root_file(options.tag, int(options.noiserun)),
                int(options.noiserun),
            )
        tmpfile = rt.TFile.Open(options.tmpname)
        if options.ped_rand==True:
            l = list(tmpfile.GetListOfKeys())
            tmphist = tmpfile.Get(l[int(random.random()*len(l))].GetName())
            bckg_array = rn.hist2array(tmphist)
        else:
            n_pics = len(
            [k.GetName() for k in tmpfile.GetListOfKeys() if "pic" in k.GetName()]
        )
            tmphist = tmpfile.Get("pic_run%05d_ev%d" % (int(options.noiserun), i % n_pics))
            bckg_array = rn.hist2array(tmphist)
    # print(bckg_array)

    return bckg_array


def SaveValues(par, out):

    out.cd()
    out.mkdir("param_dir")
    rt.gDirectory.ls()
    out.cd("param_dir")

    for k, v in par.items():
        if k != "tag" and k != "bckg_path" and k != "Vig_Map" and k != "ped_rand":
            h = rt.TH1F(k, "", 1, 0, 1)
            h.SetBinContent(1, v)
            h.Write()
    out.cd()


def round_up_to_even(f):
    return math.ceil(f / 2.0) * 2

def TrackVignetting (arrTr,xpix,ypix,VignMap):
    for i in range(0,ypix):
        for j in range(0,xpix):
            if(arrTr[i][j]!=0):

                arrTr[i][j]=round(arrTr[i][j]*VignMap.GetBinContent(VignMap.GetXaxis().FindBin(j),VignMap.GetYaxis().FindBin(i)),0)
    return arrTr


def compute_cmos_with_saturation(x_hits, y_hits, z_hits, e_hits, opt):
    # vectorized smearing
    S3D_x, S3D_y, S3D_z = cloud_smearing3D_vectorized(
        x_hits, y_hits, z_hits, e_hits, opt
    )

    # if there are no electrons on GEM3, just use empty image
    if S3D_x.size == 0:
        array2d_Nph = np.zeros((opt.x_pix, opt.y_pix))

    # if there are electrons on GEM3, apply saturation effect
    else:

        # numpy histo is faster than ROOT histo
        histo_cloud_entries = np.array([S3D_x, S3D_y, S3D_z]).transpose()

        histo_cloud_entries=histo_cloud_entries[histo_cloud_entries[:, 2].argsort()]

        xmin, xmax = dip.MaximumAndMinimum(S3D_x)
        ymin, ymax = dip.MaximumAndMinimum(S3D_y)
        zmin, zmax = dip.MaximumAndMinimum(S3D_z)

        deltaX = np.absolute(xmax-xmin)
        deltaY = np.absolute(ymax-ymin)

        # FIXME: create a function for the saturation loop
        # FIXME: find best value of maxvolume. 1e8 might not me the best one
        max_3Dhisto_volume=1*1e8      # (volume in number of voxels) that's around 0.5*1.6 GB of RAM
        deltaZ=max(2*opt.z_vox_dim , opt.z_vox_dim * max_3Dhisto_volume/(deltaX/opt.x_vox_dim)/(deltaY/opt.y_vox_dim))
        split_vals=np.arange(zmin, zmax, deltaZ)
        split_at = histo_cloud_entries[:, 2].searchsorted(split_vals)
        histo_cloud_entries_list = np.split(histo_cloud_entries, split_at)

        xbin_dim=opt.x_vox_dim #opt.x_dim/opt.x_pix
        ybin_dim=opt.y_vox_dim #opt.y_dim/opt.y_pix

        x_n_bin=round_up_to_even((xmax-xmin)/xbin_dim)
        y_n_bin=round_up_to_even((ymax-ymin)/ybin_dim)

        hout=np.zeros(shape=(x_n_bin-1, y_n_bin-1))
        i=0
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

    return array2d_Nph


def compute_cmos_without_saturation(x_hits, y_hits, z_hits, e_hits, opt):
    signal = rt.TH2I(
        "temp",
        "",
        opt.x_pix,
        0,
        opt.x_pix - 1,
        opt.y_pix,
        0,
        opt.y_pix - 1,
    )
    tot_ph_G3 = 0
    for x, y, z, e in zip(x_hits, y_hits, z_hits, e_hits):
        S2DX, S2DY = ph_smearing2D(x, y, z, e, opt)

        for sx, sy in zip(S2DX, S2DY):
            tot_ph_G3 += 1

            signal.Fill(
                (0.5 * opt.x_dim + sx) * opt.x_pix / opt.x_dim,
                (0.5 * opt.y_dim + sy) * opt.y_pix / opt.y_dim,
            )
    # print("tot num of sensor counts after GEM3 without saturation: %d"%(tot_ph_G3))
    return rn.hist2array(signal)

def is_NR(array, pdg):
    index = np.where(array > pdg)[0]
    if index.size > 0:
        return array[index[0]]
    else:
        return -999


######################################### MAIN EXECUTION ###########################################

if __name__ == "__main__":

    # CALL: python img_generator.py ConfigFile.txt -I <<INPUT_FOLDER>>

    parser = optparse.OptionParser("usage: %prog [options] arg1 arg2")

    parser.add_option(
        "-I",
        "--inputfolder",
        dest="infolder",
        default=os.getcwd() + "/src",
        help="specify the folder containing input files",
    )
    parser.add_option(
        "-O",
        "--outputfolder",
        dest="outfolder",
        default=os.getcwd() + "/out",
        help="specify the output destination folder",
    )

    parser.add_option("-f", "--inputfile", dest="inputfile", default=None, help="specify the specific file in inputfolder you want to digitize")
    parser.add_option("-R", "--GEANT4isotopes", dest="GEANT4isotopes", default=False, action='store_true', help="specify if you are digitizing tracks from radioactive isotope GEANT4 simulations")


    (opt, args) = parser.parse_args()

    config = open(args[0], "r")  # GET CONFIG FILE
    params = eval(config.read())  # READ CONFIG FILE

    for k, v in params.items():
        setattr(opt, k, v)

    ## fit from Fernando Amaro's single GEM gain measurement
    GEM1_gain = 0.0347 * np.exp((0.0209) * opt.GEM1_HV)
    GEM2_gain = 0.0347 * np.exp((0.0209) * opt.GEM2_HV)
    GEM3_gain = 0.0347 * np.exp((0.0209) * opt.GEM3_HV)
    print("GEM1_gain = %d" % GEM1_gain)
    print("GEM2_gain = %d" % GEM2_gain)
    print("GEM3_gain = %d" % GEM3_gain)

    ## dividing Fernando's to Francesco&Karolina's single GEM gain measurement
    extraction_eff_GEM1 = 0.87319885 * np.exp(-0.0020000000 * opt.GEM1_HV)
    extraction_eff_GEM2 = 0.87319885 * np.exp(-0.0020000000 * opt.GEM2_HV)
    extraction_eff_GEM3 = 0.87319885 * np.exp(-0.0020000000 * opt.GEM3_HV)
    print("extraction eff GEM1 = %f" % extraction_eff_GEM1)
    print("extraction eff GEM2 = %f" % extraction_eff_GEM2)
    print("extraction eff GEM3 = %f" % extraction_eff_GEM3)

    demag = opt.y_dim / opt.sensor_size
    a = opt.camera_aperture
    omega = 1.0 / math.pow((4 * (demag + 1) * a), 2)  # solid angle ratio
    # print(omega)

    #### CODE EXECUTION ####
    run_count = 1
    t0 = time.time()

    # check if 'fixed seed' for random distributions (for debugging purposes)
    if opt.fixed_seed == True:
        np.random.seed(seed=0)

    if not os.path.exists(opt.outfolder):  # CREATING OUTPUT FOLDER
        os.makedirs(opt.outfolder)

    dict_isotopes={}
    if opt.GEANT4isotopes:
        with open(f'{os.getcwd()}/Z_A_isotopes.csv','r') as data:
            for line in csv.reader(data):
                dict_isotopes[line[0]]=int('%05d'%(int(1000*int(line[1])+int(line[2]))))

    for infile in os.listdir(opt.infolder):  # READING INPUT FOLDER

        if infile.endswith(".root"):  # KEEPING .ROOT FILES ONLY
            if opt.inputfile:
                if not opt.inputfile==str(infile): continue       

            # FIXME

            rootfile = rt.TFile.Open(opt.infolder + "/" + infile)
            tree = rootfile.Get("nTuple")  # GETTING NTUPLES

            infilename = infile[:-5]

            #name of output file = digi_<name_of_input_file>
            #outfilename = '{}/digi_{}'.format(opt.outfolder,infile)

            #standard: name of output file = histograms_RunRRRRR.root (R run number)
            outfilename = "%s/histograms_Run%05d.root" % (opt.outfolder, run_count)

            ##for radioisotope simulation: histograms_RunZZAAANN.root (Z=atomic number, A=mass numbe$
            ##NOTE: this is a 7 digit run number, while reconstruction currently looks for 5
            #if opt.GEANT4isotopes and infilename.find('part')>0: 
            #    outfilename = '%s/histograms_Run%05d%02d.root' % (opt.outfolder,dict_isotopes[str(infile.split("_")[1].split(".")[0])],int(infilename.split('part')[1]))
            #elif opt.GEANT4isotopes and infilename.find('part')==-1: 
            #    outfilename='%s/histograms_Run%05d00.root' % (opt.outfolder,dict_isotopes[str(infile.split("_")[1].split(".")[0])])

            outfile = rt.TFile(outfilename, "RECREATE") 

            outfile.mkdir("event_info")
            SaveValues(params, outfile)  ## SAVE PARAMETERS OF THE RUN

            eventnumber = np.array([-999], dtype="int")
            particle_type = np.array([-999], dtype="int")
            energy = np.array([-999], dtype="float32")
            theta = np.array([-999], dtype="float32")
            phi = np.array([-999], dtype="float32")
            track_length_3D = np.empty((1), dtype="float32")
            x_vertex = np.empty((1), dtype="float32")
            y_vertex = np.empty((1), dtype="float32")
            z_vertex = np.empty((1), dtype="float32")
            x_vertex_end = np.empty((1), dtype="float32")
            y_vertex_end = np.empty((1), dtype="float32")
            z_vertex_end = np.empty((1), dtype="float32")
            x_min = np.empty((1), dtype="float32")
            x_max = np.empty((1), dtype="float32")
            y_min = np.empty((1), dtype="float32")
            y_max = np.empty((1), dtype="float32")
            z_min = np.empty((1), dtype="float32")
            z_max = np.empty((1), dtype="float32")
            px = np.empty((1), dtype="float32")
            py = np.empty((1), dtype="float32")
            pz = np.empty((1), dtype="float32")
            proj_track_2D = np.empty((1), dtype="float32")
            nhits_og = np.empty((1), dtype="int32")
            #xhits_og = np.empty((10000), dtype="float32")
            #yhits_og = np.empty((10000), dtype="float32")
            #zhits_og = np.empty((10000), dtype="float32")
            #EDepHit_og = np.empty((10000), dtype="float32")
            N_photons = np.empty((1), dtype="int32")

            outtree = rt.TTree("info_tree", "info_tree")
            outtree.Branch("eventnumber", eventnumber, "eventnumber/I")
            outtree.Branch("particle_type", particle_type, "particle_type/I")
            outtree.Branch("energy", energy, "energy/F")
            outtree.Branch("theta", theta, "theta/F")
            outtree.Branch("phi", phi, "phi/F")
            outtree.Branch("track_length_3D", track_length_3D, "track_length_3D/F")
            outtree.Branch("proj_track_2D", proj_track_2D, "proj_track_2D/F")
            outtree.Branch("x_vertex", x_vertex, "x_vertex/F")
            outtree.Branch("y_vertex", y_vertex, "y_vertex/F")
            outtree.Branch("z_vertex", z_vertex, "z_vertex/F")
            outtree.Branch("x_vertex_end", x_vertex_end, "x_vertex_end/F")
            outtree.Branch("y_vertex_end", y_vertex_end, "y_vertex_end/F")
            outtree.Branch("z_vertex_end", z_vertex_end, "z_vertex_end/F")
            outtree.Branch("x_min", x_min, "x_min/F")
            outtree.Branch("x_max", x_max, "x_max/F")
            outtree.Branch("y_min", y_min, "y_min/F")
            outtree.Branch("y_max", y_max, "y_max/F")
            outtree.Branch("z_min", z_min, "z_min/F")
            outtree.Branch("z_max", z_max, "z_max/F")
            outtree.Branch("N_photons", N_photons, "N_photons/I")
            outtree.Branch("px", px, "px/F")
            outtree.Branch("py", py, "py/F")
            outtree.Branch("pz", pz, "pz/F")
            outtree.Branch("nhits_og", nhits_og, "nhits_og/I")
            #outtree.Branch("xhits_og", xhits_og, "xhits_og[nhits_og]/F")
            #outtree.Branch("yhits_og", yhits_og, "yhits_og[nhits_og]/F")
            #outtree.Branch("zhits_og", zhits_og, "zhits_og[nhits_og]/F")
            #outtree.Branch("EDepHit_og", EDepHit_og, "EDepHit_og[nhits_og]/F")

            max_events = tree.GetEntries()
            totev = max_events if opt.events == -1 else opt.events
            totev = min(totev, max_events)

            VignMap=rt.TH2D()
            if(opt.Vignetting):
                VignFile=rt.TFile.Open('VignettingMap/'+opt.Vig_Map,'read')
                VignMap=VignFile.Get("normmap_lime")
                VignMap.Smooth
                print(VignMap)

            for entry in range(totev):  # RUNNING ON ENTRIES
                tree.GetEntry(entry)
                print("Entry %d of %d" % (entry, totev))#, end="\r")
                print("Energy %d keV" % (int(tree.energyDep)))
                if tree.energyDep>200: continue
                if tree.energyDep < opt.ion_pot:
                    background = AddBckg(opt, entry)
                    total = background
                    final_image = rt.TH2I(
                        "pic_run" + str(run_count) + "_ev" + str(entry),
                        "",
                        opt.x_pix,
                        0,
                        opt.x_pix - 1,
                        opt.y_pix,
                        0,
                        opt.y_pix - 1,
                    )  
                    final_image = rn.array2hist(total, final_image)

                    outfile.cd()
                    final_image.Write()
                    outtree.Fill()

                    continue


                if opt.NR==True:
                    x_hits_tr = np.array(tree.x_hits) + opt.x_offset
                    y_hits_tr = np.array(tree.y_hits) + opt.y_offset
                    z_hits_tr = np.array(tree.z_hits) + opt.z_offset
                else:
                    x_hits_tr = np.array(tree.z_hits) + opt.x_offset
                    y_hits_tr = np.array(tree.y_hits) + opt.y_offset
                    z_hits_tr = np.array(tree.x_hits) + opt.z_offset

                # add random Z to tracks
                if opt.randZ_range:
                    rand = (random.random() - 0.5) * (opt.randZ_range)
                    for ihit in range(0, tree.numhits):
                        z_hits_tr[ihit] += rand


                eventnumber[0] = tree.eventnumber
                # FIXME
                if opt.NR == True:
                    energy[0] = tree.ekin_particle
                    particle_type[0] = tree.particle_type
                else:
                    ##this would be the energy of the primary particle - equal to deposited energy only if it is completely contained in the sensitive volume
                    #energy[0] = tree.ekin_particle[0] * 1000
                    energy[0] = tree.energyDep
                    #particle_type[0] = 0
                    if tree.energyDep_NR>0: particle_type[0]=is_NR(np.array(tree.pdgID_hits),int(1e9))
                    else: particle_type[0]=tree.pdgID_hits[0] #this will tell us if the deposit was started by a photon or an electron


                proj_track_2D[0] = np.sum(
                    np.sqrt(
                        np.power(np.ediff1d(x_hits_tr), 2)
                        + np.power(np.ediff1d(y_hits_tr), 2)
                    )
                )

                # if there are at least 2 hits compute theta and phi
                if len(x_hits_tr) > 1:
                    phi[0] = np.arctan2(
                        (y_hits_tr[1] - y_hits_tr[0]), (x_hits_tr[1] - x_hits_tr[0])
                    )
                    theta[0] = np.arccos(
                        (z_hits_tr[1] - z_hits_tr[0])
                        / np.sqrt(
                            np.power((x_hits_tr[1] - x_hits_tr[0]), 2)
                            + np.power((y_hits_tr[1] - y_hits_tr[0]), 2)
                            + np.power((z_hits_tr[1] - z_hits_tr[0]), 2)
                        )
                    )
                else:
                    phi[0] = -999.0
                    theta[0] = -999.0

                track_length_3D[0] = np.sum(np.array(tree.tracklen_hits))
                #NOTE: all *_particle branches refer to the primary particle in GEANT4, not the track in the gas
                px[0] = np.array(tree.px_particle)[0]
                py[0] = np.array(tree.py_particle)[0]
                pz[0] = np.array(tree.pz_particle)[0]

                x_vertex[0] = (
                    (np.array(x_hits_tr)[0] + 0.5 * opt.x_dim)
                    * opt.x_pix
                    / opt.x_dim
                ) #in pixels
                y_vertex[0] = (
                    (np.array(y_hits_tr)[0] + 0.5 * opt.y_dim)
                    * opt.y_pix
                    / opt.y_dim
                ) #in pixels
                z_vertex[0] = abs(np.array(z_hits_tr)[0]-opt.z_gem) #distance from GEMs in mm

                x_vertex_end[0] = (
                    (np.array(x_hits_tr)[-1] + 0.5 * opt.x_dim)
                    * opt.x_pix
                    / opt.x_dim
                ) #in pixels
                y_vertex_end[0] = (
                    (np.array(y_hits_tr)[-1] + 0.5 * opt.y_dim)
                    * opt.y_pix
                    / opt.y_dim
                ) #in pixels
                z_vertex_end[0] = abs(np.array(z_hits_tr)[-1]-opt.z_gem) #distance from GEMs in mm

                x_min[0] = (np.min(x_hits_tr) + 0.5*opt.x_dim)*opt.x_pix/opt.x_dim
                x_max[0] = (np.max(x_hits_tr) + 0.5*opt.x_dim)*opt.x_pix/opt.x_dim
                y_min[0] = (np.min(y_hits_tr) + 0.5*opt.y_dim)*opt.y_pix/opt.y_dim
                y_max[0] = (np.max(y_hits_tr) + 0.5*opt.y_dim)*opt.y_pix/opt.y_dim
                z_min[0] = min(abs(np.max(z_hits_tr) - opt.z_gem),abs(np.min(z_hits_tr) - opt.z_gem))
                z_max[0] = max(abs(np.max(z_hits_tr) - opt.z_gem),abs(np.min(z_hits_tr) - opt.z_gem))

                nhits_og[0] = tree.numhits

                #outtree.Fill()

                #xhits_og = np.array(x_hits_tr) + opt.x_offset
                #yhits_og = np.array(tree.y_hits) + opt.y_offset
                #zhits_og = np.array(tree.z_hits) + opt.z_offset
                #EDepHit_og = np.array(tree.energyDep_hits)

                ## with saturation
                if opt.saturation:
                    array2d_Nph = compute_cmos_with_saturation(
                        x_hits_tr, y_hits_tr, z_hits_tr, np.array(tree.energyDep_hits), opt
                    )
                ## no saturation
                else:
                    array2d_Nph = compute_cmos_without_saturation(
                        x_hits_tr, y_hits_tr, z_hits_tr, np.array(tree.energyDep_hits), opt
                    )

                N_photons[0] = np.sum(array2d_Nph)
                background = AddBckg(opt, entry)
                if(opt.Vignetting):
                   array2d_Nph = TrackVignetting(array2d_Nph, opt.y_pix, opt.x_pix, VignMap)
                if(opt.exposure_time_effect):
                   array2d_Nph = generate_exposure_time_map(opt.exposure_time, opt) * array2d_Nph  
                total = array2d_Nph + background


                final_image = rt.TH2I(
                    "pic_run" + str(run_count) + "_ev" + str(entry),
                    "",
                    opt.x_pix,
                    0,
                    opt.x_pix - 1,
                    opt.y_pix,
                    0,
                    opt.y_pix - 1,
                )  # smeared track with background
                final_image = rn.array2hist(total, final_image)

                outtree.Fill()

                outfile.cd()
                final_image.Write()

            outfile.cd("event_info")
            outtree.Write()
            print("COMPLETED RUN %d" % (run_count))
            run_count += 1
            # outfile.Close()

    if opt.donotremove == False:
        sw.swift_rm_root_file(opt.tmpname)

    print("\nGeneration took %d seconds" % (time.time() - t0))
