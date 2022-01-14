from ROOT import TFile, TTree
from array import array
import pandas as pd


# example used to create a root tree in python:
#         https://root-forum.cern.ch/t/filling-branch-with-pyroot-limited/34744


# NOTE: this code doesn't create a ntuple for phi_ini, theta_ini and energy_ini yet.
#       So you can't use these quantities in the digitization code, if using a ROOT file converted from SRIM


# set here input file and input folder
infile="ionization_profile_He_10keV.txt"
infolder="/home/pietro/Analysis_TESI/srim2root/SRIM_files/QF/"



srim_df = pd.read_csv(infolder+infile, 
                  sep='\t', 
                  names=["ion_number",
                         "hit_number",
                         "primary_secondary",
                         "X",                          # [mm]
                         "Y",                          # [mm]
                         "Z",                          # [mm]
                         "energy_deposit",             # [keV]
                         "ionization_energy_deposit"]  # [keV]
                     )

# remove NaN rows
srim_df.dropna() 

infilename=infile[:-4]  
outfile=TFile('%s.root' % (infilename), 'RECREATE')       


outtree = TTree( 'nTuple', 'nTuple' )

max_hits = 99999   # max number of hits per event
eventnumber_data = array( 'i', [ 0 ] )
numhits_data = array( 'i', [ 0 ] )
particle_type_data =  array( 'i', [ 0 ] )
x_hits_data = array( 'f', max_hits*[ 0] )
y_hits_data = array( 'f', max_hits*[ 0] )
z_hits_data = array( 'f', max_hits*[ 0] )
energyDep_hits_data = array( 'f', max_hits*[ 0] )


# see the number of events in srim file
events= 5  ##srim_df.Ion_number.iloc[-1]

outtree.Branch( 'eventnumber', eventnumber_data, 'eventnumber/I' )
outtree.Branch( 'numhits', numhits_data, 'numhits/I' )
outtree.Branch( 'particle_type', particle_type_data, 'particle_type/I' )
outtree.Branch( 'x_hits', x_hits_data, 'x_hits[numhits]/F' )
outtree.Branch( 'y_hits', y_hits_data, 'y_hits[numhits]/F' )
outtree.Branch( 'z_hits', z_hits_data, 'z_hits[numhits]/F' )
outtree.Branch( 'energyDep_hits', energyDep_hits_data, 'energyDep_hits[numhits]/F' )


     
for ev in range(0,events):
    eventnumber_data[0] = ev                                             # event number
    numhits_data[0] = len(srim_df[srim_df["ion_number"]==ev+1].index)    # number of hits
    particle_type_data[0] = 1
    
    
    for hit in range(0,numhits_data[0]):
        x_hits_data[hit] = srim_df[srim_df.ion_number==ev+1]["X"].iloc[hit]
        y_hits_data[hit] = srim_df[srim_df.ion_number==ev+1]["Y"].iloc[hit]
        z_hits_data[hit] = srim_df[srim_df.ion_number==ev+1]["Z"].iloc[hit]
        energyDep_hits_data[hit] = srim_df[srim_df.ion_number==ev+1]["ionization_energy_deposit"].iloc[hit]
    
    outtree.Fill()
    print("Processing event n.: ", ev, end="\r")
    
outfile.Write()
outfile.Close()
