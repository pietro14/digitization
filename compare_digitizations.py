#!/usr/bin/python
import numpy as np
import sys
import ROOT as rt

# This script compares the ouput of two different digitization. It can be useful for debugging, or to do a simple check after the code has been changed.

# The script takes in input two root files with the same number of digitized images (output of MC_data_new.py) and compute the difference bin by bin of each image.

# Note: in order to compare two root files, the same seed in the probabilità distribution must be set. In numpy, for instance:       np.random.seed(seed=0)

filename1=str(sys.argv[1])
filename2=str(sys.argv[2])

rootFile1 = rt.TFile.Open(filename1)
rootFile2 = rt.TFile.Open(filename2)

num_of_events1=rootFile1.Get("event_info/info_tree").GetEntries()
num_of_events2=rootFile2.Get("event_info/info_tree").GetEntries()


if num_of_events1 != num_of_events2:
    print("\n\nThe number of events in these two root files is different! I can't compare them!\n")
    print(filename1)
    print(filename2)
    print("\n")
    sys.exit()

for ev in range(num_of_events1):
    print("\nComparing event n.", ev , " ...")
    exec('pic1=rootFile1.pic_run1_ev'+str(ev))
    exec('pic2=rootFile2.pic_run1_ev'+str(ev))

    x_bins = pic1.GetNbinsX()
    y_bins = pic1.GetNbinsY()

    bins1 = np.zeros((x_bins,y_bins))

    for y_bin in range(y_bins): 
        for x_bin in range(x_bins): 
            bins1[x_bin,y_bin] = pic1.GetBinContent(x_bin + 1,y_bin + 1)

    x_bins = pic2.GetNbinsX()
    y_bins = pic2.GetNbinsY()

    bins2 = np.zeros((x_bins,y_bins))

    for y_bin in range(y_bins): 
        for x_bin in range(x_bins): 
            bins2[x_bin,y_bin] = pic2.GetBinContent(x_bin + 1,y_bin + 1)


    if (bins1 == bins2).all():
        print("Comparison result: equal\n")
    else:
        print("Comparison result: not equal\n")

