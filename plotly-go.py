#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 14 14:15:48 2023

@author: dozeduck
"""
import getopt
import sys
import re
# import plotly
import plotly.graph_objs as go
import plotly.io as pio
# for PCA
import numpy as np
from rpy2.robjects import r

args=sys.argv[1:]
file1 = ''
file2 = ''
file3 = ''
renumber = 'fault'
ave = 'fault'
output_name = 'plotly.png'
font_size = 40
font_family = 'Arial'
font_color = 'black'
xaxis_name = 0
yaxis_name = 0
rdf_cutoff = 0
multi_files = 0
plot_name = ''
pca = 0
rscript = 0
nbin = 500
size = 400

try:
    opts, args_lala = getopt.getopt(args,"f:s:t:r:a:o:x:y:c:m:p:l:j:n:z:h",["file1=",
                                             "file2=",
                                             "file3=",
                                             "renumber=",
                                             "average=",
                                             "output_name=",
                                             "xaxis_name=",
                                             "yaxis_name=",
                                             "rdf_cutoff=",
                                             "multi_files=",
                                             "plot_name=",
                                             "x_file=",
                                             "pca",
                                             "nbin=",
                                             "size=",
                                             "help"])
except getopt.GetoptError:
    print('Version: 1.3 \
          \nCurrently it works for rmsd, rmsf, sasa_time, sasa_residue, gyration, dipole movement, rdf, distance, PCA  \
          \nIt can read one, two or three same type files \
          \n-o output_name.png, suitable file format for output: png, jpg, jpeg, webp, svg, pdf, eps, json \
          \n-r true: default is fault, renumber the residues, mainly for the duplicated residue numbers while work on rmsf and sasa per residue \
          \n-c number: rdf distance cutoff value \
          \n-a 2 or 3 or true: default is fault, output the average score for each replica \
          \n-p plot_name: the title shoed in plot \
          \n-j whether reading cluster_PCA0.xvg cluster_PCA1.xvg ..., default is 0, if you want generate PCA contour heat map, then type 1 \
          \n-n represent "nbin", mainly used for pca ploting, default value=1000, the larger the smoother, however, if there is white line in your PCA_Density plot, please change the value to get a cleaner plot \
          \n-z represent "size", mainly used for pca ploting, default value=500, the larger the higher resolution, if there is white line in your PCA Density plot, please change the value to get a cleaner plot \
          \nUsage: \
          \n ./plotly_go -f <file1> -s <file2> -t <file3> -o <output_name> -r <true/fault> -a <2/3/true> -x <xaxis_name> -y <yaxis_name> -c <rdf_cutoff> -p <plot_title> \
          \n ./plotly_go -m <file1> <file2> <file3>  -o <output_name> -r <true/fault> -a <2/3/true> -x <xaxis_name> -y <yaxis_name> -c <rdf_cutoff> -p <plot_title>  \
          \n ./plotly_go -f rmsd1.xvg -o rmsd1.png \
          \n ./plotly_go -f rmsd1.xvg -s rmsd2.xvg -o rmsd12.png \
          \n ./plotly_go -f rmsd1.xvg -s rmsd2.xvg -t rmsd3.xvg -o rmsd123.png \
          \n ./plotly_go -f rmsf1.xvg -s rmsf2.xvg -t rmsf3.xvg -o rmsf123.png -r true \
          \n ./plotly_go -f rmsf1.xvg -s rmsf2.xvg -t rmsf3.xvg -o rmsf123.png -a 3 \
          \n ./plotly_go -f rmsf1.xvg -s rmsf2.xvg -o rmsf12.png -a 2 \
          \n ./plotly_go -f sasa1.xvg -s sasa2.xvg -t sasa3.xvg -o sasa123.png -r true -x "Time (ns)" -y "SASA (nm<sup>2</sup>)" \
          \n ./plotly_go -f rdf1.xvg -s rdf2.xvg -t rdf3.xvg -o rdf123.png -c 4.7 -o rdf123.png \
          \n ./plotly_go -m sasa1.xvg sasa2.xvg sasa3.xvg sasa4.xvg sasa5.xvg sasa6.xvg sasa7.xvg sasa8.xvg sasa9.xvg -o multy_sasa_9.png -a true \
          \n ./plotly_go -m resarea1.xvg resarea2.xvg resarea3.xvg -o test_resare123.png -a true -r true \
          \n ./plotly_go -m 2dproj1.xvg -n 1000 -z 500 -o pca_myname.png')
    sys.exit(2)

for opt, arg in opts:
    if opt == '-h':
        print('Version: 1.3 \
          \nCurrently it works for rmsd, rmsf, sasa_time, sasa_residue, gyration, dipole movement, rdf, distance, PCA  \
          \nIt can read one, two or three same type files \
          \n-o output_name.png, suitable file format for output: png, jpg, jpeg, webp, svg, pdf, eps, json \
          \n-r true: default is fault, renumber the residues, mainly for the duplicated residue numbers while work on rmsf and sasa per residue\
          \n-a 2 or 3 or true: default is fault, output the average score for each replica \
          \n-c number: rdf distance cutoff value \
          \n-p plot_name: the title shoed in plot \
          \n-j whether reading cluster_PCA0.xvg cluster_PCA1.xvg ..., default is 0, if you want generate PCA contour heat map, then type 1 \
          \n-n represent "nbin", mainly used for pca ploting, default value=1000, the larger the smoother, however, if there is white line in your PCA_Density plot, please change the value(-n 500) to get a cleaner plot \
          \n-z represent "size", mainly used for pca ploting, default value=500, the larger the higher resolution, if there is white line in your PCA Density plot, please change the value to get a cleaner plot \
          \nUsage: \
          \n ./plotly_go -f <file1> -s <file2> -t <file3> -o <output_name> -r <true/fault> -a <2/3/true> -x <xaxis_name> -y <yaxis_name> -c <rdf_cutoff> -p <plot_title> \
          \n ./plotly_go -m <file1> <file2> <file3>  -o <output_name> -r <true/fault> -a <2/3/true> -x <xaxis_name> -y <yaxis_name> -c <rdf_cutoff> -p <plot_title> \
          \n ./plotly_go -f rmsd1.xvg -o rmsd1.png \
          \n ./plotly_go -f rmsd1.xvg -s rmsd2.xvg -o rmsd12.png \
          \n ./plotly_go -f rmsd1.xvg -s rmsd2.xvg -o rmsd12.png -a 2 \
          \n ./plotly_go -f rmsd1.xvg -s rmsd2.xvg -t rmsd3.xvg -o rmsd123.png \
          \n ./plotly_go -f rmsd1.xvg -s rmsd2.xvg -t rmsd3.xvg -o rmsd123.png -a 3\
          \n ./plotly_go -f sasa1.xvg -s sasa2.xvg -t sasa3.xvg -o sasa123.png -r true -x "Time (ns)" -y "SASA (nm<sup>2</sup>)" \
          \n ./plotly_go -f rdf1.xvg -s rdf2.xvg -t rdf3.xvg -o rdf123.png -c 4.7 -o rdf123.png \
          \n ./plotly_go -m sasa1.xvg sasa2.xvg sasa3.xvg sasa4.xvg sasa5.xvg sasa6.xvg sasa7.xvg sasa8.xvg sasa9.xvg -o multy_sasa_9.png -a true \
          \n ./plotly_go -m resarea1.xvg resarea2.xvg resarea3.xvg -o test_resare123.png -a true -r true \
          \n ./plotly_go -m 2dproj1.xvg -n 1000 -z 500 -o pca_myname.png')
        sys.exit()
    elif opt in ("-f", "--file1"):
        file1 = str(arg)
    elif opt in ("-s", "--file2"):
        file2 = str(arg)
    elif opt in ("-t", "--file3"):
        file3 = str(arg)
    elif opt in ("-r", "--renumber"):
        renumber = str(arg)
    elif opt in ("-a", "--average"):
        ave = arg
    elif opt in ("-o", "--output_name"):
        output_name = str(arg)
    elif opt in ("-x", "--xaxis_name"):
        xaxis_name = str(arg)
    elif opt in ("-y", "--yaxis_name"):
        yaxis_name = str(arg)
    elif opt in ("-c", "--rdf_cutoff"):
        rdf_cutoff = round(float(arg),1)
    elif opt in ("-l", "--multi_files"):
        multi_files = arg.split(',')
    elif opt in ("-p", "--plot_name"):
        plot_name = str(arg)
    # how to recognize the space seperated command line file inputs
    elif opt in ("-j", "--pca"):
        pca = int(arg)
    elif opt in ("-n", "--nbin"):
        nbin = int(arg)
    elif opt in ("-z", "--size"):
        size = int(arg)
    elif opt in ("-m", "--multi_files"):
        value = 1
        multi_files = []
        index = args.index(opt)
        for number in range(index, len(args)):
            if number+1 < len(args):
                if not args[number+1].startswith("-") and value <= len(args):
                    multi_files.append(args[number+1])
                    value += 1
                else:
                    value = value + 100000
        args.remove('-m')
        for value in multi_files:
            args.remove(value)
        try:
            opts, args_lala = getopt.getopt(args,"f:s:t:r:a:o:x:y:c:m:p:l:j:n:z:h",["file1=",
                                                     "file2=",
                                                     "file3=",
                                                     "renumber=",
                                                     "average=",
                                                     "output_name=",
                                                     "xaxis_name=",
                                                     "yaxis_name=",
                                                     "rdf_cutoff=",
                                                     "multi_files=",
                                                     "plot_name=",
                                                     "x_file=",
                                                     "pca",
                                                     "nbin=",
                                                     "size=",
                                                     "help"])
        except getopt.GetoptError:
            sys.exit(2)

        for opt, arg in opts:
            if opt in ("-f", "--file1"):
                file1 = str(arg)
            elif opt in ("-s", "--file2"):
                file2 = str(arg)
            elif opt in ("-t", "--file3"):
                file3 = str(arg)
            elif opt in ("-r", "--renumber"):
                renumber = str(arg)
            elif opt in ("-a", "--average"):
                ave = arg
            elif opt in ("-o", "--output_name"):
                output_name = str(arg)
            elif opt in ("-x", "--xaxis_name"):
                xaxis_name = str(arg)
            elif opt in ("-y", "--yaxis_name"):
                yaxis_name = str(arg)
            elif opt in ("-c", "--rdf_cutoff"):
                rdf_cutoff = round(float(arg),1)
            elif opt in ("-l", "--multi_files"):
                multi_files = arg.split(',')
            elif opt in ("-p", "--plot_name"):
                plot_name = str(arg)
            elif opt in ("-j", "--pca"):
                pca = int(arg)
            elif opt in ("-n", "--nbin"):
                nbin = int(arg)
            elif opt in ("-z", "--size"):
                size = int(arg)





class plotly_go():
    flag = ''
    sasa_flag = ''
    pca_flag = ''
    time1 = []
    values1 = []
    sd1 = []
    time2 = []
    values2 = []
    sd2 = []
    time3 = []
    values3 = []
    sd3 = []
    max_value = []
    average_value = []
    multi_flag = ''

    def __init__(self,file1, file2, file3, output_name, renumber, ave, xaxis_name, yaxis_name, rdf_cutoff, multi_files, plot_name, pca, nbin, size):
        if multi_files == 0:
            self.flag_recognizer(file1, file2, file3)
            if self.flag == 'rmsd':
                self.rmsd_averager(file1,file2,file3, ave, output_name, xaxis_name, yaxis_name)
            elif self.flag == 'rmsf':
                self.rmsf_averager(file1, file2,file3, renumber, ave, output_name, xaxis_name, yaxis_name)
            elif self.flag == 'sasa' and self.sasa_flag != '-or':
                self.sasa_averager(file1, file2, file3, ave, output_name, xaxis_name, yaxis_name)
            elif self.flag == 'sasa' and self.sasa_flag == '-or':
                self.sasa_res_averager(file1, file2, file3,renumber,ave, output_name, xaxis_name, yaxis_name)
            elif self.flag == 'gyrate':
                self.gyrate_averager(file1, file2,file3, ave, output_name, xaxis_name, yaxis_name)
            elif self.flag == 'dipoles':
                self.dipoles_averager(file1, file2,file3, ave, output_name, xaxis_name, yaxis_name)
            elif self.flag == 'distance':
                self.distance_averager(file1, file2,file3, ave, output_name, xaxis_name, yaxis_name)
            elif self.flag == 'rdf':
                self.rdf_averager(file1,file2,file3, ave, output_name, xaxis_name, yaxis_name, rdf_cutoff)

        elif len(multi_files) >=1:
            # print(multi_files)
            file1 = multi_files[0]
            self.flag_recognizer(file1, file2, file3)
            if self.pca_flag != 1:
                self.plotly_multy(self.flag, multi_files, xaxis_name, yaxis_name, renumber, ave, output_name, rdf_cutoff, plot_name)
            elif self.pca_flag == 1:
                self.plotly_pca(multi_files, xaxis_name, yaxis_name, renumber, ave, output_name, rdf_cutoff, plot_name, nbin, size)





    def flag_recognizer(self,file1, file2, file3):                                                   # first method to be called in __main__, used for creating object and charactors.
        #self.atomic_index.clear()                                                   # cleaveage the information of previous object before put new record into these charactors
        with open(file1, 'r') as f:
            lines = f.readlines()                                                 # "filename" = $PATH/crankpep_docking_results/PLDAYL_corrected_top_1.pdb
            if len(lines)  >= 3:
                try:
                    self.flag = lines[2].split()[5]
                    if self.flag == 'rms,' :
                        self.flag = 'rmsd'
                    elif self.flag == 'rmsf,':
                        self.flag = 'rmsf'
                    elif self.flag == 'sasa,':
                        self.flag = 'sasa'
                    elif self.flag == 'gyrate,':
                        self.flag = 'gyrate'
                    elif self.flag == 'dipoles,':
                        self.flag = 'dipoles'
                    elif self.flag == 'distance,':
                        self.flag = 'distance'
                    elif self.flag == 'rdf,':
                        self.flag = 'rdf'
                except:
                    pass

            if len(lines) >= 9 and '-or' in lines[8]:
                self.sasa_flag = '-or'

            if 'pca' in str(file1).lower() or '2dproj' in str(file1):
                self.pca_flag = 1




    def rmsd_averager(self, file1, file2, file3, ave, output_name, xaxis_name, yaxis_name):
        a = 18
        with open(file1, 'r') as f:
            lines = f.readlines()
            for i in range(len(lines)):
                if i >= a:
                    #print(lines[i])
                    self.time1.append(float(lines[i].split()[0]))
                    self.values1.append(float(lines[i].split()[1]))
        try:
            with open(file2, 'r') as f:
                lines = f.readlines()
                for i in range(len(lines)):
                    if i >= a:
                        #print(lines[i])
                        self.time2.append(float(lines[i].split()[0]))
                        self.values2.append(float(lines[i].split()[1]))
        except:
            self.values2 = ''
            pass
        try:
            with open(file3, 'r') as f:
                lines = f.readlines()
                for i in range(len(lines)):
                    if i >= a:
                        #print(lines[i])
                        self.time3.append(float(lines[i].split()[0]))
                        self.values3.append(float(lines[i].split()[1]))
        except:
            self.values3 = ''
            pass

        # calculate the average and max value for each replica
        try:
            for i in [self.values1, self.values2, self.values3]:
                self.max_value.append(max(i))
                self.average_value.append(sum(i)/len(i))
        except:
            pass

        # define time unit is ps or ns
        with open(file1, 'r') as f: # define time unit
            lines = f.readlines()
            if len(lines) >= 15 and '(ns)"' in lines[14]:
                timeflag = 'ns'
            elif len(lines) >= 15 and '(ps)"' in lines[14]:
                timeflag = 'ps'
        if timeflag == 'ps':
            divisor = 1000
            self.time1 = [x / divisor for x in self.time1]


        # Line chart
        plot_title = 'RMSD'
        if xaxis_name == 0:
            x_name = 'Time (ns)'
        else:
            x_name = xaxis_name
        if yaxis_name == 0:
            y_name = 'RMSD (nm)'
        else:
            y_name = yaxis_name
        self.plotly_noSD(output_name, plot_title, x_name, y_name, self.time1, self.values1, self.values2, self.values3, ave)




    def rmsf_averager(self, file1, file2, file3, renumber, ave, output_name, xaxis_name, yaxis_name):
        a = 17
        with open(file1, 'r') as f:
            lines = f.readlines()
            for i in range(len(lines)):
                if i >= a:
                    #print(lines[i])
                    self.time1.append(float(lines[i].split()[0]))
                    self.values1.append(float(lines[i].split()[1]))
        try:
            with open(file2, 'r') as f:
                lines = f.readlines()
                for i in range(len(lines)):
                    if i >= a:
                        #print(lines[i])
                        self.time2.append(float(lines[i].split()[0]))
                        self.values2.append(float(lines[i].split()[1]))
        except:
            self.values2 = ''
            pass
        try:
            with open(file3, 'r') as f:
                lines = f.readlines()
                for i in range(len(lines)):
                    if i >= a:
                        #print(lines[i])
                        self.time3.append(float(lines[i].split()[0]))
                        self.values3.append(float(lines[i].split()[1]))
        except:
            self.values3 = ''
            pass

        if renumber == 'true':
            discontinuous_positions = []
            for i in range(1, len(self.time1)):
                if self.time1[i] != self.time1[i-1]+1:
                    discontinuous_positions.append(i)
            for i in range(len(self.time1)):
                self.time1[i] = i+1

        # define axis is atom or residues
        if xaxis_name == 0:
            with open(file1, 'r') as f: # define time unit
                lines = f.readlines()
                if len(lines) >= 15 and '"Residue"' in lines[14]:
                    x_name = 'Residue Number'
                elif len(lines) >= 15 and '"Atom"' in lines[14]:
                    x_name = 'Atom Number'
        else:
            x_name = xaxis_name

        plot_title = 'RMS fluctuation'
        if yaxis_name == 0:
            y_name = 'RMSF (nm)'
        else:
            y_name = yaxis_name
        self.plotly_noSD(output_name, plot_title, x_name, y_name, self.time1, self.values1, self.values2, self.values3, ave)






    def sasa_averager(self, file1, file2, file3, ave, output_name, xaxis_name, yaxis_name):
        a = 24
        with open(file1, 'r') as f:
            lines = f.readlines()
            for i in range(len(lines)):
                if i >= a:
                    #print(lines[i])
                    self.time1.append(float(lines[i].split()[0]))
                    self.values1.append(float(lines[i].split()[1]))
        try:
            with open(file2, 'r') as f:
                lines = f.readlines()
                for i in range(len(lines)):
                    if i >= a:
                        #print(lines[i])
                        self.time2.append(float(lines[i].split()[0]))
                        self.values2.append(float(lines[i].split()[1]))
        except:
            self.values2 = ''
            pass
        try:
            with open(file3, 'r') as f:
                lines = f.readlines()
                for i in range(len(lines)):
                    if i >= a:
                        #print(lines[i])
                        self.time3.append(float(lines[i].split()[0]))
                        self.values3.append(float(lines[i].split()[1]))
        except:
            self.values3 = ''
            pass

        # define time unit is ps or ns
        with open(file1, 'r') as f: # define time unit
            lines = f.readlines()
            if len(lines) >= 15 and '(ns)"' in lines[14]:
                timeflag = 'ns'
            elif len(lines) >= 15 and '(ps)"' in lines[14]:
                timeflag = 'ps'
        if timeflag == 'ps':
            divisor = 1000
            self.time1 = [x / divisor for x in self.time1]


        # Line chart
        plot_title = 'Solvent Accessible Surface'
        if xaxis_name == 0:
            x_name = 'Time (ns)'
        else:
            x_name = xaxis_name
        if yaxis_name == 0:
            y_name = 'Area (nm<sup>2</sup>)'
        else:
            y_name = yaxis_name
        self.plotly_noSD(output_name, plot_title, x_name, y_name, self.time1, self.values1, self.values2, self.values3, ave)



    def sasa_res_averager(self, file1, file2, file3,renumber,ave, output_name, xaxis_name, yaxis_name):
        a = 25
        with open(file1, 'r') as f:
            lines = f.readlines()
            for i in range(len(lines)):
                if i >= a:
                    #print(lines[i])
                    self.time1.append(float(lines[i].split()[0]))
                    self.values1.append(float(lines[i].split()[1]))
                    self.sd1.append(float(lines[i].split()[2]))
        try:
            with open(file2, 'r') as f:
                lines = f.readlines()
                for i in range(len(lines)):
                    if i >= a:
                        #print(lines[i])
                        self.time2.append(float(lines[i].split()[0]))
                        self.values2.append(float(lines[i].split()[1]))
                        self.sd2.append(float(lines[i].split()[2]))
        except:
            self.values2 = ''
            self.sd2 = ''
            pass
        try:
            with open(file3, 'r') as f:
                lines = f.readlines()
                for i in range(len(lines)):
                    if i >= a:
                        #print(lines[i])
                        self.time3.append(float(lines[i].split()[0]))
                        self.values3.append(float(lines[i].split()[1]))
                        self.sd3.append(float(lines[i].split()[2]))
        except:
            self.values3 = ''
            self.sd3 = ''
            pass

        if renumber == 'true':
            discontinuous_positions = []
            for i in range(1, len(self.time1)):
                if self.time1[i] != self.time1[i-1]+1:
                    discontinuous_positions.append(i)
            for i in range(len(self.time1)):
                self.time1[i] = i+1


        # Line chart
        plot_title_noSD = 'Area per residue over the trajectory (noSD)'
        plot_title_SD = 'Area per residue over the trajectory (SD)'
        if xaxis_name == 0:
            x_name = 'Residue Number'
        else:
            x_name = xaxis_name
        if yaxis_name == 0:
            y_name = 'Area (nm<sup>2</sup>)'
        else:
            y_name = yaxis_name
        self.plotly_noSD(output_name, plot_title_noSD, x_name, y_name, self.time1, self.values1, self.values2, self.values3, ave)
        self.plotly_SD(output_name, plot_title_SD, x_name, y_name, self.time1, self.values1, self.values2, self.values3, self.sd1, self.sd2, self.sd3, ave)




    def gyrate_averager(self, file1, file2, file3, ave, output_name, xaxis_name, yaxis_name):
        a = 27
        with open(file1, 'r') as f:
            lines = f.readlines()
            for i in range(len(lines)):
                if i >= a:
                    #print(lines[i])
                    self.time1.append(float(lines[i].split()[0]))
                    self.values1.append(float(lines[i].split()[1]))
        try:
            with open(file2, 'r') as f:
                lines = f.readlines()
                for i in range(len(lines)):
                    if i >= a:
                        #print(lines[i])
                        self.time2.append(float(lines[i].split()[0]))
                        self.values2.append(float(lines[i].split()[1]))
        except:
            self.values2 = ''
            pass
        try:
            with open(file3, 'r') as f:
                lines = f.readlines()
                for i in range(len(lines)):
                    if i >= a:
                        #print(lines[i])
                        self.time3.append(float(lines[i].split()[0]))
                        self.values3.append(float(lines[i].split()[1]))
        except:
            self.values3 = ''
            pass

        try:
            for i in [self.values1, self.values2, self.values3]:
                self.max_value.append(max(i))
                self.average_value.append(sum(i)/len(i))
        except:
            pass

        # define time unit is ps or ns
        with open(file1, 'r') as f: # define time unit
            lines = f.readlines()
            if len(lines) >= 15 and '(ns)"' in lines[14]:
                timeflag = 'ns'
            elif len(lines) >= 15 and '(ps)"' in lines[14]:
                timeflag = 'ps'
        if timeflag == 'ps':
            divisor = 1000
            self.time1 = [x / divisor for x in self.time1]

        # Line chart
        plot_title = 'Radius of Gyration'
        if xaxis_name == 0:
            x_name = 'Time (ns)'
        else:
            x_name = xaxis_name
        if yaxis_name == 0:
            y_name = 'Rg (nm)'
        else:
            y_name = yaxis_name
        self.plotly_noSD(output_name, plot_title, x_name, y_name, self.time1, self.values1, self.values2, self.values3, ave)



    def dipoles_averager(self, file1, file2, file3, ave, output_name, xaxis_name, yaxis_name):
        a = 27
        with open(file1, 'r') as f:
            lines = f.readlines()
            for i in range(len(lines)):
                if i >= a:
                    #print(lines[i])
                    self.time1.append(float(lines[i].split()[0]))
                    self.values1.append(float(lines[i].split()[4]))
        try:
            with open(file2, 'r') as f:
                lines = f.readlines()
                for i in range(len(lines)):
                    if i >= a:
                        #print(lines[i])
                        self.time2.append(float(lines[i].split()[0]))
                        self.values2.append(float(lines[i].split()[4]))
        except:
            self.values2 = ''
            pass
        try:
            with open(file3, 'r') as f:
                lines = f.readlines()
                for i in range(len(lines)):
                    if i >= a:
                        #print(lines[i])
                        self.time3.append(float(lines[i].split()[0]))
                        self.values3.append(float(lines[i].split()[4]))
        except:
            self.values3 = ''
            pass

        try:
            for i in [self.values1, self.values2, self.values3]:
                self.max_value.append(max(i))
                self.average_value.append(sum(i)/len(i))
        except:
            pass

        # define time unit is ps or ns
        with open(file1, 'r') as f: # define time unit
            lines = f.readlines()
            if len(lines) >= 15 and '(ns)"' in lines[14]:
                timeflag = 'ns'
            elif len(lines) >= 15 and '(ps)"' in lines[14]:
                timeflag = 'ps'
        if timeflag == 'ps':
            divisor = 1000
            self.time1 = [x / divisor for x in self.time1]

        # Line chart
        plot_title = 'Total dipole moment of the simulation box vs. time'
        if xaxis_name == 0:
            x_name = 'Time (ns)'
        else:
            x_name = xaxis_name
        if yaxis_name == 0:
            y_name = 'Total Dipole Moment (Debye)'
        else:
            y_name = yaxis_name



        self.plotly_noSD(output_name, plot_title, x_name, y_name, self.time1, self.values1, self.values2, self.values3, ave)


    def distance_averager(self, file1, file2, file3, ave, output_name, xaxis_name, yaxis_name):
        a = 24
        with open(file1, 'r') as f:
            lines = f.readlines()
            for i in range(len(lines)):
                if i >= a:
                    #print(lines[i])
                    self.time1.append(float(lines[i].split()[0]))
                    self.values1.append(float(lines[i].split()[1]))
        try:
            with open(file2, 'r') as f:
                lines = f.readlines()
                for i in range(len(lines)):
                    if i >= a:
                        #print(lines[i])
                        self.time2.append(float(lines[i].split()[0]))
                        self.values2.append(float(lines[i].split()[1]))
        except:
            self.values2 = ''
            pass
        try:
            with open(file3, 'r') as f:
                lines = f.readlines()
                for i in range(len(lines)):
                    if i >= a:
                        #print(lines[i])
                        self.time3.append(float(lines[i].split()[0]))
                        self.values3.append(float(lines[i].split()[1]))
        except:
            self.values3 = ''
            pass

        # calculate the average and max value for each replica
        try:
            for i in [self.values1, self.values2, self.values3]:
                self.max_value.append(max(i))
                self.average_value.append(sum(i)/len(i))
        except:
            pass

        # define time unit is ps or ns
        with open(file1, 'r') as f: # define time unit
            lines = f.readlines()
            if len(lines) >= 15 and '(ns)"' in lines[14]:
                timeflag = 'ns'
            elif len(lines) >= 15 and '(ps)"' in lines[14]:
                timeflag = 'ps'
        if timeflag == 'ps':
            divisor = 1000
            self.time1 = [x / divisor for x in self.time1]


        # Line chart
        plot_title = 'Distance'
        if xaxis_name == 0:
            x_name = 'Time (ns)'
        else:
            x_name = xaxis_name
        if yaxis_name == 0:
            y_name = 'Distance (nm)'
        else:
            y_name = yaxis_name
        self.plotly_noSD(output_name, plot_title, x_name, y_name, self.time1, self.values1, self.values2, self.values3, ave)

    def rdf_averager(self, file1, file2, file3, ave, output_name, xaxis_name, yaxis_name, rdf_cutoff):
        a = 25
        with open(file1, 'r') as f:
            lines = f.readlines()
            for i in range(len(lines)):
                if i >= a:
                    if rdf_cutoff != 0 and float(lines[i].split()[0]) <= rdf_cutoff:
                    #print(lines[i])
                        self.time1.append(float(lines[i].split()[0]))
                        self.values1.append(float(lines[i].split()[1]))
                    elif rdf_cutoff == 0:
                        self.time1.append(float(lines[i].split()[0]))
                        self.values1.append(float(lines[i].split()[1]))
        try:
            with open(file2, 'r') as f:
                lines = f.readlines()
                for i in range(len(lines)):
                    if i >= a:
                        if rdf_cutoff != 0 and float(lines[i].split()[0]) <= rdf_cutoff:
                            #print(lines[i])
                            self.time2.append(float(lines[i].split()[0]))
                            self.values2.append(float(lines[i].split()[1]))
                        elif rdf_cutoff == 0:
                            self.time2.append(float(lines[i].split()[0]))
                            self.values2.append(float(lines[i].split()[1]))
        except:
            self.values2 = ''
            pass
        try:
            with open(file3, 'r') as f:
                lines = f.readlines()
                for i in range(len(lines)):
                    if i >= a:
                        if rdf_cutoff != 0 and float(lines[i].split()[0]) <= rdf_cutoff:
                            #print(lines[i])
                            self.time3.append(float(lines[i].split()[0]))
                            self.values3.append(float(lines[i].split()[1]))
                        elif rdf_cutoff == 0:
                            self.time3.append(float(lines[i].split()[0]))
                            self.values3.append(float(lines[i].split()[1]))
        except:
            self.values3 = ''
            pass

        # calculate the average and max value for each replica
        try:
            for i in [self.values1, self.values2, self.values3]:
                self.max_value.append(max(i))
                self.average_value.append(sum(i)/len(i))
        except:
            pass

        # define distance unit is nm



        # Line chart
        plot_title = 'RDF'
        if xaxis_name == 0:
            x_name = 'r (nm)'
        else:
            x_name = xaxis_name
        if yaxis_name == 0:
            y_name = 'g(r)'
        else:
            y_name = yaxis_name
        self.plotly_noSD(output_name, plot_title, x_name, y_name, self.time1, self.values1, self.values2, self.values3, ave)


    def plotly_noSD(self, output_file_name, plot_title, x_name, y_name, time1, values1, values2, values3, ave):
            trace1 = go.Scatter(x=time1, y=values1, name='Replica 1')
            try:
                trace2 = go.Scatter(x=time1, y=values2, name='Replica 2')
            except:
                trace2 = 0
                pass
            try:
                trace3 = go.Scatter(x=time1, y=values3, name='Replica 3')
            except:
                trace3 = 0
                pass
            # data = [trace1, trace2, trace3]
            data = []
            for i in trace1, trace2, trace3:
                if i != 0:
                    data.append(i)

            layout = go.Layout(title=plot_title, title_x=0.5, title_y=1, font=dict(size=24),
                               xaxis=dict(title=x_name, titlefont=dict(size=40, color='black', family='Arial'), zeroline=False, autorange=True,
                                          showgrid=True, gridwidth=1, gridcolor='rgba(235,240,248,100)', tickfont=dict(size=30)),
                               yaxis=dict(title=y_name, titlefont=dict(size=40, color='black', family='Arial'), zeroline=False, autorange=True,
                                          showgrid=True, gridwidth=1, gridcolor='rgba(235,240,248,100)', tickfont=dict(size=30)),
                               legend=dict(x=1, y=1, orientation='v', font=dict(size=30)), showlegend=True,
                               plot_bgcolor='rgba(255, 255, 255, 0.1)',
                               paper_bgcolor='rgba(255, 255, 255, 0.2)')
            fig = go.Figure(data=data, layout=layout)
            # fig.show()
            pio.write_image(fig, output_file_name)
            try:
                if int(ave) == 3:
                    average_value = [(values1[i] + values2[i] + values3[i]) / 3 for i in range(len(values1))]
                elif int(ave) == 2:
                    average_value = [(values1[i] + values2[i]) / 2 for i in range(len(values1))]
                trace_ave = trace1 = go.Scatter(x=time1, y=average_value, name='Mean Value')
                data_ave = [trace_ave]
                fig_ave = go.Figure(data=data_ave, layout=layout)
                pio.write_image(fig_ave, 'mean_'+output_file_name)
            except:
                pass


    def plotly_SD(self, output_file_name, plot_title, x_name, y_name, time1, values1, values2, values3, sd1, sd2, sd3, ave):
            trace1 = go.Scatter(x=time1, y=values1, name='Replica 1')
            error_y1 = dict(type='data', array=sd1, visible=True, thickness=1)
            trace1.update(error_y=error_y1)
            try:
                trace2 = go.Scatter(x=time1, y=values2, name='Replica 2')
                error_y2 = dict(type='data', array=sd2, visible=True, thickness=1)
                trace2.update(error_y=error_y2)
            except:
                trace2 = 0
                pass
            try:
                trace3 = go.Scatter(x=time1, y=values3, name='Replica 3')
                error_y3 = dict(type='data', array=sd3, visible=True, thickness=1)
                trace3.update(error_y=error_y3)
            except:
                trace3 = 0
                pass

            # data = [trace1, trace2, trace3]
            data = []
            for i in trace1, trace2, trace3:
                if i != 0:
                    data.append(i)
            layout = go.Layout(title=plot_title, title_x=0.5, title_y=1, font=dict(size=24),
                               xaxis=dict(title=x_name, titlefont=dict(size=40, color='black', family='Arial'), zeroline=False, autorange=True,
                                          showgrid=True, gridwidth=1, gridcolor='rgba(235,240,248,100)', tickfont=dict(size=30)),
                               yaxis=dict(title=y_name, titlefont=dict(size=40, color='black', family='Arial'), zeroline=False, autorange=True,
                                          showgrid=True, gridwidth=1, gridcolor='rgba(235,240,248,100)', tickfont=dict(size=30)),
                               legend=dict(x=1, y=1, orientation='v', font=dict(size=30)),
                               plot_bgcolor='rgba(255, 255, 255, 0.1)',
                               paper_bgcolor='rgba(255, 255, 255, 0.2)')
            fig = go.Figure(data=data, layout=layout)
            # fig.show()
            pio.write_image(fig, 'SD_' + output_file_name)
            try:
                if int(ave) == 3:
                    average_value = [(values1[i] + values2[i] + values3[i]) / 3 for i in range(len(values1))]
                    average_sd = [(sd1[i] + sd2[i] + sd3[i]) / 3 for i in range(len(sd1))]
                elif int(ave) == 2:
                    average_value = [(values1[i] + values2[i]) / 2 for i in range(len(values1))]
                    average_sd = [(sd1[i] + sd2[i]) / 2 for i in range(len(sd1))]
                trace_ave = go.Scatter(x=time1, y=average_value, name='Mean Value')
                error_y_ave = dict(type='data', array=average_sd, visible=True, thickness=1)
                trace_ave.update(error_y=error_y_ave)
                data_ave = [trace_ave]
                fig_ave = go.Figure(data=data_ave, layout=layout)
                pio.write_image(fig_ave, 'SD_mean_'+output_file_name)
            except:
                pass


    def plotly_multy(self, flag, multi_files, xaxis_name, yaxis_name, renumber, ave, output_file_name, rdf_cutoff, plot_name):
        a=1
        data = []

        ################## for regular expression substance ##################
        regex = r"\[|\]|'"

        ################## define plot title, x axis name and y axis name ##################
        if plot_name == '':
            with open(multi_files[0], "r") as f:
                plot_title = re.sub(regex, "", str(re.findall('"([^"]*)"', f.readlines()[13])))
        else:
            plot_title = str(plot_name)
        if xaxis_name == 0:
            with open(multi_files[0], "r") as f:
                x_name = re.sub(regex, "", str(re.findall('"([^"]*)"', f.readlines()[14])))
        else:
            x_name = xaxis_name
            pass

        if yaxis_name == 0 and plot_title not in ['Solvent Accessible Surface', 'Area per residue over the trajectory']:
            with open(multi_files[0], "r") as f:
                y_name = re.sub(regex, "", str(re.findall('"([^"]*)"', f.readlines()[15])))
        elif yaxis_name == 0 and plot_title in ['Solvent Accessible Surface', 'Area per residue over the trajectory']:
            y_name = 'Area (nm<sup>2</sup>)'
        else:
            y_name = yaxis_name
            pass
        ################## reading the datas!! ##################
        for i in multi_files:
            # create empty list
            locals()["x_" + str(a)] = []
            locals()["y_" + str(a)] = []
            locals()["sd_" + str(a)] = []
            # grab datas from input files
            with open(i, "r") as f:
                lines = f.readlines()
                for num in range(len(lines)):
                    if lines[num].startswith("#") or lines[num].startswith("@"):
                        pass
                    else:
                        if x_name == 'Time (ps)':       # set time value into ns
                            locals()["x_" + str(a)].append(float(lines[num].split()[0])/1000)
                        else:
                            locals()["x_" + str(a)].append(float(lines[num].split()[0]))
                        locals()["y_" + str(a)].append(float(lines[num].split()[1]))
                        try:
                            locals()["sd_" + str(a)].append(float(lines[num].split()[2]))
                        except:
                            pass
            if x_name == 'Residue' and renumber == 'true':
                for k in range(len(locals()["x_" + str(a)])):
                    locals()["x_" + str(a)][k] = k+1

            ################## define traces ##################
            locals()["trace" + str(a)] = go.Scatter(x=locals()["x_" + str(a)], y=locals()["y_" + str(a)], name=str(i).split('.')[0])
            data.append(locals()["trace" + str(a)])
            a += 1
        ################## test if time unit is ns, if not then change it from ps to ns ##################
        if x_name == 'Time (ps)':
            x_name = 'Time (ns)'

        ################## plot the datas ##################
        layout = go.Layout(title=plot_title, title_x=0.5, title_y=1, font=dict(size=24),
                           xaxis=dict(title=x_name, titlefont=dict(size=40, color='black', family='Arial'), zeroline=False, autorange=True,
                                      showgrid=True, gridwidth=1, gridcolor='rgba(235,240,248,100)', tickfont=dict(size=30)),
                           yaxis=dict(title=y_name, titlefont=dict(size=40, color='black', family='Arial'), zeroline=False, autorange=True,
                                      showgrid=True, gridwidth=1, gridcolor='rgba(235,240,248,100)', tickfont=dict(size=30)),
                           legend=dict(x=1, y=1, orientation='v', font=dict(size=30)), showlegend=True,
                           plot_bgcolor='rgba(255, 255, 255, 0.1)',
                           paper_bgcolor='rgba(255, 255, 255, 0.2)',
                           width=800, height=600)
        fig = go.Figure(data=data, layout=layout)
        pio.write_image(fig, output_file_name)
        ################## if user ask for average the inputs ##################
        if ave == 'true':
            number = len(multi_files)
            # average_value = [(locals()["x_" + str(a)][i] + values2[i] + values3[i]) / 3 for i in range(len(values1))]
            # average_sd = [(sd1[i] + sd2[i] + sd3[i]) / 3 for i in range(len(sd1))]
            average_value = locals()["y_" + str(1)]
            for a in range(1, number):
                average_value = [x + y for x, y in zip(average_value, locals()["y_" + str(a+1)])]
            average_value = [x/number for x in average_value]

            while len(locals()["x_" + str(a)]) != len(average_value):
                if len(locals()["x_" + str(a)]) > len(average_value):
                    locals()["x_" + str(a)].pop()
                else:
                    average_value.pop()
            data = []
            trace_ave = go.Scatter(x=locals()["x_" + str(a)], y=average_value, name='Average Values')
            data.append(trace_ave)
            layout = go.Layout(title=plot_title, title_x=0.5, title_y=1, font=dict(size=24),
                               xaxis=dict(title=x_name, titlefont=dict(size=40, color='black', family='Arial'), zeroline=False, autorange=True,
                                          showgrid=True, gridwidth=1, gridcolor='rgba(235,240,248,100)', tickfont=dict(size=30)),
                               yaxis=dict(title=y_name, titlefont=dict(size=40, color='black', family='Arial'), zeroline=False, autorange=True,
                                          showgrid=True, gridwidth=1, gridcolor='rgba(235,240,248,100)', tickfont=dict(size=30)),
                               legend=dict(x=1, y=1, orientation='v', font=dict(size=30)), showlegend=True,
                               plot_bgcolor='rgba(255, 255, 255, 0.1)',
                               paper_bgcolor='rgba(255, 255, 255, 0.2)',
                               width=800, height=600)
            fig = go.Figure(data=data, layout=layout)
            pio.write_image(fig, "Average_" + output_file_name)
            ######### output the average.xvg file ########
            with open(output_file_name[:-4]+"_average.xvg", 'w') as f:
                with open(multi_files[0], "r") as a:
                    lines = a.readlines()
                    for num in range(len(lines)):
                        if lines[num].startswith("#") or lines[num].startswith("@"):
                            f.write(lines[num])
                        else:
                            pass
                for num in range(len(average_value)):
                    average_line = "{}     {}\n".format(locals()["x_" + str(1)][num], average_value[num])
                    f.write(average_line)

    def plotly_pca(self, multi_files, xaxis_name, yaxis_name, renumber, ave, output_name, rdf_cutoff, plot_name, nbin, size):
        x_points = []
        y_points = []

        ################## define plot title, x axis name and y axis name ##################
        if plot_name == '':
            plot_title = 'PCA 2D projection of trajectory'
        else:
            plot_title = str(plot_name)
            pass

        if xaxis_name == 0:
            x_name = 'projection on eigenvector 1 (nm)'
        else:
            x_name = xaxis_name
            pass

        if yaxis_name == 0 and plot_title not in ['Solvent Accessible Surface', 'Area per residue over the trajectory']:
            y_name = 'projection on eigenvector 2 (nm)'
        else:
            y_name = yaxis_name
            pass
        ################## reading the datas!! ##################
        for i in multi_files:
            x_points=[]
            y_points=[]
            labels = []
            # create empty list
            with open(i, "r") as f:
                lines = f.readlines()
                for num in range(len(lines)):
                    if lines[num].startswith("#") or lines[num].startswith("@"):
                        pass
                    else:
                        x_points.append(float(lines[num].split()[0]))
                        y_points.append(float(lines[num].split()[1]))
                        labels.append(num)

            # # 创建散点图轨迹
            # scatter_trace = go.Scatter(
            #     x=x_points,
            #     y=y_points,
            #     mode='markers',
            #     marker=dict(
            #         size=5,
            #         color='black'
            #     )
            # )

            # layout = go.Layout(title=plot_title, title_x=0.5, title_y=1, font=dict(size=20),
            #                     xaxis=dict(title=x_name, titlefont=dict(size=20, color='black', family='Arial'), zeroline=False, autorange=True,
            #                               showgrid=False, gridwidth=1, gridcolor='rgba(0,0,0,0.1)', tickfont=dict(size=20)),
            #                     yaxis=dict(title=y_name, titlefont=dict(size=20, color='black', family='Arial'), zeroline=False, autorange=True,
            #                               showgrid=False, gridwidth=1, gridcolor='rgba(0,0,0,0.1)', tickfont=dict(size=20)),
            #                     legend=dict(x=1, y=1, orientation='v', font=dict(size=30)), showlegend=False,
            #                     plot_bgcolor='rgba(255, 255, 255, 0.1)',
            #                     paper_bgcolor='rgba(255, 255, 255, 0.2)')

            # # 创建图形对象
            # fig = go.Figure(data=scatter_trace, layout=layout)

            # # 显示图形
            # # fig.show()
            # if output_name == 'plotly.png':
            #     pio.write_image(fig, "PCA_Scatter_"+i.split('.')[0]+".png")
            # else:
            #     pio.write_image(fig, "Scatter_" + output_name)

            # 创建一个 Scatter 对象
            scatter = go.Scatter(
                x=x_points,
                y=y_points,
                mode='markers',
                marker=dict(
                    color=labels,  # 设置颜色为标签的数值
                    colorscale='rainbow',  # 颜色映射，你可以根据需要选择不同的颜色映射
                    colorbar=dict(title='Label Range'),  # 添加颜色条
                ),
            )

            # 创建数据列表
            data = [scatter]

            # 创建布局
            layout = go.Layout(
                title='PCA plot with Color Bar for frame order', title_x=0.5, title_y=1, font=dict(size=24),
                xaxis=dict(title='PC1 (nm)', titlefont=dict(size=40, color='black', family='Arial'), zeroline=False, autorange=True,
                           showgrid=True, gridwidth=1, gridcolor='rgba(235,240,248,100)', tickfont=dict(size=30)),
                yaxis=dict(title='PC2 (nm)', titlefont=dict(size=40, color='black', family='Arial'), zeroline=False, autorange=True,
                           showgrid=True, gridwidth=1, gridcolor='rgba(235,240,248,100)', tickfont=dict(size=30)),
                plot_bgcolor='rgba(255, 255, 255, 0.1)',
                paper_bgcolor='rgba(255, 255, 255, 0.2)',
                width=800, height=600
            )

            # 创建 Figure 对象
            fig = go.Figure(data=data, layout=layout)
            if output_name == 'plotly.png':
                pio.write_image(fig, "PCA_Scatter_"+i.split('.')[0]+".png")
            else:
                pio.write_image(fig, "Scatter_" + output_name)

####################################################################################################################################
        # 用R创建PCA图
        x_diff = np.max(x_points) - np.min(x_points)
        y_diff = np.max(y_points) - np.min(y_points)
        width= (x_diff / y_diff) * size
        height = 1 * size
        # 读取数据
        for i in multi_files:
            r('data <- read.table(%r, skip = 17, header = FALSE)' % (i))

            # 设置变量名
            r('names(data) <- c("PC1", "PC2")')

            # 导入所需的R包
            r('library(RColorBrewer)')

            # 执行其他指令
            r('zBot <- 1.52')
            r('zTop <- 3.42')
            r('zW <- 0.83')
            r('buylrd <- rev(brewer.pal(11,"RdYlBu"))')

            # 保存为PNG文件
            if output_name == 'plotly.png':
                # r('png(file=%r, height=600, width=450)' % ("PCA_Density_" + i.split('.')[0] + ".png"))
                # r('png(file=%r)' % ("PCA_Density_" + i.split('.')[0] + ".png"))
                r('png(file=%r, height=%r, width=%r)' % ("PCA_Density_" + i.split('.')[0] + ".png", height, width))
                r('smoothScatter(data$PC2 ~ data$PC1, nbin=%r, colramp = colorRampPalette(c(buylrd)),nrpoints=Inf, pch="", cex=.7,col="black",main=%r, xlab=%r, ylab=%r,transformation = function(x) x^.45)' % (nbin, plot_title, x_name, y_name))
                r('dev.off()')
            else:
                # r('png(file=%r, height=600, width=450)' % (output_name))
                # r('png(file=%r)' % (output_name))
                r('png(file=%r, height=%r, width=%r)' % ("Density_" + output_name, height, width))
                r('smoothScatter(data$PC2 ~ data$PC1, nbin=%r, colramp = colorRampPalette(c(buylrd)),nrpoints=Inf, pch="", cex=.7,col="black",main=%r, xlab=%r, ylab=%r,transformation = function(x) x^.45)' % (nbin, plot_title, x_name, y_name))
                r('dev.off()')

            # # 保存为PDF文件
            # r('pdf(file="PCA.pdf", height=1600, width=1600,paper = "a4")')
            # r('smoothScatter(data$PC2 ~ data$PC1, nbin=1000, colramp = colorRampPalette(c(buylrd)),nrpoints=Inf, pch="", cex=.7,col="black",main="Shangze MD simulation PCA analysis", xlab="PCA1", ylab="PCA2",transformation = function(x) x^.5)')
            # r('dev.off()')








x = plotly_go(file1,file2,file3, output_name, renumber, ave, xaxis_name, yaxis_name, rdf_cutoff, multi_files, plot_name, pca, nbin, size)
