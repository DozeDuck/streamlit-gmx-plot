#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 14 16:15:48 2024

@author: dozeduck
"""
# import getopt
# import sys
import re
import os
import pandas as pd
import mimetypes
# import plotly
import plotly.graph_objs as go
import plotly.io as pio
# for PCA
import numpy as np
# from rpy2.robjects import r
# for histogram dist plot
import plotly.figure_factory as ff
# import argparse
# for bool values
import ast
# metal restraints adding
import math
# for Streamlit usage, wide screen display
import streamlit as st
st.set_page_config(layout="wide")
from tempfile import NamedTemporaryFile
import base64
# for contact map
import MDAnalysis as mda
from MDAnalysis.analysis import contacts
import csv
import matplotlib.pyplot as plt

#################################################################################################################################################

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

    def __init__(self, multi_files, output_name, renumber, rdf_cutoff, average, plot_name, nbin, size, move_average, mean_value, histogram, xaxis_name, yaxis_name, xaxis_size, yaxis_size, xy_font, title_font, legend_show, legend_font, font_family, font_color, grid_show, uploaded_filenames, l,r,t,b, violin):

        if len(multi_files) >=1:
            # print(multi_files)
            file1 = multi_files[0]
            self.flag_recognizer(file1, plot_name)
            if self.pca_flag != 1 and self.flag != 'pca':
                self.plotly_multy(multi_files, output_name, renumber, rdf_cutoff, average, plot_name, nbin, size, move_average, mean_value, histogram, xaxis_name, yaxis_name, xaxis_size, yaxis_size, xy_font, title_font, legend_show, legend_font, font_family, font_color, grid_show, self.flag, uploaded_filenames, l,r,t,b, violin)
            elif self.pca_flag == 1:
                self.plotly_pca(multi_files, output_name, renumber, rdf_cutoff, average, plot_name, nbin, size, move_average, mean_value, histogram, xaxis_name, yaxis_name, xaxis_size, yaxis_size, xy_font, title_font, legend_show, legend_font, font_family, font_color, grid_show, self.flag, uploaded_filenames, l,r,t,b)
            elif self.flag == 'pca':
                self.plotly_pca(multi_files, output_name, renumber, rdf_cutoff, average, plot_name, nbin, size, move_average, mean_value, histogram, xaxis_name, yaxis_name, xaxis_size, yaxis_size, xy_font, title_font, legend_show, legend_font, font_family, font_color, grid_show, self.flag, uploaded_filenames, l,r,t,b)


    def flag_recognizer(self,file1, plot_name):                                                   # first method to be called in __main__, used for creating object and charactors.
        flags_map = {
            'rms,': 'rmsd',
            'rmsf,': 'rmsf',
            'sasa,': 'sasa',
            'gyrate,': 'gyrate',
            'dipoles,': 'dipoles',
            'distance,': 'distance',
            'rdf,': 'rdf',
            'convergence': 'convergence',
            'anaeig,': 'pca',
            'angle,': 'angle'
        }
                 
        if file1.endswith(".xvg"):
            with open(file1, 'r') as f:
                lines = f.readlines()
                if len(lines) >= 3:
                    try:
                        flag = lines[2].split()[5]
                        self.flag = flags_map.get(flag, flag)
                    except IndexError:
                        pass
                if len(lines) >= 9 and '-or' in lines[8]:
                    self.sasa_flag = '-or'
    
                if 'pca' in str(file1).lower() or '2dproj' in str(file1):
                    self.pca_flag = 1
                print("I know you are plotting " + self.flag + " figures!")
       
        elif file1.endswith(".csv"):
            base_name = os.path.basename(file1)
            filename = os.path.splitext(base_name)[0]
            self.flag = filename.lower()
            if plot_name.lower == 'pca':
                self.flag = 'pca'
            print("I know you are plotting " + self.flag + " figures!")
        st.text(self.flag)
        


    def consist(self,x_values):
        seq1 = [x_values[0]]
        seq2 = []
        # seq3 = []
        for i in range(1, len(x_values)):
            # find the break point
            if x_values[i] <= x_values[i-1]+2 and seq2 == []:
                seq1.append(x_values[i])
            else:
                seq2.append(x_values[i])
        return seq1, seq2


    def read_data(self, file, x_name, renumber):
        # 从文件中读取数据   
        x_data, y_data, sd_data = [], [], []
        if file.endswith(".xvg"):
            with open(file, "r") as f:
                lines = f.readlines()
                for line in lines:
                    if line.startswith("#") or line.startswith("@"):
                        continue
                    else:
                        # 解析数据行
                        split_line = line.split()
                        x_value = float(split_line[0])
                        y_value = float(split_line[1])
    
                        if x_name == 'Time (ps)':  # 将时间从ps转换为ns
                            x_value /= 1000
    
                        if x_name == 'Residue' and renumber == 'true':
                            x_value = len(x_data) + 1
    
                        x_data.append(x_value)
                        y_data.append(y_value)
    
                        # 读取标准差（如果存在）
                        try:
                            sd_data.append(float(split_line[2]))
                        except IndexError:
                            pass
      
        elif file.endswith(".csv"):
            df = pd.read_csv(file, skiprows=0)
            x_data = df[df.columns[0]].tolist()  # 第一列作为X轴数据
            if x_name == 'Time (ps)':  # 将时间从ps转换为ns
                scaled_list = [x / 1000 for x in x_data]
                x_data = scaled_list
            elif x_name in ['Residue', 'residue'] and renumber == 'true':
                x_data = list(range(1,len(x_data)+1))
            # 将除第一列外的所有列作为Y轴数据，每列一个Y数据序列
            y_data = [df[col].tolist() for col in df.columns[1:]]
            # CSV文件不包含标准差数据，因此sd_data保持为空 
                   
        return x_data, y_data, sd_data  
    
    def extract_plot_details(self, multi_files, plot_name, xaxis_name, yaxis_name, flag, histogram):
        traces_name_list = []
        ## Read XVG files
        if multi_files[0].endswith(".xvg"):
            regex = r"\[|\]|'"
            
            # 提取或设置图表标题
            if plot_name == 'auto detect':
                with open(multi_files[0], "r") as f:
                    plot_title = re.sub(regex, "", str(re.findall('"([^"]*)"', f.readlines()[13])))
            else:
                plot_title = str(plot_name)
    
            # 提取或设置X轴名称
            if xaxis_name == 'auto detect':
                with open(multi_files[0], "r") as f:
                    x_name = re.sub(regex, "", str(re.findall('"([^"]*)"', f.readlines()[14])))
            else:
                x_name = xaxis_name
    
            # 提取或设置Y轴名称
            if yaxis_name == 'auto detect':
                with open(multi_files[0], "r") as f:
                    y_name = re.sub(regex, "", str(re.findall('"([^"]*)"', f.readlines()[15])))
                if plot_title in ['Solvent Accessible Surface', 'Area per residue over the trajectory']:
                    y_name = 'Area (nm<sup>2</sup>)'
                elif flag == 'dihedral_distribution' and histogram == 'true':
                    y_name = 'Probability'
            else:
                y_name = yaxis_name
    

        ## Read CSV files
        elif multi_files[0].endswith(".csv"):
            df = pd.read_csv(multi_files[0])
            # 提取或设置图表标题
            if plot_name == 'auto detect':
                base_name = os.path.basename(multi_files[0])
                filename = os.path.splitext(base_name)[0]
                plot_title = str(filename)
            else:
                plot_title = str(plot_name)
    
            # 提取或设置X轴名称
            if xaxis_name == 'auto detect':
                x_name = df.columns[0]
            else:
                x_name = xaxis_name
    
            # 提取或设置Y轴名称
            if yaxis_name == 'auto detect':
                y_name = ''
                traces_name_list.extend(df.columns[1:])
            else:
                y_name = yaxis_name
                traces_name_list.extend(df.columns[1:])
    
        return plot_title, x_name, y_name, traces_name_list

    def define_trace(self, x_data, y_data, file_name, colour, violine='False', flag=0, labels=0):
        # 创建并返回迹线
        if flag == 'pca':
            trace = go.Scatter(
                x=x_data,
                y=y_data,
                mode='markers',
                marker=dict(
                    color=labels,  # 设置颜色为标签的数值
                    colorscale=colour,  # 颜色映射，你可以根据需要选择不同的颜色映射
                    colorbar=dict(title='Label Range'),  # 添加颜色条
                ),
            )
        elif violine != 'False':
            trace = go.Violin(x0=str(file_name).split('.')[0], y=y_data, line=dict(color='black'), fillcolor=colour, name=str(file_name).split('.')[0], box_visible=True, meanline_visible=True, opacity=0.6)            
        else:
            trace = go.Scatter(x=x_data, y=y_data, line=dict(color=colour), name=str(file_name).split('.')[0])
        return trace

    def setup_layout(self, plot_title, title_font, x_name, y_name, xy_font, xaxis_size, yaxis_size, font_color, legend_show, legend_font, font_family, grid_show, l,r,t,b, violine='False', flag=0):
        # 设置布局
        if flag == 'pca':
            legend_show = False
        if violine != 'False':
            x_name = ''

        layout = go.Layout(
            title=plot_title, title_x=0.5, title_y=0.99, font=dict(size=title_font, color=font_color),
            xaxis=dict(title=x_name, titlefont=dict(size=xy_font, color=font_color, family=font_family), zeroline=False, autorange=True,
                       showgrid=grid_show, gridwidth=1, gridcolor='rgba(235,240,248,100)', tickfont=dict(size=30)),
            yaxis=dict(title=y_name, titlefont=dict(size=xy_font, color=font_color, family=font_family), zeroline=False, autorange=True,
                       showgrid=grid_show, gridwidth=1, gridcolor='rgba(235,240,248,100)', tickfont=dict(size=30)),
            legend=dict(x=1, y=1, orientation='v', font=dict(size=legend_font, color=font_color)), showlegend=legend_show,
            plot_bgcolor='rgba(255, 255, 255, 0.1)',
            paper_bgcolor='rgba(255, 255, 255, 0.2)',
            margin=dict(l=int(l), r=int(r), t=int(t), b=int(b)),
            width=xaxis_size, height=yaxis_size
        )
        return layout

    def streamlit_download_file_plotly(self, download_name, content_file):
        # 读取文件内容
        with open(content_file, "rb") as file:
            file_content = file.read()
        
        # 获取文件的 MIME 类型
        mime_type, _ = mimetypes.guess_type(content_file)
        
        # 创建下载按钮
        st.download_button(
            label=f"Download {download_name}",
            data=file_content,
            file_name=download_name,
            mime=mime_type)


    def plot_graph(self, data, layout, output_file_name):
        # 使用数据和布局绘制图形
        fig = go.Figure(data=data, layout=layout)
        pio.write_image(fig, "/tmp/" + output_file_name)
        self.streamlit_download_file_plotly(output_file_name, "/tmp/" + output_file_name)

    def plot_histogram(self, histogram_data, group_labels, plot_title, output_file_name, colors, nbin):
        # 处理直方图
        fig_hist = ff.create_distplot(histogram_data, group_labels, colors=colors, bin_size=nbin, curve_type='normal')
        fig_hist.update_layout(title_text=plot_title)
        pio.write_image(fig_hist, "/tmp/" + output_file_name)
        self.streamlit_download_file_plotly("hist_" + output_file_name, "/tmp/" + output_file_name)

    def calculate_average(self, multi_files, xaxis_name, renumber):
        # 计算平均值
        sum_data = None
        for file in multi_files:
            x_data, y_data, _ = self.read_data(file, xaxis_name, renumber)
            if sum_data is None:
                sum_data = np.array(y_data)
            else:
                sum_data += np.array(y_data)
        return sum_data / len(multi_files)
    
    def output_average_file(self, output_file_name, average_value, multi_files, xaxis_name, renumber, x_data):
        with open(output_file_name[:-4]+"_average.xvg", 'w') as f:
            with open(multi_files[0], "r") as a:
                lines = a.readlines()
                for num in range(len(lines)):
                    if lines[num].startswith("#") or lines[num].startswith("@"):
                        f.write(lines[num])
                    else:
                        pass
            for num in range(len(average_value)):
                average_line = "{}     {}\n".format(x_data[num], average_value[num])
                f.write(average_line)

    def moving_average(self, y_data, window_size):
        # 计算移动平均
        return np.convolve(y_data, np.ones(window_size) / window_size, mode='valid')


    def plotly_multy(self, multi_files, output_name, renumber, rdf_cutoff, average, plot_name, nbin, size, move_average, mean_value, histogram, xaxis_name, yaxis_name, xaxis_size, yaxis_size, xy_font, title_font, legend_show, legend_font, font_family, font_color, grid_show, flag, uploaded_filenames, l,r,t,b, violin):
        Plotly = ['#636EFA', '#EF553B', '#00CC96', '#AB63FA', '#FFA15A', '#19D3F3', '#FF6692', '#B6E880', '#FF97FF', '#FECB52']
        data, histogram_data, group_labels = [], [], []

        # 读取plot_title, x_name, y_name
        plot_title, x_name, y_name, traces_name_list = self.extract_plot_details(multi_files, plot_name, xaxis_name, yaxis_name, flag, histogram)
        # 读取数据并创建迹线
        if multi_files[0].endswith(".xvg"):
            for i, file in enumerate(multi_files):
                x_data, y_data, _ = self.read_data(file, x_name, renumber)
                trace = self.define_trace(x_data, y_data, uploaded_filenames[i], Plotly[i % len(Plotly)], violine=violin)
                data.append(trace)
    
                # 添加直方图数据
                if histogram == 'true':
                    histogram_data.append(y_data)
                    group_labels.append(str(uploaded_filenames[i]).split('.')[0])
        elif multi_files[0].endswith(".csv"):
            for i, trace in enumerate(traces_name_list):
                x_data, y_data, _ = self.read_data(multi_files[0], x_name, renumber)
                trace = self.define_trace(x_data, y_data[i], trace, Plotly[i % len(Plotly)], violine=violin)
                data.append(trace)
    
                # 添加直方图数据
                if histogram == 'true':
                    histogram_data.append(y_data)
                    group_labels.append(str(file).split('.')[0])
    
        # change Time (ps) to Time (ns)
        if x_name == 'Time (ps)':
            x_name = 'Time (ns)'
        
        # 设置布局
        layout = self.setup_layout(plot_title, title_font, x_name, y_name, xy_font, xaxis_size, yaxis_size, font_color, legend_show, legend_font, font_family, grid_show, l, r, t , b, violine=violin)

        # 绘制图形
        self.plot_graph(data, layout, output_name)

        # 处理直方图
        if histogram == 'true':
            self.plot_histogram(histogram_data, group_labels, plot_title, output_name, Plotly, nbin)

        # 处理平均值
        if average == 'true':
            average_data = self.calculate_average(multi_files, xaxis_name, renumber)
            average_trace = self.define_trace(x_data, average_data, "Average", 'black')
            # data.append(average_trace)
            data = average_trace
            self.plot_graph(data, layout, "Average_" + output_name)
            self.output_average_file("/tmp/" + output_name, average_data, multi_files, xaxis_name, renumber, x_data)

        # 处理移动平均
        if move_average != 0:
            ma_data = []
            for file in multi_files:
                _, y_data, _ = self.read_data(file, xaxis_name, renumber)
                ma_y_data = self.moving_average(y_data, move_average)
                ma_trace = self.define_trace(x_data[move_average - 1:], ma_y_data, str(file).split('.')[0], Plotly[0])
                ma_data.append(ma_trace)
            self.plot_graph(ma_data, layout, "MovingAverage_" + output_name)


        # # 调用上述方法
        # for i, file in enumerate(multi_files):
        #     x_data, y_data, sd_data = self.read_data(file, xaxis_name, renumber)
        #     trace = self.define_trace(x_data, y_data, file, Plotly[i % len(Plotly)])
        #     data.append(trace)

        # layout = self.setup_layout(plot_title, title_font, x_name, y_name, xy_font, xaxis_size, yaxis_size, font_color, legend_show, legend_font, font_family, grid_show)
        # self.plot_graph(data, layout, output_name)


         
            
    def plotly_pca(self, multi_files, output_name, renumber, rdf_cutoff, average, plot_name, nbin, size, move_average, mean_value, histogram, xaxis_name, yaxis_name, xaxis_size, yaxis_size, xy_font, title_font, legend_show, legend_font, font_family, font_color, grid_show, flag, uploaded_filenames, l, r, t ,b):
        data = []
        color = ['rainbow']
        # labels = []
        # 使用 extract_plot_details 方法获取图表标题和轴标签
        plot_title, x_name, y_name, traces_name_list = self.extract_plot_details(multi_files, plot_name, xaxis_name, yaxis_name, flag, histogram)

        # 处理 PCA 数据
        for i, file in enumerate(multi_files):          
            x_data, y_data, _ = self.read_data(file, "PC1", renumber)  # 假设 "PC1" 和 "PC2" 是合适的轴名称
            labels = [x for x in range(len(y_data))]
            
            # 使用 define_trace 创建迹线
            trace = self.define_trace(x_data, y_data, file, 'rainbow', flag=flag, labels=labels)  # 假设使用 'rainbow' 作为颜色
            data.append(trace)

        # 使用 setup_layout 设置布局
        layout = self.setup_layout(plot_title, title_font, 'PC1 (nm)', 'PC2 (nm)', xy_font, xaxis_size, yaxis_size, font_color, legend_show, legend_font, font_family, grid_show, l, r, t ,b, flag=flag)

        # 使用 plot_graph 绘制图形
        self.plot_graph(data, layout, "Scatter_" + output_name)

##########################################################################################################################################################################################
class mr(): # read content from the uploaded file directly.
    head = ''
    total_atom = 0
    resid = []
    resname= []
    atomname = []
    index = []
    x = []
    y = []
    z = []
    xyz = []
    last = ''
    metals = []
    coordinators = []
    metal1 = 0
    metal2 = 0
    metal3 = 0
    metal4 = 0
    metal5 = 0
    metal6 = 0
    metal7 = 0
    atom1 = []
    atom2 = []
    atom3 = []
    atom4 = []
    atom5 = []
    atom6 = []
    atom7 = []    
    
    def __init__(self, gro,num_neighbours, distance_value, atom_list, metal_list, residue_list, bond_strength, angle_strength):
        self.output = ""
        self.GROreader(gro)
        self.MetalMiner(metal_list)
        self.coordinator(num_neighbours, distance_value, atom_list, metal_list, residue_list)
        # self.bond_cal(atom6,bond_strength)
        # self.pair_cal()
        # self.angle_cal(angle_strength)
        self.bond_cal( self.metal1, self.atom1, self.metal2, self.atom2, self.metal3, self.atom3, self.metal4, self.atom4, self.metal5, self.atom5, self.metal6, self.atom6, bond_strength)
        self.pair_cal( self.metal1, self.atom1, self.metal2, self.atom2, self.metal3, self.atom3, self.metal4, self.atom4, self.metal5, self.atom5, self.metal6, self.atom6)
        self.angle_cal( self.metal1, self.atom1, self.metal2, self.atom2, self.metal3, self.atom3, self.metal4, self.atom4, self.metal5, self.atom5, self.metal6, self.atom6, angle_strength)

    def GROreader(self,gro): 
        lines = gro.splitlines()  # for streamlit 如果 'gro' 是一个二进制文件，使用 gro.decode().splitlines() 

            
        # extra lines
        self.head = lines[0].strip()
        self.total_atom = int(lines[1])
        self.last = lines[-1]
        
        # 忽略前两行和最后一行
        lines = lines[2:-1]
        
        # 逐行解析内容
        for line in lines:
            line = line.strip()  # 去除首尾空格和换行符
            match = re.match(r'(\d+)([A-Za-z]{2,})', line)
            if match:
                self.resid.append(int(match.group(1)))
                self.resname.append(str(match.group(2)))
            self.atomname.append(str(line.split()[1]))                        # The 3rd column is the atom name C CA CD1 CD2 and so on
            self.index.append(int(line.split()[2]))                   # Column 4 is the residue name TYR ALA etc.
            self.x.append(float(line.split()[3]))                         # The 5th column is the name of the chain it is on
            self.y.append(float(line.split()[4]))               # The sixth column is the residue number
            self.z.append(float(line.split()[5]))                   # Column 7 is the x-coordinate of the atom
            self.xyz.append([float(line.split()[3]),float(line.split()[4]),float(line.split()[5])])
    
    def MetalMiner(self, metal_list):
        print(metal_list)
        for i in range(len(self.atomname)):
            if self.atomname[i] in metal_list and self.resname[i] in metal_list:
                self.metals.append(self.index[i])
        # the index in list should -1   
        for i in range(len(self.metals)):
            try:
                setattr(self, f'metal{i+1}', self.metals[i] - 1)
            except IndexError:
                pass

        metals_name = [self.atomname[i-1] for i in self.metals]
        sentence = "The metals atom index are: {}".format(list(zip(metals_name, self.metals)))
        print(sentence)
        # print(self.metal1,self.metal2,self.metal3)
            
        
        
        
    def coordinator(self, num_neighbours, distance_value, atom_list, metal_list, residue_list):
        # find the atom index
        if self.metal1 != 0:
            atom_distances = {}
            for i in range(len(self.index)):
                if self.resname[i] in residue_list and self.atomname[i] in atom_list:
                    dist = self.distance(self.metal1, i)
                    if dist <= distance_value:
                        atom_distances[i] = dist
            
            # 对满足条件的原子按照距离 metal 的距离进行排序
            sorted_atoms = sorted(atom_distances.items(), key=lambda x: x[1])
            
            # 仅保留前 num_neighbours 个邻居
            self.atom1 = [atom_index for atom_index, _ in sorted_atoms[:num_neighbours]]
                
                
        if self.metal2 != 0:
            atom_distances = {}
            for i in range(len(self.index)):
                if self.resname[i] in residue_list and self.atomname[i] in atom_list:
                    dist = self.distance(self.metal2, i)
                    if dist <= distance_value:
                        atom_distances[i] = dist
            
            # 对满足条件的原子按照距离 metal 的距离进行排序
            sorted_atoms = sorted(atom_distances.items(), key=lambda x: x[1])
            
            # 仅保留前 num_neighbours 个邻居
            self.atom2 = [atom_index for atom_index, _ in sorted_atoms[:num_neighbours]]
                
        if self.metal3 != 0:
            atom_distances = {}
            for i in range(len(self.index)):
                if self.resname[i] in residue_list and self.atomname[i] in atom_list:
                    dist = self.distance(self.metal3, i)
                    if dist <= distance_value:
                        atom_distances[i] = dist
            
            # 对满足条件的原子按照距离 metal 的距离进行排序
            sorted_atoms = sorted(atom_distances.items(), key=lambda x: x[1])
            
            # 仅保留前 num_neighbours 个邻居
            self.atom3 = [atom_index for atom_index, _ in sorted_atoms[:num_neighbours]]

        if self.metal4 != 0:
            atom_distances = {}
            for i in range(len(self.index)):
                if self.resname[i] in residue_list and self.atomname[i] in atom_list:
                    dist = self.distance(self.metal4, i)
                    if dist <= distance_value:
                        atom_distances[i] = dist
            
            # 对满足条件的原子按照距离 metal 的距离进行排序
            sorted_atoms = sorted(atom_distances.items(), key=lambda x: x[1])
            
            # 仅保留前 num_neighbours 个邻居
            self.atom4 = [atom_index for atom_index, _ in sorted_atoms[:num_neighbours]]

        if self.metal5 != 0:
            atom_distances = {}
            for i in range(len(self.index)):
                if self.resname[i] in residue_list and self.atomname[i] in atom_list:
                    dist = self.distance(self.metal5, i)
                    if dist <= distance_value:
                        atom_distances[i] = dist
            
            # 对满足条件的原子按照距离 metal 的距离进行排序
            sorted_atoms = sorted(atom_distances.items(), key=lambda x: x[1])
            
            # 仅保留前 num_neighbours 个邻居
            self.atom5 = [atom_index for atom_index, _ in sorted_atoms[:num_neighbours]]

        if self.metal6 != 0:
            atom_distances = {}
            for i in range(len(self.index)):
                if self.resname[i] in residue_list and self.atomname[i] in atom_list:
                    dist = self.distance(self.metal6, i)
                    if dist <= distance_value:
                        atom_distances[i] = dist
            
            # 对满足条件的原子按照距离 metal 的距离进行排序
            sorted_atoms = sorted(atom_distances.items(), key=lambda x: x[1])
            
            # 仅保留前 num_neighbours 个邻居
            self.atom6 = [atom_index for atom_index, _ in sorted_atoms[:num_neighbours]]

        if self.metal7 != 0:
            atom_distances = {}
            for i in range(len(self.index)):
                if self.resname[i] in residue_list and self.atomname[i] in atom_list:
                    dist = self.distance(self.metal7, i)
                    if dist <= distance_value:
                        atom_distances[i] = dist
            
            # 对满足条件的原子按照距离 metal 的距离进行排序
            sorted_atoms = sorted(atom_distances.items(), key=lambda x: x[1])
            
            # 仅保留前 num_neighbours 个邻居
            self.atom7 = [atom_index for atom_index, _ in sorted_atoms[:num_neighbours]]

        
    def distance(self, index1, index2):
        distance = math.sqrt((self.x[index2] - self.x[index1])**2 + (self.y[index2] - self.y[index1])**2 + (self.z[index2] - self.z[index1])**2)
        return distance
    
    def calculate_distance(self, point1, point2):
        x1, y1, z1 = point1
        x2, y2, z2 = point2
        distance = math.sqrt((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2)
        return distance
    
    def calculate_angle(self, point1, point2, point3):
        vector_ab = np.array(point2) - np.array(point1)
        vector_bc = np.array(point2) - np.array(point3)
    
        dot_product = np.dot(vector_ab, vector_bc)
        norm_ab = np.linalg.norm(vector_ab)
        norm_bc = np.linalg.norm(vector_bc)
    
        cos_angle = dot_product / (norm_ab * norm_bc)
        angle_rad = np.arccos(cos_angle)
        angle_deg = np.degrees(angle_rad)
    
        return angle_deg

    def bond_cal(self,  metal1, atom1, metal2, atom2, metal3, atom3, metal4, atom4, metal5, atom5, metal6, atom6, bond_strength):
        for i in [1,2,3,4,5,6]:
            if locals()["metal" + str(i)] !=0:
                metal = locals()["metal" + str(i)]+1
                atom = [i + 1 for i in locals()["atom" + str(i)]]
                metal_point = [self.x[metal-1],self.y[metal-1],self.z[metal-1]]
                print("; please add below to topol.top's distance part")
                self.output += "; please add below to topol.top's distance part\n" # for streamlit
                for i in atom:
                    locals()["atom" + str(i)] = [self.x[i-1],self.y[i-1],self.z[i-1]]
                    print("%5d%6d%6d%7.2f%9d" % (metal, i, 6, self.calculate_distance(metal_point, locals()["atom" + str(i)]), bond_strength))
                    self.output += "%5d%6d%6d%7.2f%9d\n" % (metal, i, 6, self.calculate_distance(metal_point, locals()["atom" + str(i)]), bond_strength) # for streamlit
        
    def pair_cal(self,  metal1, atom1, metal2, atom2, metal3, atom3, metal4, atom4, metal5, atom5, metal6, atom6):
        for i in [1,2,3,4,5,6]:
            if locals()["metal" + str(i)] !=0:
                metal = locals()["metal" + str(i)]+1
                atom = [i + 1 for i in locals()["atom" + str(i)]]
                target_resid = [self.resid[x-1] for x in atom]
        
                print("; I added the pairs - so the zn will not nonbonded interact with the cyx residues")
                self.output += "; add the pairs \n" # for streamlit
                for i in range(len(self.atomname)):
                    if self.atomname[i] == 'CA' and self.resid[i] in target_resid:
                        print("%5d%6d%6d" % (metal,i+1,1))
                        self.output += "%5d%6d%6d\n" % (metal,i+1,1) # for streamlit
           
    
    def angle_cal(self,  metal1, atom1, metal2, atom2, metal3, atom3, metal4, atom4, metal5, atom5, metal6, atom6, angle_strength):
        for i in [1,2,3,4,5,6]:
            if locals()["metal" + str(i)] !=0:
                metal = locals()["metal" + str(i)]+1
                atom = [i + 1 for i in locals()["atom" + str(i)]]
                # define how many neighbour atoms
                neighbour = len(atom)
                if neighbour >= 2:
                    metal_point = [self.x[metal-1],self.y[metal-1],self.z[metal-1]]
                    for i in atom:
                        locals()["atom" + str(i)] = [self.x[i-1],self.y[i-1],self.z[i-1]]
                    print("[ angle_restraints ]")
                    self.output += "[ angle_restraints ]\n" # for streamlit
                    if neighbour == 2:
                        print("%5d%6d%6d%6d%5d%9.2f%9d%9d" % (metal, atom[0], metal, atom[1],1,self.calculate_angle(locals()["atom" + str(atom[0])],metal_point, locals()["atom" + str(atom[1])]),angle_strength, 1))
                        self.output += "%5d%6d%6d%6d%5d%9.2f%9d%9d\n" % (metal, atom[0], metal, atom[1],1,self.calculate_angle(locals()["atom" + str(atom[0])],metal_point, locals()["atom" + str(atom[1])]),angle_strength, 1) # for streamlit
                    elif neighbour == 3:
                        print("%5d%6d%6d%6d%5d%9.2f%9d%9d" % (metal, atom[0], metal, atom[1],1,self.calculate_angle(locals()["atom" + str(atom[0])],metal_point, locals()["atom" + str(atom[1])]),angle_strength, 1))
                        print("%5d%6d%6d%6d%5d%9.2f%9d%9d" % (metal, atom[0], metal, atom[2],1,self.calculate_angle(locals()["atom" + str(atom[0])],metal_point, locals()["atom" + str(atom[2])]),angle_strength, 1))
                        print("%5d%6d%6d%6d%5d%9.2f%9d%9d" % (metal, atom[1], metal, atom[2],1,self.calculate_angle(locals()["atom" + str(atom[1])],metal_point, locals()["atom" + str(atom[2])]),angle_strength, 1))
                        self.output += "%5d%6d%6d%6d%5d%9.2f%9d%9d\n" % (metal, atom[0], metal, atom[1],1,self.calculate_angle(locals()["atom" + str(atom[0])],metal_point, locals()["atom" + str(atom[1])]),angle_strength, 1) # for streamlit
                        self.output += "%5d%6d%6d%6d%5d%9.2f%9d%9d\n" % (metal, atom[0], metal, atom[2],1,self.calculate_angle(locals()["atom" + str(atom[0])],metal_point, locals()["atom" + str(atom[2])]),angle_strength, 1) # for streamlit
                        self.output += "%5d%6d%6d%6d%5d%9.2f%9d%9d\n" % (metal, atom[1], metal, atom[2],1,self.calculate_angle(locals()["atom" + str(atom[1])],metal_point, locals()["atom" + str(atom[2])]),angle_strength, 1) # for streamlit
                    elif neighbour == 4:
                        print("%5d%6d%6d%6d%5d%9.2f%9d%9d" % (metal, atom[0], metal, atom[1],1,self.calculate_angle(locals()["atom" + str(atom[0])],metal_point, locals()["atom" + str(atom[1])]),angle_strength, 1))
                        print("%5d%6d%6d%6d%5d%9.2f%9d%9d" % (metal, atom[1], metal, atom[2],1,self.calculate_angle(locals()["atom" + str(atom[1])],metal_point, locals()["atom" + str(atom[2])]),angle_strength, 1))
                        print("%5d%6d%6d%6d%5d%9.2f%9d%9d" % (metal, atom[2], metal, atom[3],1,self.calculate_angle(locals()["atom" + str(atom[2])],metal_point, locals()["atom" + str(atom[3])]),angle_strength, 1))
                        print("%5d%6d%6d%6d%5d%9.2f%9d%9d" % (metal, atom[0], metal, atom[3],1,self.calculate_angle(locals()["atom" + str(atom[0])],metal_point, locals()["atom" + str(atom[3])]),angle_strength, 1))
                        self.output += "%5d%6d%6d%6d%5d%9.2f%9d%9d\n" % (metal, atom[0], metal, atom[1],1,self.calculate_angle(locals()["atom" + str(atom[0])],metal_point, locals()["atom" + str(atom[1])]),angle_strength, 1) # for streamlit
                        self.output += "%5d%6d%6d%6d%5d%9.2f%9d%9d\n" % (metal, atom[1], metal, atom[2],1,self.calculate_angle(locals()["atom" + str(atom[1])],metal_point, locals()["atom" + str(atom[2])]),angle_strength, 1) # for streamlit
                        self.output += "%5d%6d%6d%6d%5d%9.2f%9d%9d\n" % (metal, atom[2], metal, atom[3],1,self.calculate_angle(locals()["atom" + str(atom[2])],metal_point, locals()["atom" + str(atom[3])]),angle_strength, 1) # for streamlit
                        self.output += "%5d%6d%6d%6d%5d%9.2f%9d%9d\n" % (metal, atom[0], metal, atom[3],1,self.calculate_angle(locals()["atom" + str(atom[0])],metal_point, locals()["atom" + str(atom[3])]),angle_strength, 1) # for streamlit
            
            
    def GROwriter(self, gro):
        print(self.head)
        print("%5d"  % (self.total_atom))
        for i in range(len(self.resid)):    
            print("%5d%-3s%7s%5d%8.3f%8.3f%8.3f" %  (self.resid[i], self.resname[i], self.atom[i], self.index[i], self.x[i], self.y[i], self.z[i]))
        print(self.last)

##########################################################################################################################################################################################

class gromerger(): # read uploaded files
    def __init__(self, receptor_gro, ligand_gro, ligand_itp, receptor_top, rec_name, lig_name, ligitp_name, rectop_name):
        self.merge_gro_files(receptor_gro, ligand_gro, ligand_itp, receptor_top, rec_name, lig_name, ligitp_name, rectop_name)
        
    def streamlit_download_file(self, download_name, content_file):
        # Download topol.top file  #      
        # 打开 content_file 文件并读取其内容
        with open(content_file, 'r') as top_file:
            content = top_file.read()
        
        # 添加一个下载按钮，传递 receptor_top_content 作为文件内容
        st.download_button(
            label = "Download " +  download_name,
            data = content,
            key = download_name,
            file_name = download_name
            )
    def merge_gro_files(self, receptor_gro, ligand_gro, ligand_itp, receptor_top, rec_name, lig_name, ligitp_name, rectop_name):
        # Merge the two gro files
        with open(ligand_gro, 'r') as ligand_file, open(receptor_gro, 'r') as receptor_file:
            ligand_lines = ligand_file.readlines()[1:-1]
            receptor_lines = receptor_file.readlines()
            a = int(receptor_lines[1].split()[0])
            b = int(ligand_lines[0].split()[0])
            c = a + b
            receptor_lines[1] = f"{c}\n"
            
            with open(ligand_gro, 'w') as complex_file:
                complex_file.writelines(receptor_lines[0:-1])
                complex_file.writelines(ligand_lines[1:])
                complex_file.writelines(receptor_lines[-1])
        
        # Edit the topol.top file         
        # with open('include.dat', 'w') as include_file:
        #     include_file.write("; Include ligand topology\n")
        #     include_file.write(f"#include \"{ligitp_name}\"\n")
        #     include_file.write("#ifdef POSRES_LIG\n")
        #     include_file.write("#include \"posre_lig.itp\"\n")
        #     include_file.write("#endif\n")
        
        # adding content to the end of topol file
        with open(ligand_itp, 'r') as ligand_itp_file:
            # find ligand molecule name
            nrexcl_line = None
            lines = ligand_itp_file.readlines()
            for i in range(len(lines)):
                if 'nrexcl' in lines[i]:
                    nrexcl_line = lines[i+1].strip()
                    break
            
            if nrexcl_line:
                B = nrexcl_line.split()[0]
                with open(receptor_top, 'r') as topol_file:
                    topol_lines = topol_file.readlines()
                # add lig_GMX.itp record to topol.top
                for i, line in enumerate(topol_lines):
                    if 'moleculetype' in line:
                        topol_lines.insert(i, f"#endif\n")
                        topol_lines.insert(i, f"#include \"posre_lig.itp\"\n")
                        topol_lines.insert(i, f"#ifdef POSRES_LIG\n")
                        topol_lines.insert(i, f"#include \"{ligitp_name}\"\n")
                        topol_lines.insert(i, f"; Include ligand topology\n")
                        break
                # add ligand molecule type to topol.top
                for i, line in enumerate(topol_lines):
                    if 'molecules' in line:
                        topol_lines.insert(-1, f"{B}                1\n")
                        break
                
                with open(receptor_top, 'w') as topol_file:
                    topol_file.writelines(topol_lines)
        
        
        # download topol.top
        self.streamlit_download_file(rectop_name, receptor_top)
        self.streamlit_download_file("complex.gro", ligand_gro)
##########################################################################################################################################################################################    

class contact_map_detect(): # read uploaded files
    protein = ''
    ligand = ''
    def __init__(self, topol, traj, lig, output, distance):
        contact_map = self.calculate_contact(topol, traj, lig, output, distance)
        self.plot(contact_map, output_name, distance)
        self.csv_writer(contact_map, output_name)
        
    def calculate_contact(self, topol, traj, lig, output, distance):
        # 加载蛋白质和配体的拓扑和轨迹文件
        u = mda.Universe(topol, traj)
        
        # 选择蛋白质和配体
        self.protein = u.select_atoms('protein')
        self.ligand = u.select_atoms('resname ' + lig)  
        
        # 初始化接触图矩阵
        contact_map = np.zeros((len(u.trajectory), len(self.protein.residues)))
    
        # 计算每帧的接触情况
        for ts in u.trajectory:
            y_ticks = [] 
            y_labels= []
            frame_index = ts.frame
            for i, residue in enumerate(self.protein.residues):
                min_dist = np.min(contacts.distance_array(residue.atoms.positions, self.ligand.positions))
                y_ticks.append(i)
                y_labels.append(f'{residue.resid} {residue.resname}')
                if min_dist < distance:
                    contact_map[frame_index, i] = 1
        
        return contact_map
    
    def plot(self, contact_map, output_name, distance):
        resid_list = [i for i in range(self.protein.residues.resids[0], len(self.protein.residues.resids)+self.protein.residues.resids[0], 5)]
        resname_list =[]
        for i, residue in enumerate(self.protein.residues):
            if residue.resid in resid_list:
                resname_list.append(residue.resname)
        i_list = range(len(resid_list))
        # print(resid_list)
        # print(resname_list)
        plt.imshow(contact_map.T, aspect='auto', origin='lower', cmap='Greys')
        plt.xlabel('Time (ns)')
        plt.yticks(ticks=resid_list, labels=['' if (resid_list[i]-self.protein.residues.resids[0]) % 20 != 0 else f'{resid_list[i]} {resname_list[i]}' for i in i_list])
        plt.ylabel('Residue Index')
        plt.title('Protein-Ligand Contact Map')
        plt.suptitle('Distance < ' + str(distance))
        plt.colorbar(label='Contact', ticks=[0, 1])
        # plt.show()
        plt.savefig("/tmp/" + output_name)
        self.streamlit_download_file_plotly(output_name, "/tmp/" + output_name)


    def csv_writer(self, contact_map, output_name):
        # 将contact_map数据写入CSV文件
        with open('/tmp/ligand_contact.csv', 'w', newline='') as f:
            writer = csv.writer(f)
            # 写入标题行，假设每列代表一个残基，每行代表一个时间帧
        #     header = ['Frame'] + [f'Residue_{i}' for i in range(1, len(protein.residues) + 1)]
            header = ['Time (ns)'] + [f'{residue.resid}_{residue.resname}' for i, residue in enumerate(self.protein.residues)]
            writer.writerow(header)
        
            # 写入每一行的数据
            for frame_index, contacts_ in enumerate(contact_map):
                row = [frame_index] + contacts_.tolist()
                writer.writerow(row)
        
        self.streamlit_download_file("ligand_contact.csv", '/tmp/ligand_contact.csv')
        # print("Contact map data has been written to 'ligand_contact.csv'.")
        
    def streamlit_download_file(self, download_name, content_file):
        # Download topol.top file  #      
        # 打开 content_file 文件并读取其内容
        with open(content_file, 'r') as top_file:
            content = top_file.read()
        
        # 添加一个下载按钮，传递 receptor_top_content 作为文件内容
        st.download_button(
            label = "Download " +  download_name,
            data = content,
            key = download_name,
            file_name = download_name
            )
    
    def streamlit_download_file_plotly(self, download_name, content_file):
        # 读取文件内容
        with open(content_file, "rb") as file:
            file_content = file.read()
        
        # 获取文件的 MIME 类型
        mime_type, _ = mimetypes.guess_type(content_file)
        
        # 创建下载按钮
        st.download_button(
            label=f"Download {download_name}",
            data=file_content,
            file_name=download_name,
            mime=mime_type)
    
    # download topol.top
    
##########################################################################################################################################################################################

class pep2lig(): 
    pdb     = ''
    pepname = ''
    chain   = 'A'
    resnum  = 1
    atomic_index        = []
    atomic_name         = []
    residue_name        = []
    chain_name          = []
    residue_index       = []
    X_peratom           = []
    Y_peratom           = []
    Z_peratom           = []
    bfactor_per_factor  = []
    temp_factor         = []
    Atomtype_per_atom   = []
    def __init__(self, pdb, pepname):
        self.pdb = pdb
        self.pepname = pepname
        self.converter(pdb, pepname, self.chain, self.resnum)
        self.PDBwriter("/tmp/modified_PEP.pdb")
        
    def converter(self, pdb, pepname, chain, resnum):
            with open(pdb, 'r') as infile:
                for line in infile:                                                              # iterate each line in file "f"                     
                        if(line.split()[0] in["ATOM","HETATM"]):       # Judgment Sentence，Used to split each row and then determine whether the first column of the row == ATOM or HETATM
                            self.atomic_index.append(int(line[6:11].strip()))                # The second column is the atomic number
                            self.atomic_name.append(line[12:16].strip())                        # The 3rd column is the atom name C CA CD1 CD2 and so on
                            # self.residue_name.append(line[17:20].strip())                       # Column 4 is the residue name TYR ALA etc.
                            self.residue_name.append(pepname)
                            # self.chain_name.append(line[21].strip())                         # The 5th column is the name of the chain it is on
                            self.chain_name.append(chain)
                            # self.residue_index.append(int(line[22:26].strip()))               # The sixth column is the residue number
                            self.residue_index.append(resnum)
                            self.X_peratom.append(float(line[30:38].strip()))
                            self.Y_peratom.append(float(line[38:46].strip()))
                            self.Z_peratom.append(float(line[46:54].strip()))
                            self.bfactor_per_factor.append(float(line[54:60].strip()) if line[54:60].strip() else 0.0)
                            self.temp_factor.append(float(line[60:66].strip()) if line[60:66].strip() else 0.0 )
                            try:
                                self.Atomtype_per_atom.append(line[76:78].strip())
                            except:
                                self.Atomtype_per_atom.append(" ")              
        
    def PDBwriter(self,filename):
        f = open(filename, "w")                                                             # e.g: f = linesplit[0]+"_PO3.pdb"
        for i in range (0 ,len(self.atomic_index)):                                         # Create a loop, i is a sequence starting from 0, and the number of atoms is the length  
            print("%4s%7d  %-4s%1s%2s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f%12s" %  ("ATOM" ,     # Formatted output, %4s, right-aligned, the output occupies 4 columns in total. If the length is less than 4 columns, the left end will be filled with spaces. If it is greater than 4 columns, the actual length will be output as a string
                                             self.atomic_index[i],                          # %7d, right-aligned, the output occupies a total of 7 columns, if the length is less than 7 columns, the left end is filled with spaces, signed decimal certificate integer
                                             self.atomic_name[i],                           # %-4s, left-aligned, the output occupies a total of 4 columns, if the length is less than 4 columns, the right end is filled with spaces, if it is greater than 4 columns, the actual length is output as a string
                                             self.residue_name[i],                          # %1s, right-aligned, the output occupies a total of 1 column. If it is less than 1 column, it will be filled with spaces from the left end. If it is greater than 1 column, the actual length will be output as a string
                                             "A",                            # %2s, right-aligned, the output occupies 2 columns in total. If it is less than 2 columns, it will be filled with spaces from the left end. If it is greater than 2 columns, the actual length will be output as a string
                                             self.residue_index[i],                         # %4d, right-aligned, the output occupies a total of 4 columns, if the length is less than 4 columns, the left end is filled with spaces, a signed decimal certificate integer
                                             self.X_peratom[i],                             # %8.3f, right-aligned, the output occupies a total of 8 columns, including 3 decimal places, if the width of the value is less than 8, fill in a space at the left end, decimal
                                             self.Y_peratom[i],                             # %8.3f, right-aligned, the output occupies a total of 8 columns, including 3 decimal places, if the width of the value is less than 8, fill in a space at the left end, decimal
                                             self.Z_peratom[i],                             # %8.3f, right-aligned, the output occupies a total of 8 columns, including 3 decimal places, if the width of the value is less than 8, fill in a space at the left end, decimal
                                             self.bfactor_per_factor[i],                    # %6.2f, right-aligned, the output occupies a total of 6 columns, including 2 decimal places, if the width of the value is less than 6, fill in a space at the left end, decimal
                                             self.temp_factor[i],                     # %6.2f, right-aligned, the output occupies a total of 6 columns, including 2 decimal places, if the width of the value is less than 6, fill in a space at the left end, decimal
                                             self.Atomtype_per_atom[i]), file = f )         # %12s, right-aligned, the output occupies a total of 12 columns, if it is less than 12 columns, it will be filled with spaces from the left end
        print("END", file = f)
        f.close()
        
        self.streamlit_download_file("modified_peptide.pdb", filename)

        
    def streamlit_download_file(self, download_name, content_file):
        # Download topol.top file  #      
        # 打开 content_file 文件并读取其内容
        with open(content_file, 'r') as top_file:
            content = top_file.read()
        
        # 添加一个下载按钮，传递 receptor_top_content 作为文件内容
        st.download_button(
            label = "Download " +  download_name,
            data = content,
            key = download_name,
            file_name = download_name
            )
##########################################################################################################################################################################################

# Title
st.title("Welcome to gmx tool box")

# 创建4栏布局
plot, mradder, gromerge, contact_map = st.columns(4)

# 保存文件并获取临时文件路径
def save_uploaded_file(uploaded_file):
    try:
        with NamedTemporaryFile(delete=False, suffix='.' + uploaded_file.name.split('.')[-1]) as temp_file:
            temp_file.write(uploaded_file.getvalue())
            return temp_file.name
    except Exception as e:
        st.error(f"Error saving file: {e}")
        return None

# 在第一栏中添加内容
with plot:
    st.header("Plot figures")
    # Collecting user inputs
    multi_files = st.file_uploader("Upload files", accept_multiple_files=True, type=['xvg', 'gro', 'itp', 'top', 'csv'])
    # 保存上传文件的文件名
    uploaded_filenames = [uploaded_file.name for uploaded_file in multi_files]
    # 保存上传的文件到临时位置
    tmp_path = [save_uploaded_file(multi_files[i]) for i in range(len(multi_files))]
    output_name = st.text_input("Output file name", 'output.png')
    renumber = st.selectbox("Renumber residues", ['false', 'true'])
    rdf_cutoff = st.number_input("RDF cutoff value", min_value=0.0, step=0.1, value=0.0)
    average = st.selectbox("Calculate average", ['false', 'true'])
    plot_name = st.text_input("Plot title", value="auto detect")
    nbin = st.number_input("Number of bins (for PCA)", min_value=1, step=1, value=1)
    size = st.number_input("Size (for PCA)", min_value=1, step=1, value=500)
    move_average = st.number_input("Window size for moving average", min_value=0, step=1, value=0)
    mean_value = st.selectbox("Draw mean value line", ['false', 'true'])
    histogram = st.selectbox("Generate histogram plot", ['false', 'true'])
    xaxis_name = st.text_input("X-axis name", value="auto detect")
    yaxis_name = st.text_input("Y-axis name", value="auto detect")
    xaxis_size = st.number_input("Width of plot", min_value=0, step=100, value=800)
    yaxis_size = st.number_input("Height of plot", min_value=0, step=100, value=600)
    xy_font = st.number_input("XY-axis font size", min_value=0, step=1, value=40)
    title_font = st.number_input("Title font size", min_value=0, step=1, value=24)
    legend_show = st.selectbox("Show legend", ['True', 'False'])
    legend_show = ast.literal_eval(legend_show)
    legend_font = st.number_input("Legend font size", min_value=0, step=1, value=30)
    font_family = st.text_input("Font family", 'Arial')
    grid_show = st.selectbox("Show grid", ['True', 'False'])
    grid_show = ast.literal_eval(grid_show)
    margin_l  = st.number_input("left margin", min_value=0, step=10, value=100)
    margin_r  = st.number_input("right margin", min_value=0, step=10, value=60)
    margin_t  = st.number_input("top margin", min_value=0, step=10, value=60)
    margin_b  = st.number_input("bottom margin", min_value=0, step=10, value=100)
    font_color = st.color_picker("Font color", '#000000')
    violin = st.selectbox("violin style", ['False', 'True'])

    if st.button('Plotting') and multi_files[0] != 0:
        x = plotly_go(tmp_path, output_name, renumber, rdf_cutoff, average, plot_name, nbin, size, move_average, mean_value, histogram, xaxis_name, yaxis_name, xaxis_size, yaxis_size, xy_font, title_font, legend_show, legend_font, font_family, font_color, grid_show, uploaded_filenames, margin_l, margin_r, margin_t, margin_b, violin)

# 在第二栏中添加内容
with mradder:
    st.header("Add restraints")
    mr_file = st.file_uploader("Upload gro file for metal restraints", type='gro')
    mr_num_neighbours = st.number_input("Number of neighbours for metal restraints", min_value=1, step=1, value=3)
    mr_distance_value = st.number_input("Distance value for metal restraints", min_value=0.0, step=0.1, value=0.4)
    mr_atom_list = st.text_area("Atom list for metal restraints", value="ND1 OE1 OE2 OD1 OD2 ND2 SG NZ OH")
    mr_metal_list = st.text_area("Metal list", value="MG MN ZN CA")
    mr_residue_list = st.text_area("Residue list", value="HIS GLU ASP ASN CYS LYS TYR")
    mr_bond_strength = st.number_input("Bond strength for metal restraints", min_value=0, step=1000, value=200000)
    mr_angle_strength = st.number_input("Angle strength for metal restraints", min_value=0, step=1000, value=10000)
    if st.button('Add') and mr_file != 0:
        file_content = mr_file.getvalue().decode("utf-8")
        x = mr(file_content, mr_num_neighbours, mr_distance_value, mr_atom_list, mr_metal_list, mr_residue_list, mr_bond_strength, mr_angle_strength)
        st.text(x.output)
    
# 在第三栏中添加内容
with gromerge:
    st.header("Merge gro files")
    gm_receptor_gro = st.file_uploader("Upload the receptor GRO file", type=['gro'])
    gm_ligand_gro = st.file_uploader("Upload the ligand GRO file", type=['gro'])
    gm_ligand_itp = st.file_uploader("Upload the ligand ITP file", type=['itp'])
    gm_receptor_top = st.file_uploader("Upload the receptor TOP file", type=['top'])
    # save the files name:
    if gm_receptor_gro and gm_ligand_gro and gm_ligand_itp and gm_receptor_top:
        receptor_gro_filename = gm_receptor_gro.name
        ligand_gro_filename = gm_ligand_gro.name
        ligand_itp_filename = gm_ligand_itp.name
        receptor_top_filename = gm_receptor_top.name
    # 保存上传的文件到临时位置
    if st.button('Merge') and gm_receptor_gro and gm_ligand_gro and gm_ligand_itp and gm_receptor_top:
        # 保存文件并获取临时文件路径
        receptor_gro_path = save_uploaded_file(gm_receptor_gro)
        ligand_gro_path = save_uploaded_file(gm_ligand_gro)
        ligand_itp_path = save_uploaded_file(gm_ligand_itp)
        receptor_top_path = save_uploaded_file(gm_receptor_top)

        # 如果文件保存成功，则进行合并
        if receptor_gro_path and ligand_gro_path and ligand_itp_path and receptor_top_path:
            with open(receptor_gro_path, 'r') as file:
                file_content = file.read()
                # st.text(file_content)
            x = gromerger(receptor_gro_path, ligand_gro_path, ligand_itp_path, receptor_top_path, receptor_gro_filename, ligand_gro_filename, ligand_itp_filename, receptor_top_filename)
            st.success("Files merged successfully!")
        else:
            st.error("Failed to save files paths.")
    else:
        pass

    # if st.button('Merge') and gm_receptor_gro != 0:
    #     x = gromerger(gm_receptor_gro, gm_ligand_gro, gm_ligand_itp, gm_receptor_top )

# 在第4栏中添加内容
with contact_map:
    st.header("Ligand contact map")
    cm_topol_top = st.file_uploader("Upload the topology file", type=['gro', 'top', 'pdb', 'tpr'])
    cm_trajectory_pdb = st.file_uploader("Upload the trajectory pdb file", type=['xtc', 'trr', 'pdb'])
    cm_ligand_name = st.text_input("ligand res name", 'LIG')
    cm_output_name = st.text_input("Output file name", 'contact.png')
    cm_distance = st.number_input("Threshold value for the distance (Angstrom)", min_value=1.0, step=0.1, value=3.5)
    # save the files name:
    if cm_topol_top and cm_trajectory_pdb and cm_ligand_name and cm_output_name:
        topol_top_filename = cm_topol_top.name
        trajectory_pdb_filename = cm_trajectory_pdb.name

    # 保存上传的文件到临时位置
    if st.button('detact') and cm_topol_top and cm_trajectory_pdb:
        # 保存文件并获取临时文件路径
        topol_top_path = save_uploaded_file(cm_topol_top)
        trajectory_pdb_path = save_uploaded_file(cm_trajectory_pdb)
        # 如果文件保存成功，则进行合并
        if topol_top_path and trajectory_pdb_path:
            # with open(receptor_gro_path, 'r') as file:
                # file_content = file.read()
                # st.text(file_content)
            x = contact_map_detect(topol_top_path, trajectory_pdb_path, cm_ligand_name, cm_output_name, cm_distance)
            st.success("Ligand contact map generated successfully!")
        else:
            st.error("Failed to generate Ligand contact map.")
    else:
        pass

    
    ### subcolum in column 4, for peptide convert to ligand ###
    contact_map.subheader("Peptide to Ligand")
    p2l_pdb = st.file_uploader("Upload the Peptide pdb file", type=['pdb'])
    p2l_name =  st.text_input("Give the peptide a name", 'PEP') 
    # save the files name:
    if p2l_pdb and p2l_name:
        peptide_pdb = p2l_pdb.name
    # 保存上传的文件到临时位置
    if st.button('convert') and p2l_pdb and p2l_name:
        # 保存文件并获取临时文件路径
        peptide_pdb_path = save_uploaded_file(p2l_pdb)
        # 如果文件保存成功，则进行合并
        if peptide_pdb_path:
            # with open(receptor_gro_path, 'r') as file:
                # file_content = file.read()
                # st.text(file_content)
            x = pep2lig(peptide_pdb_path, p2l_name)
            # st.success("Peptide converted successfully!")
        else:
            st.error("Failed to convert the Peptide, do you have the last column (atom type column) in you peptide PDB?")
    else:
        pass           
#################################################################################################################################################
