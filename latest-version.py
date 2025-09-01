#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Version 2.0 Add x_scale and y_scale for plotly;continus error band added, DSSP(+percentage plot), Renumber MODELS in PDB; modified the restraints adder function; Fit to PLotly 6.0.0
Created on Tue 16:15:48 2024-05-14

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
# import rpy2.robjects as ro
# from rpy2.robjects import pandas2ri
# from rpy2.robjects.packages import importr
# from rpy2.robjects import conversion, default_converter

# for free energy
import io

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
# for plot dssp
import json
# for RMSD per Residue
from Bio.PDB import PDBParser
import matplotlib.colors as mcolors
# for peptide to ligand
import zipfile


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

    def __init__(self, multi_files, output_name, renumber, rdf_cutoff, average, ls
                 , nbin, size, move_average, mean_value, histogram, xaxis_name, yaxis_name, xaxis_size, yaxis_size, xy_font, title_font, legend_show, legend_font, font_family, font_color, grid_show, uploaded_filenames, l,r,t,b, violin, smooth, error_bar, replica_number, axis_show, line_width, transparency
                 , x_low, x_high, y_low, y_high, traces_color_schemes, pca_color_by_replicas):
        trace_color_scheme = [
            c.strip()
            for c in re.split(r"[,\s]+", traces_color_schemes.strip())  # 逗号或任意空白都能分隔
            if c                                             # 去掉可能出现的空元素
        ]

        if len(multi_files) >=1:
            pca_color_by_replica = int(pca_color_by_replicas)
            # print(multi_files)
            file1 = multi_files[0]
            self.flag_recognizer(file1, plot_name)
            if self.pca_flag != 1 and self.flag != 'pca' and self.flag != 'free energy':
                self.plotly_multy(multi_files, output_name, renumber, rdf_cutoff, average, plot_name, nbin, size, move_average, mean_value, histogram, xaxis_name, yaxis_name, xaxis_size, yaxis_size, xy_font, title_font, legend_show, legend_font, font_family, font_color, grid_show, self.flag, uploaded_filenames, l,r,t,b, violin, smooth, error_bar, replica_number, axis_show, line_width, transparency, x_low, x_high, y_low, y_high, trace_color_scheme)
            elif self.pca_flag == 1:
                self.plotly_pca(multi_files, output_name, renumber, rdf_cutoff, average, plot_name, nbin, size, move_average, mean_value, histogram, xaxis_name, yaxis_name, xaxis_size, yaxis_size, xy_font, title_font, legend_show, legend_font, font_family, font_color, grid_show, self.flag, uploaded_filenames, l,r,t,b, smooth, axis_show, line_width, x_low, x_high, y_low, y_high, trace_color_scheme, pca_color_by_replica)
            elif self.flag == 'pca':
                self.plotly_pca(multi_files, output_name, renumber, rdf_cutoff, average, plot_name, nbin, size, move_average, mean_value, histogram, xaxis_name, yaxis_name, xaxis_size, yaxis_size, xy_font, title_font, legend_show, legend_font, font_family, font_color, grid_show, self.flag, uploaded_filenames, l,r,t,b, smooth, axis_show, line_width, x_low, x_high, y_low, y_high, trace_color_scheme, pca_color_by_replica)
            elif self.flag == 'free energy':
                self.plotly_free_energy(multi_files, output_name, plot_name, nbin, size, xaxis_name, yaxis_name, xaxis_size, yaxis_size, xy_font, title_font, legend_show, legend_font, font_family, font_color, grid_show, self.flag, uploaded_filenames, l,r,t,b, violin, smooth, error_bar, replica_number, axis_show, line_width, transparency, x_low, x_high, y_low, y_high, trace_color_scheme)

    def flag_recognizer(self,file1, plot_name):                                                   # first method to be called in __main__, used for creating object and charactors.
        flags_map = {
            'rms,': 'rmsd',
            'rmsd' : 'rmsd',
            'rmsf,': 'rmsf',
            'rmsf' : 'rmsf',
            'sasa,': 'sasa',
            'sasa' : 'sasa',
            'gyrate,': 'gyrate',
            'gyrate' : 'gyrate',
            'dipoles,': 'dipoles',
            'dipoles' : 'dipoles',
            'distance,': 'distance',
            'distance' : 'distance',
            'rdf,': 'rdf',
            'rdf' : 'rdf',
            'convergence': 'convergence',
            'anaeig,': 'pca',
            'pca' : 'pca',
            'angle,': 'angle',
            'angle' : 'angle',
            'free': 'free energy'
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
            found = False
            for key in flags_map:
                if key.strip(',') in plot_name.lower():  
                    found = True
                    self.flag = flags_map[key]  
                    break  
        
        elif file1.endswith(".dat"):
            with open(file1, 'r') as f:
                lines = f.readlines()            
                first_line = lines[0]
                flag = 'free' if 'free' in first_line else None   
                self.flag = flags_map.get(flag, flag)   
            print("I know you are plotting " + self.flag + " figures!")


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
     
    def read_data_xvg(self, file, x_name, renumber):
        # 从文件中读取数据   
        x_data, y_data, sd_data = [], [], []
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
        return x_data, y_data, sd_data  
    
    def read_data_csv(self, file, x_name, renumber):
        x_data, y_data, sd_data = [], [], []
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
    
    def read_data_dat(self, file_name):
        with open(file_name, 'r') as file:
            lines = file.readlines()
        # Identifying the line with column names
        columns_line = [line for line in lines if line.startswith('#! FIELDS')][0]
        # Extracting column names
        column_names = columns_line.strip().split()[2:]  # Skip '#! FIELDS'
        # 判断free_energy data 在第几列
        index_of_free_energy = [index for index, name in enumerate(column_names) if 'free' in name][0]
        # Extracting data lines (those not starting with '#')
        data_lines = [line for line in lines if not line.startswith('#')]
        # Converting data lines to a pandas DataFrame
        df = pd.read_csv(io.StringIO('\n'.join(data_lines)), delim_whitespace=True, names=column_names)
        
        if index_of_free_energy == 2:
            x_data = df[df.columns[0]].tolist()  # 第一列作为X轴数据
            y_data = df[df.columns[1]].tolist()  # 第2列作为y轴数据
            z_data = df[df.columns[2]].tolist()  # 第3列作为z轴数据
        elif index_of_free_energy == 1:
            x_data = df[df.columns[0]].tolist()  # 第一列作为X轴数据
            y_data = df[df.columns[1]].tolist()  # 第2列作为y轴数据
            z_data = []

        # CSV文件不包含标准差数据，因此sd_data保持为空 
                   
        return x_data, y_data, z_data, df, index_of_free_energy, column_names
    
    
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


    def define_trace(self, x_data, y_data, file_name, colour, violine='False', flag=0, labels=0, smooth=0):
        # 创建并返回迹线
        if flag == 'pca' and smooth != 'true':
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
        elif flag =='pca' and smooth == 'true':
            trace = go.Heatmap(z=x_data, colorscale='Viridis', showscale=True, connectgaps=True, zsmooth='best')
        elif flag =='angle' and smooth == 'true':
            trace = go.Heatmap(z=x_data, colorscale='Viridis', showscale=True, connectgaps=True, zsmooth='best', x=[-180, -120, -60, 60, 120,180], y=[-180, -120, -60, 60, 120,180])   
        elif flag not in ['pca', 'angle'] and smooth == 'true':
            trace = go.Heatmap(z=x_data, colorscale='Viridis', showscale=True, connectgaps=True, zsmooth='best')   
        elif violine != 'False':
            trace = go.Violin(x0=str(file_name).split('.')[0], y=y_data, line=dict(color='black'), fillcolor=colour, name=str(file_name).split('.')[0], box_visible=True, meanline_visible=True, opacity=0.6)            
        else:
            trace = go.Scatter(x=x_data, y=y_data, line=dict(color=colour), name=str(file_name).split('.')[0])
        return trace

    def calculate_for_error_bar_or_band(self, multi_files, x_name, replica_number, uploaded_filenames):
        df_data = pd.DataFrame()
        df_average = pd.DataFrame()
        df_sd   = pd.DataFrame()
        count = 1
        if multi_files[0].endswith(".xvg"):
            for i, file in enumerate(multi_files):
                x_data, y_data, _ = self.read_data_xvg(file, x_name, renumber)
                if i == 0:
                    x_datas = x_data
                df_data[f'y_data_{i+1}'] = y_data
                if i == (count * replica_number) - 1:  # 检查是否达到组内文件数量
                    # 计算当前df_data的所有列的平均值和标准差
                    mean_vals = df_data.mean(axis=1)
                    std_vals = df_data.std(axis=1)
                    # 将计算得到的平均值和标准差添加到相应的DataFrame中
                    raw_name  = uploaded_filenames[(count-1)*3]
                    legend    = os.path.splitext(os.path.basename(raw_name))[0]
                    legend   = legend[:-2]
                    df_average[legend] = mean_vals
                    df_sd[legend] = std_vals
                    # 重置df_data以便下一组的使用，并更新计数器
                    df_data = pd.DataFrame()
                    count += 1
        elif multi_files[0].endswith(".csv"):
            for i, file in enumerate(multi_files):
                x_data, y_data, _ = self.read_data_csv(file, x_name, renumber) 
                if i == 0:
                    x_datas = x_data
                df_data[f'y_data_{i+1}'] = y_data
                if i == (count * replica_number) - 1:  # 检查是否达到组内文件数量
                    # 计算当前df_data的所有列的平均值和标准差
                    mean_vals = df_data.mean(axis=1)
                    std_vals = df_data.std(axis=1)
                    # 将计算得到的平均值和标准差添加到相应的DataFrame中
                    df_average[uploaded_filenames[(count-1)*3]] = mean_vals
                    df_sd[uploaded_filenames[(count-1)*3]] = std_vals
                    # 重置df_data以便下一组的使用，并更新计数器
                    df_data = pd.DataFrame()
                    count += 1
        elif multi_files[0].endswith(".dat"):
            for i, file in enumerate(multi_files):
                x_data, y_data, z_data, df, index_of_free_energy, column_names = self.read_data_dat(file, x_name, renumber) 
                if i == 0:
                    x_datas = x_data
                df_data[f'y_data_{i+1}'] = y_data
                if i == (count * replica_number) - 1:  # 检查是否达到组内文件数量
                    # 计算当前df_data的所有列的平均值和标准差
                    mean_vals = df_data.mean(axis=1)
                    std_vals = df_data.std(axis=1)
                    # 将计算得到的平均值和标准差添加到相应的DataFrame中
                    df_average[uploaded_filenames[(count-1)*3]] = mean_vals
                    df_sd[uploaded_filenames[(count-1)*3]] = std_vals
                    # 重置df_data以便下一组的使用，并更新计数器
                    df_data = pd.DataFrame()
                    count += 1
        return df_average, df_sd, x_datas
    def hex_to_rgba(self, hex_color, alpha=0.2):
        # 移除可能的 "#" 符号
        hex_color = hex_color.lstrip('#')
        # 通过列表推导从十六进制字符串中提取并转换为RGB整数值
        rgb = tuple(int(hex_color[i:i+2], 16) for i in (0, 2, 4))
        # 将RGB整数值和透明度alpha组合成RGBA字符串
        return f"rgba({rgb[0]}, {rgb[1]}, {rgb[2]}, {alpha})"
   
    def define_trace_for_error_bands(self, error_bar, df_average, df_sd, x_data, transparency, color_scheme):
        Plotly = ['#636EFA', '#EF553B', '#00CC96', '#AB63FA', '#FFA15A', '#19D3F3', '#FF6692', '#B6E880', '#FF97FF', '#FECB52']
        Plotly = color_scheme
        Dark24 = ['#2E91E5', '#E15F99', '#1CA71C', '#FB0D0D', '#DA16FF', '#222A2A', '#B68100', '#750D86', '#EB663B', '#511CFB', '#00A08B', '#FB00D1', '#FC0080', '#B2828D', '#6C7C32', '#778AAE', '#862A16', '#A777F1', '#620042', '#1616A7', '#DA60CA', '#6C4516', '#0D2A63', '#AF0038']
        AMPK_color = ['#222A2A', '#FB0D0D', '#2E91E5']
        # AMPK_color = ['#2E91E5']
        colors = ['rgb(0,100,80)', 'rgb(255,0,0)']  # 不同组使用不同颜色
        traces = []
        if error_bar =='error band':
            for idx, (col_name_avg, col_name_sd) in enumerate(zip(df_average.columns, df_sd.columns)):
                y = df_average[col_name_avg]
                y_std = df_sd[col_name_sd]
                y_upper = y + y_std
                y_lower = y - y_std
                fill_color = self.hex_to_rgba(Plotly[idx], alpha=transparency)
                
                traces.append(go.Scatter(
                    x=x_data,
                    y=y,
                    line=dict(color=Plotly[idx]),
                    mode='lines',
                    name=col_name_avg  # 使用列名作为轨迹名称
                ))
                traces.append(go.Scatter(
                    x=list(x_data) + list(x_data)[::-1],
                    y=list(y_upper) + list(y_lower)[::-1],
                    fill='toself',
                    fillcolor=fill_color,
                    line=dict(color='rgba(255,255,255,0)'),
                    hoverinfo="skip",
                    showlegend=False
                ))
        elif error_bar == 'error bar':
            for idx, (col_name_avg, col_name_sd) in enumerate(zip(df_average.columns, df_sd.columns)):
                y = df_average[col_name_avg]
                y_std = df_sd[col_name_sd]
                traces.append(go.Scatter(
                    x=x_data,
                    y=df_average[col_name_avg],
                    error_y=dict(type='data', array=df_sd[col_name_sd], visible=True)
                ))
        return traces

    def setup_layout(self, plot_title, title_font, x_name, y_name, xy_font, xaxis_size, yaxis_size, font_color, legend_show, legend_font, font_family, grid_show, l,r,t,b, x_low, x_high, y_low, y_high, violine='False', flag=0, axis_shows='True', line_width=2):
        # 设置布局
        if flag == 'pca':
            legend_show = False
        if violine != 'False':
            x_name = ''
        if y_low == y_high == 0:
            y_autorange = True
        else:
            y_autorange = False
        if x_low == x_high == 0:
            x_autorange = True
        else:
            x_autorange = False

        layout = go.Layout(
            title=plot_title, title_x=0.5, title_y=0.99, font=dict(size=title_font, color=font_color),
            xaxis=dict(title=dict(text=x_name, font=dict(size=xy_font, color=font_color, family=font_family)), zeroline=False, autorange=x_autorange, range=[x_low,x_high],
                       showgrid=grid_show, gridwidth=1, gridcolor='rgba(235,240,248,100)', tickfont=dict(size=30), showline=axis_shows, linewidth=line_width, linecolor=font_color),
            yaxis=dict(title=dict(text=y_name, font=dict(size=xy_font, color=font_color, family=font_family)), zeroline=False, autorange=y_autorange, range=[y_low,y_high],
                       showgrid=grid_show, gridwidth=1, gridcolor='rgba(235,240,248,100)', tickfont=dict(size=30), showline=axis_shows, linewidth=line_width, linecolor=font_color),
            legend=dict(x=1, y=1, orientation='v', font=dict(size=legend_font, color=font_color)), showlegend=legend_show,
            plot_bgcolor='rgba(255, 255, 255, 0.1)',
            paper_bgcolor='rgba(255, 255, 255, 0.2)',
            margin=dict(l=int(l), r=int(r), t=int(t), b=int(b)),
            width=xaxis_size, height=yaxis_size
        )
        return layout
    def pca_bins_density_define(self, nbins, data):
        # 确定边界
        x_min, y_min = np.min(data, axis=0)
        x_max, y_max = np.max(data, axis=0)

        # 创建网格
        n_bins = nbins
        x_bins = np.linspace(x_min, x_max, n_bins + 1)
        y_bins = np.linspace(y_min, y_max, n_bins + 1)

        # 计算每个格子内的点的数量（密度）
        density_matrix = np.zeros((n_bins, n_bins))
        for x, y in data:
            x_idx = np.digitize(x, x_bins) - 1
            y_idx = np.digitize(y, y_bins) - 1

            # 确保索引不超出density_matrix的范围
            x_idx = min(x_idx, n_bins - 1)
            y_idx = min(y_idx, n_bins - 1)

            density_matrix[y_idx, x_idx] += 1  # 注意矩阵的索引和坐标系的方向

        return density_matrix
    
    def plot_heatmap(self, x_points, y_points, plot_title, x_name, y_name, size, output_name):
        pandas2ri.activate()
        
        # 用Python计算高度和宽度
        x_diff = np.max(x_points) - np.min(x_points)
        y_diff = np.max(y_points) - np.min(y_points)
        width= (x_diff / y_diff) * size
        height = 1 * size
        with conversion.localconverter(default_converter):
        # 将Python列表转换为R向量，并创建data.frame
            data = ro.DataFrame({x_name: ro.FloatVector(x_points), y_name: ro.FloatVector(y_points)})
            ro.globalenv['data'] = data
        
            # 手动定义RdYlBu颜色渐变
            rdylbu_colors = [
                "#d73027",  # 红色
                "#f46d43",  # 红橙色
                "#fdae61",  # 橙黄色
                "#fee090",  # 浅黄色
                "#ffffbf",  # 最浅的黄色
                "#e0f3f8",  # 浅蓝色
                "#abd9e9",  # 天蓝色
                "#74add1",  # 亮蓝色
                "#4575b4",  # 蓝色
            ]
            
        
            
            # 执行其他R指令
            ro.globalenv['zBot'] = 1.52
            ro.globalenv['zTop'] = 3.42
            ro.globalenv['zW'] = 0.83
            # 在R全局环境中创建颜色向量
            # ro.globalenv['buylrd'] = ro.StrVector(rdylbu_colors)
            # color_ramp_function = ro.r('colorRamp')(ro.globalenv['buylrd'])  # 从颜色向量创建颜色渐变函数
            # ro.globalenv['color_ramp'] = color_ramp_function  # 将颜色渐变函数保存到全局环境中
            ro.r('library(RColorBrewer)')
            ro.globalenv['buylrd'] = ro.r('rev(brewer.pal(11,"RdYlBu"))')
            
            # 设置文件名并保存为PNG
            if output_name == 'output.png':
                output_filename = "DensityMap_output.png"  # 示例文件名，你可能需要根据multi_files调整
            else:
                output_filename = "DensityMap_" + output_name
            
            # 绘制图形并保存
            grDevices = importr('grDevices')
            grDevices.png(file="/tmp/" + output_filename, height=height, width=width)
            # ro.r('smoothScatter(data$%s ~ data$%s, colramp = color_ramp, nrpoints=Inf, pch="", cex=.7, col="black", main=%r, xlab=%r, ylab=%r, transformation = function(x) x^.45)' % (y_name, x_name, plot_title, x_name, y_name))
            ro.r('smoothScatter(data$y_name ~ data$x_name, colramp=colorRampPalette(c(buylrd)), nrpoints=Inf, pch="", cex=.7, col="black", main="Title")')
            # ro.r('''
    #  smoothScatter(data$%s ~ data$%s, colramp=color_ramp, nrpoints=Inf, pch="", cex=.7, col="black", main=%r, xlab=%r, ylab=%r, transformation=function(x) x^.45)
    # ''' % (y_name, x_name, plot_title, x_name, y_name))
            grDevices.dev_off()
        
        pandas2ri.deactivate()
        # Download the file
        self.streamlit_download_file_plotly(output_filename, "/tmp/" + output_filename)
    
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

    
    def split_even_ranges(self, n_points: int, k: int):
        """
        将 n_points 个点平均切分为 k 组，返回 [(start, end), ...] 半开区间。
        前 remainder 个分组多 1 个点，保证尽可能均分。
        """
        k = max(1, min(k, n_points))  # 防越界
        base = n_points // k
        rem  = n_points %  k
        ranges = []
        start = 0
        for i in range(k):
            size = base + (1 if i < rem else 0)
            end  = start + size
            ranges.append((start, end))
            start = end
        return ranges


    def plotly_multy(self, multi_files, output_name, renumber, rdf_cutoff, average, plot_name, nbin, size, move_average, mean_value, histogram, xaxis_name, yaxis_name, xaxis_size, yaxis_size, xy_font, title_font, legend_show, legend_font, font_family, font_color, grid_show, flag, uploaded_filenames, l,r,t,b, violin, smooth, error_bar, replica_number, axis_show, linewidth, transparency, x_low, x_high, y_low, y_high, trace_color):
        # Plotly = ['#636EFA', '#EF553B', '#00CC96', '#AB63FA', '#FFA15A', '#19D3F3', '#FF6692', '#B6E880', '#FF97FF', '#FECB52']
        Plotly = trace_color
        data, histogram_data, group_labels = [], [], []

        # 读取plot_title, x_name, y_name
        plot_title, x_name, y_name, traces_name_list = self.extract_plot_details(multi_files, plot_name, xaxis_name, yaxis_name, flag, histogram)
        # 读取数据并创建迹线
        if multi_files[0].endswith(".xvg"):
            for i, file in enumerate(multi_files):
                x_data, y_data, _ = self.read_data_xvg(file, x_name, renumber)
                points = [list(pair) for pair in zip(x_data, y_data)]
                density_matrix = self.pca_bins_density_define(nbin, points)
                # 使用 define_trace 创建迹线
                if smooth == 'true':
                    trace = self.define_trace(density_matrix, density_matrix, file, 'rainbow', flag=flag, smooth=smooth)  # 假设使用 'rainbow' 作为颜色
                else:
                    trace = self.define_trace(x_data, y_data, uploaded_filenames[i], Plotly[i % len(Plotly)], violine=violin)
                data.append(trace)
                # 添加直方图数据
                if histogram == 'true':
                    histogram_data.append(y_data)
                    group_labels.append(str(uploaded_filenames[i]).split('.')[0])
        elif multi_files[0].endswith(".csv"):
            for i, trace in enumerate(traces_name_list):
                x_data, y_data, _ = self.read_data_csv(multi_files[0], x_name, renumber)
                points = [list(pair) for pair in zip(x_data, y_data)]
                density_matrix = self.pca_bins_density_define(nbin, points)
                # 使用 define_trace 创建迹线
                if smooth == 'true':
                    trace = self.define_trace(density_matrix, density_matrix, file, 'rainbow', flag=flag, smooth=smooth)  # 假设使用 'rainbow' 作为颜色
                else:
                    trace = self.define_trace(x_data, y_data, uploaded_filenames[i], Plotly[i % len(Plotly)], violine=violin)
                data.append(trace)
    
                # 添加直方图数据
                if histogram == 'true':
                    histogram_data.append(y_data)
                    group_labels.append(str(file).split('.')[0])
    
        # change Time (ps) to Time (ns)
        if x_name == 'Time (ps)':
            x_name = 'Time (ns)'
        
        # 设置布局
        layout = self.setup_layout(plot_title, title_font, x_name, y_name, xy_font, xaxis_size, yaxis_size, font_color, legend_show, legend_font, font_family, grid_show, l, r, t , b, x_low, x_high, y_low, y_high, violine=violin, axis_shows=axis_show, line_width=linewidth)

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
        
        # 处理error band or error bar

        if error_bar != 'false':
            plot_title, x_name, y_name, traces_name_list = self.extract_plot_details(multi_files, plot_name, xaxis_name, yaxis_name, flag, histogram)
            df_average, df_sd, x_data = self.calculate_for_error_bar_or_band(multi_files, x_name, replica_number, uploaded_filenames)
            error_data = self.define_trace_for_error_bands(error_bar, df_average, df_sd, x_data, transparency, Plotly)
            # change Time (ps) to Time (ns)
            if x_name == 'Time (ps)':
                x_name = 'Time (ns)'
            error_layout = self.setup_layout(plot_title, title_font, x_name, y_name, xy_font, xaxis_size, yaxis_size, font_color, legend_show, legend_font, font_family, grid_show, l, r, t , b, x_low, x_high, y_low, y_high, violine=violin, axis_shows=axis_show, line_width=linewidth)
            self.plot_graph(error_data, error_layout, "error_bar_" + output_name)

    def plotly_pca(self, multi_files, output_name, renumber, rdf_cutoff, average, plot_name, nbin, size, move_average, mean_value, histogram, xaxis_name, yaxis_name, xaxis_size, yaxis_size, xy_font, title_font, legend_show, legend_font, font_family, font_color, grid_show, flag, uploaded_filenames, l, r, t ,b, smooth, axis_show, linewidth, x_low, x_high, y_low, y_high, trace_color, pca_color_by_replica):
        data = []
        if pca_color_by_replica == 0:
            color = ['rainbow']
        elif pca_color_by_replica !=0:
            color = trace_color[:pca_color_by_replica]
        # labels = []
        # 使用 extract_plot_details 方法获取图表标题和轴标签
        plot_title, x_name, y_name, traces_name_list = self.extract_plot_details(multi_files, plot_name, xaxis_name, yaxis_name, flag, histogram)
        if xaxis_name == 'auto detect':
            x_name = 'PC1 (nm)'
        if yaxis_name == 'auto detect':
            y_name = 'PC2 (nm)'
    # 处理 PCA 数据
        if smooth != 'true' and pca_color_by_replica == 0:
            if multi_files[0].endswith(".xvg"):
                for i, file in enumerate(multi_files):          
                    x_data, y_data, _ = self.read_data_xvg(file, "PC1", renumber)  # 假设 "PC1" 和 "PC2" 是合适的轴名称
                    labels = [x for x in range(len(y_data))]
                    
                    # 使用 define_trace 创建迹线
                    trace = self.define_trace(x_data, y_data, file, 'rainbow', flag=flag, labels=labels)  # 假设使用 'rainbow' 作为颜色
                    data.append(trace)
            elif multi_files[0].endswith(".csv"):
                for i, trace in enumerate(traces_name_list):
                    x_data, y_data, _ = self.read_data_csv(multi_files[0], "PC1", renumber)
                    labels = [x for x in range(len(y_data[i]))]
                    trace = self.define_trace(x_data, y_data[i], multi_files[0], 'rainbow', flag=flag, labels=labels)
                    data.append(trace)
            layout = self.setup_layout(plot_title, title_font, x_name, y_name, xy_font, xaxis_size, yaxis_size, font_color, legend_show, legend_font, font_family, grid_show, l, r, t ,b, x_low, x_high, y_low, y_high, flag=flag, axis_shows=axis_show, line_width=linewidth)

        elif smooth != 'true' and pca_color_by_replica != 0:
            number_replicas = pca_color_by_replica
            if multi_files[0].endswith(".xvg"):
                for i, file in enumerate(multi_files):          
                    x_data, y_data, _ = self.read_data_xvg(file, "PC1", renumber)  # 假设 "PC1" 和 "PC2" 是合适的轴名称
                    labels = ["replica-" + str(x+1) for x in range(number_replicas)]
                    n=len(y_data)
                    ranges = self.split_even_ranges(n, number_replicas)
                    for j, (a, bb) in enumerate(ranges):
                        labels = list(range(a, bb))
                        # legend 名称：文件名 + 分段编号
                        trace_name = f"{file} · replica-{j+1}"
                        trace = self.define_trace(
                            x_data[a:bb], y_data[a:bb], trace_name,
                            colors[j], flag=flag, labels=labels
                        )
                        data.append(trace)

        
        elif smooth == 'true':
            if multi_files[0].endswith(".xvg"):
                for i, file in enumerate(multi_files):          
                    x_data, y_data, _ = self.read_data_xvg(file, "PC1", renumber)  # 假设 "PC1" 和 "PC2" 是合适的轴名称
                    points = [list(pair) for pair in zip(x_data, y_data)]
                    density_matrix = self.pca_bins_density_define(nbin, points)
                    # 使用 define_trace 创建迹线
                    trace = self.define_trace(density_matrix, density_matrix, file, 'rainbow', flag=flag, smooth=smooth)  # 假设使用 'rainbow' 作为颜色
                    data.append(trace)
                    #     def plot_heatmap(self, x_points, y_points, plot_title, x_name, y_name, size, height, width, output_name):
            elif multi_files[0].endswith(".csv"):
                for i, trace in enumerate(traces_name_list):
                    x_data, y_data, _ = self.read_data_csv(multi_files[0], "PC1", renumber)
                    points = [list(pair) for pair in zip(x_data, y_data[0])]
                    density_matrix = self.pca_bins_density_define(nbin, points)
                    trace = self.define_trace(density_matrix, density_matrix, multi_files[0], 'rainbow', flag=flag, smooth=smooth)
                    data.append(trace)
            # 使用 setup_layout 设置布局
            layout = self.setup_layout(plot_title, title_font, x_name, y_name, xy_font, xaxis_size, yaxis_size, font_color, legend_show, legend_font, font_family, grid_show, l, r, t ,b, x_low, x_high, y_low, y_high, flag=flag, axis_shows=axis_show, line_width=linewidth)
        # 使用 rpy2 绘制图形
        # self.plot_heatmap(x_data, y_data, plot_title, x_name, y_name, size, output_name)
        # 使用 plot_graph 绘制图形
        self.plot_graph(data, layout, "Scatter_" + output_name)

    def plotly_free_energy(self, multi_files, output_name, plot_name, nbin, size, xaxis_name, yaxis_name, xaxis_size, yaxis_size, xy_font, title_font, legend_show, legend_font, font_family, font_color, grid_show, flag, uploaded_filenames, l,r,t,b, violin, smooth, error_bar, replica_number, axis_show, linewidth, transparency, x_low, x_high, y_low, y_high, trace_color):
        # Plotly = ['#636EFA', '#EF553B', '#00CC96', '#AB63FA', '#FFA15A', '#19D3F3', '#FF6692', '#B6E880', '#FF97FF', '#FECB52']
        Plotly = trace_color
        # AMPK_color = ['#222A2A', '#FB0D0D', '#2E91E5'] # black red blue
        # AMPK_color = ['#FB0D0D', '#2E91E5'] # black red blue
        # AMPK_color = ['#2E91E5'] # black red blue

        data, histogram_data, group_labels = [], [], []
        plot_title = 'Free Energy Surface'
        for i, file in enumerate(multi_files):
            x_data, y_data, z_data, df, index_of_free_energy, column_names = self.read_data_dat(file)
            # 如果有3列，则为phi psi 自由能
            if index_of_free_energy == 2:
                if 'phi' in column_names or 'psi' in column_names:
                    x_values = np.degrees(np.unique(x_data))
                    y_values = np.degrees(np.unique(y_data))
                else:
                    x_values = np.unique(x_data)
                    y_values = np.unique(y_data)                   
                x_grid, y_grid = np.meshgrid(x_values, y_values)
                z_data_array = np.array(z_data)
                free_energy_grid = z_data_array.reshape(len(y_values), len(x_values))
                # Plot the FES
                plt.contourf(x_grid, y_grid, free_energy_grid, levels=100)
                if xaxis_name == 'auto detect' and yaxis_name == 'auto detect':
                    xaxis_name = column_names[0]
                    yaxis_name = column_names[1]
                    plt.xlabel(xaxis_name)
                    plt.ylabel(yaxis_name)
                else:
                    plt.xlabel(xaxis_name)
                    plt.ylabel(yaxis_name)                
                plt.colorbar(label='Free energy / kJ mol^-1')
                plt.title('Free Energy Surface')
                plt.savefig("/tmp/" + str(i) + "_" + output_name)
                self.streamlit_download_file_plotly(str(i) + "_" + output_name, "/tmp/" + str(i) + "_" + output_name)
                # 3D plot
                fig=plt.figure(figsize=(24,14), dpi=600)
                ax = plt.axes(projection='3d')
                surf = ax.plot_surface(x_grid, y_grid, free_energy_grid, cmap = 'jet', rstride=1, cstride=1, alpha=None, linewidth=0, antialiased=True)
                # Set axes label
                ax.set_title('Free Energy Surface')
                ax.set_xlabel(xaxis_name, labelpad=5)
                ax.set_ylabel(yaxis_name, labelpad=5)
                ax.set_zlabel('Free energy / kJ mol^-1', labelpad=5)
                fig.colorbar(surf, shrink=0.7, aspect=15)
                plt.savefig("/tmp/" + str(i) + "_3D_" + output_name)
                self.streamlit_download_file_plotly(str(i) + "_3D_" + output_name, "/tmp/" + str(i) + "_3D_" + output_name)


            # 如果有2列，则为distance 自由能
            elif index_of_free_energy == 1:
                if xaxis_name == 'auto detect':
                    x_name = column_names[0]
                else:
                    x_name = xaxis_name
                if yaxis_name == 'auto detect':
                    y_name = 'Free energy / kJ mol<sup>-1</sup>'
                else:
                    y_name = yaxis_name
                # 使用 define_trace 创建迹线
                trace = self.define_trace(x_data, y_data, uploaded_filenames[i], Plotly[i % len(Plotly)], violine=violin)
                # trace = self.define_trace(x_data, y_data, uploaded_filenames[i], AMPK_color[i % len(AMPK_color)], violine=violin) # 曾经用于AMPK 论文的设置
                data.append(trace)
        if data != []:
           # 设置布局
            layout = self.setup_layout(plot_title, title_font, x_name, y_name, xy_font, xaxis_size, yaxis_size, font_color, legend_show, legend_font, font_family, grid_show, l, r, t , b, x_low, x_high, y_low, y_high, violine=violin, axis_shows=axis_show, line_width=linewidth)
            # 绘制图形
            self.plot_graph(data, layout, output_name)

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

    def GROreader_not_good(self,gro): 
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
    def GROreader(self, gro): 
        lines = gro.splitlines()  # for streamlit 如果 'gro' 是一个二进制文件，使用 gro.decode().splitlines() 
    
        self.head = lines[0].strip()
        self.total_atom = int(lines[1])
        self.last = lines[-1].strip()
        
        # 忽略前两行和最后一行
        lines = lines[2:-1]
        
        # 初始化属性列表
        self.resid = []
        self.resname = []
        self.atomname = []
        self.index = []
        self.x = []
        self.y = []
        self.z = []
        self.xyz = []
        
        # 逐行解析内容
        for line in lines:
            resid = int(line[0:5].strip())
            resname = line[5:10].strip()
            atomname = line[10:15].strip()
            atomindex = int(line[15:20].strip())
            x = float(line[20:28].strip())
            y = float(line[28:36].strip())
            z = float(line[36:44].strip())
            
            self.resid.append(resid)
            self.resname.append(resname)
            self.atomname.append(atomname)
            self.index.append(atomindex)
            self.x.append(x)
            self.y.append(y)
            self.z.append(z)
            self.xyz.append([x, y, z])

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
        st.text(sentence)
        # print(self.metal1,self.metal2,self.metal3)
            
        
        
        
    def coordinator(self, num_neighbours, distance_value, atom_list, metal_list, residue_list):
        # find the atom index
        # st.text(atom_list)
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
                        topol_lines.append(f"{B}                1\n")
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
    pdb     = []
    pdb_file_name = []
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
    def __init__(self, pdb, pdb_file_name, pepname):
        self.pdb = pdb
        self.pdb_file_name = pdb_file_name
        self.pepname = pepname
        for i in range(len(pdb)):
            self.converter(pdb[i], pepname, self.chain, self.resnum)
            self.PDBwriter("ligand-" + pdb_file_name[i])
        zip_buffer = self.create_zip_file(pdb_file_name)
        self.streamlit_download_zip_file(zip_buffer)    
        
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
        f = open("/tmp/" + filename, "w")                                                             # e.g: f = linesplit[0]+"_PO3.pdb"
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

        self.streamlit_download_file(filename, "/tmp/" + filename)
        
    def create_zip_file(self, ligand_files):
        # 创建一个BytesIO对象来保存zip文件内容
        zip_buffer = io.BytesIO()
        
        with zipfile.ZipFile(zip_buffer, "w", zipfile.ZIP_DEFLATED) as zip_file:
            for ligand_file in ligand_files:
                # 将每个ligand文件写入zip包中
                zip_file.write("/tmp/ligand-" + ligand_file, os.path.basename("/tmp/ligand-" + ligand_file))
        
        zip_buffer.seek(0)  # 重置指针到文件的开头
        return zip_buffer        

        
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
    def streamlit_download_zip_file(self, zip_buffer):
        st.download_button(
            label="Download All Ligands as ZIP",
            data=zip_buffer,
            file_name="ligands.zip",
            mime="application/zip"
        )
##########################################################################################################################################################################################
class gmx_dssp():
    def __init__(self, ds_data, ds_traj, ds_output_name, ds_original, ds_color, original_map, original_colorbar, original_colorscale, simple_map, simple_colorbar, simple_colorscale):
        self.plot_figure(ds_data, ds_traj, ds_output_name, ds_original, ds_color, original_map, original_colorbar, original_colorscale, simple_map, simple_colorbar, simple_colorscale)
    
    # def read_time_and_residue(self, traj):
        # 载入PDB文件
      #  u = mda.Universe(traj)
        
        # 提取时间戳和残基列表 # 2024-07-05
      #  times = [ts.time for ts in u.trajectory]  # 假设时间单位是ns  # 2024-07-05
      #  residues = [res.resname + str(res.resid) for res in u.residues] # 2024-07-05
      #  return times, residues # 2024-07-05
    def read_time_and_residue(self, traj):
        # 载入PDB文件
        u = mda.Universe(traj)
        
        # 提取时间戳
        times = [ts.time for ts in u.trajectory]  # 假设时间单位是ns 
        
        # 提取残基列表并自动添加链ID
        residues = []
        previous_resid = None
        chain_id = 'A'
        
        for res in u.residues:
            if previous_resid is not None and res.resid < previous_resid:
                # Residue ID decreased, indicating a new chain
                chain_id = chr(ord(chain_id) + 1)
            residues.append(res.resname + str(res.resid) + "_" + chain_id)
            previous_resid = res.resid
        
        return times, residues
    
    def detect_break(self, first_line, residue_list):
        # Find the positions of '=' in the first line
        equal_positions = [pos for pos, char in enumerate(first_line) if char == "="]
        print("Here is the break points ", equal_positions)
        # Sort the positions in reverse order
        equal_positions.sort(reverse=True)
        # Insert 'break' into the residue_list at the corresponding positions from the end
        for pos in equal_positions:
            if pos < len(residue_list):  # Only insert if the position is within the bounds of the residue_list
                residue_list.insert(pos, 'break')
        return residue_list
        
    def read_data(self, data, traj, original, unique_color_bar): 
        times, residues = self.read_time_and_residue(traj)
        # 读取DSSP数据（dssp.dat）
        with open(data, "r") as file:
            dssp_lines = file.readlines()
            first_line = dssp_lines[0]
            residues = self.detect_break(first_line, residues)
    
        # 转换DSSP数据为列表
        dssp_data = [list(line.strip()) for line in dssp_lines]
        
        # 创建DataFrame
        df = pd.DataFrame(data=dssp_data, index=times, columns=residues)
        
        if original == 'false' and unique_color_bar == 'true':
            # times = times.append(times[-1]+1) # for scale the color
            # Calculate the length of each third
            third_length = len(df.columns) // 3
            
            # Create a new row
            new_row = [1] * third_length + [2] * third_length + [3] * (len(df.columns) - 2 * third_length)
            
            # Append the new row at the bottom of the DataFrame
            df.loc[len(df)] = new_row
        elif original != 'false' and unique_color_bar == 'true':
            tenth_length = len(df.columns) // 10
            # Initialize an empty list for the new row
            new_row = []
            # Create new row by appending numbers from 1 to 10, each repeated 'tenth_length' times
            # Note that the last segment fills the remainder of the row if it's not evenly divisible by 10
            for i in range(10):
                if i < 9:
                    new_row += [i] * tenth_length
                else:
                    # The last segment includes any extra columns
                    new_row += [i] * (len(df.columns) - len(new_row))

            # Append the new row at the bottom of the DataFrame
            df.loc[len(df)] = new_row
            pass
        
        return df

    def replace_letters(self, original):
        # 定义转换字典
        if original != 'false':
            structure_values = {
                'H': 9,   # 'alpha-helix',
                'B': 8,   # 'residue in isolated beta-bridge',
                'E': 7,   # 'extended strand that participates in beta-ladder',
                'G': 6,   # '3_10-helix',
                'I': 5,   # 'pi-helix',
                'P': 4,   # 'kappa-helix',  # Assuming kappa-helix is synonymous 
                'S': 3,   # 'bend',
                'T': 2,   # 'hydrogen-bonded turn',
                '=': 1,   # 'break',
                '~': 0    # 'loop'  # No special secondary structure designation
            }
        else:
            structure_values = {
                'H': 2,   # 'alpha-helix',
                'B': 1,   # 'residue in isolated beta-bridge',
                'E': 3,   # 'beta-sheet',
                'G': 1,   # '3_10-helix',
                'I': 1,   # 'pi-helix',
                'P': 1,   # 'kappa-helix',  # Assuming kappa-helix is synonymous 
                'S': 1,   # 'beta-bend',
                'T': 1,   # 'beta-turn',
                '=': 1,   # 'break',
                '~': 1    # 'loop'  # No special secondary structure designation
            }
        return structure_values
        
    def calculate_ss_percentage(self, df, structure_values, unique_color):
        if unique_color == "true":
            # 删除最后一行
            df = df.drop(df.index[-1])
        else:
            pass
        inverse_structure_values = {v: k for k, v in structure_values.items()}
        ss_percentage = pd.DataFrame(index=df.index)
        for ss_type in set(structure_values.values()):
            ss_char = inverse_structure_values[ss_type]
            ss_percentage[ss_char] = (df == ss_type).sum(axis=1) / df.shape[1] * 100
        return ss_percentage

    def plot_ss_percentage(self, df, structure_values, original, original_colorbar, simple_colorbar, outputname, unique_color):
        # 计算二级结构百分比
        ss_percentage = self.calculate_ss_percentage(df, structure_values, unique_color)

        # 定义映射字典
        if original != 'false':
            legend_mapping = {v: k for k, v in original_colorbar.items()}
        else:
            legend_mapping = {v: k for k, v in simple_colorbar.items()}

        # Plotly颜色列表
        plotly_colors = ['#636EFA', '#EF553B', '#00CC96', '#AB63FA', '#FFA15A', '#19D3F3', '#FF6692', '#B6E880', '#FF97FF', '#FECB52']

        # 绘制百分比图
        fig_ss_percentage = go.Figure()
        for idx, ss_type in enumerate(ss_percentage.columns):
            ss_type_name = legend_mapping[structure_values[ss_type]]
            fig_ss_percentage.add_trace(go.Scatter(
                x=ss_percentage.index,
                y=ss_percentage[ss_type],
                mode='lines',
                name=ss_type_name,
                line=dict(color=plotly_colors[idx % len(plotly_colors)])
            ))

        fig_ss_percentage.update_layout(
            title='DSSP Percentage',
            xaxis_title='Time (ns)',
            yaxis_title='Percentage (%)',
            width=800,
            height=600
        )

        # 显示百分比图
        pio.write_image(fig_ss_percentage, "/tmp/percentage-" + outputname)
        self.streamlit_download_file_plotly("percentage-" + outputname, "/tmp/percentage-" + outputname)        
    
    def plot_figure(self, data, traj, outputname, original, unique_color, original_map, original_colorbar, original_colorscale, simple_map, simple_colorbar, simple_colorscale):
        df = self.read_data(data, traj, original,unique_color)   
        # 使用字典转换DataFrame中的值
        if original == 'false':
            structure_values = simple_map
        else:
            structure_values = original_map
        df.replace(structure_values, inplace=True)

        # 绘制percentage 图
        self.plot_ss_percentage(df, structure_values, original, original_colorbar, simple_colorbar, outputname, unique_color)
        
        # Define color scale and color bar settings
        simple_colorscale = simple_colorscale
        original_colorscale = original_colorscale

        # colorscale = [[0.00, "red"],   [0.33, "red"], [0.33, "green"], [0.66, "green"], [0.66, "blue"],  [1.00, "blue"]]
        # 将颜色条的刻度设置为固定的值
        # 从字典中提取键和值
        simple_colorbar_ticktext = list(simple_colorbar.keys())   # 提取键，作为标签文本
        simple_colorbar_ticks = list(simple_colorbar.values())    # 提取值，作为颜色条的刻度
        colorbar_ticks = list(original_colorbar.values())
        colorbar_ticktext = list(original_colorbar.keys())
        if original != 'false':
            # 创建热图
            fig = go.Figure(data=go.Heatmap(
                z=df.T.values,
                x=df.index,
                y=df.columns,
                colorscale=original_colorscale,
                colorbar=dict(tickmode = 'array', tickvals=colorbar_ticks, ticktext=colorbar_ticktext),
                hoverongaps=False))
        elif original == 'false':
            colorbar_ticktext = simple_colorbar_ticktext
            colorbar_ticks = simple_colorbar_ticks
            # 创建热图
            fig = go.Figure(data=go.Heatmap(
                z=df.T.values,
                x=df.index,
                y=df.columns,
                colorscale=simple_colorscale,
                colorbar=dict(tickmode = 'array', tickvals=colorbar_ticks, ticktext=colorbar_ticktext),
                hoverongaps=False))
           
        fig.update_layout(
            title='Secondary Structure Analysis Over Time',
            yaxis_title='Residue',
            xaxis_title='Time (ns)',
            width=800,
            height=600)

        # 显示热图
        pio.write_image(fig, "/tmp/" + outputname)

        self.streamlit_download_file_plotly(outputname, "/tmp/" + outputname)
    
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
##########################################################################################################################################################################################
class renumber_MODEL():
    def __init__(self, files, name):
        # 初始化计数器
        count = 1
        lines = 1
        # 打开输入文件和输出文件
        with open(files, 'r') as file, open(f'/tmp/renumbered.pdb', 'w') as output:
            # 遍历文件中的每一行
            for line in file:
                # 使用正则表达式检查行是否包含"MODEL"和后面的数字
                if re.match(r"^\s*MODEL\s+\d+", line):
                    # 替换匹配的行为"MODEL"后面接计数器的值，并将计数器加一
                    # new_line = re.sub(r"^\s*MODEL\s+\d+", f"MODEL        {count}", line)
                    new_line = re.sub(r"^\s*MODEL\s+\d+", f"MODEL{count:9d}", line)
                    count += 1
                    output.write(new_line)
                    lines += 1
                else:
                    output.write(line)
                    lines += 1
        #with open('/tmp/' + name +'_renumbered.pdb', 'r') as file:
        #    lines = file.readlines()
        #    st.write(f"before download, the file includes: {len(lines)}")
            
        self.streamlit_download_file(str(name[:-4]) + "_renumbered.pdb", "/tmp/renumbered.pdb")            
        
    def streamlit_download_file(self, download_name, content_file):
        # Download topol.top file  #      
        # 打开 content_file 文件并读取其内容
        with open(content_file, 'r') as top_file:
            content = top_file.read()   

        #with open(content_file, 'r') as file:
        #    lines = file.readlines()
        #    st.write(f"after download, the file includes: {len(lines)}")

            
        # 添加一个下载按钮，传递 receptor_top_content 作为文件内容
        st.download_button(
            label = "Download " +  download_name,
            data = content,
            key = download_name,
            file_name = download_name
            )
##########################################################################################################################################################################################
class ff_res_adder():
    def __init__(self, lig_itp, lig_gro, ff_rtp, ff_hdb, ff_bonded, ff_nonbonded, ff_aomtypes, ff_restype, atom_types, isamino_acid, output_name):
        self.add_res_to_ff(lig_itp, lig_gro, ff_rtp, ff_hdb, ff_bonded, ff_nonbonded, ff_aomtypes, ff_restype, atom_types, isamino_acid, output_name)


    def add_res_to_ff(self, lig_itp, lig_gro, ff_rtp, ff_hdb, ff_bonded, ff_nonbonded, ff_aomtypes, ff_restype, atom_types, isamino_acid, output_name):
        # parse the lig_GMX.itp file
        itp_content_dicts   =   self.parse_itp(lig_itp, atom_types)
        # parse the lig_GMX.gro file
        if lig_gro != 0:
            residue_numbers, residue_names, atom_names, atom_numbers, positions = self.parse_gro(lig_gro)
        else:
            pass
        # replace the digit number to atom_types for bonds, angles, dihedrals, impropers
        bond_records_dict   =   self.generate_bond_records(itp_content_dicts)
        # the notes for residuetypes.dat
        ff_restype_dict     =   self.parse_residuetypes(itp_content_dicts, output_name)
        # the content to be added to aminoacids.rtp
        ff_rtp_dict         =   self.parse_rtp(itp_content_dicts, isamino_acid, output_name)
        # the content to be added to aminoacids.hdb          
        ff_hdb_dict         =   self.parse_hdb(itp_content_dicts, output_name) 
        # parse the ffbonded.itp file
        ff_bonded_dict      =   self.parse_bonded(ff_bonded)
        # find the missed bonds records
        missing_bonds_dict  =   self.find_missing_bonds(bond_records_dict, ff_bonded_dict, output_name)
        # 打印嵌套字典
        # self.print_nested_dict(bond_records_dict) 
        # self.print_nested_dict(ff_bonded_dict) 
        # Get the list of atom types
        # ff_nonbonded_dict   =   self.parse_nonbonded(ff_nonbonded)  
        if ff_aomtypes != 0 and ff_nonbonded != 0:
            ff_aomtypes_dict    =   self.parse_atomtypes(ff_aomtypes, ff_nonbonded) 
        # self.print_nested_dict(ff_aomtypes_dict)  
        # Compare ff_atomtype and user provided atom_types
        try:
            missing_atom_type   =   self.find_missing_atomtypes(ff_aomtypes_dict, atom_types, output_name)
        except:
            pass
        # bond_values_dict    =   self.calculate_bond_values(positions, itp_content_dicts)
        self.streamlit_download_file(output_name, "/tmp/" + output_name)
    # def parse_itp(self, lig_itp, atom_types):
    def parse_itp(self, lig_itp, atom_types):
        header_atomtypes = ['name', 'bond_type', 'mass', 'charge', 'ptype', 'sigma', 'epsilon', 'Amb_sigma', 'Amb_epsilon']      
        header_moleculetype = ['name', 'nrexcl']        
        header_atoms = ['atom_number', 'atom_type', 'resid', 'resname', 'atom_name', 'cgnr', 'partial_charge', 'mass', 'qtot', 'total_charge']        
        header_bonds = ['a1', 'a2', 'funct', 'bond_length', 'bond_strength', 'a1_name', 'a2_name']        
        header_pairs = ['a1', 'a2', 'funct', 'a1_name', 'a2_name']        
        header_angles = ['a1', 'a2', 'a3', 'funct', 'angle_degree', 'angle_strength', 'a1_name', 'a2_name', 'a3_name']        
        header_dihedrals = ['a1', 'a2', 'a3', 'a4', 'funct', 'C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'a1_name', 'a2_name', 'a3_name', 'a4_name']        
        header_impropers = ['a1', 'a2', 'a3', 'a4', 'funct', 'angle_degree', 'angle_strength', 'multiplicity', 'a1_name', 'a2_name', 'a3_name', 'a4_name']    
        header_list = [header_atomtypes, header_moleculetype,  header_atoms ,  header_bonds ,  header_pairs ,  header_angles ,  header_dihedrals ,  header_impropers ]
    
        # Define dictionaries to store different sections
        sections = {
            'atomtypes': {},
            'moleculetype': {},
            'atoms': {},
            'bonds': {},
            'pairs': {},
            'angles': {},
            'dihedrals': {},
            'impropers': {}
        }
    
        # Open and read the file
        with open(lig_itp, 'r') as file:
            lines = file.readlines()

        # Helper variables to track sections and headers
        current_section = None
        headers = []
        data_rows = []
        count = 0

        for line in lines:
            line = line.strip()
            if line.startswith('['):  # Detect new section
                if current_section and headers:
                    # save the datas to target section
                    transposed_data = list(zip(*data_rows))
                    sections[current_section] = {key: tuple(values) for key, values in zip(headers, transposed_data)}          
                # Extract section name from format [ section_name ] and reset the previous section
                section_name = line.split()[1]
                current_section = section_name
                # headers = []  # Reset headers for the new section
                headers = header_list[count]
                count += 1 
                data_rows = []
            # define whether it's improper dihedral part
            elif current_section == 'dihedrals' and line.startswith(';') and 'propers' in line:
                current_section = 'impropers'
            # read contents
            elif line and not line.startswith(';'):  # Valid data lines
                if current_section and headers:  # We already have headers, parse data
                    a = line.split()
                    # remove the values = ; or -
                    a = [item for item in a if item not in [';', '-']]
                    # remove the "-" for dihedrals or improper dihedrals
                    if current_section in ['dihedrals', 'impropers']:
                        a = a[:-4] + [item.replace('-', '') for item in a[-4:]]
                    # a = [item for item in a if item != ';']
                    data_rows.append(a)
        # 最后一个节的数据保存
        if current_section and headers:
            transposed_data = list(zip(*data_rows))
            sections[current_section] = {key: tuple(values) for key, values in zip(headers, transposed_data)}
           
       
        # change the atom_type to user's input
        if len(atom_types) >= len(sections['atoms']['atom_type']):
            sections['atoms']['atom_type'] = tuple(atom_types)
        else:
            st.text("Error: atom_types list does not have enough elements to replace all atom_types.")

        return sections
    
    def parse_gro(self, lig_gro):
        
        
        # 定义空列表以存储每列的数据
        residue_numbers = []
        residue_names = []
        atom_names = []
        atom_numbers = []
        positions    = []
        
        # 读取文件并提取所需行
        with open(lig_gro, 'r') as file:
            lines = file.readlines()
            # 获取第3行到倒数第二行的数据行
            data_lines = lines[2:-1]
        
        # 解析每行数据
        for line in data_lines:
            # 假设每个元素之间由多个空格分隔，这是一个典型的固定列宽格式
            parts = line.split()
            if len(parts) >= 7:  # 确保每行数据都是完整的
                residue_numbers.append(parts[0])
                residue_names.append(parts[1])
                atom_names.append(parts[2])
                atom_numbers.append(parts[3])
                positions.append([float(parts[4]), float(parts[5]), float(parts[6])])


        return residue_numbers, residue_names, atom_names, atom_numbers, positions 
    def parse_rtp(self, itp_content_dict, isamino_acids, output_name):  
        with open("/tmp/" + output_name, "a") as file:
            file.write("######################## Please add below to aminoacids.rtp! ########################\n")
            file.write(f"[ {itp_content_dict['atoms']['resname'][0]} ]\n")
            file.write(" [ atoms ]\n")
            for atomname, atomtype, charge, indices in zip(itp_content_dict['atoms']['atom_name'], itp_content_dict['atoms']['atom_type'], itp_content_dict['atoms']['partial_charge'], itp_content_dict['atoms']['atom_number']):
                file.write("%6s    %-12s%8.5f%6d\n" % (atomname, atomtype, float(charge), int(indices)))
            file.write(" [ bonds ]\n")
            for atom1, atom2 in zip(itp_content_dict['bonds']['a1_name'], itp_content_dict['bonds']['a2_name']):
                file.write("%6s%6s\n" % (atom1,atom2))
            if isamino_acids == "true":
                file.write("%6s%6s\n" % ("-C","N"))
            file.write(" [ impropers ]\n")
            if isamino_acids == "true":
                file.write("%6s%6s%6s%6s\n" % ("-C","CA","N","H"))
                file.write("%6s%6s%6s%6s\n" % ("CA","+N","C", "O"))
            for atom1,atom2,atom3,atom4 in zip(itp_content_dict['impropers']['a1_name'], itp_content_dict['impropers']['a2_name'], itp_content_dict['impropers']['a3_name'], itp_content_dict['impropers']['a4_name']):
                file.write("%6s%6s%6s%6s\n" % (atom1,atom2,atom3,atom4))
        

         
    def parse_hdb(self, itp_content_dict, output_name):  
        hdb_bond_types = """
1. one planar hydrogen, e.g. rings or peptide bond
2. one single hydrogen, e.g. hydroxyl
3. two planar hydrogens, e.g. ethylene -C=CH2, or amide -C(=O)NH2
4. two or three tetrahedral hydrogens, e.g. -CH3
5. one tetrahedral hydrogen, e.g. C3* CH*
6. two tetrahedral hydrogens, e.g. C-CH2 *-C*
"""
        def find_connections(atom, bonds):
            connections = {'H': [], 'non_H': []}
            for a1, a2 in zip(bonds['a1_name'], bonds['a2_name']):
                if a1 == atom:
                    (connections['H'] if 'H' in a2 else connections['non_H']).append(a2)
                if a2 == atom:
                    (connections['H'] if 'H' in a1 else connections['non_H']).append(a1)
            return connections
        def define_bond_type(atoms, counts):
            bond_type = "5"
            if atoms == "N":
                bond_type = str(1)
            elif atoms == "CA" and counts == 2:
                bond_type = str(6)
            elif atoms == "CA" and counts == 1:
                bond_type = str(5)
            elif atoms == "CB" and counts == 2:
                bond_type = str(6)
            elif atoms == "CB" and counts == 1:
                bond_type = str(5)
            elif counts == 1:
                bond_type = str(5)
            elif counts == 2:
                bond_type = str(6)
            elif counts == 3:
                bond_type = str(4)
            return bond_type
        # Main processing
        output = []
        for atom in itp_content_dict['atoms']['atom_name']:
            if 'H' not in atom:
                conns = find_connections(atom, itp_content_dict['bonds'])
                if conns['H']:
                    # Only consider atoms directly bonded to H
                    count_hydrogens = len(conns['H'])
                    bond_type = define_bond_type(atom, count_hydrogens)
                    hydrogens = "H" + str(atom[1:] if len(atom) > 1 else "")
                    others = '       '.join(conns['non_H'])
                    # output.append(f"{count_hydrogens}\t6\t{hydrogens}\t{atom}\t{', '.join(conns['non_H'])}")
                    output.append(f"{count_hydrogens}\t{bond_type}\t{hydrogens}\t{atom}\t{others}")
        # Print output as per the format discussed
        with open("/tmp/" + output_name, "a") as file:
            file.write("######################## Please add below to aminoacids.hdb! ########################\n")
            # file.write("CR1\t" + str(len(output)) + "\n")
            file.write(f"{itp_content_dict['atoms']['resname'][0]}\t" + str(len(output)) + "\n")
            for line in output:
                file.write(line)
                file.write("\n")
            file.write("##Rules of the Hydrogen adding!##")        
            file.write(hdb_bond_types)
        
        
       
               
    def parse_residuetypes(self, itp_content_dict, output_name):  
        with open("/tmp/" + output_name, "w") as file:
            file.write("######################## Please add below to residuetypes.dat! ########################\n")
            file.write("%-8s%-7s\n" % (itp_content_dict['atoms']['resname'][0],"Protein"))
        

    def generate_bond_records(self, itp_content_dict):
        # 建立atom_number到atom_type的映射
        atom_type_mapping = {num: typ for num, typ in zip(itp_content_dict['atoms']['atom_number'], itp_content_dict['atoms']['atom_type'])}
        
        # 替换bonds字典中的a1和a2为对应的atom_type
        updated_a1 = [atom_type_mapping[num] for num in itp_content_dict['bonds']['a1']]
        updated_a2 = [atom_type_mapping[num] for num in itp_content_dict['bonds']['a2']]
        
        # 更新bonds字典
        itp_content_dict['bonds']['a1'] = tuple(updated_a1)
        itp_content_dict['bonds']['a2'] = tuple(updated_a2)

        # 替换angles字典中的a1, a2, a3
        updated_a1 = [atom_type_mapping[num] for num in itp_content_dict['angles']['a1']]
        updated_a2 = [atom_type_mapping[num] for num in itp_content_dict['angles']['a2']]
        updated_a3 = [atom_type_mapping[num] for num in itp_content_dict['angles']['a3']]

        # 更新angles字典
        itp_content_dict['angles']['a1'] = tuple(updated_a1)
        itp_content_dict['angles']['a2'] = tuple(updated_a2)
        itp_content_dict['angles']['a3'] = tuple(updated_a3)

        # 替换dihedrals字典中的a1, a2, a3, a4
        updated_a1 = [atom_type_mapping[num] for num in itp_content_dict['dihedrals']['a1']]
        updated_a2 = [atom_type_mapping[num] for num in itp_content_dict['dihedrals']['a2']]
        updated_a3 = [atom_type_mapping[num] for num in itp_content_dict['dihedrals']['a3']]
        updated_a4 = [atom_type_mapping[num] for num in itp_content_dict['dihedrals']['a4']]

        # 更新dihedrals字典
        itp_content_dict['dihedrals']['a1'] = tuple(updated_a1)
        itp_content_dict['dihedrals']['a2'] = tuple(updated_a2)
        itp_content_dict['dihedrals']['a3'] = tuple(updated_a3)
        itp_content_dict['dihedrals']['a4'] = tuple(updated_a4)

        # 替换impropers字典中的a1, a2, a3, a4
        updated_a1 = [atom_type_mapping[num] for num in itp_content_dict['impropers']['a1']]
        updated_a2 = [atom_type_mapping[num] for num in itp_content_dict['impropers']['a2']]
        updated_a3 = [atom_type_mapping[num] for num in itp_content_dict['impropers']['a3']]
        updated_a4 = [atom_type_mapping[num] for num in itp_content_dict['impropers']['a4']]

        # 更新impropers字典
        itp_content_dict['impropers']['a1'] = tuple(updated_a1)
        itp_content_dict['impropers']['a2'] = tuple(updated_a2)
        itp_content_dict['impropers']['a3'] = tuple(updated_a3)
        itp_content_dict['impropers']['a4'] = tuple(updated_a4)

        return itp_content_dict

    def parse_bonded(self, ff_bonded): 
        header_bondtypes = ['a1', 'a2', 'funct', 'bond_length', 'bond_strength']      
        header_constrainttypes = ['a1', 'a2', 'funct', 'bond_length']        
        header_angletypes = ['a1', 'a2', 'a3', 'funct', 'angle_degree', 'angle_strength']        
        header_dihedral_f1 = ['a1', 'a2', 'a3', 'a4', 'funct', 'angle_degree', 'angle_strength', 'multiplicity']        # acpype for non-periodic improper dihedrals, to fix the dihedral
        header_dihedral_f3 = ['a1', 'a2', 'a3', 'a4', 'funct', 'c0', 'c1', 'c2', 'c3', 'c4', 'c5']                      # acpype for proper dihedrals
        header_dihedral_f4 = ['a1', 'a2', 'a3', 'a4', 'funct', 'angle_degree', 'angle_strength', 'multiplicity']        # Gromacs for periodic improper dihedral, have some flexibility to spin
        header_dihedral_f9 = ['a1', 'a2', 'a3', 'a4', 'funct', 'angle_degree', 'angle_strength', 'multiplicity']        # Gromacs for proper dihedral (multiple), multiple is the number of energy peaks or valleys
        header_list = [header_bondtypes, header_constrainttypes,  header_angletypes ,  header_dihedral_f1 ,  header_dihedral_f3 ,  header_dihedral_f4 ,  header_dihedral_f9]
    
        # Define dictionaries to store different sections
        sections = {
            'bondtypes': {},
            'constrainttypes': {},
            'angletypes': {},
            'dihedral_f1': {},
            'dihedral_f3': {},
            'dihedral_f4': {},
            'dihedral_f9': {}
        }
    
        # Open and read the file
        with open(ff_bonded, 'r') as file:
            lines = file.readlines()

        # Helper variables to track sections and headers
        current_section = None
        headers = []
        data_rows = []
        count = 0

        for line in lines:
            line = line.strip()
            if line.startswith('['):  # Detect new section
                # step-5. write the data to dictionary!
                if current_section and headers:
                    # save the datas to target section
                    transposed_data = list(zip(*data_rows))
                    sections[current_section] = {key: tuple(values) for key, values in zip(headers, transposed_data)}          
                # step-1. Extract section name from format [ section_name ] and reset the previous section
                section_name = line.split()[1]
                current_section = section_name
                # step-2. Assign headers for the columns
                if current_section == "bondtypes":
                    count = 0
                elif current_section == "constrainttypes":
                    count = 1
                elif current_section == "angletypes":
                    count = 2
                elif current_section == "dihedraltypes":
                    count = 4
                elif current_section == "dihedraltypes":
                    count = 4
                elif current_section == "dihedraltypes":
                    count = 4
                elif current_section == "dihedraltypes":
                    count = 4                    
                headers = header_list[count]

                data_rows = []
            # step-3. Define the type of dihedrals, also write the 1st row to data_rows
            elif current_section == 'dihedraltypes' and not line.startswith((';', '[', '#')):
                a = line.split()
                # remove the values = ; or -
                a = [item for item in a if item not in [';', '-']]
                data_rows.append(a)                
                if str(line.split()[4]) == '1':
                    current_section = 'dihedral_f1'
                    count = 3
                if str(line.split()[4]) == '3':
                    current_section = 'dihedral_f3'
                    count = 4
                if str(line.split()[4]) == '4':
                    current_section = 'dihedral_f4'
                    count = 5
                if str(line.split()[4]) == '9':
                    current_section = 'dihedral_f9' 
                    count = 6
                headers = header_list[count]
            # step-4. Read contents
            elif line and not line.startswith((';', '[', '#')):  # Valid data lines
                if current_section and headers:  # We already have headers, parse data
                    a = line.split()
                    # remove the values = ; or -
                    a = [item for item in a if item not in [';', '-']]
                    data_rows.append(a)
        # step-6. 最后一个节的数据保存
        if current_section and headers:
            transposed_data = list(zip(*data_rows))
            sections[current_section] = {key: tuple(values) for key, values in zip(headers, transposed_data)}

        return sections
    
    def print_nested_dict(self, dictionary, indent=0):
        for key, value in dictionary.items():
            if isinstance(value, dict):
                print(f"{' ' * indent}{key}:")
                self.print_nested_dict(value, indent + 4)
            else:
                print(f"{' ' * indent}{key}: {value}")
       
    # Used in find_missing_bonds method
    def matches_with_x(self, ff_entry, lig_entry):
        # ff_entry 和 lig_entry 都是包含四个元素的列表，例如 dihedral ['C', 'N', 'CA', 'X']
        for ff_atom, lig_atom in zip(ff_entry, lig_entry):
            if ff_atom != 'X' and ff_atom != lig_atom:
                return False
        return True
    
    def find_missing_bonds(self, lig_dict, ffbonded_dict, output_name):
        # 提取 ligand 的键合信息
        lig_bonds_list = [list(pair) for pair in zip(lig_dict['bonds']['a1'], lig_dict['bonds']['a2'])]  
        # lig_angles
        lig_angles_list = [list(pair) for pair in zip(lig_dict['angles']['a1'], lig_dict['angles']['a2'], lig_dict['angles']['a3'])]  
        # lig_dihedrals
        lig_dihedrals_list = [list(pair) for pair in zip(lig_dict['dihedrals']['a1'], lig_dict['dihedrals']['a2'], lig_dict['dihedrals']['a3'], lig_dict['dihedrals']['a4'])]   
        # lig_impropers
        lig_impropers_list = [list(pair) for pair in zip(lig_dict['impropers']['a1'], lig_dict['impropers']['a2'], lig_dict['impropers']['a3'], lig_dict['impropers']['a4'])]          
        # 提取 forcefield 中已有的键合信息
        ff_bonds_list = [list(pair) for pair in zip(ffbonded_dict['bondtypes']['a1'], ffbonded_dict['bondtypes']['a2'])]  
        # ff_angles
        ff_angles_list = [list(pair) for pair in zip(ffbonded_dict['angletypes']['a1'], ffbonded_dict['angletypes']['a2'], ffbonded_dict['angletypes']['a3'])]  
        # ff_dihedrals
        ff_dihedrals_list = []
        for dihedral_dict in [ffbonded_dict['dihedral_f3'], ffbonded_dict['dihedral_f9']]:
            if dihedral_dict:
                dihedrals_list = [list(pair) for pair in zip(dihedral_dict['a1'], dihedral_dict['a2'], dihedral_dict['a3'], dihedral_dict['a4'])]
                ff_dihedrals_list.extend(dihedrals_list)
        ff_impropers_list = []
        for dihedral_dict in [ffbonded_dict['dihedral_f1'], ffbonded_dict['dihedral_f4']]:
            if dihedral_dict:
                dihedrals_list = [list(pair) for pair in zip(dihedral_dict['a1'], dihedral_dict['a2'], dihedral_dict['a3'], dihedral_dict['a4'])]
                ff_impropers_list.extend(dihedrals_list)
        # 找出 forcefield 中缺少的键合信息
        missing_bonds = [bond for bond in lig_bonds_list if bond not in ff_bonds_list and bond[::-1] not in ff_bonds_list]
        missing_angles = [angle for angle in lig_angles_list if angle not in ff_angles_list and angle[::-1] not in ff_angles_list]
        # missing_dihedrals = [dihedral for dihedral in lig_dihedrals_list if dihedral not in ff_dihedrals_list and dihedral[::-1] not in ff_dihedrals_list]
        missing_dihedrals = []
        for lig_dihedral in lig_dihedrals_list:
            if not any(self.matches_with_x(ff_dihedral, lig_dihedral) or self.matches_with_x(ff_dihedral[::-1], lig_dihedral) for ff_dihedral in ff_dihedrals_list):
                missing_dihedrals.append(lig_dihedral)
        # missing_impropers = [dihedral for dihedral in lig_impropers_list if dihedral not in ff_impropers_list and dihedral[::-1] not in ff_impropers_list]
        missing_impropers = []
        for lig_improper in lig_impropers_list:
            if not any(self.matches_with_x(ff_dihedral, lig_improper) or self.matches_with_x(ff_dihedral[::-1], lig_improper) for ff_dihedral in ff_impropers_list):
                missing_impropers.append(lig_improper)
        # 删除重复
        missing_bonds_dedu = []
        missing_angles_dedu =[]
        missing_dihedrals_dedu = []
        missing_impropers_dedu = []
        for bond in missing_bonds:
            if bond not in missing_bonds_dedu:
                missing_bonds_dedu.append(bond)
        # [missing_bonds_dedu.append(bonds) for bonds in missing_bonds if bonds not in missing_bonds_dedu]
        # 如果missing_angles中的元素不在missing_angles_dedu中，则添加到missing_angles_dedu中
        for angle in missing_angles:
            if angle not in missing_angles_dedu:
                missing_angles_dedu.append(angle)
        # 如果missing_dihedrals中的元素不在missing_dihedrals_dedu中，则添加到missing_dihedrals_dedu中
        for dihedral in missing_dihedrals:
            if dihedral not in missing_dihedrals_dedu:
                missing_dihedrals_dedu.append(dihedral)
        # 如果missing_impropers中的元素不在missing_impropers_dedu中，则添加到missing_impropers_dedu中
        for improper in missing_impropers:
            if improper not in missing_impropers_dedu:
                missing_impropers_dedu.append(improper)

        with open("/tmp/" + output_name, "a") as file:
            file.write("######################## Please add below to ffbonded.itp! ########################\n") 
            file.write(f"There are in total {len(missing_bonds_dedu)} missed bonds: {missing_bonds_dedu}\n")
            file.write(f"There are in total {len(missing_angles_dedu)} missed angles: {missing_angles_dedu}\n")
            file.write(f"There are in total {len(missing_dihedrals_dedu)} missed proper dihedrals: {missing_dihedrals_dedu}\n")
            file.write(f"There are in total {len(missing_impropers_dedu)} missed improper dihedrals: {missing_impropers_dedu}\n")
            file.write("########### Add below to [ bondtypes ] segment ###########\n")
            count = 0
            for bond in missing_bonds_dedu:
                for a1, a2, funct, bond_length, bond_strength in zip(lig_dict['bonds']['a1'], lig_dict['bonds']['a2'], lig_dict['bonds']['funct'], lig_dict['bonds']['bond_length'], lig_dict['bonds']['bond_strength']):
                    if bond[0] == a1 and bond[1] == a2:
                        file.write("   %-3s%-11s%1d%11.5f%11.1f\n" % (a1, a2, int(funct), float(bond_length), float(bond_strength)))
                        break
            file.write("########### Add below to [ angletypes ] segment ###########\n")
            for angle in missing_angles_dedu:
                for a1, a2, a3, funct, angle_degree, angle_strength in zip(lig_dict['angles']['a1'], lig_dict['angles']['a2'], lig_dict['angles']['a3'], lig_dict['angles']['funct'], lig_dict['angles']['angle_degree'], lig_dict['angles']['angle_strength']):
                    if angle[0] == a1 and angle[1] == a2 and angle[2] == a3:
                        file.write("%-4s%-4s%-13s%1d%10.3f%11.3f\n" % (a1, a2, a3, int(funct), float(angle_degree), float(angle_strength)))
                        break
            file.write("########### Add below to [ dihedraltypes ] function = 1 segment ###########\n")
            for dihedral in missing_impropers_dedu:
                for a1, a2, a3, a4, funct, angle_degree, angle_strength, multiplicity in zip(lig_dict['impropers']['a1'], lig_dict['impropers']['a2'], lig_dict['impropers']['a3'], lig_dict['impropers']['a4'], lig_dict['impropers']['funct'], lig_dict['impropers']['angle_degree'], lig_dict['impropers']['angle_strength'], lig_dict['impropers']['multiplicity']):
                    if dihedral[0] == a1 and dihedral[1] == a2 and dihedral[2] == a3 and dihedral[3] == a4:
                        file.write("%-4s%-4s%-4s%-9s%1d%12.2f%12.5f%6d\n" % (a1, a2, a3, a4, int(funct), float(angle_degree), float(angle_strength), int(multiplicity)))
                        break
            file.write("########### Add below to [ dihedraltypes ] function = 3 segment ###########\n")
            for dihedral in missing_dihedrals_dedu:
                for a1, a2, a3, a4, funct, c0, c1, c2, c3, c4, c5 in zip(lig_dict['dihedrals']['a1'], lig_dict['dihedrals']['a2'], lig_dict['dihedrals']['a3'], lig_dict['dihedrals']['a4'], lig_dict['dihedrals']['funct'], lig_dict['dihedrals']['C0'], lig_dict['dihedrals']['C1'], lig_dict['dihedrals']['C2'], lig_dict['dihedrals']['C3'], lig_dict['dihedrals']['C4'], lig_dict['dihedrals']['C5']):
                    if dihedral[0] == a1 and dihedral[1] == a2 and dihedral[2] == a3 and dihedral[3] == a4:
                        file.write("%-4s%-4s%-4s%-9s%1d%12.6f%11.6f%11.6f%11.6f%11.6f%11.6f\n" % (a1, a2, a3, a4, int(funct), float(c0), float(c1), float(c2), float(c3), float(c4), float(c5)))
                        break

    # def parse_nonbonded(self, ff_nonbonded):  
        
    def parse_atomtypes(self, ff_aomtypes, ff_nonbonded):  # check if the user provided atom_types has any new atom type, if yes, then ask user the new atomtype's at.num(atomtypes), mass(atomtypes), sigma(ffnonbonded) and epsilon(ffnonbonded), give a guess value based on the database
        # Dictionary to store combined atom data
        atom_data = {}
        
        # Read data from atomtypes.atp
        with open(ff_aomtypes, 'r') as file:
            for line in file:
                if line.strip() and not line.startswith(';'):
                    parts = line.split()
                    atom_type = parts[0]
                    mass_atp = float(parts[1])
                    # Initialize or update the dictionary for this atom type
                    if atom_type not in atom_data:
                        atom_data[atom_type] = {'mass_atp': mass_atp}
                    else:
                        atom_data[atom_type]['mass_atp'] = mass_atp
        
        # Read data from ffnonbonded.itp
        with open(ff_nonbonded, 'r') as file:
            for line in file:
                if line.strip() and not line.startswith(';') and not line.startswith('['):
                    parts = line.split()
                    atom_type = parts[0]
                    at_num = int(parts[1])
                    mass_ffn = float(parts[2])
                    charge = float(parts[3])
                    ptype = parts[4]
                    sigma = float(parts[5])
                    epsilon = float(parts[6])
                    # Initialize or update the dictionary for this atom type
                    if atom_type not in atom_data:
                        atom_data[atom_type] = {
                            'at_num': at_num, 'mass_ffn': mass_ffn, 'charge': charge,
                            'ptype': ptype, 'sigma': sigma, 'epsilon': epsilon
                        }
                    else:
                        atom_data[atom_type].update({
                            'at_num': at_num, 'mass_ffn': mass_ffn, 'charge': charge,
                            'ptype': ptype, 'sigma': sigma, 'epsilon': epsilon
                        })
        
        return atom_data

    def find_missing_atomtypes(self, ff_atomtype, atom_types, output_name):
        # Check for missing atom types from the user-provided list
        missing_atom_types = []
        missing_atom_types = [atype for atype in atom_types if atype not in ff_atomtype]
        missing_atom_type = list(set(missing_atom_types))
        if missing_atom_types:
            with open("/tmp/" + output_name, "a") as file:
                file.write("######################## Please follow the steps below to edit atomtypes.atp and ffnonbonded.itp separately! ########################\n")
                file.write(f"There are missed atom types: {missing_atom_type}\n")
                file.write(f"Add the new atom type: {missing_atom_type} parameters to atomtypes.atp and ffnonbonded.itp\n")
        
        return missing_atom_type
        
    def streamlit_download_file(self, download_name, content_file):
        # Download topol.top file  #      
        # 打开 content_file 文件并读取其内容
        with open(content_file, 'r') as top_file:
            content = top_file.read()   

        #with open(content_file, 'r') as file:
        #    lines = file.readlines()
        #    st.write(f"after download, the file includes: {len(lines)}")

            
        # 添加一个下载按钮，传递 receptor_top_content 作为文件内容
        st.download_button(
            label = "Download " +  download_name,
            data = content,
            key = download_name,
            file_name = download_name
            )    
##########################################################################################################################################################################################
class RMSD_per_Residue():
    def __init__(self, pdb_file, output_name, name_list=None):
        self.pdb_file = pdb_file
        self.name_list = name_list if name_list is not None else []
        self.pdb_parser = PDBParser(QUIET=True)
        self.structure = self.pdb_parser.get_structure("Models", self.pdb_file)
        self.models = list(self.structure.get_models())

    def get_ca_coordinates_and_residues(self, model):
        coords, residues_info, residues_id = [], [], []
        for chain in model.get_chains():
            for residue in chain.get_residues():
                if "CA" in residue:
                    ca_atom = residue["CA"]
                    coords.append(ca_atom.get_coord())
                    residues_info.append(residue.resname + str(residue.id[1]))
                    residues_id.append(residue.id[1])
        return np.array(coords), residues_info, residues_id

    def calculate_distances(self):
        first_model = self.models[0]
        first_model_coords, first_model_residues, first_model_ids = self.get_ca_coordinates_and_residues(first_model)
        distances_per_model = []

        for model in self.models[1:]:  # Exclude the first model
            model_coords, model_residues, model_ids = self.get_ca_coordinates_and_residues(model)
            distances = [None] * len(first_model_ids)
            model_id_to_coord = {res_id: coord for res_id, coord in zip(model_ids, model_coords)}
            for i, res_id in enumerate(first_model_ids):
                if res_id in model_id_to_coord:
                    distances[i] = np.linalg.norm(first_model_coords[i] - model_id_to_coord[res_id])
            distances_per_model.append(distances)
        return first_model_residues, distances_per_model

    def bck_plot_distances(self, first_model_residues, output_name, distances_per_model):
        colors = list(mcolors.TABLEAU_COLORS) * ((len(distances_per_model) // len(mcolors.TABLEAU_COLORS)) + 1)
        plt.figure(figsize=(12, 8))
        for i, distances in enumerate(distances_per_model):
            plt.scatter(first_model_residues, distances, color=colors[i], label=self.name_list[i] if i < len(self.name_list) else f'Model {i+2}')
        plt.title('Distance from Model 1 to Other Models')
        plt.xlabel('Residue Name_Index')
        plt.ylabel('Distance (Å)')
        plt.legend(title='Model Comparison', bbox_to_anchor=(1.05, 1), loc='upper left')
        ticks = range(0, len(first_model_residues), 10)
        plt.xticks(ticks, [first_model_residues[i] for i in ticks], rotation=45)
        plt.grid(True)
        plt.tight_layout()
        plt.savefig("/tmp/" + output_name)
        self.streamlit_download_file(output_name, "/tmp/" + output_name)
        # plt.show()
    def plot_distances(self, first_model_residues, output_name, distances_per_model):
        colors = list(mcolors.TABLEAU_COLORS) * ((len(distances_per_model) // len(mcolors.TABLEAU_COLORS)) + 1)
        plt.figure(figsize=(12, 8))  # Size of the figure in inches
        for i, distances in enumerate(distances_per_model):
            plt.scatter(first_model_residues, distances, color=colors[i], label=self.name_list[i] if i < len(self.name_list) else f'Model {i+2}')
        plt.title('Distance from Model 1 to Other Models')
        plt.xlabel('Residue Name_Index')
        plt.ylabel('Distance (Å)')
        plt.legend(title='Model Comparison', bbox_to_anchor=(1.05, 1), loc='upper left')
        ticks = range(0, len(first_model_residues), 10)
        plt.xticks(ticks, [first_model_residues[i] for i in ticks], rotation=45)
        plt.grid(True)
        plt.tight_layout()

        # Save the figure with a higher resolution
        plt.savefig("/tmp/" + output_name, dpi=600)  # Increase the DPI for higher resolution
        self.streamlit_download_file(output_name, "/tmp/" + output_name)
    def streamlit_download_file(self, download_name, content_file):
        # Download topol.top file  #      
        # 打开 content_file 文件并读取其内容
        with open(content_file, 'rb') as top_file:
            content = top_file.read()   

        #with open(content_file, 'r') as file:
        #    lines = file.readlines()
        #    st.write(f"after download, the file includes: {len(lines)}")

            
        # 添加一个下载按钮，传递 receptor_top_content 作为文件内容
        st.download_button(
            label = "Download " +  download_name,
            data = content,
            key = download_name,
            file_name = download_name
            )   
##########################################################################################################################################################################################
class PDBModifier:
    def __init__(self, input_filename, output_filename, resnumber=None):
        self.headers = []  # 存储每个帧的头部信息（如REMARK, TITLE, CRYST1等）
        self.models = []  # 存储每个MODEL块的数据
        self.atomic_index = []
        self.atomic_name = []
        self.residue_name = []
        self.chain_name = []
        self.residue_index = []
        self.X_peratom = []
        self.Y_peratom = []
        self.Z_peratom = []
        self.bfactor_per_factor = []
        self.charge_per_factor = []
        self.Atomtype_per_atom = []
        self.chain_id_map = {}
        self.number_of_models = 0
        self.number_of_atoms_per_model = 0
        self.number_of_residues = 1
        ########################################
        self.input_filename = input_filename
        self.output_filename = output_filename
        self.resnumber = resnumber

        # 读取PDB文件，检查是否有链ID信息
        self.has_chain_id = self.PDBreader(self.input_filename)

        # 如果没有链ID，则添加链ID
        if not self.has_chain_id:
            chain_info = self.parse_chain_info(resnumber)
            self.add_chain_id(chain_info=chain_info)
            self.write_pdb(output_filename)

    def PDBreader(self, filename):
        # 清除之前的数据
        self.headers.clear()
        self.models.clear()
        self.atomic_index.clear()
        self.atomic_name.clear()
        self.residue_name.clear()
        self.chain_name.clear()
        self.residue_index.clear()
        self.X_peratom.clear()
        self.Y_peratom.clear()
        self.Z_peratom.clear()
        self.bfactor_per_factor.clear()
        self.charge_per_factor.clear()
        self.Atomtype_per_atom.clear()
        self.chain_id_map.clear()
        self.number_of_models = 0
        self.number_of_atoms_per_model = 0
        self.number_of_residues = 1

        has_chain_id = False
        current_model = []
        current_header = []
        in_model = False

        with open(filename, "r") as f:
            for line in f:
                if line.startswith(("REMARK", "TITLE", "CRYST1")):
                    current_header.append(line)
                elif line.startswith("MODEL"):
                    self.number_of_atoms_per_model = 0
                    self.number_of_models += 1
                    if current_header:  # 保存上一个头部信息
                        self.headers.append(current_header)
                        current_header = []
                    in_model = True
                    current_model = [line]

                elif in_model and not line.startswith("END"):
                    if line.startswith("ATOM") or line.startswith("HETATM"):
                        self.number_of_atoms_per_model += 1
                        self.atomic_index.append(float(line[6:11].strip()))
                        self.atomic_name.append(line[12:16].strip())
                        self.residue_name.append(line[17:20].strip())
                        chain_id = line[21].strip()
                        self.chain_name.append(chain_id)
                        self.residue_index.append(float(line[22:26].strip()))
                        if len(self.residue_index) > 1 and self.residue_index[-1] != self.residue_index[-2]:
                            self.number_of_residues += 1
                        self.X_peratom.append(float(line[30:38].strip()))
                        self.Y_peratom.append(float(line[38:46].strip()))
                        self.Z_peratom.append(float(line[46:54].strip()))
                        self.bfactor_per_factor.append(float(line[60:66].strip()))
                        self.charge_per_factor.append(float(line[54:60].strip()))
                        atom_type = line[76:78].strip() if len(line) > 76 else "null"
                        self.Atomtype_per_atom.append(atom_type)

                        if chain_id:
                            has_chain_id = True
                    current_model.append(line)
                elif line.startswith("END"):
                    current_model.append(line)
                    self.models.append(current_model)
                    in_model = False
            st.text(f"Number of atoms is: {self.number_of_atoms_per_model}")
            st.text(f"Number of models is: {self.number_of_models}")
            st.text(f"Number of residues is: {self.number_of_residues/self.number_of_models}")

        return has_chain_id

    def parse_chain_info(self, resnumber_list):
        if resnumber_list != "None":
            chain_info = {}
            for item in resnumber_list:
                chain_id, range_str = item.split(':')
                start, end = map(int, range_str.split('-'))
                chain_info[chain_id] = (start, end)
            # print(chain_info)
        else:
            chain_info = None
        
        return chain_info

    def add_chain_id(self, chain_info=None):
        chain_id_list = ['A', 'B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z']
        count = 0
        res_count = 1
        current_chain_id = chain_id_list[count]
        last_residue_index = None
        cumulative_atoms = 0  # 用于累计遍历过的原子数
        # max_res = 0
        max_res = 0
    
        if chain_info:
            # 提取所有范围的最大值
            max_res = max([end for _, (_, end) in chain_info.items()])
            # print(max_res)
        
        for i in range(len(self.residue_index)):
            
            if cumulative_atoms == self.number_of_atoms_per_model:
                cumulative_atoms = 0
                count = 0
                current_chain_id = chain_id_list[count]
            # if chain_info:
            #     # 根据用户提供的信息添加链ID
            #     if i > 0 and self.residue_index[i] != self.residue_index[i-1]:
            #         res_count += 1
            #     for chain_id, (start, end) in chain_info.items():
            #         if res_count == 1 or start <= res_count <= end:
            #             self.chain_name[i] = chain_id
            #             break
            if chain_info:
                # 判断是否进入下一个残基
                if i > 0 and self.residue_index[i] != self.residue_index[i-1]:
                    res_count += 1
                
                # 判断当前res_count在哪个范围内，并分配对应的chain ID
                if res_count <= max_res:
                    for chain_id, (start, end) in chain_info.items():
                        if start <= res_count <= end:
                            current_chain_id = chain_id
                            break
                # 如果res_count超过最大范围值，继续使用最后的chain ID
                self.chain_name[i] = current_chain_id
            else:
                # 自动检测并添加链ID
                if last_residue_index is not None and (self.residue_index[i] < last_residue_index or self.residue_index[i] > last_residue_index + 1):
                    count += 1
                    current_chain_id = chain_id_list[count]
                
                self.chain_name[i] = current_chain_id
    
            last_residue_index = self.residue_index[i]
            # 保存到映射中，以便应用于所有MODEL
            self.chain_id_map[self.residue_index[i]] = self.chain_name[i]
            cumulative_atoms += 1


    def write_pdb(self, output_filename):
        with open("/tmp/" + output_filename, "w") as f:
            for idx, model in enumerate(self.models):
                # 写入每个MODEL的头部信息
                for line in self.headers[idx]:
                    f.write(line)

                # 写入每个MODEL的数据
                for line in model:
                    if line.startswith("ATOM") or line.startswith("HETATM"):
                        atom_idx = self.atomic_index.index(float(line[6:11].strip()))
                        chain_id = self.chain_name[atom_idx]
                        atom_type = self.Atomtype_per_atom[atom_idx].replace("null", "  ")
                        f.write(f"{line[:21]}{chain_id:<1}{line[22:76]}{atom_type}\n")
                    else:
                        f.write(line)
        
        self.streamlit_download_file(output_filename, "/tmp/" + output_filename)            
        
    def streamlit_download_file(self, download_name, content_file):
        # Download topol.top file  #      
        # 打开 content_file 文件并读取其内容
        with open(content_file, 'r') as top_file:
            content = top_file.read()   

        #with open(content_file, 'r') as file:
        #    lines = file.readlines()
        #    st.write(f"after download, the file includes: {len(lines)}")

            
        # 添加一个下载按钮，传递 receptor_top_content 作为文件内容
        st.download_button(
            label = "Download " +  download_name,
            data = content,
            key = download_name,
            file_name = download_name
            )


##########################################################################################################################################################################################

class acpype4ligand(): 
    mol2_files      = []
    charge_type     = "user"
    ligand_names    = []

    def __init__(self, mol2_files, ligand_names, charge_type):
        self.mol2_files = mol2_files
        self.charge_type = charge_type
        
        # 1) 移动/重命名上传的文件到 /tmp
        for i in range(len(self.mol2_files)):
            old_path = mol2_files[i]
            ligand_name = ligand_names[i]+".mol2"
            new_path = os.path.join("/tmp", ligand_name)
            # st.text(new_path)
            os.rename(old_path, new_path)

            # 2) 调用 acpype 转换
            self.converter(new_path, charge_type)
            # files = os.listdir(ligand_names[i] + ".acpype")
            # st.text(files) 

        # 3) 打包 ACPYPE 输出并提供下载
        zip_buffer = self.create_zip_file_of_acpype_dirs()
        self.streamlit_download_zip_file(zip_buffer)    

        # 4) 删除生成的文件夹
        os.system("rm -r *.acpype")

    def converter(self, ligand, charge_type):
        # 注意: ACPYPE 通常会在同目录下生成一个 ligand_name.acpype/ 文件夹
        os.system(f"acpype -i {ligand} -c {charge_type}")

    def create_zip_file_of_acpype_dirs(self):
        """只打包当前文件夹中所有 *.acpype 的目录，并返回一个 ZIP 的 BytesIO。"""
        zip_buffer = io.BytesIO()

        with zipfile.ZipFile(zip_buffer, "w", zipfile.ZIP_DEFLATED) as zip_file:
            # 遍历当前目录，查找后缀为 .acpype 的目录
            for item in os.listdir('.'):
                if item.endswith('.acpype') and os.path.isdir(item):
                    # 递归打包整个目录
                    for root, dirs, files in os.walk(item):
                        for file in files:
                            full_path = os.path.join(root, file)
                            # arcname要使用相对路径，以保持目录层级结构
                            arcname = os.path.relpath(full_path, start='.')
                            zip_file.write(full_path, arcname=arcname)

        zip_buffer.seek(0)  # 将指针重置到开头，方便后续读取/下载
        return zip_buffer      

    def streamlit_download_zip_file(self, zip_buffer):
        st.download_button(
            label="Download Ligand Acpype Output as ZIP",
            data=zip_buffer,
            file_name="ligands_acpype.zip",
            mime="application/zip"
        )   
##########################################################################################################################################################################################



# Title
st.title("Welcome to gmx tool box 2.0")

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
    multi_files = st.file_uploader("Upload files", accept_multiple_files=True, type=['xvg', 'gro', 'itp', 'top', 'csv', 'dat'])
    # 保存上传文件的文件名
    uploaded_filenames = [uploaded_file.name for uploaded_file in multi_files]
    # 保存上传的文件到临时位置
    tmp_path = [save_uploaded_file(multi_files[i]) for i in range(len(multi_files))]
    output_name = st.text_input("Output file name", 'output.png')
    renumber = st.selectbox("Renumber residues", ['false', 'true'])
    rdf_cutoff = st.number_input("RDF cutoff value", min_value=0.0, step=0.1, value=0.0)
    error_bar = st.selectbox("Error bar or Error band", ['false', 'error bar', 'error band'])
    replica_number = st.number_input("Number of replicas, name format: rmsd-1.xvg rmsd-2.xvg rmsd-3.xvg; you may change 'rmsd' to any name you like.", min_value=1, step=1, value=3)
    transparency = st.number_input("Number of transparent (for error band)", min_value=0.0, max_value=1.0, step=0.1, value=0.2)
    average = st.selectbox("Calculate average", ['false', 'true'])
    plot_name = st.text_input("Plot title", value="auto detect")
    smooth = st.selectbox("Heatmap style", ['false', 'true'])
    pca_color_by_replicas = st.text_input("Whether color the pca dots by the order of replicas", value="0")
    nbin = st.number_input("Number of bins (for PCA set to 6)", min_value=1, step=1, value=6)
    size = st.number_input("Size (for PCA)", min_value=1, step=1, value=500)
    move_average = st.number_input("Window size for moving average", min_value=0, step=1, value=0)
    mean_value = st.selectbox("Draw mean value line", ['false', 'true'])
    histogram = st.selectbox("Generate histogram plot", ['false', 'true'])
    xaxis_name = st.text_input("X-axis name", value="auto detect")
    yaxis_name = st.text_input("Y-axis name", value="auto detect")
    subcolX, subcolY = st.columns(2)
    with subcolX:
        x_low = st.number_input("Lower boundary of X", min_value=0.0, step=0.1, value=0.0)
        y_low = st.number_input("Lower boundary of Y", min_value=0.0, step=0.1, value=0.0)
        width_size = st.number_input("Width of plot", min_value=0, step=100, value=800)
    with subcolY:
        x_high = st.number_input("Higher boundary of X", min_value=0.0, step=0.1, value=0.0)
        y_high = st.number_input("Higher boundary of Y", min_value=0.0, step=0.1, value=0.0)
        height_size = st.number_input("Height of plot", min_value=0, step=100, value=600)
    xy_font = st.number_input("XY-axis font size", min_value=0, step=1, value=40)
    title_font = st.number_input("Title font size", min_value=0, step=1, value=24)
    legend_show = st.selectbox("Show legend", ['True', 'False'])
    legend_show = ast.literal_eval(legend_show)
    legend_font = st.number_input("Legend font size", min_value=0, step=1, value=30)
    font_family = st.text_input("Font family", 'Arial')
    grid_show = st.selectbox("Show grid", ['True', 'False'])
    grid_show = ast.literal_eval(grid_show)
    axis_show = st.selectbox("Show axis", ['True', 'False'])
    axis_show = ast.literal_eval(axis_show)
    line_width = st.number_input("Width of axis", min_value=1, step=1, value=2)
    margin_l  = st.number_input("left margin", min_value=0, step=10, value=120)
    margin_r  = st.number_input("right margin", min_value=0, step=10, value=60)
    margin_t  = st.number_input("top margin", min_value=0, step=10, value=60)
    margin_b  = st.number_input("bottom margin", min_value=0, step=10, value=100)
    font_color = st.color_picker("Font color", '#000000')
    trace_color_scheme = st.text_area("Traces Color Scheme", value="#636EFA #EF553B #00CC96 #AB63FA #FFA15A #19D3F3 #FF6692 #B6E880 #FF97FF #FECB52")
    violin = st.selectbox("violin style", ['False', 'True'])

    if st.button('Plotting') and multi_files[0] != 0:
        x = plotly_go(tmp_path, output_name, renumber, rdf_cutoff, average, plot_name, nbin, size, move_average, mean_value, histogram, xaxis_name, yaxis_name, width_size, height_size, xy_font, title_font, legend_show, legend_font, font_family, font_color, grid_show, uploaded_filenames, margin_l, margin_r, margin_t, margin_b, violin, smooth, error_bar, replica_number, axis_show, line_width, transparency, x_low, x_high, y_low, y_high, trace_color_scheme, pca_color_by_replicas)

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

    ### 1st subcolum in column 2, for renumber the MODEL number in extreme1.pdb ###
    mradder.subheader("Renumber the MODEL number for pdb files")
    extreme_pdb = st.file_uploader("Upload the extreme.pdb file", type=['pdb'])
    # save the files name:
    if extreme_pdb:
        extreme_pdb_name = extreme_pdb.name
    # 保存上传的文件到临时位置
    if st.button('Renumber') and extreme_pdb:
        # 保存文件并获取临时文件路径
        extreme_pdb_path = save_uploaded_file(extreme_pdb)
        # 如果文件保存成功，则进行合并
        if extreme_pdb_path:
            x = renumber_MODEL(extreme_pdb_path,extreme_pdb_name)
        else:
            st.error("Failed to renumber the pdb file, Do you have a MODEL line in your pdb file?")
    else:
        pass  

    ### 2nd subcolum in column 2, for adding new residues to gromacs forcefield ###
    mradder.subheader("Generate forcefield parameters for new residues")
    st.text("Ps: you can upload files in below combination:")
    st.text("1. lig_GMX.itp + ffbonded.itp")
    st.text("2. lig_GMX.itp + ffbonded.itp + ffnonbonded.itp + atomtypes.atp")
    # Collecting user inputs
    tmp_lig_gro_path = 0
    tmp_ff_rtp_path = 0
    tmp_ff_hdb_path = 0
    tmp_ff_restype_path = 0
    tmp_ff_nonbonded_path = 0
    tmp_ff_aomtypes_path = 0
    ###
    ar_lig_itp = st.file_uploader("Upload lig_GMX.itp", accept_multiple_files=True, type=['itp'])
    uploaded_lig_itp_name = [uploaded_file.name for uploaded_file in ar_lig_itp]
    try:
        tmp_lig_itp_path = save_uploaded_file(ar_lig_itp[0])
    except:
        pass
    ###
    ar_ff_bonded = st.file_uploader("Upload ffbonded.itp", accept_multiple_files=True, type=['itp'])
    uploaded_ff_bonded_name = [uploaded_file.name for uploaded_file in ar_ff_bonded]
    try:
        tmp_ff_bonded_path = save_uploaded_file(ar_ff_bonded[0])
    except:
        pass
    ###
    ar_ff_nonbonded = st.file_uploader("Upload ffnonbonded.itp", accept_multiple_files=True, type=['itp'])
    uploaded_ff_nonbonded_name = [uploaded_file.name for uploaded_file in ar_ff_nonbonded]
    try:
        tmp_ff_nonbonded_path = save_uploaded_file(ar_ff_nonbonded[0])
    except:
        pass
    ###
    ar_ff_aomtypes = st.file_uploader("Upload atomtypes.atp", accept_multiple_files=True, type=['atp'])
    uploaded_ff_aomtypes_name = [uploaded_file.name for uploaded_file in ar_ff_aomtypes]
    try:
        tmp_ff_aomtypes_path = save_uploaded_file(ar_ff_aomtypes[0])
    except:
        pass
    ###
    ar_lig_gro = st.file_uploader("Upload lig_GMX.gro", accept_multiple_files=True, type=['gro'])
    uploaded_lig_gro_name = [uploaded_file.name for uploaded_file in ar_lig_gro]
    try:
        tmp_lig_gro_path = save_uploaded_file(ar_lig_gro[0])
    except:
        pass
    ###
    ar_ff_rtp = st.file_uploader("Upload aminoacids.rtp", accept_multiple_files=True, type=['rtp'])
    uploaded_ff_rtp_name = [uploaded_file.name for uploaded_file in ar_ff_rtp]
    try:
        tmp_ff_rtp_path = save_uploaded_file(ar_ff_rtp[0])
    except:
        pass
    ###
    ar_ff_hdb = st.file_uploader("Upload aminoacids.hdb", accept_multiple_files=True, type=['hdb'])
    uploaded_ff_hdb_name = [uploaded_file.name for uploaded_file in ar_ff_hdb]
    try:
        tmp_ff_hdb_path = save_uploaded_file(ar_ff_hdb[0])
    except:
        pass

    ###
    ar_ff_restype = st.file_uploader("Upload residuetypes.dat", accept_multiple_files=True, type=['dat'])
    uploaded_ff_restype_name = [uploaded_file.name for uploaded_file in ar_ff_restype]
    try:
        tmp_ff_restype_path = save_uploaded_file(ar_ff_restype[0])
    except:
        pass
    ###
    isamino_acid = st.selectbox("Whethe treat it as an amino acids?", ['true', 'false'])
    ###
    ar_atom_types_input = st.text_area("Atom Types ordered as the atoms in mol2 file, Example: N CT C O CT S CT CT N3 C O CT CT NX NX O O O O H1 H1 H1 H1 H1 H1 H1 H1 H1 H")
    atom_types = ar_atom_types_input.split()
    ###
    ar_output_name = st.text_input("Output file name", 'ff_parameter.txt')
    ###
    if st.button('Generate') and ar_lig_itp[0] != 0:
        x = ff_res_adder(tmp_lig_itp_path, tmp_lig_gro_path, tmp_ff_rtp_path, tmp_ff_hdb_path, tmp_ff_bonded_path, tmp_ff_nonbonded_path, tmp_ff_aomtypes_path, tmp_ff_restype_path, atom_types, isamino_acid, ar_output_name)

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
    
    ### subcolum in column 3, for dssp plotting ###
    gromerge.subheader("Plot DSSP")
    ds_data = st.file_uploader("Upload the dssp.dat file", type=['dat'])
    ds_traj = st.file_uploader("Upload the trajectory file for time & residues", type=['pdb', 'xtc', 'trr']) 
    ds_output_name = st.text_input("Output file name", 'dssp.png')
    ds_origianl = st.selectbox("Whether use full structure types", ['true', 'false'])
    ds_unique_color = st.selectbox("Whether include a unique colorbar", ['true', 'false'])
    ds_original_map = st.text_area("Enter the full mapping in JSON format:", '{"H": 9, "B": 8, "E": 7, "G": 6, "I": 5, "P": 4, "S": 3, "T": 2, "=": 1, "~": 0}')
    ds_original_colorbar = st.text_area("Enter the full color bar:", '{"loop": 0, "break": 1, "h-bond turn": 2, "bend": 3, "kappa helix": 4, "pi helix": 5, "3_10 helix": 6, "strand": 7, "beta bridge": 8, "a-helix": 9}')
    ds_original_color_scale = st.text_area("Enter original colorscale in JSON format:", '[[0.0, "#636EFA"], [0.1, "#636EFA"], [0.1, "#EF553B"], [0.2, "#EF553B"], [0.2, "#00CC96"], [0.3, "#00CC96"], [0.3, "#AB63FA"], [0.4, "#AB63FA"], [0.4, "#FFA15A"], [0.5, "#FFA15A"], [0.5, "#19D3F3"], [0.6, "#19D3F3"], [0.6, "#FF6692"], [0.7, "#FF6692"], [0.7, "#B6E880"], [0.8, "#B6E880"], [0.8, "#FF97FF"], [0.9, "#FF97FF"], [0.9, "#FECB52"], [1.0, "#FECB52"]]')
    ds_simple_map = st.text_area("Enter a simple mapping in JSON format:", '{"H": 2, "B": 1, "E": 3, "G": 1, "I": 1, "P": 1, "S": 1, "T": 1, "=": 1, "~": 1}')
    ds_simple_colorbar = st.text_area("Enter the simple color bar:", '{ "loop": 1, "helix": 2, "beta sheet": 3 }')
    ds_simple_color_scale = st.text_area("Enter simple colorscale in JSON format:", '[ [0.00, "gold"], [0.33, "gold"], [0.33, "mediumturquoise"], [0.66, "mediumturquoise"], [0.66, "lightsalmon"], [1.00, "lightsalmon"] ]')
    # transfer JSON inputs to list or dictionary:
    original_map = json.loads(ds_original_map)
    original_colorbar = json.loads(ds_original_colorbar)
    original_colorscale = json.loads(ds_original_color_scale)
    simple_map = json.loads(ds_simple_map)
    simple_colorbar = json.loads(ds_simple_colorbar)
    simple_colorscale = json.loads(ds_simple_color_scale)


    # save the files name:
    if ds_data and ds_traj:
        dssp_data = ds_data.name
        dssp_traj = ds_traj.name
    # 保存上传的文件到临时位置
    if st.button('plot') and ds_data and ds_traj:
        # 保存文件并获取临时文件路径
        dssp_data_path = save_uploaded_file(ds_data)
        dssp_traj_path = save_uploaded_file(ds_traj)
        # 如果文件保存成功，则进行合并
        if dssp_data_path and dssp_traj_path:
            # with open(receptor_gro_path, 'r') as file:
                # file_content = file.read()
                # st.text(file_content)
            x = gmx_dssp(dssp_data_path, dssp_traj_path, ds_output_name, ds_origianl, ds_unique_color, original_map, original_colorbar, original_colorscale, simple_map, simple_colorbar, simple_colorscale)
            # st.success("Peptide converted successfully!")
        else:
            st.error("Failed to plot the dssp, Is the dssp.dat generated by: \n  \
                     gmx dssp -f md_noPBC_dt1000.pdb -s md_noPBC_dt1000.pdb -o dssp.dat -xvg xmgrace")
    else:
        pass 

    ### subcolum in column 3, for acpype converting ligand ###
    gromerge.subheader("Acpype Converting Ligand")
    acl_ligand = st.file_uploader("Upload the ligand mol2 file", accept_multiple_files=True, type=['mol2'])
    acl_charge_type = st.selectbox("Charge Types", ['user', 'bcc']) 



    # 保存上传的文件到临时位置
    if st.button('convert') and acl_ligand:
        # 保存上传文件的文件名
        acl_ligand_names = [uploaded_file.name for uploaded_file in acl_ligand]
        # Get rid of the subaddress
        acl_ligand_names = [os.path.splitext(a)[0] for a in acl_ligand_names]
        # st.text(acl_ligand_names)
        # 保存上传的文件到临时位置
        acl_ligand_path = [save_uploaded_file(acl_ligand[i]) for i in range(len(acl_ligand))]
        # 如果文件保存成功，则进行合并
        if acl_ligand_path:
            x = acpype4ligand(acl_ligand_path, acl_ligand_names, acl_charge_type)
            # st.success("Peptide converted successfully!")
        else:
            st.error("Failed to convert the Peptide, do you have the last column (atom type column) in you peptide PDB?")
    else:
        pass   
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

    


    ### subcolum in column 4, for Calculate the superimposed models distances per residue ###
    contact_map.subheader("Distances per residues' alpha carbon")
    rmsd_pdb = st.file_uploader("Upload the Sperimposed multi-model pdb file", type=['pdb']) 
    rmsd_name = st.text_area("Name list for the legend, e.g swiss itasser modeller", value="Model-1 Model-2 Model-3 Model-4 Model-5") 
    # 将输入的字符串分割成列表
    name_list = rmsd_name.split()  # 默认以空格分割
    rmsd_output_name =  st.text_input("Output name", 'figure.png') 
    # save the files name:
    if rmsd_pdb and rmsd_name:
        multimodel_pdb = rmsd_pdb.name
    # 保存上传的文件到临时位置
    if st.button('Plot') and rmsd_pdb and rmsd_name:
        # 保存文件并获取临时文件路径
        peptide_pdb_path = save_uploaded_file(rmsd_pdb)
        # 如果文件保存成功，则进行合并
        if peptide_pdb_path:
            # with open(receptor_gro_path, 'r') as file:
                # file_content = file.read()
                # st.text(file_content)
            rmsd = RMSD_per_Residue(peptide_pdb_path, rmsd_output_name, name_list)
            first_model_residues, distances_per_model = rmsd.calculate_distances()
            rmsd.plot_distances(first_model_residues, rmsd_output_name, distances_per_model)
            # st.success("Peptide converted successfully!")
        else:
            st.error("Failed to plot the RMSD per residues, please delete the TER and CONNECT informations in the PDB file")
    else:
        pass      

    ### subcolum in column 4, for adding ChainID to the pdb trajectories ###
    contact_map.subheader("ChainID Adder to Gromacs pdb trajectory")
    ca_traj_files    = st.file_uploader("Upload the PDB trajectory file", type=['pdb']) 
    ca_output_name   = st.text_area("Output Name", value="traj_with_chain.pdb") 
    ca_resnumber     = st.text_area("The chainID and number of residues, e.g: A:1-100 C:101-300 D:301-660 optional", value="None")

    # save the files name:
    if ca_traj_files and ca_output_name:
        traj_name = ca_traj_files.name
    # 保存上传的文件到临时位置
    if st.button('Add ChainID') and ca_traj_files:
        # 保存文件并获取临时文件路径
        traj_path = save_uploaded_file(ca_traj_files)
        # 如果文件保存成功，则进行合并
        if traj_path:
            PDBModifier(traj_path, ca_output_name, resnumber=ca_resnumber)
            # st.success("Peptide converted successfully!")
        else:
            st.error("Failed to add ChainID to the trajectory file")
    else:
        pass      
#################################################################################################################################################
    ### subcolum in column 4, for peptide convert to ligand ###
    contact_map.subheader("Peptide to Ligand")
    p2l_pdb = st.file_uploader("Upload the Peptide pdb file", accept_multiple_files=True, type=['pdb'])

    p2l_name =  st.text_input("Give the peptide a name", 'PEP') 

    # 保存上传的文件到临时位置
    if st.button('convert peptide') and p2l_pdb and p2l_name:
        # 保存上传文件的文件名
        uploaded_filenames = [uploaded_file.name for uploaded_file in p2l_pdb]
        # 保存上传的文件到临时位置
        peptide_pdb_path = [save_uploaded_file(p2l_pdb[i]) for i in range(len(p2l_pdb))]
        # 如果文件保存成功，则进行合并
        if peptide_pdb_path:
            # with open(receptor_gro_path, 'r') as file:
                # file_content = file.read()
                # st.text(file_content)
            x = pep2lig(peptide_pdb_path, uploaded_filenames, p2l_name)
            # st.success("Peptide converted successfully!")
        else:
            st.error("Failed to convert the Peptide, do you have the last column (atom type column) in you peptide PDB?")
    else:
        pass 
#################################################################################################################################################
