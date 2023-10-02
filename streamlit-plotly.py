import streamlit as st
import plotly.graph_objects as go
import plotly.io as pio
import re
import tempfile
import os
import io
import mimetypes
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

# buffer = io.BytesIO()
# print(type(buffer))
# print(buffer)


class gmxplotly:
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
    
    def __init__(self, file1, file2, file3, output_name, renumber, ave, xaxis_name, yaxis_name, rdf_cutoff, multi_files, plot_name, pca, nbin, size):
            file1 = multi_files[0]
            self.flag_recognizer(file1)
            if self.pca_flag != 1:
                self.plotly_multy(self.flag, multi_files, xaxis_name, yaxis_name, renumber, ave, output_name, rdf_cutoff, plot_name)
            elif self.pca_flag == 1:
                self.plotly_pca(multi_files, xaxis_name, yaxis_name, renumber, ave, output_name, rdf_cutoff, plot_name, nbin, size)
    
    def flag_recognizer(self,file1):                                                   # first method to be called in __main__, used for creating object and charactors.
        #self.atomic_index.clear()                                                   # cleaveage the information of previous object before put new record into these charactors

        # 使用临时文件路径打开文件进行处理
        with open(file1, 'r') as f:
            lines = f.readlines()
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
                    elif self.flag == 'anaeig,':
                        self.flag = 'pca'
                except: pass

            if len(lines) >= 9 and '-or' in lines[8]:
                self.sasa_flag = '-or'
                
            # if 'pca' in str(file1).lower() or '2dproj' in str(file1):
            if self.flag == 'pca':
                self.pca_flag = 1

  
    def plotly_multy(self, flag, multi_files, xaxis_name, yaxis_name, renumber, ave, output_file_name, rdf_cutoff, plot_name):
        # 将您的原始函数放在这里
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
        # pio.write_image(fig, buffer, format="png") # 将图像写入 BytesIO 对象
        # print(type(buffer))
        # print(buffer)
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
            # pio.write_image(fig, buffer, format="png") # 将图像写入 BytesIO 对象
            # print(type(buffer))
            # print(buffer)
    def plotly_pca(self, multi_files, xaxis_name, yaxis_name, renumber, ave, output_name, rdf_cutoff, plot_name, nbin, size):
        # 将您的原始函数放在这里
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
            # create empty list
            with open(i, "r") as f:
                lines = f.readlines()
                for num in range(len(lines)):
                    if lines[num].startswith("#") or lines[num].startswith("@"):
                        pass
                    else:
                        x_points.append(float(lines[num].split()[0]))
                        y_points.append(float(lines[num].split()[1])) 
                
            # 创建散点图轨迹
            # print(mode)
            scatter_trace = go.Scatter(
                x=x_points,
                y=y_points,
                mode='markers',
                marker=dict(
                    size=5,
                    color='black'
                )
            )
            #print("###################################################")
            #print(scatter_trace.mode)
            layout = go.Layout(title=plot_title, title_x=0.5, title_y=1, font=dict(size=20),
                                xaxis=dict(title=x_name, titlefont=dict(size=20, color='black', family='Arial'), zeroline=False, autorange=True,
                                          showgrid=False, gridwidth=1, gridcolor='rgba(0,0,0,0.1)', tickfont=dict(size=20)),
                                yaxis=dict(title=y_name, titlefont=dict(size=20, color='black', family='Arial'), zeroline=False, autorange=True,
                                          showgrid=False, gridwidth=1, gridcolor='rgba(0,0,0,0.1)', tickfont=dict(size=20)),
                                legend=dict(x=1, y=1, orientation='v', font=dict(size=30)), showlegend=False,
                                plot_bgcolor='rgba(255, 255, 255, 0.1)',
                                paper_bgcolor='rgba(255, 255, 255, 0.2)')
            
            # 创建图形对象
            fig = go.Figure(data=scatter_trace, layout=layout)
                    
            # 显示图形
            # fig.show()
            if output_name == 'plotly.png':
                pio.write_image(fig, "PCA_Scatter_"+i.split('.')[0]+".png")
                #pio.write_image(fig, "plot.png")
            else:
                # pio.write_image(fig, "PCA_Scatter_" + output_name)
                pio.write_image(fig, "PCA_Scatter_" + i.split('.')[0]+'_'+output_name.split('.')[0] + ".png")
                #pio.write_image(fig, "plot.png")
                
            # 将图像写入 BytesIO 对象
            # pio.write_image(fig, buffer, format="png") # 将图像写入 BytesIO 对象
            # print(type(buffer))
            # print(buffer)
            # buffer.seek(0) # 移动到缓冲区的开始位置，以便读取图像数据

            # # 清除 buffer 以便重复使用（如果需要）
            # buffer.truncate(0)
            # buffer.seek(0)
# 创建Streamlit应用程序
def main():
    # st.cache(allow_output_mutation=True)
    st.title("Shang's first streamlit app!")  # 设置应用程序标题
    # 添加用户界面元素，例如文本框、文件上传、选择框等
    st.write("Welcome to the plotly for gromacs!!")

    # 创建多个文件上传元素，并将它们添加到uploaded_files列表中
    multi_files = []
    upload_files = st.file_uploader(f"upload the xvg files：", type=["xvg"], accept_multiple_files=True)
    # 检查哪些文件已经上传
    if upload_files is not None:
        for upload_file in upload_files:
            file_name = upload_file.name
            # 使用tempfile模块创建临时文件
            # temp_file = tempfile.NamedTemporaryFile(prefix=file_name,delete=False, suffix=".xvg")

            # 将上传文件的内容写入临时文件
            # temp_file.write(upload_file.getvalue())

            # 关闭临时文件以确保内容被刷新到磁盘
            # temp_file.close()

            # 获取临时文件的路径，并将其添加到multiple_files列表中
            # temp_file_path = temp_file.name
            with open(file_name, 'w') as f:
                 f.write(upload_file.getvalue().decode("utf-8"))
            # multi_files.append(temp_file_path)
            multi_files.append(file_name)
            # 打印临时文件的路径（可选）
            st.write(f"Temporary files have been created：{file_name}")    
    else:
        st.write("Haven't upload your data")

    # 添加用户界面元素

    output_name = st.text_input("Please give the filename for output:")
    # 当用户输入文件名并点击提交按钮时执行以下操作

    ave = st.selectbox('Do you want to have the average?', ['true', 'fault'])
    renumber = st.selectbox("Do you want to renumber the sequence from 1?", ['true', 'fault'])
    
    st.button("提交")
    # # 获取文件列表
    # file_list = st.selectbox('Do you want to show the current folder?', ['true', 'fault'])
    # if file_list == 'true':
    #     folder_path='./'
    #     files = os.listdir(folder_path)
        
    #     # 在 Streamlit 应用中展示文件列表
    #     for file in files:
    #         st.write(file)
    
    
    # # Download files
    # outputname = st.text_input("Please input the png file name, e.g plot.png ")
    # if st.button("Submit"):
    #     if output_name:
    #         st.write(f"Your output filename is: {outputname}")
    #     else:
    #         st.write("Please input the filename for output")
    # with open("plot.png", "rb") as file:
    #     st.download_button(
    #         label="Download Plot as PNG",
    #         data=file,
    #         file_name=outputname,
    #         mime="image/png",
    #     )
    # with open("Average_plot.png", "rb") as file:
    #     if file:
    #         st.download_button(
    #             label="Download Average Plot as PNG",
    #             data=file,
    #             file_name="Average_" + outputname,
    #             mime="image/png",
    #         )
    
    # 获取当前文件夹中的文件
    folder_path = '.'  # 使用 '.' 代表当前文件夹       
    files = [f for f in os.listdir(folder_path) if os.path.isfile(os.path.join(folder_path, f))]
    
    # 让用户在 UI 中选择文件
    selected_file = st.selectbox('Pick the file you want to download', files)
    
    # 获取文件的绝对路径
    file_path = os.path.join(folder_path, selected_file)
    
    # 读取文件内容
    with open(file_path, "rb") as file:
        file_content = file.read()
    
    # 获取文件的 MIME 类型
    mime_type, _ = mimetypes.guess_type(file_path)
    
    # 创建下载按钮
    download_button = st.download_button(
        label=f"Download {selected_file}",
        data=file_content,
        file_name=selected_file,
        mime=mime_type,
    )

    # 如果用户点击下载按钮，显示一条消息
    if download_button:
        st.write(f"You have donwloaded {selected_file}")


    # 检查按钮是否被按下
    if st.button('Clear Files'):
        # 定义要删除的文件路径（这里可以根据需要修改）
        file_path = '\*png'
        
        # 检查文件是否存在
        if os.path.exists(file_path):
            # 删除文件
            os.remove(file_path)
            st.write(f"File {file_path} has been removed!")
        else:
            st.write(f"File {file_path} does not exist!")
            
    # 创建一个实例
    app = gmxplotly(file1,file2,file3, output_name, renumber, ave, xaxis_name, yaxis_name, rdf_cutoff, multi_files, plot_name, pca, nbin, size)
    # 删除临时文件
    for i in multi_files:
        os.remove(i)
        




if __name__ == '__main__':
    main()

