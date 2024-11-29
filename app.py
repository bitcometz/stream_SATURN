import streamlit as st
import pandas as pd
from streamlit import session_state as ss
import scanpy as sc
import json
from helpers.file_handling import *
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import plotly.express as px
import time
import pickle

# inputs 展示两个物种的umap图，zebrafish为细胞注释，frog为cluster id
# run SATURN，sleep 20秒
# 整合umap，展示物种、细胞标签的混合情况
# 另外一个为frog本身的，可以选择cluster和预测的
# 小于1%的注释结果不要，然后先用3w距离，然后subsample


def showSGR():
    singleron_base64 = read_image("src/singleron.png")
    with st.sidebar:
        st.markdown("---")

        st.markdown(
            f'<div style="margin-top: 0.25em; text-align: left;"><a href="https://singleron.bio/" target="_blank"><img src="data:image/png;base64,{singleron_base64}" alt="Homepage" style="object-fit: contain; max-width: 174px; max-height: 41px;"></a></div>',
            unsafe_allow_html=True,
        )
        st.markdown(
            '''
            sessupport@singleronbio.com
            '''
        )

def binaryswitch(session_state, keys):
    for key in keys:
        if session_state[key] is True:
            session_state[key] = False
        else:
            session_state[key] = True

### 画图函数
def plot4umap(df, label, symbol, title, pos, y_b=-0.6):
    fig   = go.Figure()
    x_min = min( df["umap_x"] )
    x_max = max( df["umap_x"] )
    y_min = min( df["umap_y"] )
    y_max = max( df["umap_y"] )
    # 将字符串标签映射到内置颜色序列中的索引
    label_to_index = {label: i for i, label in enumerate(df[label].unique())}
    color_indices = df[label].map(label_to_index)
    # 使用 Plotly Express 默认的颜色循环进行映射
    colors = px.colors.qualitative.Plotly

    lab2color = {}
    tmp = -1
    for k in df[label].unique():
        tmp = tmp+1 if tmp < len(colors)-1 else 0
        lab2color[k] = colors[tmp]


    for cell_type in df[label].unique():
        df_tmp = df[df[label] == cell_type]
        cx     = df_tmp["umap_x"]
        cy     = df_tmp["umap_y"]
        ct     = df_tmp[label]

        fig.add_trace(
            go.Scatter(
                x=cx, 
                y=cy,
                mode='markers',
                marker=dict(color=lab2color[cell_type], size=4, symbol=symbol),  # 使用内置颜色序列进行映射
                name=str(cell_type),
                text=ct,  # 使用 'text' 列的值作为提示信息
                hoverinfo='text'  # 显示 text 参数中的信息
            ),
        )

    # 图框
    m0 = 1.18
    fig.add_shape(
        type="rect",
        x0=x_min * m0,
        y0=y_min * m0,
        x1=x_max * m0,
        y1=y_max * m0,
        line=dict(color="black", width=1),
        opacity=0.5,
    )
    #fig.update_layout(height=fh, width=fw)
    #fig.update_layout( legend_title_text=label, title="Umap of " + label)
    # 创建一个空白的布局
    layout = go.Layout(
        plot_bgcolor  = 'rgba(0,0,0,0)', 
        paper_bgcolor = 'rgba(0,0,0,0)', 
        xaxis         = dict(visible=False),  # 隐藏 x 轴
        yaxis         = dict(visible=False),  # 隐藏 y 轴
    )
    fig.update_layout(layout)
    if pos == "bottom":
        fig.update_layout(
            go.Layout(
                legend=dict(orientation="h", yanchor="bottom", y=y_b, xanchor="right", x=0.9)
            )
        )
    
    # 添加标题
    fig.update_layout(
        title=title,
        font=dict(
            size=23  # 设置字体大小为18
        )

    )
    if label == "old_labels":
        fig.update_layout(
            width=700,  # 图片宽度为整个显示区域宽度的 80%
            height=700  # 图片高度为整个显示区域高度的 60%
        )   

    return fig    


@st.cache_data
def readCSV(file):
    df = pd.read_csv(file)
    return df

@st.cache_data
def readH5AD(file):
    adata = sc.read_h5ad(file)
    return adata

@st.cache_data
def readPKL(file):
    with open(file, "rb") as f:
        macrogene_weights = pickle.load(f)
    return macrogene_weights

def main():
    # initial variables
    default_variables = {
        "dataOK": False,
        "adata_ref": None,
        "adata_que": None,
        "run_flag":False,
    }
    if "dataOK" not in ss:
        for key, value in default_variables.items():
            if key not in st.session_state:
                ss[key] = value

    ## 设置页面宽一些
    st.set_page_config(layout="wide")
    st.sidebar.markdown("# SATURN")
    st.sidebar.markdown("## 1.Inputs")

    input_selectbox = st.sidebar.selectbox(
        "What data would you like to use for analysis?",
        ("Demo", "Upload new"),
        index=None, placeholder="Please select ..."
    )
    ## Inputs
    if input_selectbox is None:
        st.title("Welcome !!!")
        st.write("Please select the inputs from the left slidebar!")
    else:
        if input_selectbox == "Demo":
            ss["adata_ref"] = pd.read_csv("./demo/frog_zebrafish/zebrafish.ori.csv")
            ss["adata_que"] = pd.read_csv("./demo/frog_zebrafish/frog.ori.csv")
            ss["dataOK"]    = True
        else:
            st.subheader("Upload your RNA-Seq data (h5ad format)")
            exp_file  = st.file_uploader("Choose a h5ad file for reference",   type="csv", disabled=True)
            #st.write("Note: The first column of the matrix contains gene names, followed by the expression matrix of each sample, with column names representing the sample names.")
            spa_file = st.file_uploader("Choose a h5ad file for query", type="csv", disabled=True)
            #st.write("Note: This file contains grouping information. The first column is the sample name, with column name \"sample\". The values in the first column correspond to the expression matrix of the samples mentioned earlier. The second column contains grouping information with column name \"group\".")

        if ss["dataOK"] == True:
            st.markdown("## 1.Inputs")
            nlist       = ['barcode', 'umap_x', 'umap_y']
            col1, col2  = st.columns([1, 1])
            ref_cols    = [i for i in ss["adata_ref"].columns if i not in nlist]

            ref_label   = st.sidebar.selectbox( "Reference cell label:", ref_cols, index=0, placeholder="Please select ...")
            umap_ref    = plot4umap(ss["adata_ref"], ref_label, "circle", "Zebrafish umap of "+ref_label, "bottom")
            col1.plotly_chart(umap_ref, use_container_width=True)

            
            que_cols    = [i for i in ss["adata_que"].columns if i not in nlist]
            que_label   = st.sidebar.selectbox( "Query cell label:", que_cols, index=0, placeholder="Please select ...")
            umap_que    = plot4umap(ss["adata_que"], que_label, "triangle-up", "Frog umap of "+que_label, "bottom")
            col2.plotly_chart(umap_que, use_container_width=True)

        else:
            st.write("Currently, this app only supports demo data demonstration!")
            st.stop()

        ## Running
        if st.button("Run SATURN!"):
            with st.spinner('Wait for running ...'):
                time.sleep(1)

            ss["run_flag"] = True
            st.write("Done!")

        if ss["run_flag"]:
            st.markdown("## 2.Results")
            st.sidebar.markdown("## 2.Results")
            #st.markdown("#### 物种整合umap结果展示")
            df = readCSV("./demo/frog_zebrafish/saturn.6k.csv")
            nlist       = ['barcode', 'umap_x', 'umap_y']
            all_cols    = [i for i in df.columns if i not in nlist]

            all_label   = st.sidebar.selectbox( "2.1 Results cell label:", all_cols, index=0, placeholder="Please select ...")
            
            if all_label == "old_labels":
                umap_all    = plot4umap(df, all_label, "circle", "2.1 Integration results umap of "+all_label, "bottom", y_b=-0.6)
            else:
                umap_all    = plot4umap(df, all_label, "circle", "2.1 Integration results umap of "+all_label, "right")
            st.plotly_chart(umap_all, use_container_width=True)

            ## 如果选择的是细胞类型，则进行macrogene的weight分析
            if all_label == "cell_type":

                mac_cols = df["cell_type"].unique()
                st.markdown("### 2.2 Macrogene Differential Expression")
                adata_m = readH5AD("./demo/frog_zebrafish/macrogene_6k.h5ad")
                m_label   = st.sidebar.selectbox( "2.2 Celltype to vs rest:", mac_cols, index=0, placeholder="Please select ...")

                sc.tl.rank_genes_groups(adata_m, groupby="labels2", groups=[m_label], method="wilcoxon")
                dot_plot = sc.pl.rank_genes_groups_dotplot(adata_m, return_fig=True, n_genes=7)
                dot_plot.swap_axes()
                st.pyplot(dot_plot)

                
                de_df = sc.get.rank_genes_groups_df(adata_m, group=m_label).head(10)
                

                mgenes            = de_df['names'].tolist()
                macrogene         = st.sidebar.selectbox( "2.3 Top genes of macrogene ", mgenes, index=0, placeholder="Please select ...")
                weight_file       = "./demo/frog_zebrafish/test256_data_frog.3w_zebrafish.3w_org_saturn_seed_0_genes_to_macrogenes.pkl"
                macrogene_weights = readPKL(weight_file)

                ## 根据 gene x macrogene 获得 每个macrogene对应的gene scores
                def get_scores(macrogene):
                    '''
                    Given the index of a macrogene, return the scores by gene for that centroid
                    '''
                    scores = {}
                    for (gene), score in macrogene_weights.items():
                        scores[gene] = score[int(macrogene)]
                    return scores

                df_gene = pd.DataFrame(get_scores(macrogene).items(), columns=["gene", "weight"])\
                    .sort_values("weight", ascending=False)\
                    .head(10)

                col3, col4  = st.columns([1, 1])

                col3.markdown( "Top macrogenes for **{}**".format(m_label) )
                col3.dataframe(de_df, hide_index=True)
                col4.markdown( "Top genes of macrogene **{}** ".format(macrogene) )
                col4.dataframe(df_gene, hide_index=True)




if __name__ == "__main__":
    main()
