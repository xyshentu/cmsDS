################################################################
# Version : 0.1
# Author: Xinyi Shentu
# Date: 2024-11-20
################################################################

import pandas as pd
import plotly.express as px
import scanpy as sc
import squidpy as sq
import warnings
import ipywidgets as widgets
from IPython.display import display
from matplotlib import rc_context
import matplotlib.pyplot as plt


class InteractH5ad:
    def __init__(self, adata):
        self.data = adata.copy()

        if "counts" not in self.data.layers:
            self.data.layers["counts"] = self.data.X.copy()

        self.is_pp = False
        self.is_dr = False

    def copy(self):
        new_instance = InteractH5ad(self.data)
        new_instance.is_dr = self.is_dr
        new_instance.is_pp = self.is_pp
        return new_instance

    def fun_pp(
        self,
        use_hvg_for_pca=False,
        n_top_genes=8000,
        flavor="seurat",
        layer=None,
        batch_key=None,
    ):
        """

        Step 1.
        Preprocessing with your raw adata (:)
        Return to the raw count by `adata.X = adata.layers['counts'].copy()`, for whenever you want ! (^ - ^)

        """

        data = self.data.copy()
        data.X = data.layers["counts"].copy()

        # sc.pp.calculate_qc_metrics(data, log1p = False, inplace=True)

        sc.pp.normalize_total(data)
        sc.pp.log1p(data)

        sc.pp.highly_variable_genes(
            data,
            n_top_genes=n_top_genes,
            flavor=flavor,
            layer=None,
            batch_key=batch_key,
        )

        sc.pp.pca(data, mask_var="highly_variable" if use_hvg_for_pca else None)

        self.data = data
        self.is_pp = True

    # def fun_remove_batch(self, batch_key=None, obsm_key="X_pca"):
    #     """

    #     (Optional): Remove Batch Effect

    #     """

    #     import harmonypy as hm

    #     data = self.data

    #     out = hm.run_harmony(data.obsm[obsm_key], adata.obs, batch_key)
    #     data.obsm[f"{obsm_key}_hm"] = out.Z_corr.T

    #     self.data = data

    def fun_reduce_dimension(self, n_neighbors=15, n_pcs=None, use_rep=None):
        """

        Step 2.
        After data preprocessing,
        Reduce the dimensions !

        """
        if not self.is_pp:
            print(
                "It seems the input data hasn't been prepreprocess? Please run `vis_pp(your_data)` to do it ! "
            )
            return

        data = self.data

        if use_rep is None:
            use_rep_val = "X_pca"
        else:
            use_rep_val = use_rep

        # Find neighbor
        warnings.warn(f"Find neighbor based on : {use_rep_val} (Step 1/2)")

        sc.pp.neighbors(data, n_neighbors=n_neighbors, n_pcs=n_pcs, use_rep=use_rep_val)

        warnings.warn(f"Calculate UMAP (Step 2/2)")

        sc.tl.umap(data)

        print("Done Dimension Reduction !")

        self.data = data
        self.is_dr = True

    def fun_cluster(
        self, method, resolution=None, restrict_to=None, flavor=None, **kwargs
    ):
        """

        Step 3.
        After dimension reduction,
        Get the variable cluster !

        method : [leiden , louvain]

        """
        if method not in ["leiden", "louvain"]:
            print("Please choose method within [leiden , louvain] ")
            return

        data = self.data

        if method == "leiden":
            res = 1 if resolution is None else resolution
            sc.tl.leiden(
                data,
                resolution=res,
                restrict_to=restrict_to,
                flavor="leidenalg" if flavor is None else flavor,
                key_added=f"{method}_res{res}",
            )
        else:
            res = resolution
            sc.tl.louvain(
                data,
                resolution=res,
                restrict_to=restrict_to,
                flavor="vtraag" if flavor is None else flavor,
                key_added=f"{method}_res{res}",
            )

        self.data = data

    def fun_extract_cells(self, query):

        qdf = self.data.obs.query(query)

        qdata = self.copy()
        qdata.data = qdata.data[qdf.index, :].copy()

        return qdata

    def fun_extract_genes(self, genelist):

        qdata = self.copy()
        qdata.data = qdata.data[:, qdata.data.var_names.isin(genelist)].copy()

        return qdata

    def fun_get_deg(
        self,
        groupby,
        mask_var=None,
        groups="all",
        reference="rest",
        pts=True,
        method=None,
        pval_cutoff=None,
        log2fc_min=None,
        log2fc_max=None,
    ):
        """

        Rank DEGs by group.
        This function is based on `scanpy.tl.rank_genes_groups() scanpy.get.rank_genes_groups_df()`

        method : The default method is `'t-test'`,
        `'t-test_overestim_var'` overestimates variance of each group,
        `'wilcoxon'` uses Wilcoxon rank-sum,
        `'logreg'` uses logistic regression. See :cite:t:`Ntranos2019`

        """
        data = self.data
        sc.tl.rank_genes_groups(
            data,
            groupby=groupby,
            mask_var=mask_var,
            groups=groups,
            reference=reference,
            pts=pts,
            method=method,
            use_raw = False
        )
        df = sc.get.rank_genes_groups_df(
            data,
            group=None,
            pval_cutoff=pval_cutoff,
            log2fc_min=log2fc_min,
            log2fc_max=log2fc_max,
        )
        return df

    def fun_output(self):
        """
        Here is your final data! May you have a good day (^o^)
        """
        return self.data





def vis_pp(interactor):
    """
    Use plotly to interactly do preprocessing!
    """

    candi_batchkeys = interactor.data.obs.columns.tolist().copy()
    candi_batchkeys.append(None)
    
    #############################################
    # Set the widgets
    #############################################
    # Create interactive widgets
    n_top_genes_slider = widgets.IntSlider(
        value=2000, min=1000, max=20000, step=100, description="n_top_genes"
    )
    use_hvg_for_pca_checkbox = widgets.Checkbox(
        value=False, description="Use HVG for PCA"
    )
    batchkey_dropdown = widgets.Dropdown(
        options=candi_batchkeys,
        value=None,
        description="BatchKey HVG",
    )
    flavor_dropdown = widgets.Dropdown(
        options=["seurat", "cell_ranger", "seurat_v3", "seurat_v3_paper"],
        value="seurat",
        description="Flavor",
    )
    run_button = widgets.Button(description="Run Preprocessing")
    progress_message = widgets.HTML(value="")

    #############################################
    # Bind function to widgets
    #############################################
    def on_button_clicked(b):
        progress_message.value = '<h3>Start</h3>'
        interactor.fun_pp(
            use_hvg_for_pca=use_hvg_for_pca_checkbox.value,
            n_top_genes=n_top_genes_slider.value,
            flavor=flavor_dropdown.value,
            batch_key=batchkey_dropdown.value,
        )
        progress_message.value = '<h3 style="color: blue;">Done</h3>'
    run_button.on_click(on_button_clicked)
    
    #############################################
    # UI
    #############################################
    display(
        widgets.HTML(
            value="<h1>Step 1: Preprocessing your anndata</h1>\n<h2>Normalization, Logrize, PCA</h2>"
        ),
        widgets.HBox(
            children=[
                widgets.HTML("<h3>HVG</h3>"),
                widgets.Box(layout=widgets.Layout(width="50px")),
                widgets.VBox(
                    [
                        n_top_genes_slider,
                        batchkey_dropdown,
                        flavor_dropdown,
                    ]
                ),
            ]
        ),
        widgets.HBox(
            children=[
                widgets.HTML("<h3>PCA</h3>"),
                widgets.Box(layout=widgets.Layout(width="50px")),
                widgets.VBox(
                    [widgets.HTML(value="<br>"),use_hvg_for_pca_checkbox],
                ),
            ]
        ),
        widgets.HTML(value="<br>"),
        run_button,
        progress_message,
    )





def vis_plotSpatial(interactor):

    options = list(
        set(interactor.data.var_names).union(set(interactor.data.obs.columns))
    )

    #############################################
    # Set the widgets
    #############################################
    output = widgets.Output()
    progress_message = widgets.HTML(value="")
    plot_spatial_button = widgets.Button(description="Plot on Spatial")
    color_inputbox = widgets.Text(
        value="total_counts",
        placeholder="Input a value to show on spatial",
        description="Input:",
    )
    color_inputdropdown = widgets.Dropdown(
        options=["total_counts"],
        value=None,
        description="Suggestions:",
    )
    vmin_checkbox = widgets.Checkbox(value=False, description="Enable vmin")
    vmax_checkbox = widgets.Checkbox(value=False, description="Enable vmax")
    vmin_inputbox = widgets.FloatText(value=0, description="vmin")
    vmax_inputbox = widgets.FloatText(value=0, description="vmax")

    cmap_dropdown = widgets.Dropdown(
        options=plt.colormaps(),
        value=None,
        description="cmap",
    )

    #############################################
    # Bind function to widgets
    #############################################
    def on_input_change(change):
        input_value = change["new"]
        if input_value:
            filtered_options = [
                option for option in options if input_value.lower() in option.lower()
            ]
            color_inputdropdown.options = filtered_options[:10]  # 只显示前10个匹配项
        else:
            color_inputdropdown.options = []

    color_inputbox.observe(on_input_change, names="value")

    def on_dropdown_select(change):
        if change["new"]:
            color_inputbox.value = change["new"]

    color_inputdropdown.observe(on_dropdown_select, names="value")

    def on_plot_spatial_button_click(b):
        with output:
            output.clear_output()

            transparent_color = (1.0, 1.0, 1.0, 0.0)
            with rc_context(
                {
                    "figure.facecolor": transparent_color,  # 图形背景颜色，透明度 0%
                    "axes.facecolor": transparent_color,  # 坐标轴背景颜色，透明度 0%
                    "savefig.facecolor": transparent_color,  # 保存图形时的背景颜色，透明度 0%
                }
            ):
                sq.pl.spatial_scatter(
                    interactor.data,
                    color=color_inputbox.value,
                    shape=None,
                    library_id="spatial",
                    cmap=cmap_dropdown.value,
                    edges_color=None,
                    vmin=vmin_inputbox.value if vmin_checkbox.value else None,
                    vmax=vmax_inputbox.value if vmax_checkbox.value else None,
                )
                plt.axis("off")
                plt.show()
                progress_message.value = "<h3>Done</h3>"

    plot_spatial_button.on_click(on_plot_spatial_button_click)

    #############################################
    # UI
    #############################################
    display(
        widgets.HTML(value="<h1>Plot: Spatial</h1>\n"),
        widgets.HBox(
            children=[
                widgets.VBox(
                    [
                        color_inputbox,
                        color_inputdropdown,
                        widgets.VBox(
                            children=[vmin_checkbox, vmin_inputbox],
                        ),
                        widgets.VBox(
                            children=[vmax_checkbox, vmax_inputbox],
                        ),
                    ],
                ),
                cmap_dropdown,
                plot_spatial_button,
            ]
        ),
        output,
    )

def vis_dr(interactor):

    candi_repkeys = list(interactor.data.obsm.keys())
    candi_repkeys = [x for x in candi_repkeys if x.startswith("X_")]

    #############################################
    # Set the widgets
    #############################################
    n_neighbor_slider = widgets.IntSlider(
        value=15, min=3, max=50, step=1, description="n_neighbors"
    )
    n_pcs_slider = widgets.IntSlider(
        value=50, min=1, max=100, step=1, description="n_pcs"
    )
    repkey_dropdown = widgets.Dropdown(
        options=candi_repkeys,
        value=candi_repkeys[0] if candi_repkeys else None,  # 设置初始值
        description="RepKey",
    )
    dr_button = widgets.Button(description="Reducing Dimension")
    progress_message = widgets.HTML(value="")

    #############################################
    # Bind function to widgets
    #############################################

    def on_dr_button_clicked(b):
        progress_message.value = "<h3>Start</h3>"
        print("Button clicked!")  # 调试信息
        try:
            interactor.fun_reduce_dimension(
                n_neighbors=n_neighbor_slider.value,
                n_pcs=n_pcs_slider.value,
                use_rep=repkey_dropdown.value,
            )
            progress_message.value = '<h3 style="color: blue;">Done</h3>'
        except Exception as e:
            progress_message.value = f"<h3 style='color: red;'>Error: {e}</h3>"
            print(e)  # 打印异常信息
    dr_button.on_click(on_dr_button_clicked)

    #############################################
    # UI
    #############################################
    display(
        widgets.HTML(
            value="<h1>Step 2: Dimension Reduction</h1>\n<h2>Find neighbor, UMAP</h2>"
        ),
        widgets.HBox(
            children=[
                widgets.HTML("<h3>Neighbor</h3>"),
                widgets.Box(layout=widgets.Layout(width="50px")),
                widgets.VBox(
                    children=[n_neighbor_slider, n_pcs_slider, repkey_dropdown]
                ),
            ]
        ),
        dr_button,
        progress_message,
    )

def vis_cluster(interactor):

    method_flavors = {
        "leiden": ["leidenalg", "igraph"],
        "louvain": ["vtraag", "igraph", "rapids"],
    }

    #############################################
    # Set the widgets
    #############################################
    method_dropdown = widgets.Dropdown(
        options=["leiden", "louvain"],
        value="leiden",
        description="Method",
    )
    flavor_dropdown = widgets.Dropdown(
        options=method_flavors[method_dropdown.value],
        value=method_flavors[method_dropdown.value][0],
        description="Flavor",
    )
    res_slider = widgets.FloatSlider(
        value=1.0, min=0.0, max=10.0, step=0.1, description="Resolution"
    )
    cluster_button = widgets.Button(description="Clustering")
    progress_message = widgets.HTML(value="")

    #############################################
    # Bind function to widgets
    #############################################
    def on_method_select(change):
        dp_value = change["new"]
        if dp_value:
            flavor_dropdown.options = method_flavors[dp_value]
            flavor_dropdown.values = flavor_dropdown.options[0]

    method_dropdown.observe(on_method_select, names="value")

    def on_cluster_button_clicked(b):
        progress_message.value = "<h3>Start</h3>"
        print("Button clicked!")  # 调试信息
        try:
            interactor.fun_cluster(
                method=method_dropdown.value,
                flavor=flavor_dropdown.value,
                resolution=res_slider.value,
            )
            progress_message.value = '<h3 style="color: blue;">Done</h3>'
        except Exception as e:
            progress_message.value = f"<h3 style='color: red;'>Error: {e}</h3>"
            print(e)  # 打印异常信息

    cluster_button.on_click(on_cluster_button_clicked)

    #############################################
    # UI
    #############################################
    display(
        widgets.HTML(
            value="<h1>Step 3: Get clustering</h1>\n<h2>Method: Leiden or Louvain</h2>"
        ),
        method_dropdown,
        flavor_dropdown,
        res_slider,
        cluster_button,
        progress_message,
    )

def vis_operate_data(interactor, global_ns):

    #############################################
    # Set the widgets
    #############################################
    query_inputbox = widgets.Text(
        value="total_counts > 3000",
        placeholder="How to select cells",
        description="Cells:",
    )
    gene_inputbox = widgets.Text(
        value=None, placeholder="Give me the name of genelist", description="Genes:"
    )
    output_inputbox = widgets.Text(
        value="sub_interactor", placeholder="newname", description="Output variable:"
    )
    operate_button = widgets.Button(description="Operating")
    progress_message = widgets.HTML(value="")

    #############################################
    # Bind function to widgets
    #############################################
    
    def on_operate_button_clicked(b):
        progress_message.value = "<h3>Start</h3>"
        try:

            tmp = interactor

            if query_inputbox.value != "":
                tmp = tmp.fun_extract_cells(query_inputbox.value)
            
            if gene_inputbox.value != "":
                if gene_inputbox.value in global_ns.keys():
                    tmp = tmp.fun_extract_genes(global_ns[gene_inputbox.value])
                    progress_message.value = '<h3 style="color: blue;">Done</h3>'
                else:
                    warnings.warn(f"{gene_inputbox.value} is not exist")
            global_ns[output_inputbox.value] = tmp
            

        except Exception as e:
            progress_message.value = f"<h3 style='color: red;'>Error: {e}</h3>"
            print(e)  # 打印异常信息

        if output_inputbox.value in global_ns:
            progress_message.value = f"The extracted cells x genes shape is : {(global_ns[output_inputbox.value]).data.shape}"
        else:
            progress_message.value = f"<h3 style='color: red;'>Error: {output_inputbox.value} is not defined in the global namespace.</h3>"

    operate_button.on_click(on_operate_button_clicked)

    #############################################
    # UI
    #############################################
    display(
        widgets.HTML(value="<h1>Operating data</h1>"),
        query_inputbox, gene_inputbox, output_inputbox, operate_button, progress_message
    )


def vis_deg(interactor):

    groupby_list = interactor.data.obs.columns.tolist()

    #############################################
    # Set the widgets
    #############################################
    groupby_dropdown = widgets.Dropdown(
        options=groupby_list,
        value=None,
        description="Groupby:",
    )
    method_dropdown = widgets.Dropdown(
        options=["t-test", 't-test_overestim_var', 'wilcoxon', 'logreg'],
        value="t-test",
        description="Method:",
    )
    groups_inputbox = widgets.Text(
        value="all",
        placeholder="Group to find marker",
        description="Groups:",
    )
    reference_inputbox = widgets.Text(
        value="rest",
        placeholder="Contrast to",
        description="Reference:",
    )
    pval_inputbox = widgets.FloatText(value=0.05, description="p.adj cutoff:")
    log2fcmin_inputbox = widgets.FloatText(value=1, description="log2FC_min:")
    savefile_inputbox = widgets.Text(
        value="markers.csv",
        placeholder="Path to outout.csv",
        description="Save to:",
        layout=widgets.Layout(flex="flex-grow"),
    )

    deg_button = widgets.Button(description="Get DEG")
    progress_message = widgets.HTML(value="")

    #############################################
    # Bind function to widgets
    #############################################
    def on_deg_button_clicked(b):
        progress_message.value = "<h3>Start</h3>"
        print("Button clicked!")  # 调试信息
        try:
            df = interactor.fun_get_deg(
                groupby=groupby_dropdown.value,
                groups=groups_inputbox.value,
                reference=reference_inputbox.value,
                method=method_dropdown.value,
                pval_cutoff=pval_inputbox.value,
                log2fc_min=log2fcmin_inputbox.value,
            )
            df.to_csv(savefile_inputbox.value)
            progress_message.value = '<h3 style="color: blue;">Done</h3>'
        except Exception as e:
            progress_message.value = f"<h3 style='color: red;'>Error: {e}</h3>"
            print(e)  # 打印异常信息

    deg_button.on_click(on_deg_button_clicked)

    #############################################
    # UI
    #############################################
    display(
        widgets.HTML(value="<h1>Get DEG result</h1>"),
        widgets.HBox(
            [
                widgets.HTML(
                    value="<h3>Rank genes</h3>Please check <em>groupby</em> carefully !"
                ),
                widgets.Box(layout=widgets.Layout(width="50px")),
                widgets.VBox(
                    [
                        groupby_dropdown,
                        groups_inputbox,
                        reference_inputbox,
                        method_dropdown,
                    ],
                    layout=widgets.Layout(flex="flex-shrink"),
                ),
            ]
        ),
        widgets.HTML(value="<hr>"),
        widgets.HBox(
            [
                widgets.HTML(value="<h3>DEG cutoff</h3>Save the deg result to .csv"),
                widgets.Box(layout=widgets.Layout(width="90px")),
                widgets.VBox(
                    [pval_inputbox, log2fcmin_inputbox],
                    layout=widgets.Layout(flex="flex-shrink"),
                ),
            ]
        ),
        widgets.HTML(value="<hr>"),
        savefile_inputbox,
        deg_button,
        progress_message,
    )






def vis_plotUMAP(interactor):

    options = list(
        set(interactor.data.var_names).union(set(interactor.data.obs.columns))
    )

    output = widgets.Output()

    # Create a button to plot UMAP
    plot_umap_button = widgets.Button(description="Plot UMAP")

    color_inputbox = widgets.Text(value = 'total_counts', description = "Input")
    
    color_dropdown = widgets.Dropdown(
        options=['total_counts'],
        value=None,
        description="Color",
    )

    def on_input_change(change):
        input_value = change["new"]
        if input_value:
            filtered_options = [
                option for option in options if input_value.lower() in option.lower()
            ]
            color_dropdown.options = filtered_options[:10]  # 只显示前10个匹配项
        else:
            color_dropdown.options = []

    color_inputbox.observe(on_input_change, names="value")

    def on_dropdown_select(change):
        if change["new"]:
            color_inputbox.value = change["new"]

    color_dropdown.observe(on_dropdown_select, names="value")

    def on_plot_umap_button_click(b):
        with output:
            output.clear_output()

            transparent_color = (1.0, 1.0, 1.0, 0.0)
            with rc_context(
                {
                    "figure.facecolor": transparent_color,  # 图形背景颜色，透明度 0%
                    "axes.facecolor": transparent_color,  # 坐标轴背景颜色，透明度 0%
                    "savefig.facecolor": transparent_color,  # 保存图形时的背景颜色，透明度 0%
                }
            ):
                sc.pl.umap(
                    interactor.data,
                    color=color_inputbox.value,
                    frameon=False,
                    legend_loc="on data",
                    title=color_dropdown.value,
                )

    plot_umap_button.on_click(on_plot_umap_button_click)

    display(
        widgets.HTML(value="<h1>Plot: UMAP</h1>\n"),
        color_inputbox,
        color_dropdown,
        plot_umap_button,
        output,
    )

def vis_dot(interactor, global_ns):

    output = widgets.Output()

    #############################################
    # Set the widgets
    #############################################

    plot_button = widgets.Button(description="Plot dotplot")

    gene_inputbox = widgets.Text(
        value=None, placeholder="Give me the name of genelist", description="Genes:"
    )
    
    groupby_dropdown = widgets.Dropdown(
        options=interactor.data.obs.columns.tolist(),
        value=None,
        description="Groupby:",
    )
    scale_dropdown = widgets.Dropdown(
        options=[None, 'var', 'group'],
        value=None,
        description="Scale:",
    )

    #############################################
    # Binding
    #############################################

    def on_plot_button_click(b):
        with output:
            output.clear_output()

            if gene_inputbox.value != "":
                if gene_inputbox.value in global_ns.keys():
                    gl = interactor.data.var_names[ interactor.data.var_names.isin(global_ns[gene_inputbox.value]) ]


                    transparent_color = (1.0, 1.0, 1.0, 0.0)
                    with rc_context(
                        {
                            "figure.facecolor": transparent_color,  # 图形背景颜色，透明度 0%
                            "axes.facecolor": transparent_color,  # 坐标轴背景颜色，透明度 0%
                            "savefig.facecolor": transparent_color,  # 保存图形时的背景颜色，透明度 0%
                        }
                    ):
                        sc.pl.dotplot(
                                    interactor.data,
                                    var_names=gl,
                                    groupby=groupby_dropdown.value,
                                    standard_scale=scale_dropdown.value,
                                    var_group_rotation=90,
                                    cmap="bwr",
                                )

    plot_button.on_click(on_plot_button_click)

    #############################################
    # UI
    #############################################
    display(
        widgets.HTML(value="<h1>Plot: Dotplot</h1>\n"),
        gene_inputbox,
        groupby_dropdown,
        scale_dropdown,
        plot_button,
        output,
    )

def vis_violin(interactor):

    options = list(
        set(interactor.data.var_names).union(set(interactor.data.obs.columns))
    )

    output = widgets.Output()

    #############################################
    # Set the widgets
    #############################################

    plot_button = widgets.Button(description="Plot violin")
    color_inputbox = widgets.Text(value="total_counts", description="Input:")
    color_dropdown = widgets.Dropdown(
        options=["total_counts"],
        value=None,
        description="Suggestion:",
    )
    groupby_dropdown = widgets.Dropdown(
        options=interactor.data.obs.columns.tolist(),
        value=None,
        description="Groupby:",
    )

    #############################################
    # Binding
    #############################################

    def on_input_change(change):
        input_value = change["new"]
        if input_value:
            filtered_options = [
                option for option in options if input_value.lower() in option.lower()
            ]
            color_dropdown.options = filtered_options[:10]  # 只显示前10个匹配项
        else:
            color_dropdown.options = []

    color_inputbox.observe(on_input_change, names="value")

    def on_dropdown_select(change):
        if change["new"]:
            color_inputbox.value = change["new"]

    color_dropdown.observe(on_dropdown_select, names="value")

    def on_plot_button_click(b):
        with output:
            output.clear_output()

            transparent_color = (1.0, 1.0, 1.0, 0.0)
            with rc_context(
                {
                    "figure.facecolor": transparent_color,  # 图形背景颜色，透明度 0%
                    "axes.facecolor": transparent_color,  # 坐标轴背景颜色，透明度 0%
                    "savefig.facecolor": transparent_color,  # 保存图形时的背景颜色，透明度 0%
                }
            ):
                sc.pl.violin(
                    interactor.data,
                    groupby=groupby_dropdown.value,
                    keys=color_inputbox.value,
                    rotation=90,
                    size=0
                )

    plot_button.on_click(on_plot_button_click)

    #############################################
    # UI
    #############################################
    display(
        widgets.HTML(value="<h1>Plot: Violin</h1>\n"),
        groupby_dropdown,
        color_inputbox,
        color_dropdown,
        plot_button,
        output,
    )