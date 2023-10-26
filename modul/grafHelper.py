import networkx as nx
import matplotlib.pyplot as plt

def _set_networkx_graph(df_node, df_edge):
    # parameter --------------------------
    # df_node: dataframe dari node
    # df_edge: dataframe dari edge
    # output -----------------------------
    # return: networkx graph
    # ====================================

    #konversi graph 
    gnx = nx.MultiDiGraph()
    #node
    for i,a in df_node.iterrows():
        #mulai disini akan digunakan taksonomi bahasa indonesia pada data.
        gnx.add_node(
            a['taxon_id'],
            label=a['taxon_name'],
            superkingdom=a['superkingdom'],
            kingdom=a['kingdom'],
            filum=a['phylum'],
            kelas=a['class'],
            ordo=a['order'],
            famili=a['family'],
            genus=a['genus'],
            spesies=a['species'],
            group=a['group'],
            color=a['color'],
        )
    #edge
    for i,a in df_edge.iterrows():
        gnx.add_edge(
            a['source_taxon_id'],
            a['target_taxon_id'],
            label=a['interaction_type'],
        )
    return gnx



def _plot_nx_by_matplotlib(G):
    # parameter input nx graph
    # no return cuma tampilan, visualisasi graf


    fig, ax = plt.subplots(figsize=(20, 20))

    # Generate layout for visualization
    # pos = nx.kamada_kawai_layout(G)
    # pos = nx.spring_layout(G)
    pos = nx.nx_agraph.graphviz_layout(G, prog="neato", args="")

    # Visualize graph components
    nx.draw_networkx_edges(G, pos, alpha=0.3, edge_color='g')
    nx.draw_networkx_nodes(G, pos, node_color=list(nx.get_node_attributes(G, "color").values()), alpha=0.9)

    #node label
    # for i in ['#b22222','#671f92','#1f922b','#EADDCA']: # filtering dengan bedakan warna node
    #     label_options = {"ec": i, "fc": 'white', "alpha": 0.7}
    #     nx.draw_networkx_labels(
    #         nx.subgraph_view(G, filter_node=lambda n1: G.nodes(data=True)[n1].get("color", True) == i),
    #         pos, 
    #         font_size=10, 
    #         bbox=label_options
    #     )

    #edge labels
    edge_labels={x:i for i,x in zip(nx.get_edge_attributes(G, "label").values(),G.edges())}
    nx.draw_networkx_edge_labels(G, pos=pos, edge_labels=edge_labels)


    # Title/legend
    font = {"fontname": "Helvetica", "color": "k", "fontweight": "bold", "fontsize": 14}
    ax.set_title("Graf Interaksi Serangga<->virus<->tanaman", font)
    # Change font color for legend
    font["color"] = "r"

    ax.text(
        0.80,
        0.10,
        "hijau = Tanaman",
        horizontalalignment="center",
        transform=ax.transAxes,
        fontdict=font,
    )
    ax.text(
        0.80,
        0.08,
        "merah = Serangga",
        horizontalalignment="center",
        transform=ax.transAxes,
        fontdict=font,
    )

    ax.text(
        0.80,
        0.06,
        "ungu = Virus",
        horizontalalignment="center",
        transform=ax.transAxes,
        fontdict=font,
    )

    ax.text(
        0.80,
        0.04,
        "abu-abu = Nogroup",
        horizontalalignment="center",
        transform=ax.transAxes,
        fontdict=font,
    )

    # Resize figure for label readibility
    ax.margins(0.1, 0.05)
    fig.tight_layout()
    plt.axis("off")
    plt.show()