import plotly.graph_objects as go
import networkx as nx
import matplotlib.pyplot as plt

# cuma tampilan
def pie_proporsi(df_node):
    data = df_node.groupby(['group','color']).agg({'group': ['count'], }).reset_index().sort_values(
        ('group', 'count'),ascending=False
    ).reset_index(drop=True).values
    labels = [i[0] for i in data]
    colors = [i[1] for i in data]
    slices = [i[2] for i in data]

    fig = go.Figure(data=[go.Pie(labels=labels,values=slices)])
    fig.update_traces(hoverinfo='label+percent', textinfo='value+percent', textfont_size=20, marker=dict(colors=colors, line=dict(color='#000000', width=0.1)))
    return fig.to_html(full_html=False) #, include_plotlyjs='cdn'


def networkxGraph(df_node, df_edge):
    #3
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

    # cuma tampilan, visualisasi graf
    G=gnx

    fig, ax = plt.subplots(figsize=(20, 20))

    # Generate layout for visualization
    # pos = nx.kamada_kawai_layout(G)
    # pos = nx.spring_layout(G)
    pos = nx.nx_agraph.graphviz_layout(G, prog="neato", args="")

    # Visualize graph components
    nx.draw_networkx_edges(G, pos, alpha=0.3, edge_color='g')
    nx.draw_networkx_nodes(G, pos, node_color=list(nx.get_node_attributes(G, "color").values()), alpha=0.9)

    #edge labels
    edge_labels={x:i for i,x in zip(nx.get_edge_attributes(G, "label").values(),G.edges())}
    nx.draw_networkx_edge_labels(G, pos=pos, edge_labels=edge_labels)


    # Title/legend
    font = {"fontname": "Helvetica", "color": "k", "fontweight": "bold", "fontsize": 14}
    ax.set_title("Interaksi Tanaman-Serangga-Virus", font)
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