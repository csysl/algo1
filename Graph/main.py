import numpy as np
import matplotlib.pyplot as plt
import networkx as nx

if __name__ == "__main__":
    G = nx.DiGraph()
    u, v, w = np.loadtxt("C:/Users/92582/Desktop/Project1/Project1/GstFlow.txt",
                         delimiter=' ', unpack=True)

    k = int(u[len(u) - 1])
    q = int(v[len(u) - 1])

    # get locations of notes
    locations = []
    for i in range(k):
        for j in range(q - 1, -1, -1):
            locations.append([(2 * (i + 1) - 1), (j + 1)])
            locations.append([(2 * (i + 1)), (j + 1)])

    # notes list
    notes = []
    for i in range(2 * k * q + 2):
        notes.append(i + 1)

    # edges list
    edges = []
    for i in range(len(u) - 1):
        edges.append((u[i], v[i], w[i]))

    # zip the note with its location
    pos = dict(zip(notes, locations))
    st_dict = {(2 * k * q + 1): [0, q / 2], (2 * k * q + 2): [2 * k + 1, q / 2]}
    pos.update(st_dict)

    # draw notes
    nx.draw_networkx_nodes(G, pos, notes, node_size=1, node_color='black', with_labels=True)

    # draw edges
    nx.draw_networkx_edges(G, pos, edges, edge_color='b', arrows=False)

    # get edges of a path
    tag = np.loadtxt("C:/Users/92582/Desktop/Project1/Project1/f.txt", unpack=True)
    f_edges = []
    for i in range(len(tag)):
        f_edges.append((2 * tag[i] - 1, 2 * tag[i]))

    # draw the path
    nx.draw_networkx_edges(G, pos, f_edges, edge_color='r', arrows=False)

    # draw the graph
    nx.draw(G)
    plt.title("stFlow", fontsize=18)
    plt.savefig("stFlow.jpg", dpi=1000)
    #plt.subplots_adjust(left=0.09, right=1, wspace=0.25, hspace=0.25, bottom=0.13, top=0.91)
    plt.show()
    print("stFlow set!")
