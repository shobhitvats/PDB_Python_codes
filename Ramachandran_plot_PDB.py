import matplotlib.pyplot as plt
import numpy as np
import imports


def graph(pdb_angles_file , protein_id , plot_path):

    phi_psi_list = np.loadtxt(pdb_angles_file)

    plt.figure(1, figsize=(3, 3))

    # plot function
    plt.scatter(
        phi_psi_list[:, 0],
        phi_psi_list[:, 1],
        s=2,
        color="red",
        marker=".",
    )

    # set the plot title
    plt.title("Ramachandran Plot of " + protein_id)

    # set the limitation of axes
    plt.xlim((-np.pi, np.pi))
    plt.ylim((-np.pi, np.pi))

    # set the name of axes
    plt.xlabel("phi")
    plt.ylabel("psi")

    plt.savefig(plot_path)
    plt.clf()



for protein_id in imports.pdb_list:
    graph('./angles/' + protein_id + '.angles', protein_id , './plots/' + protein_id)
