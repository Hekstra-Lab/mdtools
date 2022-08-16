from simtk.unit import *
import mdtraj
import matplotlib.pyplot as plt
import numpy as np


def getFieldStrength(e):
    """
    Convert from given unit of electric field strength to a value
    in OpenMM's standard unit
    """
    def convert(v):
        return (v*AVOGADRO_CONSTANT_NA).value_in_unit(kilojoule_per_mole / elementary_charge / nanometer)

    if isinstance(e, list):
        return [convert(e_c) for e_c in e]
    else:
        return convert(e)


def compute_RMSD_by_chain(target_traj, top, rule, chain_id_list = None):
    chains_rmsd = []
    n_frames = target_traj.n_frames
    if chain_id_list == None:
        chain_id_list = range(top.n_chains)
    for chain_id in chain_id_list:
        atoms_selection = top.select("chainid " + str(chain_id) + " and " + rule)
        if(len(atoms_selection) == 0): # no valid atoms bound to a.a. residues
            continue
        chain_rmsd = mdtraj.rmsd(target_traj, target_traj, 0, atoms_selection)
        chains_rmsd.append(chain_rmsd)

    chains_rmsd = np.array(chains_rmsd)
    return chains_rmsd, np.mean(chains_rmsd, axis = 0), np.std(chains_rmsd, axis = 0)

def plot_RMSD_by_chain(chains_rmsd, rule, n_frames):
    plt.figure(figsize=(18, 12))
    plt.rc('font', size=16)
    plt.title("PDZ domain RMSD by chain (rule: {})".format(rule))
    plt.xlabel("Time (ns)")
    plt.xticks(np.arange(0, n_frames, 100), np.arange(0, n_frames,100) / 10)
    plt.ylabel("RMSD (nm)")
    for chid, chain_rmsd in enumerate(chains_rmsd):
        plt.plot(chain_rmsd, label = "chain " + str(chid))

    plt.legend(loc = 'upper right')


def plot_RMSD_by_chain_stat(avg_rmsd, std_rmsd, rule, n_frames):
    plt.figure(figsize=(18, 12))
    plt.rc('font', size=16)
    plt.title("PDZ domain RMSD by chain (rule: {})".format(rule))
    plt.xlabel("Time (ns)")
    plt.xticks(np.arange(0, n_frames, 100), np.arange(0, n_frames,100) / 10)
    plt.ylabel("RMSD (nm)")
    plt.plot(avg_rmsd, label = "averaged")
    plt.plot(avg_rmsd + std_rmsd, label = "+1 sd")
    plt.plot(avg_rmsd - std_rmsd, label = "-1 sd")
    plt.legend(loc = 'upper right')