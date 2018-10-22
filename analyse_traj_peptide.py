import os
import sys
import mdtraj
import numpy as np
from  itertools import *
import subprocess
import argparse
from  matplotlib import pyplot as  plt



def process_traj(topol, traj):

    print( '.... processing traj .... ')

    process_traj='\
    echo 1 1 | gmx trjconv -f %s -pbc cluster -o traj_cluster.xtc -s md.tpr ; \n \
    echo 1 1 | gmx trjconv -f traj_cluster.xtc  -o traj_cluster_alg.xtc -s md.tpr -fit rot+trans  ;  \n \
    echo 1 1 | gmx trjconv -f %s -pbc cluster -o traj_cluster.pdb -s md.tpr ; \n'  %(traj, topol)

    process = subprocess.check_output(process_traj,shell=True )
    print ('.... processing traj .... DONE')

def helicity(traj,peptide_chain):
     print( '.... computing helicity ....')
     dssp=mdtraj.compute_dssp(traj, simplified=True)
     residues=[residue.index for residue in traj.topology.chain(peptide_chain).residues ]

     unique, counts = np.unique(dssp[:,residues[0]:residues[-1]], return_counts=True)
     hel = dict(zip(unique, counts)).get('H',0)/np.sum(counts)*10
     print('helicity ' + str(hel))

     print( '.... computing helicity .... DONE')
     return hel

def find_centroid(traj):
    print( '.... finding centroid ....')
    atom_indices = [a.index for a in traj.topology.atoms if a.name == 'CA']
    distances = np.empty((traj.n_frames, traj.n_frames))
    for i in range(traj.n_frames):
      distances[i] = mdtraj.rmsd(traj, traj, i, atom_indices=atom_indices)
    beta = 1
    index = np.exp(-beta*distances / distances.std()).sum(axis=1).argmax()
    print( '.... finding centroid .... DONE')
    return index

def compute_rmsd(traj, peptide_chain, centroid):
    print( '.... computing RMSD ....')
#  atoms_indices = traj.topology.chain("%s" %(peptide_chain))
    atoms_indices = traj.topology.select("chainid  %s" %(peptide_chain))

    rmsd = mdtraj.rmsd(traj, traj, frame=centroid, atom_indices=atoms_indices, parallel=True, precentered=True)

    print( '.... computing RMSD .... DONE')

    return rmsd

def compute_rmsd_sidechains(traj, centroid):
    print( '.... computing RMSD ....')
#  atoms_indices = traj.topology.chain("%s" %(peptide_chain))
    #atoms_indices = traj.topology.select("chainid  3 to 5 and sidechain" )
    atoms_indices = traj.topology.select("chainid  3 and sidechain" )
    rmsd = mdtraj.rmsd(traj, traj, frame=centroid, atom_indices=atoms_indices, parallel=True, precentered=True)

    print( '.... computing RMSD .... DONE')

    return rmsd
def contacts_bonds(traj, peptide_chain ):

  group_1 = [residue.index for residue in traj.topology.chain(peptide_chain).residues ]
  group_2 = [residue.index for residue in traj.topology.chain(0).residues or traj.topology.chain(1).residues or   traj.topology.chain(2).residues ]
  pairs = list(product(group_1, group_2))

  contacts_bonds= mdtraj.compute_contacts(traj,pairs , scheme='closest-heavy', ignore_nonprotein=True, periodic=True, soft_min=False, soft_min_beta=20)


def hydrogen_bonds(traj):
  print( '.... computing hbonds (can take few minutes) ....')

  list_peptide_hbonds=[]
  hbonds = mdtraj.baker_hubbard(traj, freq =0.4)

  for hbond in hbonds:
      #if ( hbond[0] in  traj.topology.select( 'chainid  %s to %s' % (3,5) ) and   hbond[2] not in  traj.topology.select( 'chainid  %s to %s' % (3,5)) ) :
      if ( hbond[0] in  traj.topology.select( 'chainid  %s ' % (3) ) and   hbond[2] not in  traj.topology.select( 'chainid  %s' % (3)) ) :
          list_peptide_hbonds.append( hbond[[0,2]] )
          print(' hbond : %s -- %s' % (traj.topology.atom(hbond[0]), traj.topology.atom(hbond[2])))
          print( 'hbond : %s -- %s' % (hbond[0], hbond[2]))
      #elif ( hbond[0] not in  traj.topology.select( 'chainid  %s to %s' % (3,5) ) and   hbond[2]  in  traj.topology.select( 'chainid  %s to %s' % (3,5)) ) :
      elif ( hbond[0] not in  traj.topology.select( 'chainid  %s ' % (3) ) and   hbond[2]  in  traj.topology.select( 'chainid  %s ' % (3)) ) :
          list_peptide_hbonds.append( hbond[[0,2]] )
          print(' hbond : %s -- %s' % (traj.topology.atom(hbond[0]), traj.topology.atom(hbond[2])))
          print('hbond : %s -- %s' % (hbond[0], hbond[2]))
  '''
  da_distances = mdtraj.compute_distances(traj,  np.array(list_peptide_hbonds), periodic=False)

  color = cycle(['r', 'b', 'gold'])
  print(len(np.array(list_peptide_hbonds)))
  for i in  range(len(np.array(list_peptide_hbonds))):
      plt.hist(da_distances[:, i], color=next(color), label=label(hbonds[i]), alpha=0.5)
  plt.legend()
  plt.ylabel('Freq');
  plt.xlabel('Donor-acceptor distance [nm]')
  plt.savefig('hbonds.png', bbox_inches='tight')
  '''
  return np.array(list_peptide_hbonds)

def str2bool(v):
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

#------------------------------------------------------------




parser = argparse.ArgumentParser(description='Process MDtraj')
parser.add_argument('--peptide_chains', metavar='N', type=str, nargs='+',
                    help='peptide_chains')
parser.add_argument('--topol', default='nvt.gro', type=str,
                    help='topology file (.gro or .pdb)')
parser.add_argument('--traj', default='md.xtc', type=str,
                    help='md.xtc')
parser.add_argument('--process_traj',  type=str2bool, nargs='?',
                        const=True, default=True, help='cluster and align protein')




args = parser.parse_args()


if  args.process_traj ==True:
    process_traj(args.topol, args.traj )
    traj = mdtraj.load('traj_cluster_alg.xtc', top='traj_cluster.pdb')

else :     traj = mdtraj.load(args.traj, top=args.topol)

traj.center_coordinates()
results=[]


for i in [3]:
    contact = contacts_bonds(traj, i )

for i in [3]:
    hel = helicity(traj,i)
    results.append(hel)



plt.figure()
centroid=find_centroid(traj)
for i in [3]:
        rmsd=compute_rmsd(traj, i, centroid)
        results.append(1/np.average(rmsd))
        plt.plot(range(traj.n_frames), rmsd, '-')

plt.title('rmsd')
plt.xlabel('time')
plt.ylabel('rmsd')
plt.savefig('rmsd.png', bbox_inches='tight')
plt.figure()

#hydrogen_bonds= hydrogen_bonds(traj)
#results.append(len(hydrogen_bonds)/3)
print(hydrogen_bonds)
side = compute_rmsd_sidechains(traj, centroid)
results.append(side)



def draw_radar(results):
    N = len(results)
    theta = radar_factory(N, frame='polygon')
    r= np.array(results)
    #   theta = np.arange(0, 2 * np.pi , 2 * np.pi /len(results) ,dtype=float )
    area = 200
    fig, axes = plt.subplots(figsize=(N, N), nrows=2, ncols=2,
                             subplot_kw=dict(projection='radar'))
    fig.subplots_adjust(wspace=0.25, hspace=0.20, top=0.85, bottom=0.05)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='polar')
    ax.set_xticklabels(['helicity', 'rmsd' , 'hbond'])
    c = ax.scatter(theta, r, c=colors, s=area, cmap='hsv', alpha=0.75)
    fig.savefig('score.png', bbox_inches='tight')
