import subprocess
import zipfile
import os
import pandas as pd
from cell_population import *
from expression_normalization import *
from differential_expression import *

# get directory where data will reside
sceptre2_dir = subprocess.run('source ~/.research_config; echo $LOCAL_SCEPTRE2_DATA_DIR', shell = True, capture_output = True).stdout.decode('utf-8').rstrip()
weissman_check_dir = sceptre2_dir + 'results/weissman_check'

# create data directory if it is not already created
if not os.path.isdir(weissman_check_dir):
  os.mkdir(weissman_check_dir)

# download data if it is not already downloaded
if not os.path.isdir(weissman_check_dir + "/sequencing"):
  print('Downloading data...')
  data_url = 'https://www.dropbox.com/s/22yi75rw21bsvdg/perturb_seq_demo_sequencing.zip?dl=1'
  zipfile_path = weissman_check_dir + '/perturb_seq_demo_sequencing.zip'
  wget.download(data_url, zipfile_path)
  zipfile.ZipFile(zipfile_path, 'r').extractall(weissman_check_dir)
  os.remove(zipfile_path)

# extract Thapsigargin subpopulation and save to disk if not already done
if not os.path.exists(weissman_check_dir + "/thaps_pop.hdf"):
  print('Extracting and saving Thapsigargin cells...')
  # read sequencing data into CellPopulation object
  pop = CellPopulation.from_file(weissman_check_dir + '/sequencing/10X005_new/',
                               genome='hg19')
  
  # restrict attention to droplets with single cell                             
  pop = pop.subpopulation(cells='single_cell')

  # list of targeting gRNAs
  perturbations = ['ATF6_only_pMJ145', 'PERK_only_pMJ146',  'IRE1_only_pMJ148', 
                 'ATF6_PERK_pMJ150', 'PERK_IRE1_pMJ154', 'ATF6_IRE1_pMJ152', 
                 'ATF6_PERK_IRE1_pMJ158']

  # list of control gRNAs
  controls = ['3x_neg_ctrl_pMJ144-1', '3x_neg_ctrl_pMJ144-2']

  # mapping from gRNA name to more informative label
  guide_renamer = {'3x_neg_ctrl_pMJ144-1': 'control',
                 '3x_neg_ctrl_pMJ144-2': 'control',
                 'ATF6_only_pMJ145': 'ATF6',
                 'PERK_only_pMJ146': 'PERK',
                 'IRE1_only_pMJ148': 'IRE1',
                 'ATF6_PERK_pMJ150': 'ATF6-PERK',
                 'PERK_IRE1_pMJ154': 'IRE1-PERK',
                 'ATF6_IRE1_pMJ152': 'ATF6-IRE1',
                 'ATF6_PERK_IRE1_pMJ158': 'ATF6-IRE1-PERK'}

  # mapping from gRNAs to targets
  guide_targeter = {'3x_neg_ctrl_pMJ144-1': [],
                 '3x_neg_ctrl_pMJ144-2': [],
                 'ATF6_only_pMJ145': ['ATF6'],
                 'PERK_only_pMJ146': ['EIF2AK3'],
                 'IRE1_only_pMJ148': ['ERN1'],
                 'ATF6_PERK_pMJ150': ['ATF6', 'EIF2AK3'],
                 'PERK_IRE1_pMJ154': ['ERN1', 'EIF2AK3'],
                 'ATF6_IRE1_pMJ152': ['ATF6', 'ERN1'],
                 'ATF6_PERK_IRE1_pMJ158': ['ATF6', 'EIF2AK3', 'ERN1']}

  # mapping from sample number to label
  sample_renamer = {1: 'Tm',
                   2: 'Thaps',
                   3: 'DMSO'}

  # read gem group information from cell barcodes
  pop.cells['gem_group'] = pop.cells.index.map(lambda x: int(x.split('-')[-1]))
  
  # extract list of guides
  pop.guides = sorted(pop.cells[pop.cells['single_cell']]['guide_identity'].unique())

  # clean up the cell metadata
  experiment = pop.cells['gem_group'].map(lambda x: sample_renamer[x])
  experiment.name = 'experiment'
  perturbation = pop.cells['guide_identity'].map(lambda x: guide_renamer[x])
  perturbation.name = 'perturbation'
  long_experiment = experiment + '_' + perturbation
  long_experiment.name = 'long_experiment'
  pop.add_property(cells=pd.DataFrame([experiment, perturbation, long_experiment]).T)

  # remove lowly expressed genes
  strip_low_expression(pop) 

  # find subpopulation treated with Thapsigargin (as well as DMSO control cells)
  thaps_experiments = ['DMSO_control', 'Thaps_control', 'Thaps_ATF6-IRE1', 
                       'Thaps_ATF6-PERK', 'Thaps_IRE1-PERK', 'Thaps_ATF6-IRE1-PERK',
                       'Thaps_ATF6', 'Thaps_IRE1', 'Thaps_PERK']
  thaps_pop = pop.subpopulation(cells='long_experiment in @thaps_experiments',
                             thaps_experiments=thaps_experiments)

  # save Thapsigargin population to file
  thaps_pop.to_hdf(weissman_check_dir + "/thaps_pop.hdf")

# run differential expression test if results are not already present
if not os.path.exists(weissman_check_dir + "/original_pvalues.csv"):
  print('Running differential expression analysis...')
  thaps_pop.normalized_matrix = normalize_to_control(thaps_pop,
                                                     'long_experiment == "DMSO_control"')
  ks, ps, adj_ps = ks_de(thaps_pop,
                         key='long_experiment',
                         control_cells='long_experiment == "DMSO_control"',
                         genes='mean > 0.25',
                         normalized=True,
                         alpha=0.001,
                         multi_method='fdr_by',
                         n_jobs=1)

  ps.to_csv(weissman_check_dir + "/original_pvalues.csv")
