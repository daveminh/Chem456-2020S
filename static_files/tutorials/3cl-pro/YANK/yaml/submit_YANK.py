#!/usr/bin/python

######################################
# Functions for managing GPUS on CCB #
######################################

# Determine the nodes and devices for jobs on the queue
def get_onq():
  import subprocess
  p = subprocess.Popen(['qstat','-f'],
    stdout=subprocess.PIPE, stderr=subprocess.PIPE)
  (stdoutdata, stderrdata) = p.communicate()
  p.wait()

  stdoutdata = stdoutdata.decode()
  stderrdata = stderrdata.decode()
  
  jobs = [job.strip().split('    ') for job in stdoutdata.split('Job Id:')]
  jobs = [dict([[f.replace('\n\t','').strip() for f in field.split(' = ')]
    for field in job if field.find(' = ')>-1])
      for job in jobs if job!=['']]
  cuda_jobs = [job for job in jobs if job['queue']=='cuda']

  import os
  onq = []
  for job in cuda_jobs:
    submit_script = '.'.join(job['Output_Path'].split(':')[-1].split('.')[:-1])
    if os.path.isfile(submit_script):
      F = open(submit_script,'r')
      dat = F.read()
      F.close()
      if dat.find('export CUDA_VISIBLE_DEVICES=')>-1:
        device = dat[dat.find('export CUDA_VISIBLE_DEVICES=')+28:]
        device = int(device[:device.find('\n')])
      else:
        device = None
      if dat.find('#PBS -l nodes=compute-1-')>-1:
        node = dat[dat.find('#PBS -l nodes=compute-1-')+24:]
        node = int(node[:node.find(':')])
      else:
        node = None
      if (device is not None) and (node is not None):
        if job['job_state']=='R':
          onq.append((node, device, int(job['Walltime.Remaining'])))
        elif job['job_state']=='Q':
          onq.append((node, device, max_remaining_time))
    else:
      print('Could not find submission script '+submit_script)
  return onq

def node_gpu_status(node):
  from threading import Timer
  kill = lambda process: process.kill()

  import subprocess
  p = subprocess.Popen(['ssh','compute-1-%d'%node,'"nvidia-smi"'],
    stdout=subprocess.PIPE, stderr=subprocess.PIPE)

  # Waits up to 30 seconds
  my_timer = Timer(30, kill, [p])
  try:
    my_timer.start()
    (stdoutdata, stderrdata) = p.communicate()
    p.wait()
  finally:
    my_timer.cancel()

  stdoutdata = stdoutdata.decode()
  stderrdata = stderrdata.decode()
  
  if stdoutdata=="":
    return []

  # Keep track of gpu IDs and memory used
  gpu_IDs = []
  memory_used = []

  lines = [l.strip() for l in stdoutdata.split('\n')]
  lines = lines[7:]
  while lines!=[] and lines[0]!="":
    gpu_IDs.append(int(lines[0].split()[1]))
    memory_used.append(int(lines[1][33:].split()[0][:-3]))
    lines = lines[3:]

  # Keep track of processes (if supported)
  gpu_used = []
  gpu_used_supported = True
  lines = lines[5:]
  while lines!=[] and lines[0][0]!='+':
    if lines[0].find('No running processes found')!=-1:
      pass
    elif lines[0].find('Not Supported')==-1:
      gpu_used.append(int(lines[0].split()[1]))
    else:
      gpu_used_supported = False
    lines = lines[1:]

  # Determine which GPUs are free
  gpu_free = []
  if gpu_used_supported:
    for id in gpu_IDs:
      if not id in gpu_used:
        gpu_free.append(id)
  else:
    for id in gpu_IDs:
      if memory_used[id]<100:
        gpu_free.append(id)

  return [(node, gpu_IDs[id], id in gpu_free, memory_used[id]) for id in gpu_IDs]

def pickGPU():
  # Prioritize GPUs
  gpu_status = []
  for node in range(1,9):
    gpu_status.extend(node_gpu_status(node))

  print('GPU status', gpu_status)

  onq = get_onq()
  jobs = [(job[0],job[1]) for job in onq]

  gpu_status.sort(key=lambda x:x[0],reverse=True) # 2nd priority: higher nodes
  gpu_status.sort(key=lambda x:x[2],reverse=True) # Free nodes according to nvidia-smi

  gpu_available = [g for g in gpu_status if (g[0],g[1]) not in jobs]
  if len(gpu_available)>0:
    print('Available GPUs:', gpu_available)
    (node,device,free,mem_used) = gpu_available[0]
  else:
    gpu_status.sort(key=lambda x:x[3]) # Prioritize by memory usage
    (node,device,free,mem_usedwalltime) = gpu_status[0]
  return (node, device)

###############
# Main script #
###############

import argparse
parser = argparse.ArgumentParser(\
  description = 'Submits YANK calculation to Bridges or CCB.')
parser.add_argument('--type', \
  choices=['shared','small','test','CCB-CPU','CCB-GPU'], \
  default='test', \
  help="The type of job." + \
    "'shared' is 48 hrs on Bridges GPU-shared, " + \
    "'small' is 8 hrs on Bridges GPU-small, " + \
    "'test' is 0.5 hrs on Bridges GPU-small." + \
    "'CCB-CPU' is a CPU calculation on CCB." + \
    "'CCB-GPU' is a GPU calculation on CCB.")
parser.add_argument('--yaml', \
  default='1ubq', \
  help='The yaml script for the YANK job.')
parser.add_argument('--dry', 
  action='store_true', \
  help='Creates script but does not submit it.')
parser.add_argument('--setup_only', \
  action='store_true', \
  help='Only runs setup')
args = parser.parse_args()

job_name = args.yaml.split('.')[0]
setup_only = '--setup-only' if args.setup_only else ''

import os.path
script_path = os.path.join(os.getcwd())
os.chdir(script_path)

# Submission scripts
if os.path.exists('/pylon5'): # Bridges
  cluster = 'Bridges'
  if args.type=='shared':
    partition = 'GPU-shared'
    wall_time = '48:00:00'
  elif args.type=='small':
    partition = 'GPU-small'
    wall_time = '08:00:00'
  elif args.type=='test':
    partition = 'GPU-small'
    wall_time = '00:30:00'
  else:
    raise Exception('Job type incompatible with cluster')

  submit_script = f'''#!/bin/bash
#SBATCH -N 1
#SBATCH -p {partition}
#SBATCH --ntasks-per-node 1
#SBATCH --gres=gpu:p100:1
#SBATCH -t {wall_time}

#echo commands to stdout
set -x

# Set up the environment
source ~/.bashrc
module load anaconda3'''
if os.path.exists('/share/apps/miniconda3'): # CCB Cluster
  cluster = 'CCB'
  wall_time = '168:00:00'

  if args.type=='CCB-CPU':
    submit_script = f'''#!/bin/bash
#
#PBS -S /bin/bash
#PBS -l mem=2GB,nodes=1:ppn=1,walltime={wall_time}
#PBS -q default

source ~/.bashrc
module load miniconda/3'''
  elif args.type=='CCB-GPU':
    (node, device) = pickGPU()
    submit_script = f'''#PBS -S /bin/bash
#PBS -j oe
#PBS -q cuda
#PBS -W x=GRES:gpus@1
#PBS -l nodes=compute-1-{node}:ppn=1:gpus,walltime={wall_time}
#PBS -r y

export CUDA_VISIBLE_DEVICES={device}

echo Running on compute-1-{0}, targeting device {1}

source ~/.bashrc
module load miniconda/3'''
  else:
    raise Exception('Job type incompatible with cluster')

submit_script = submit_script + \
f'''
conda activate yank

# Run YANK based on YAML script
cd {script_path}
yank script -y {args.yaml} {setup_only}
'''

# Write and submit scripts
submit_F = open(job_name+'.job','w')
submit_F.write(submit_script)
submit_F.close()

if not args.dry:
  if cluster in ['Bridges','Comet']:
    os.system(f'sbatch {job_name}.job')
  elif cluster=='CCB':
    os.system(f'qsub {job_name}.job')
