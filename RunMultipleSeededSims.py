#!/usr/bin/env python

# This script generates and submits jobs for running multiple seeded simulations
import multiprocessing
import subprocess
import os
import time
import smtplib # Sending email
from email.mime.text import MIMEText # Email module

# This is a helper function that allows to run a bash command in a separate process
def execute_command(cmd):
	    return subprocess.call(cmd, shell=True)

# Record the start time
start_time = time.time()

# Import the command line flags
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("filename",type=str,default=None,help="the file to run")
parser.add_argument("-s","--maxseed",type=int,default=10,dest="maxseed",help="the number of seeds to use")
parser.add_argument("-s0","--minseed",type=int,default=0,dest="minseed",help="The starting value of the seed range")
parser.add_argument("-d", "--dir",type=str,default="/data/clairem/LogFiles",dest="outputdir",help="the output directory")
parser.add_argument("-o","--buildoptions",type=str,default="",dest="buildopts",help="A string with any additional flags to pass to scons")
parser.add_argument("-p","--parallel",type=str,default=None,dest="parallel",help="Number of processors for each simulation")
parser.add_argument("-e","--extra_input",nargs='*',dest="inputs",help="any extra indices to run in the code, must use -e flag before -s flag")
args = vars(parser.parse_args())

# Determine the build type
build_folder = "optimised"
build_flag = "b=GccOpt"
runner_preamble = ""
if args["parallel"] is not None:
	build_folder = build_folder + "_native"
	build_flag = build_flag + "Native_" + args["parallel"]
	runner_preamble = "mpirun -np " + args["parallel"] + " "

# Get the runner file name
runnerfile = runner_preamble + ( args['filename'].replace("test/","./projects/ClaireM/build/" + build_folder + "/").replace(".hpp","Runner") )

# Edit the log directory
log_dir = args['outputdir'] + '/' + args['filename'].rsplit('/',1)[-1].replace(".hpp","") + '/'

# Log what is happening
print "\nBuilding file %s to runner %s.\nExecuting %d simulations to output directory %s\n" % (args['filename'],runnerfile,args['maxseed'],log_dir)

# Change to the correct directory
import os
os.chdir("../..")

# Make the output directory
if not os.path.isdir(log_dir):
	execute_command("mkdir %s" % log_dir)

# Run the build
build_command = "scons %s %s co=1 %s > %s/build_output.log" % (args['buildopts'],build_flag,"projects/ClaireM/"+args['filename'],log_dir)
print "Running build command '%s'\n" % build_command
execute_command(build_command)

# Build the run commands
command_list = ["nice -20 " + runnerfile]
extra_inputs = args['inputs']
# Add in any extra inputs
if extra_inputs is not None:
	n_inputs = len(extra_inputs)/2
	for i in range(n_inputs):
		n_indices = int(extra_inputs[2*i+1])
		n_previous = len(command_list)
		command_list = [item for item in command_list for k in range(n_indices)]
		command_addition = [(" -" + extra_inputs[2*i] + " ") + s for s in map(str,range(n_indices))*n_previous]
		command_list = map(str.__add__,command_list,command_addition)
# Repeat for each seed
n_seeds = args['maxseed']
min_seed = args['minseed']
n_previous = len(command_list)
command_list = [item for item in command_list for k in range(n_seeds)]
command_addition = [" -seed " + s for s in map(str,range(min_seed, min_seed + n_seeds))*n_previous]
command_list =  map(str.__add__,command_list,command_addition)
# Add output
log_directory = args['filename'].rsplit('/',1)
log_directory = log_directory[1].replace('.hpp','')
for command_i in command_list:
	log_file = command_i.split("Runner",1)[1]
	log_file = log_file.replace("-","").replace(" ","")
	command_i = command_i + (" &> %s/%s/%s.log" % (args['outputdir'],log_directory,log_file))
	print command_i 

# generate a pool of workers
pool = multiprocessing.Pool(len(command_list))

# Check the runner is not an old version (if the code has been run previously)
runner_creation_time = os.path.getmtime(runnerfile)
if (runner_creation_time < start_time):
	print 'ERROR: Build failed. \n'
else:
	# Execute
	print "\nRunning executables %s\n" % runnerfile
	pool.map(execute_command,command_list)

