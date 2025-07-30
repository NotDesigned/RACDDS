import subprocess
import sys
import os
import multiprocessing
import time
import json

current_time = time.strftime("%Y-%m-%d-%H-%M-%S", time.localtime())
dataset_path = "/mnt/data/yingli/DSD_Dataset/directed/"
output_path = "/mnt/data/yingli/DSD_Dataset/outputs_llc/"
multiprocess = 4

# this python script is for running all the experiments for the paper
# The program will accoridng to the parameters to run the experiments, and save the results to the corresponding folder
# these parameter will be passed to the DensestSubgraph program by command line arguments or a configuration file

# read the configuration file
parameters = list(map(lambda x: str(x), sys.argv[1:]))
if len(parameters) == 0:
    print("Please provide the configuration file path")
    exit(1)
config_file = parameters[0]
with open(config_file, "r") as file:
    config = json.load(file)
name = config["name"]
common = config["common"]
datasets = config["datasets"]
epss = config["epss"]
expiration = config["expiration"]
if config["current_time"] is not None:
    current_time = config["current_time"]
print("current_time:", current_time)
if "multiprocess" in config and config["multiprocess"] is not None:
    multiprocess = config["multiprocess"]

path = os.path.join(output_path, name+'_'+current_time)
if not os.path.exists(path):
    os.mkdir(path)

def run_command(common, dataset, eps, expiration, name):
    command = common + ["-path", os.path.join(dataset_path, dataset+".txt"), "-eps", str(eps)]
    # print(" ".join(command))
    filename = os.path.join(path, name) + \
        "_".join(["-eps", str(eps)]) + \
        "_" + dataset + ".txt"
    # if already exist, do not run
    if os.path.exists(filename):
        print("Skip "+ filename)
        return 
    print("Begin"," ".join(command))
    try:
        result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, timeout= expiration)
    except subprocess.TimeoutExpired:
        print("Timeout ", " ".join(command))
        with open(filename, "w") as file:
            file.write("Timeout")
        return
    except subprocess.CalledProcessError as e:
        print("Error ", " ".join(command))
        print("Program output:", e.stdout)
        print("Error running program:", e.stderr)
    except Exception as e:
        print("Unknown error ", " ".join(command))
        print(e)
    print("Finish "," ".join(command))
    with open(filename, "w") as file:
        file.write(result.stdout)

def main():
    commands = []
    for eps in epss:
        for dataset in datasets:
            commands.append((common, dataset, eps, expiration, name))
    print("Number of commands:", len(commands))
    pool = multiprocessing.Pool(processes=4)
    pool.starmap(run_command, commands)

if __name__ == "__main__":
    main()
