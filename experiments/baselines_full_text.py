import json
from tqdm import tqdm
import os


dataset_index_file = "../data/external/ec2s/data_index.json"

with open(dataset_index_file, "r") as f:
    data_index = json.load(f)


datasets = [x for x in os.listdir("../data/external/ec2s/") if os.path.isdir(x)]
