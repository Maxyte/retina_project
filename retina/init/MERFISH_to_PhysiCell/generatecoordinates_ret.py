# *************** LOAD LIBRARIES AND DATA *************************************************
# Script to convert MERSCOPE spatially-resolved cell type data
# (Choi et al.;;https://www.biorxiv.org/content/10.1101/2022.12.04.518972v1; https://zenodo.org/records/8144355)
# into PhysiCell-interpretable initial conditions.
import scanpy as sc
import math
import pandas as pd
import numpy as np
from physicool.config import ConfigFileParser # pip install physicool
import matplotlib.pyplot as plt
import seaborn as sns
import csv

def load_coords_from_mer(coord_file, pc_config_file, cell_type_vec, num_batches, compress, resolution, batch_index_to_compress):
    # LOAD IN MERSCOPE/CONFIGURATION DATA
    # Load in MERFISH data from coord_file and initialize spatial data/cell count arrays
    # cell_type_vec for our data set is ['Rod', 'Cone', 'AC', 'BC', 'RGC', 'HC', 'MG']
    adataret = sc.read_h5ad(coord_file)
    preamble = ["x", "y", "z", "type", "volume", "cycle entry", "custom:GFP", "custom:sample"] # header for PhysiCell CSV-initialization parser (v2.0)
    batches = num_batches # tissues to analyze-- this may most often be 1
    # SPATIAL HELPER FUNCTIONS
    def update_batch_data(adataindex, batch_index, cell_type, batch_vector, data_vector):
        data_vector[batch_index].append([adataret.obs.center_x[adataindex], adataret.obs.center_y[adataindex], 0, adataret.obs.majorclass[adataindex]])
        for j in range(0, len(cell_type_vec)):
            if (cell_type == cell_type_vec[j]):
                batch_vector[batch_index][j] = batch_vector[batch_index][j] + 1
        return
    def in_window(x, y, x_min, x_max, y_min, y_max):
        in_range = False
        if (x > x_min and x < x_max and y > y_min and y < y_max):
            in_range = True
        return in_range
    def index_to_coord(i, x_mn, x_mx, y_mn, y_mx):
        cofactor = (math.floor(i/n) % 2)
        x = resolution * ( n * (i/n - math.floor(i/n)) + ( cofactor * math.cos(math.pi/3))) + x_mn
        y = resolution * math.floor(i/n) * math.sin(math.pi/3) + y_mn
        return (x, y)
    # get microenvironment voxel index of cell at coord_pair
    def coord_to_index(coord_pair, x_mn, x_mx, y_mn, y_mx):
        cofactor = math.floor((coord_pair[1]-y_mn)/(resolution * math.sin(math.pi/3))) % 2
        ind = int(n * math.floor((coord_pair[1]-y_mn)/(resolution * math.sin(math.pi/3))) + math.floor(( coord_pair[0] - (x_mn + resolution * cofactor * math.cos(math.pi/3))) / resolution))
        return ind
    # DATA INTERPRETATION
    data = [] # array to contain PhysiCell-parseable spatial data when interpreted as CSV
    cell_count_by_batch = []
    cell_proportion_by_batch = []
    for i in range(0, batches):
        data.append([])
        cell_count_by_batch.append(np.zeros(len(cell_type_vec)))
        cell_proportion_by_batch.append(0.0)
        data[i].append(preamble)
    # rectangular domain based on sample observation
    x_min_adhoc = -100
    x_max_adhoc = -37.5
    y_min_adhoc = -50
    y_max_adhoc = 50
    for i in range(1, len(adataret.obs.volume)):
        batch_index = int(adataret.obs.batch[i])
        cell_type = adataret.obs.majorclass[i]
        update_batch_data(i, batch_index, cell_type, cell_count_by_batch, data)
    for i in range(0, batches):
        if (sum(cell_count_by_batch[i]) != 0):
            cell_proportion_by_batch[i] = cell_count_by_batch[i]/sum(cell_count_by_batch[i])
    # now, "compress" tissue image by mapping PhysiCell voxel indices to cell types as a function of local aggregates and proportionality bias
    if (compress):
        data = data[batch_index_to_compress]
        xml_data = ConfigFileParser(pc_config_file)
        domain = xml_data.read_domain_params()
        xdomain_min = domain.x_min
        xdomain_max = domain.x_max
        ydomain_min = domain.y_min
        ydomain_max = domain.y_max
#        scale = 1
#        scale = ((xdomain_max - xdomain_min) * (ydomain_max - ydomain_min)) / ((x_max_adhoc - x_min_adhoc) * (y_max_adhoc - y_min_adhoc))
#        x_center = 0
#        y_center = 0
#        for i in range(1, len(data)):
#            x_center += float(data[i][0]) / (len(data) - 1)
#            y_center += float(data[i][1]) / (len(data) - 1)
#        for i in range(1, len(data)):
#            data[i][0] = scale * (data[i][0] - x_center)
#            data[i][1] = scale * (data[i][1] - y_center)
        n = int(math.floor((xdomain_max - xdomain_min)/resolution) + 1)
        m = int(math.floor((ydomain_max - ydomain_min)/(resolution * math.sin(math.pi/3))) + 1)
        num_voxels = int(m * n)
        # get 2D coordinates of the voxel at index i
        indices = []
        for i in range(1, len(data)):
            cell_type_id = None
            for j in range(0, len(cell_type_vec)):
                if (data[i][3] == cell_type_vec[j]):
                    cell_type_id = j
            if ( in_window(data[i][0], data[i][1], xdomain_min, xdomain_max, ydomain_min, ydomain_max) ):
                indices.append([coord_to_index([float(data[i][0]), float(data[i][1])], xdomain_min, xdomain_max, ydomain_min, ydomain_max), cell_type_id])
        # enumerate cells in each voxel by type
        cell_bucket = []
        for i in range(0, num_voxels):
            cell_bucket.append(np.zeros(len(cell_type_vec)))
        # carefully (violently) pour cell vectors into buckets
        for i in range(0, len(indices)):
            cell_bucket[indices[i][0]][indices[i][1]] += 1
        # assign type character based on local cell type info
        cell_char = []
        local_proportion = []
        for i in range(0, num_voxels):
            cell_char.append(np.zeros(len(cell_type_vec)))
            local_proportion.append(np.zeros(len(cell_type_vec)))
        for i in range(0, num_voxels):
            type_vec = cell_bucket[i]
            neighbour_nodes = 0
            has_left = (np.floor((i-1)/n) == np.floor(i/n))
            has_right = (np.floor((i+1)/n) == np.floor(i/n))
            has_above = (i + n < num_voxels)
            has_below = (i - n > 0)
            if (has_left):
                neighbour_nodes += 1
                type_vec = np.add(type_vec, cell_bucket[i-1])
                if (has_below):
                    neighbour_nodes += 1
                    type_vec = np.add(type_vec, cell_bucket[i-n-1])
                if (has_above):
                    neighbour_nodes += 1
                    type_vec = np.add(type_vec, cell_bucket[i+n-1])
            if (has_below):
                neighbour_nodes += 1
                type_vec = np.add(type_vec, cell_bucket[i-n])
            if (has_above):
                neighbour_nodes += 1
                type_vec = np.add(type_vec, cell_bucket[i+n])
            if (has_right):
                neighbour_nodes += 1
                type_vec = np.add(type_vec, cell_bucket[i+1])
                if (has_below):
                    neighbour_nodes += 1
                    type_vec = np.add(type_vec, cell_bucket[i-n+1])
                if (has_above):
                    neighbour_nodes += 1
                    type_vec = np.add(type_vec, cell_bucket[i+n+1])
            cell_char[i] = np.add(cell_char[i], cell_bucket[i])
            local_proportion[i] = np.add(cell_char[i], type_vec)
            if (np.sum(cell_char[i] != 0)):
                local_proportion[i] = local_proportion[i]/np.sum(local_proportion[i])
        # now, normalize cell counts by local cell proportions to assign a cell type to each voxel
        type_char = cell_bucket
        for i in range(0, num_voxels):
            for j in range(0, len(cell_type_vec)):
                if ( local_proportion[i][j] != 0 ):
                    type_char[i][j] = cell_bucket[i][j]/local_proportion[i][j]
        data = []
        data.append(preamble)
        for i in range(0, num_voxels):
            pos = index_to_coord(i, xdomain_min, xdomain_max, ydomain_min, ydomain_max)
            if (np.max(type_char[i]) != 0):
                data.append([str(pos[0]),str(pos[1]),str(0.0), str(cell_type_vec[np.argmax(type_char[i])])])
                #data.append([str(pos[0]), str(pos[1]), str(0.0), np.argmax(type_char[i])])
    return data

coord_file = "VA45_integrated.h5ad"
pc_config_file = "/Users/morgana/Documents/GitHub/retina_rules/PhysiCell/user_projects/retina_dev/config/PhysiCell_settings.xml"
cell_type_vec = ["Rod", "Cone", "AC", "BC", "RGC", "HC", "MG"]
num_batches = 6
compress = True
resolution = 16.825 # standard cell diameter in PhysiCell
batch_index_to_compress = 1
data = load_coords_from_mer(coord_file, pc_config_file, cell_type_vec, num_batches, compress, resolution, batch_index_to_compress)
x_c = 0
y_c = 0
for i in range(1, len(data)):
    x_c += float(data[i][0])/(len(data)-1)
    y_c += float(data[i][1])/(len(data)-1)
for i in range(1, len(data)):
    data[i][0] = str((float(data[i][0]) - x_c))
    data[i][1] = str((float(data[i][1]) - y_c))
target_file_name = []
target_file_name.append(r"VA45_compressed_coords" + str(batch_index_to_compress) + ".csv")
with open(target_file_name[0], 'w') as csvfile:
    csvfile.write('x,y,z,type,volume,cycle entry,custom:GFP,custom:sample\n')
    for i in range(1, len(data)):
        csvfile.write(f'{data[i][0]},{data[i][1]},{data[i][2]},{data[i][3]}\n')
