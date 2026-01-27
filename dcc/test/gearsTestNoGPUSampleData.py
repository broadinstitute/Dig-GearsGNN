

# working dir
DIR_USER = "/home/javaprog"
DIR_CODE = "/home/javaprog/Code/MlWorkspace"

DIR_USER = "/home/ubuntu"
DIR_CODE = "/home/ubuntu/Code"


# imports
import sys
sys.path.append('{}/Dig-GearsGNN/'.format(DIR_CODE))
from gears import PertData, GEARS, data_utils
import logging

# constants
DIR_GEARS = "{}/Data/Broad/GeneticsML/Gears202512".format(DIR_USER)
print("got GEAR directory: {}".format(DIR_GEARS))

FILE_LOG = "gears_obesity.log"
# FILE_LOG = "{}/Logs/gears_obesity.log".format(DIR_GEARS)
print("got log file: {}".format(FILE_LOG))
# DIR_DCC_OBESITY_DATA = "{}/Data/broad_obesity".format(DIR_GEARS)

logging.basicConfig(
    filename=FILE_LOG,
    filemode="a",               # append; use "w" to overwrite
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
)


pert_data = PertData('{}/Data'.format(DIR_GEARS)) # specific saved folder

# load data
# MODIFY HERE FOR NEW DATASETS
# pert_data.load(data_name = 'norman') # specific dataset name
pert_data.load(data_name = 'norman') # specific dataset name

pert_data.new_data_process(dataset_name='norman_test', adata=pert_data.adata)

# compute uns
data_utils.get_DE_genes(adata=pert_data.adata, skip_calc_de=False)
data_utils.get_dropout_non_zero_genes(adata=pert_data.adata)


pert_data.prepare_split(split = 'simulation', seed = 1) # get data split with seed
pert_data.get_dataloader(batch_size = 32, test_batch_size = 128) # prepare data loader

# set up and train a model
# gears_model = GEARS(pert_data, device = 'cuda:8')
gears_model = GEARS(pert_data, device = 'cpu')
# gears_model.new_data_process(dataset_name='ob_chal')

gears_model.model_initialize(hidden_size = 64)
gears_model.train(epochs = 20)

# save/load model
gears_model.save_model('gears')
gears_model.load_pretrained('gears')

# predict
gears_model.predict([['CBL', 'CNN1'], ['FEV']])
gears_model.GI_predict(['CBL', 'CNN1'], GI_genes_file=None)




