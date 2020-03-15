####
#### Echo State Network (Reservoir computing) implemented with python
####

# Import modules
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt # For visualization only
from scipy import linalg
from sklearn import preprocessing
from sklearn.model_selection import train_test_split # For data selection

# Define functions
## Generate a random internal weight matrix W0 -----------------------------#
def generate_random_W(num_reservoir_nodes, alpha, sparsity):
    np.random.seed(0)
    W0 = np.random.uniform(-1, 1, (num_reservoir_nodes * num_reservoir_nodes, 1))
    np.random.seed(0)
    rand_id = np.random.randint(0, W0.shape[0], int(num_reservoir_nodes * num_reservoir_nodes * sparsity))
    W0[rand_id] = 0
    W1 = W0.reshape(num_reservoir_nodes, num_reservoir_nodes)
    ## Normalize W0 using spectral radius
    W2 = W1 / abs(np.linalg.eig(W1)[0]).max()
    ## Scale W1 to W using alpha
    W = W2 * alpha
    return(W)
#---------------------------------------------------------------------------#

## Calculate the next state (reservoir dynamics) ---------------------------#
def calculate_next_state(input_value, reservoir_nodes_t0, input_W, reservoir_W, leak, bias = 0):
    # Set the current state
    t0 = np.array(reservoir_nodes_t0)
    # Calculate the next state
    t1 = leak @ t0 + np.tanh([input_value] @ input_W + (I - leak) @ t0 @ reservoir_W + bias)
    return(t1)
#---------------------------------------------------------------------------#

## Adjust output using regularized regression  -----------------------------#
def adjust_output(teacher_signal, reservoir_node_log, ridge_lambda):
    # Ridge Regression
    E_lambda = np.identity(reservoir_node_log.shape[1]) * ridge_lambda
    inv_x = np.linalg.inv(reservoir_node_log.T @ reservoir_node_log + E_lambda)
    # update weights of output layer
    output_W = (inv_x @ reservoir_node_log.T) @ teacher_signal
    return(output_W)
#---------------------------------------------------------------------------#


# -------------------- 1. Set parameters -------------------- #
 # Prepare output object
all_predict_outputs = []
target_var, target_ts = 'C1', 'eres_ts'
network_name = 'Random Network'

# Weight parameters
alpha = 0.9 # Adjust spectral radius (should be adjusted), smaller the faster!
w_in_strength = 0.6 # Adjust the absolute values of W_in
Win_sparsity = 0.9
W_sparsity = 0.9
leak_rate = 0 # How the previous state influence the next state (0 = no influence of the previous state)
lambda0 = 0.05 # Specify lambda in a ridge regression

# Netrowk parameters
num_reservoir_nodes = 100 # Equals to the number of species in the ecological reservoir
num_input_nodes = 1 # Input is a scalar value
num_output_nodes = 1 # Output is a scalar vlue

# Specify the test fraction
test_fraction = 0.2

# Import data to be predicted
eres_ts = pd.read_csv('./data/ESM5_Data_5spModel.csv')


# -------------------- 2. DATA PREPARATION -------------------- #
# Standardize data
data  = eres_ts[target_var]
std_data = (data - data.min()) / (data.max() - data.min())  ## Data standardization

# Split data into training and test data set
data_size = data.shape[0]
test_size = int(round(data_size*test_fraction))
train_size = int(data_size - test_size)
x_train = std_data[0:train_size].reset_index()[target_var]
t_train = std_data[1:(train_size + 1)].reset_index()[target_var]
x_test = std_data[train_size:(data_size - 1)].reset_index()[target_var]
t_test = std_data[(train_size + 1):data_size].reset_index()[target_var]

# Set leak and weight matrices
I = np.identity(num_reservoir_nodes)
R = I * leak_rate
np.random.seed(0); W_in0 = np.random.uniform(-1, 1, (num_input_nodes * num_reservoir_nodes, 1))
np.random.seed(0); rand_id = np.random.randint(0, W_in0.shape[0], int(num_input_nodes * num_reservoir_nodes * Win_sparsity))
W_in0[rand_id] = 0
W_in = W_in0.reshape(num_input_nodes, num_reservoir_nodes) * w_in_strength


# -------------------- 3. Implementation of Echo state netowrk -------------------- #
# Step.1: Prepare reservoir matrix
W = generate_random_W(num_reservoir_nodes, alpha, W_sparsity)

# Step.2: Training using reservoir
record_reservoir_train_nrow = int(x_train.shape[0] + 1) # size of the training data
record_reservoir_nodes = np.zeros((record_reservoir_train_nrow, num_reservoir_nodes))

for data_i, input_train in enumerate(x_train):
    x_n1 = calculate_next_state(input_train, record_reservoir_nodes[data_i,:], W_in, W, R)
    record_reservoir_nodes[data_i + 1,:] = x_n1

record_reservoir_nodes = record_reservoir_nodes[1:,]

# Step.3: Compute output weights using ridge regression
weights_output = adjust_output(t_train, record_reservoir_nodes, lambda0)
outputs = record_reservoir_nodes @ weights_output

# Step.4: Predict test data
record_reservoir_test_nrow = int(x_test.shape[0] + 1)
record_test_reservoir_nodes = np.zeros((record_reservoir_test_nrow, num_reservoir_nodes))
record_test_reservoir_nodes[0,:] = record_reservoir_nodes[-1]

for data_i, input_test in enumerate(x_test):
    x_n1 = calculate_next_state(input_test, record_test_reservoir_nodes[data_i,:], W_in, W, R)
    record_test_reservoir_nodes[data_i + 1,:] = x_n1

record_test_reservoir_nodes = record_test_reservoir_nodes[1:,]
predicted_outputs = record_test_reservoir_nodes @ weights_output

# Step.5: Compute summary statistics
train_cor = np.corrcoef([t_train], [outputs])[1,0]
test_cor = np.corrcoef([t_test], [predicted_outputs])[1,0]


# -------------------- 4. Summarize and visualize results -------------------- #
# Summarize results
all_predict_outputs.append(np.array([network_name, target_var,
                                     round(train_cor, 4), round(test_cor, 4),
                                     num_reservoir_nodes, lambda0,
                                     alpha, w_in_strength, round(leak_rate, 2),
                                     Win_sparsity, W_sparsity, data_size]))
all_predict_df = pd.DataFrame(all_predict_outputs)
all_predict_df = all_predict_df.rename(columns = {0:"network_name", 1:"target_variable", 2:"train_pred",
                                                  3:"test_pred", 4:"num_nodes", 5:"ridge_lambda",
                                                  6:"alpha", 7:"Win_strength", 8:"leak_rate",
                                                  9:"Win_sparsity", 10:"W_sparsity", 11:"data_size"})
pd.set_option('display.max_columns', 12)
all_predict_df

# Visualize results
plt.figure(); f1 = plt.subplot()
f1.plot(t_test); f1.plot(predicted_outputs)
f1.set_ylabel("Value", fontsize = 16)
f1.set_xlabel("Time", fontsize = 16)
plt.show()
