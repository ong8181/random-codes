####
#### MNIST classification by ESN implemented with python
####

# Import essential modules
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt # For visualization only
from scipy import linalg
from sklearn import preprocessing
from sklearn.model_selection import train_test_split # For data selection
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score

# Define functions
## Generate a random internal weight matrix W0 ------------------------------#
def generate_random_W(num_reservoir_nodes, alpha, sparsity):
    np.random.seed(1235); W0 = np.random.uniform(-1, 1, (num_reservoir_nodes * num_reservoir_nodes, 1))
    np.random.seed(1235); rand_id = np.random.choice(W0.shape[0], int(num_reservoir_nodes * num_reservoir_nodes * sparsity), replace = False)
    W0[rand_id] = 0; W1 = W0.reshape(num_reservoir_nodes, num_reservoir_nodes)
    ## Normalize W0 using spectral radius
    W2 = W1 / abs(np.linalg.eig(W1)[0]).max()
    ## Scale W1 to W using alpha
    W = W2 * alpha
    return(W)
#---------------------------------------------------------------------------#

## Calculate the next state (reservoir dynamics) ---------------------------#
def calculate_next_state(input_value, reservoir_nodes_t0, input_W, reservoir_W, leak):
    # Set the current state
    t0 = np.array(reservoir_nodes_t0)
    t1 = leak @ t0 + np.tanh([input_value] @ input_W + (I - leak) @ t0 @ reservoir_W + 1)
    return(t1)
#---------------------------------------------------------------------------#


# -------------------- 1. Set parameters -------------------- #
# Prepare output object
all_predict_outputs = []

# Weight parameters
alpha = 0.9 # Adjust spectral radius (should be adjusted), smaller the faster!
w_in_strength = 0.6 # Adjust the absolute values of W_in
Win_sparsity = 0.9
W_sparsity = 0.9
leak_rate = 0.5 # How the previous state influence the next state
lambda0 = 0.5 # Specify lambda in a ridge regression

# Netrowk parameters
num_input_nodes = 28 # No. of input pixels in a column
num_reservoir_nodes = 200 # No. of reservoir nodes
num_output_nodes = 10 # No. of output nodes

# Specify the data size
subset_size = 20000 # max = 70000
test_fraction = 0.2

# -------------------- 2. DATA PREPARATION -------------------- #
# Down MNIST image data if you do not have it
#mnist = datasets.fetch_openml('mnist_784', version=1,)
#pd.to_pickle(mnist, "data/mnist_downloaded.pkl") # Save MNIST object
# Load MNIST data
mnist = pd.read_pickle("data/mnist_downloaded.pkl")

# Check data structure
data  = mnist.data[range(subset_size)]
std_data = (data - data.min()) / data.max()  ## Data standardization
# Total 70,000 images
label = mnist.target[range(subset_size)]
label = [ int(x) for x in label ] # Convert the label into integer
#one_hot_label = np.identity(10)[label] # Convert the label into "label matrix"

# Show MNIST data examples
#plt.imshow(data[0].reshape(28, 28), cmap='gray'); plt.show()
#plt.imshow(x_train[1].reshape(28, 28), cmap='gray'); plt.show()

# Split data into training data set and test data set
#x_train, x_test, t_train, t_test = train_test_split(std_data, one_hot_label, test_size = test_fraction)
x_train, x_test, t_train, t_test = train_test_split(std_data, label, test_size = test_fraction)

# Secondary parameters
I = np.identity(num_reservoir_nodes)
R = I * leak_rate
np.random.seed(1234); W_in0 = np.random.uniform(-1, 1, (num_input_nodes * num_reservoir_nodes, 1))
np.random.seed(1234); rand_id = np.random.choice(W_in0.shape[0], int(num_input_nodes * num_reservoir_nodes * Win_sparsity), replace = False)
W_in0[rand_id] = 0; W_in = W_in0.reshape(num_input_nodes, num_reservoir_nodes) * w_in_strength

# -------------------- 3. Implementation of Echo state netowrk -------------------- #
# Step.1: Procure an untrained dynamical reservoir
W = generate_random_W(num_reservoir_nodes, alpha, W_sparsity)

# Setup label vectors for training data
t_train_long = np.zeros((t_train.shape[0]*28,10))
for i in np.arange(0, t_train.shape[0]):
    for j in range(28):
        t_train_long[i*28+j,] = t_train[i]

# Step.2: Sample network training dynamics
## Drive the network by the training data  
## Column-wise input according to Schaetti et al. (2016)
record_reservoir_train_nrow = int((subset_size * (1 - test_fraction))*28 + 1)
record_reservoir_nodes = np.zeros((record_reservoir_train_nrow, num_reservoir_nodes))

# Calculate the next state
for image_i, input_train_image in enumerate(x_train):
    # Divide image into columns
    split_train_image = np.split(input_train_image, 28)
    for col_i, input_train in enumerate(split_train_image):
        x_n1 = calculate_next_state(input_train, record_reservoir_nodes[int(image_i*28+col_i),:], W_in, W, R)
        record_reservoir_nodes[int(image_i*28+col_i + 1),:] = x_n1

record_reservoir_nodes = record_reservoir_nodes[1:]
record_res_reshape = record_reservoir_nodes.reshape(int(subset_size * (1-test_fraction)), int(num_reservoir_nodes * 28))

# Step.4: Compute output weights using softmax regression
logistic_model = LogisticRegression(max_iter = 500, solver = 'newton-cg', multi_class = 'auto', penalty = "l2", C = 1)
logistic_model.fit(record_res_reshape, t_train)
learned_model = logistic_model
# Calculate prediction accuracy for training data
train_predict_value = round(accuracy_score(t_train, learned_model.predict(record_res_reshape)), 4)

# Step.5: Exploitation
record_reservoir_test_nrow = int((subset_size * test_fraction)*28 + 1)
record_test_reservoir_nodes = np.zeros((record_reservoir_test_nrow, num_reservoir_nodes))
record_test_reservoir_nodes[0,:] = record_reservoir_nodes[-1]

for image_i, input_test_image in enumerate(x_test):
    # Divide image into columns
    split_test_image = np.split(input_test_image, 28)
    for col_i, input_test in enumerate(split_test_image):
        x_n1 = calculate_next_state(input_test, record_test_reservoir_nodes[int(image_i*28+col_i),:], W_in, W, R)
        record_test_reservoir_nodes[int(image_i*28+col_i + 1),:] = x_n1

record_test_reservoir_nodes = record_test_reservoir_nodes[1:]
record_test_reshape = record_test_reservoir_nodes.reshape(int(subset_size*test_fraction), int(num_reservoir_nodes * 28))
test_predict_value = round(accuracy_score(t_test, learned_model.predict(record_test_reshape)), 4)

# -------------------- 4. Summarize results -------------------- #
all_predict_outputs.append(np.array([train_predict_value, test_predict_value,
                                     num_reservoir_nodes, alpha, w_in_strength,
                                     round(leak_rate, 2), Win_sparsity,
                                     W_sparsity, subset_size]))

all_predict_df = pd.DataFrame(all_predict_outputs)
all_predict_df = all_predict_df.rename(columns = {0:"train_pred",
                                                  1:"test_pred",
                                                  2:"num_nodes",
                                                  3:"alpha",
                                                  4:"Win_strength",
                                                  5:"leak_rate",
                                                  6:"Win_sparsity",
                                                  7:"W_sparsity",
                                                  8:"subset_size"})
pd.set_option('display.max_columns', 10)
all_predict_df

