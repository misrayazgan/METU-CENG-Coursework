# Feel free to change / extend / adapt this source code as needed to complete the homework, based on its requirements.
# This code is given as a starting point.
#
# REFEFERENCES
# The code is partly adapted from pytorch tutorials, including https://pytorch.org/tutorials/beginner/blitz/cifar10_tutorial.html

# ---- hyper-parameters ----
# You should tune these hyper-parameters using:
# (i) your reasoning and observations,
# (ii) by tuning it on the validation set, using the techniques discussed in class.
# You definitely can add more hyper-parameters here.
batch_size = 16
max_num_epoch = 100
hps = {'lr': 0.001, 'n_conv_layers': 1, 'kernel_size': 3, 'n_kernels': 2}
# n_conv_layers: 1,2,4
# kernel_size: 3,5
# n_kernels: 2,4,8
# learning_rate: between 0.0001 and 0.1

# ---- options ----
DEVICE_ID = 'cuda' # set to 'cpu' for cpu, 'cuda' / 'cuda:0' or similar for gpu.
LOG_DIR = 'checkpoints'
VISUALIZE = False # set True to visualize input, prediction and the output from the last batch
LOAD_CHKPT = False

# --- imports ---
import torch
import os
import matplotlib.pyplot as plt
import numpy as np
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
import torchvision.transforms as transforms
import hw3utils

torch.multiprocessing.set_start_method('spawn', force=True)


# ---- utility functions -----
def get_loaders(batch_size,device):
    data_root = 'ceng483-s19-hw3-dataset'
    train_set = hw3utils.HW3ImageFolder(root=os.path.join(data_root,'train'),device=device)
    train_loader = torch.utils.data.DataLoader(train_set, batch_size=batch_size, shuffle=True, num_workers=0)
    val_set = hw3utils.HW3ImageFolder(root=os.path.join(data_root,'val'),device=device)
    val_loader = torch.utils.data.DataLoader(val_set, batch_size=batch_size, shuffle=False, num_workers=0)
    # Note: you may later add test_loader to here.

    # test_set = hw3utils.HW3ImageFolder(root=os.path.join(data_root, 'test'), device=device)
    # test_loader = torch.utils.data.DataLoader(test_set, batch_size=batch_size, shuffle=False, num_workers=0)
    return train_loader, val_loader # , test_loader

# ---- ConvNet -----
class Net(nn.Module):
    def __init__(self):
        super(Net, self).__init__()
        self.n_conv_layers = hps['n_conv_layers']
        self.n_kernels = hps['n_kernels']
        self.kernel_size = hps['kernel_size']
        self.padding = (self.kernel_size - 1) // 2
        print(self.padding)

        # nn.conv2d(in_channels, out_channels, kernel_size, stride, padding, ...)
        self.four_layers = nn.Sequential(
                            nn.Conv2d(1, self.n_kernels, self.kernel_size, padding=(self.padding, self.padding)),
                            nn.ReLU(),
                            nn.Conv2d(self.n_kernels, self.n_kernels, self.kernel_size, padding=(self.padding, self.padding)),
                            nn.ReLU(),
                            nn.Conv2d(self.n_kernels, self.n_kernels, self.kernel_size, padding=(self.padding, self.padding)),
                            nn.ReLU(),
                            nn.Conv2d(self.n_kernels, 3, self.kernel_size, padding=(self.padding, self.padding)))
        self.two_layers = nn.Sequential(
                            nn.Conv2d(1, self.n_kernels, self.kernel_size, padding=(self.padding, self.padding)),
                            nn.ReLU(),
                            nn.Conv2d(self.n_kernels, 3, self.kernel_size, padding=(self.padding, self.padding)))
        self.one_layer = nn.Sequential(
                            nn.Conv2d(1, 3, self.kernel_size, padding=(self.padding, self.padding)))

    def forward(self, grayscale_image):
        # apply your network's layers in the following lines:
        if self.n_conv_layers == 4:
            x = self.four_layers(grayscale_image)
        elif self.n_conv_layers == 2:
            x = self.two_layers(grayscale_image)
        elif self.n_conv_layers == 1:
            x = self.one_layer(grayscale_image)

        return x

# ---- training code -----
device = torch.device(DEVICE_ID)
print('device: ' + str(device))
net = Net().to(device=device)
criterion = nn.MSELoss()
optimizer = optim.SGD(net.parameters(), lr=hps['lr'])
train_loader, val_loader = get_loaders(batch_size,device)       # test_loader ekle

if LOAD_CHKPT:
    print('loading the model from the checkpoint')
    model.load_state_dict(os.path.join(LOG_DIR,'checkpoint.pt'))

prev_val_loss = float("inf")
optimal_epoch = 0
stop = False

print('training begins')
for epoch in range(20):
    running_loss = 0.0 # training loss of the network
    for iteri, data in enumerate(train_loader, 0):
        #print(" train iteri:", iteri)
        inputs, targets = data # inputs: low-resolution images(grayscale), targets: high-resolution images(rgb).
        optimizer.zero_grad() # zero the parameter gradients

        # do forward, backward, SGD step
        preds = net(inputs)                     # get train outputs
        loss = criterion(preds, targets)        # get loss for each mini-batch(16 images)
        loss.backward()
        optimizer.step()

        # print loss
        running_loss += loss.item()
        print_n = 100 # feel free to change this constant
        if iteri % print_n == (print_n-1):    # print every print_n mini-batches(16 images)
            print('[%d, %5d] network-loss: %.3f' %
                  (epoch + 1, iteri + 1, running_loss / 100))
            running_loss = 0.0
            # note: you most probably want to track the progress on the validation set as well (needs to be implemented)

        if (iteri==0) and VISUALIZE:
            hw3utils.visualize_batch(inputs,preds,targets)

    if epoch % 5 == 0:
        # Compute average validation loss every 5 epochs by a full pass over the validation set.
        val_running_loss = 0.0
        for i, val_data in enumerate(val_loader, 0):
            #print("validation i:", i)
            val_inputs, val_targets = val_data
            val_preds = net(val_inputs)                         # get validation outputs
            val_loss = criterion(val_preds, val_targets)        # get loss for each mini-batch(16 images)
            val_running_loss += val_loss.item()

        print("Epoch", epoch + 1, "is over. Validation set loss:", val_running_loss / 125)             #????????????

        # If loss has increased, apply early stopping.
        if prev_val_loss < val_running_loss:
            stop = True
            optimal_epoch = epoch - 5
            break
        else:
            # If current loss < prev loss, then save the model
            print('Saving the model, end of epoch %d' % (epoch+1))
            if not os.path.exists(LOG_DIR):
                os.makedirs(LOG_DIR)
            torch.save(net.state_dict(), os.path.join(LOG_DIR,'checkpoint.pt'))
            hw3utils.visualize_batch(inputs,preds,targets,os.path.join(LOG_DIR,'example.png'))

        prev_val_loss = val_running_loss

    if stop == True:
        print("optimal number of epochs is:", epoch)
        break;


print('Finished Training')

# number of validation images: 2000, size: 80x80x3
# number of test images: 2000, size: 80x80x3
validation_estimations = np.zeros((2000, 80, 80, 3))
test_estimations = np.zeros((2000, 80, 80, 3))


# One full pass over the validation set
with torch.no_grad():           # Run model without backpropagation
    for i, data in enumerate(val_loader):
        inputs, targets = data
        preds = net(inputs)
        print("i:", i, "preds shape:", preds.shape)

        for j, pred in enumerate(preds):    # preds.data:
            pred = pred.to(torch.float64)
            print("j:", j, "pred shape:", pred.shape)               # 3x80x80
            colored_pred = (pred/2 + 0.5) * 255             # First take to range [0, 1], then to [0, 255]
            colored_pred = colored_pred.permute(1,2,0)      # Convert too 80x80x3
            validation_estimations[i * batch_size + j] = colored_pred.cpu().numpy()
            print("colored pred shape:", colored_pred.shape)

print("val est", validation_estimations)


np.save("estimations_validation.npy", validation_estimations)
#np.save("estimations_test.npy", test_estimations)
test = np.load("estimations_validation.npy")
print(test.shape)
