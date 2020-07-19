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
hps = {'lr': 0.01, 'n_conv_layers': 2, 'kernel_size': 3, 'n_kernels': 8}
# n_conv_layers: 1,2,4
# kernel_size: 3,5
# n_kernels: 2,4,8
# learning_rate: between 0.0001 and 0.1

# ---- options ----
DEVICE_ID = 'cpu' # set to 'cpu' for cpu, 'cuda' / 'cuda:0' or similar for gpu.
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
from utils import read_image

torch.multiprocessing.set_start_method('spawn', force=True)


# ---- utility functions -----
def get_loaders(batch_size,device):
    data_root = 'ceng483-s19-hw3-dataset'
    train_set = hw3utils.HW3ImageFolder(root=os.path.join(data_root,'train'),device=device)
    train_loader = torch.utils.data.DataLoader(train_set, batch_size=batch_size, shuffle=True, num_workers=0)
    val_set = hw3utils.HW3ImageFolder(root=os.path.join(data_root,'val'),device=device)
    val_loader = torch.utils.data.DataLoader(val_set, batch_size=batch_size, shuffle=False, num_workers=0)
    test_set = hw3utils.HW3ImageFolder(root=os.path.join(data_root, 'test'), device=device)
    test_loader = torch.utils.data.DataLoader(test_set, batch_size=batch_size, shuffle=False, num_workers=0)
    return train_loader, val_loader, test_loader

# ---- ConvNet -----
class Net(nn.Module):
    def __init__(self):
        super(Net, self).__init__()
        self.n_conv_layers = hps['n_conv_layers']
        self.n_kernels = hps['n_kernels']
        self.kernel_size = hps['kernel_size']
        self.padding = (self.kernel_size - 1) // 2

        # nn.conv2d(in_channels, out_channels, kernel_size, stride, padding, ...)
        self.four_layers = nn.Sequential(
                            nn.Conv2d(1, self.n_kernels, self.kernel_size, padding=(self.padding, self.padding)),
                            #nn.BatchNorm2d(self.n_kernels),
                            nn.ReLU(),
                            nn.Conv2d(self.n_kernels, self.n_kernels, self.kernel_size, padding=(self.padding, self.padding)),
                            #nn.BatchNorm2d(self.n_kernels),
                            nn.ReLU(),
                            nn.Conv2d(self.n_kernels, self.n_kernels, self.kernel_size, padding=(self.padding, self.padding)),
                            #nn.BatchNorm2d(self.n_kernels),
                            nn.ReLU(),
                            nn.Conv2d(self.n_kernels, 3, self.kernel_size, padding=(self.padding, self.padding)))
        self.two_layers = nn.Sequential(
                            nn.Conv2d(1, self.n_kernels, self.kernel_size, padding=(self.padding, self.padding)),
                            #nn.BatchNorm2d(self.n_kernels),
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
#print('device: ' + str(device))
net = Net().to(device=device)
criterion = nn.MSELoss()
optimizer = optim.SGD(net.parameters(), lr=hps['lr'])
train_loader, val_loader, test_loader = get_loaders(batch_size,device)

if LOAD_CHKPT:
    #print('loading the model from the checkpoint')
    model.load_state_dict(os.path.join(LOG_DIR,'checkpoint.pt'))

prev_val_loss = float("inf")
optimal_epoch = 0


#print('training begins')
for epoch in range(max_num_epoch):
    running_loss = 0.0 # training loss of the network
    n_training = 0
    for iteri, data in enumerate(train_loader, 0):
        n_training += 1
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

        if (iteri==0) and VISUALIZE:
            hw3utils.visualize_batch(inputs,preds,targets)

    #print(epoch + 1, running_loss / n_training)
    '''acc = 0
    with torch.no_grad():
      for i, val_data in enumerate(val_loader, 0):
          val_inputs, val_targets = val_data
          val_preds = net(val_inputs)

          for j, pred in enumerate(val_preds):
              pred = pred.to(torch.float64)
              val_target = val_targets[j]
              val_target = val_target.to(torch.float64)
              colored_target = (val_target/2 + 0.5) * 255
              colored_target = colored_target.permute(1,2,0)
              colored_pred = (pred/2 + 0.5) * 255             # First take to range [0, 1], then to [0, 255]
              colored_pred = colored_pred.permute(1,2,0)      # Convert too 80x80x3
              #print(colored_pred.shape, colored_target.shape)
              est = colored_pred.cpu().numpy().astype(np.int64)
              cur = colored_target.cpu().numpy().astype(np.int64)

              cur_acc = (np.abs(cur - est) < 12).sum() / cur.shape[0]
              acc += cur_acc

    acc /= 5001
    print(epoch + 1, acc)'''


    if epoch % 5 == 4:
        # Compute average validation loss every 5 epochs by a full pass over the validation set.
        val_running_loss = 0.0
        n_validation = 0
        for i, val_data in enumerate(val_loader, 0):
            n_validation += 1
            val_inputs, val_targets = val_data
            val_preds = net(val_inputs)                         # get validation outputs
            val_loss = criterion(val_preds, val_targets)        # get loss for each mini-batch(16 images)
            val_running_loss += val_loss.item()

#        print("Epoch", epoch + 1, "is over. Validation set loss:", val_running_loss / 125)
#        print(epoch + 1, val_running_loss / n_validation)
        # If loss has increased, apply early stopping.
        if prev_val_loss < val_running_loss:
            optimal_epoch = epoch - 5
#            print("optimal number of epochs is:", optimal_epoch)
#            print("prev_val_loss is:", prev_val_loss / 125)
#            print("current loss is:", val_running_loss / 125)
            break
        else:
            # If current loss < prev loss, then save the model
            #print('Saving the model, end of epoch %d' % (epoch+1))
            if not os.path.exists(LOG_DIR):
                os.makedirs(LOG_DIR)
            torch.save(net.state_dict(), os.path.join(LOG_DIR,'checkpoint.pt'))
            hw3utils.visualize_batch(inputs,preds,targets,os.path.join(LOG_DIR,'example.png'))

        prev_val_loss = val_running_loss

#print('Finished Training')

# number of validation images: 2000, size: 80x80x3
# number of test images: 2000, size: 80x80x3
validation_estimations = np.zeros((2000, 80, 80, 3))
test_estimations = np.zeros((2000, 80, 80, 3))


# One full pass over the validation set
with torch.no_grad():           # Run model without backpropagation
    for i, data in enumerate(val_loader):
        inputs, targets = data
        preds = net(inputs)

        for j, pred in enumerate(preds):
            pred = pred.to(torch.float64)
            colored_pred = (pred/2 + 0.5) * 255             # First take to range [0, 1], then to [0, 255]
            colored_pred = colored_pred.permute(1,2,0)      # Convert too 80x80x3
            validation_estimations[i * batch_size + j] = colored_pred.cpu().numpy()

print("val est", validation_estimations)

# One full pass over the test set
with torch.no_grad():           # Run model without backpropagation
    for i, data in enumerate(test_loader):
        inputs, targets = data
        preds = net(inputs)

        for j, pred in enumerate(preds):
            pred = pred.to(torch.float64)
            colored_pred = (pred/2 + 0.5) * 255             # First take to range [0, 1], then to [0, 255]
            colored_pred = colored_pred.permute(1,2,0)      # Convert too 80x80x3
            test_estimations[i * batch_size + j] = colored_pred.cpu().numpy()
            print("test shape:", colored_pred.shape)


np.save("estimations_validation.npy", validation_estimations)
np.save("estimations_test.npy", test_estimations)
