import xarray as xr
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import DataLoader, Dataset
from torch.optim import Adam
from tqdm import tqdm
import matplotlib.pyplot as plt

# Load data
file_path = './current_data.nc'
data = xr.open_dataset(file_path)

# Parameter settings
input_window = 80
forecast_window = 40
epochs = 20
batch_size = 16
hidden_size = 512
learning_rate = 1e-5

# Data normalization
uo_mean, uo_std = np.mean(data['uo'].values.flatten()), np.std(data['uo'].values.flatten())
vo_mean, vo_std = np.mean(data['vo'].values.flatten()), np.std(data['vo'].values.flatten())
uo_scaled = (data['uo'] - uo_mean) / uo_std
vo_scaled = (data['vo'] - vo_mean) / vo_std

# Split dataset
raw_size = uo_scaled.shape[1] * uo_scaled.shape[2] * uo_scaled.shape[3] * 2
l = len(uo_scaled)
train_uo = uo_scaled[:int(l * 0.8)]
train_vo = vo_scaled[:int(l * 0.8)]
test_uo = uo_scaled[int(l * 0.8):]
test_vo = vo_scaled[int(l * 0.8):]

# Create datasets for model training
class MyDataset(Dataset):
    def __init__(self, uo_data, vo_data, input_window, forecast_window, step=1):
        self.depth = uo_data.depth
        self.longitude = uo_data.longitude
        self.latitude = uo_data.latitude
        self.input_window = input_window
        self.forecast_window = forecast_window
        self.step = step
        self.data = None
        self.preprocess(uo_data, vo_data)

    def preprocess(self, uo_data, vo_data):
        data = np.column_stack((uo_data.values.reshape(uo_data.shape[0], -1),
                                vo_data.values.reshape(vo_data.shape[0], -1)))
        l = self.forecast_window + self.input_window
        self.data = np.empty((len(data) - l + 1,
                              l,
                              data.shape[-1]),
                             dtype=np.float32)

        for i in range(len(data) - l + 1):
            self.data[i] = data[i * self.step:i * self.step + l]

        self.data = torch.tensor(self.data)
        print("Data size:", self.data.size())

    def __getitem__(self, index):
        return self.data[index][:self.input_window], self.data[index][self.input_window:]

    def __len__(self):
        return len(self.data)

# Define model
class CNNTransformer(nn.Module):
    def __init__(self):
        super(CNNTransformer, self).__init__()
        # Define convolutional layers
        self.conv_layers = nn.Sequential(
            nn.Conv1d(in_channels=raw_size, out_channels=hidden_size, kernel_size=3, stride=1, padding=1),
            nn.ReLU(),
            nn.Conv1d(in_channels=hidden_size, out_channels=raw_size, kernel_size=5, stride=1, padding=2),
            nn.ReLU()
        )
        self.fc1 = nn.Linear(raw_size, hidden_size)
        self.dropout = nn.Dropout(0.2)
        # Define Transformer model layer
        self.transformer = nn.Transformer(d_model=hidden_size, nhead=4, num_encoder_layers=2, num_decoder_layers=2)
        self.fc2 = nn.Linear(hidden_size, raw_size)

    def forward(self, x, t):
        x = self.conv_layers(x.permute(0, 2, 1))
        x = self.fc1(x.permute(0, 2, 1))
        y = t.clone()
        y = self.conv_layers(y.permute(0, 2, 1))
        y = self.fc1(y.permute(0, 2, 1))
        x = self.dropout(x)
        y = self.dropout(y)
        x = self.transformer(x.permute(1, 0, 2), y.permute(1, 0, 2))
        x = self.fc2(x.permute(1, 0, 2))
        return x

train_set = MyDataset(train_uo, train_vo, input_window, forecast_window)
test_set = MyDataset(test_uo, test_vo, input_window, forecast_window)

# Specify device
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
print("Use device:", device)

# Initialize model
model = CNNTransformer().to(device)

# Define loss function and optimizer
criterion = nn.MSELoss()
optimizer = Adam(model.parameters(), lr=learning_rate, weight_decay= 5e-6)

# Convert data to PyTorch DataLoader
train_loader = DataLoader(train_set, batch_size=batch_size, shuffle=True)
test_loader = DataLoader(test_set, batch_size=batch_size, shuffle=False)

# Store loss values for plotting
train_loss_list = []
test_loss_list = []

# Train model
model.train()
for epoch in range(epochs):
    model.train()
    model.to("cuda")
    train_loss_epoch = 0
    for inputs, targets in tqdm(train_loader, total=len(train_loader)):
        inputs = inputs.to(device)
        targets = targets.to(device)
        outputs = model(inputs, targets)
        train_loss = criterion(outputs, targets)
        optimizer.zero_grad()
        train_loss.backward()
        optimizer.step()
        train_loss_epoch += train_loss.item()
    avg_train_loss = train_loss_epoch / len(train_loader)
    train_loss_list.append(avg_train_loss)
    # Calculate validation set loss
    model.eval()
    model.to("cpu")
    with torch.no_grad():
        test_loss_epoch = 0
        for inputs, targets in tqdm(test_loader, total=len(test_loader)):
            outputs = model(inputs, targets)
            test_loss = F.mse_loss(outputs, targets)
            test_loss_epoch += test_loss.item()
    avg_test_loss = test_loss_epoch / len(test_loader)
    test_loss_list.append(avg_test_loss)
    print((f'Epoch [{epoch + 1}/{epochs}], Train_loss: {avg_train_loss:.4f}, Test_loss: {avg_test_loss:.4f}'))

# Save model
torch.save({'model_state_dict': model.state_dict()}, "./CNNTransformer.pth")
# Plot training and validation loss
plt.figure(figsize=(10, 5))
plt.plot(train_loss_list, label='Train Loss')
plt.plot(test_loss_list, label='Test Loss')
plt.xlabel('Epoch')
plt.ylabel('Loss')
plt.legend()
plt.savefig("./loss.png")
plt.show()
