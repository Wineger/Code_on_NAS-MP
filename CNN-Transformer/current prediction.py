import xarray as xr
import numpy as np
import torch
import torch.nn as nn
from torch.utils.data import Dataset

# load data
file_path = './current_data.nc'
data = xr.open_dataset(file_path)

input_window = 80  
forecast_window = 40  
hidden_size = 512  

# standardization
uo_mean, uo_std = np.mean(data['uo'].values.flatten()), np.std(data['uo'].values.flatten())
vo_mean, vo_std = np.mean(data['vo'].values.flatten()), np.std(data['vo'].values.flatten())
uo_scaled = (data['uo'] - uo_mean) / uo_std
vo_scaled = (data['vo'] - vo_mean) / vo_std

# Partition data set
raw_size = uo_scaled.shape[1] * uo_scaled.shape[2] * uo_scaled.shape[3] * 2

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
                             dtype=np.float16)
        for i in range(len(data) - l + 1):
            self.data[i] = data[i * self.step:i * self.step + l]
        self.data = torch.tensor(self.data)

    def __getitem__(self, index):
        return self.data[index][:self.input_window].float(), self.data[index][self.input_window:].float()

    def __len__(self):
        return len(self.data)

# define model 
class CNNTransformer(nn.Module):
    def __init__(self):
        super(CNNTransformer, self).__init__()
       
        self.conv_layers = nn.Sequential(
            nn.Conv1d(in_channels=raw_size, out_channels=hidden_size, kernel_size=3, stride=1, padding=1),
            nn.ReLU(),
            nn.Conv1d(in_channels=hidden_size, out_channels=raw_size, kernel_size=5, stride=1, padding=2),
            nn.ReLU()
        )
        self.fc1 = nn.Linear(raw_size, hidden_size)
        self.dropout = nn.Dropout(0.2)
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

test_set = MyDataset(uo_scaled, vo_scaled, input_window, forecast_window)
t_placeholder = torch.zeros(1, forecast_window, raw_size)
model = CNNTransformer()
checkpoint = torch.load('./CNNTransformer.pth')
model.load_state_dict(checkpoint['model_state_dict'])
model.to('cpu')
model.eval()

# current prediction
predict = model(test_set[-1][0].unsqueeze(0), t_placeholder)
s = predict.size()
predict = predict.squeeze(0)
predict = predict.view(forecast_window,
                       2,
                       uo_scaled.shape[1],
                       uo_scaled.shape[2],
                       uo_scaled.shape[3])
predict = predict.detach().numpy()
predict[:,0,:,:,:] = predict[:,0,:,:,:]*uo_std + uo_mean
predict[:,1,:,:,:] = predict[:,1,:,:,:]*vo_std + vo_mean
data_array = xr.DataArray(predict, dims=['time', 'uo-vo', 'depth', 'latitude', 'longitude'])
data_array.to_netcdf("current_prediction.nc")
print(data_array)
