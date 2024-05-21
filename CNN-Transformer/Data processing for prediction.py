import xarray as xr

# open file
predict_file = './current_prediction.nc'
predict_data = xr.open_dataset(predict_file)

# Extract uo and vo data
uo_vo_data = predict_data['__xarray_dataarray_variable__']
uo = uo_vo_data.sel({'uo-vo': 0})
vo = uo_vo_data.sel({'uo-vo': 1})

# Create coordinate information
time_coord = xr.DataArray(uo.coords['time'], dims='time')
depth_coord = xr.DataArray(uo.coords['depth'], dims='depth')
latitude_coord = xr.DataArray(uo.coords['latitude'], dims='latitude')
longitude_coord = xr.DataArray(uo.coords['longitude'], dims='longitude')

# Add uo and vo data variables and coordinate information
combined_dataset = xr.Dataset({'uo': uo, 'vo': vo},
                              coords={'time': time_coord,
                                      'depth': depth_coord,
                                      'latitude': latitude_coord,
                                      'longitude': longitude_coord})

# save data
combined_dataset.to_netcdf("predict_2232_8040with_coords.nc")
print("uo and vo shape：", combined_dataset.uo.shape, combined_dataset.vo.shape)
print("Finished saving uo and vo as separate variables in the same file with coordinates.")
