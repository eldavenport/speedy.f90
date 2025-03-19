import os
import sys
import xarray as xr
import re

def concatenate_netcdfs(folder_path, output_file="merged_output.nc"):
    """
    Concatenates all NetCDF files in the specified folder along the time dimension.
    
    Parameters:
        folder_path (str): Path to the folder containing NetCDF files.
        output_file (str): Name of the output merged NetCDF file (default: "merged_output.nc").
    """
    # Find all .nc files with numeric names in the directory
    nc_files = sorted([os.path.join(folder_path, f) for f in os.listdir(folder_path) 
                       if f.endswith(".nc") and re.fullmatch(r"\d+\.nc", f)])
    
    if not nc_files:
        print("No matching NetCDF files found in the directory.")
        return
    
    print(f"Found {len(nc_files)} NetCDF files. Concatenating...")
    
    # Open and concatenate along time dimension
    try:
        datasets = [xr.open_dataset(f) for f in nc_files]
        combined = xr.concat(datasets, dim="time")
        
        # Save to new NetCDF file in the same folder
        output_path = os.path.join(folder_path, output_file)
        combined.to_netcdf(output_path)
        print(f"Successfully saved concatenated NetCDF as {output_path}")
        
        # Close datasets
        for ds in datasets:
            ds.close()
    except Exception as e:
        print(f"Error during concatenation: {e}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python concatenate_netcdf.py <folder_path>")
        sys.exit(1)
    
    folder_path = sys.argv[1]
    concatenate_netcdfs(folder_path)
