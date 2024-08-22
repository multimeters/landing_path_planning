import cupy as cp
import numpy as np
import scipy.io as sio
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

def load_data(file_idx):
    file_name = f'/980pro/pixel_data/e_output/data_{file_idx * 256}-{(file_idx + 1) * 256 - 1}.mat'
    data = sio.loadmat(file_name)['data']
    return data

# Load data and convert to CuPy arrays
encountered_frequency_coefficient = cp.array(sio.loadmat('encountered_frequency_coefficient_4096_2048.mat')['encountered_frequency_coefficient'])
motionRAO_data = sio.loadmat('heave_motionRAO_at_0.mat')
motionRAO_w = cp.array(motionRAO_data['motionRAO_w'].flatten())
motionRAO_amp_at_0 = cp.array(motionRAO_data['motionRAO_amp_at_0'].flatten())

# Set parameters
N = 1000
fs = 10  # Sampling frequency

# Initialize the total sum matrix
total_sum_matrix = cp.full((4096, 2048), cp.nan)  # Initialize with NaN

start_time = time.time()  # Start timing

# Use ThreadPoolExecutor to load files in the background
with ThreadPoolExecutor(max_workers=2) as executor:
    future_to_file_idx = {executor.submit(load_data, file_idx): file_idx for file_idx in range(2)}
    
    for file_idx in range(16):
        current_data_future = next(as_completed(future_to_file_idx))  # Wait for the next completed future
        data = current_data_future.result()
        current_idx = future_to_file_idx.pop(current_data_future)
        print(file_idx)
        # Submit the next file loading task, if applicable
        if current_idx + 2 < 16:
            future_to_file_idx[executor.submit(load_data, current_idx + 2)] = current_idx + 2
        
        # Process data in chunks
        for row_chunk_start in range(0, 256, 256):  # Process 256 rows at a time
            row_chunk_end = min(row_chunk_start + 256, 256)
            
            for col_chunk_start in range(0, 2048, 256):  # Process 256 columns at a time
                col_chunk_end = min(col_chunk_start + 256, 2048)
                
                chunk_start_time = time.time()  # Start timing for the chunk

                data_chunk = data[row_chunk_start:row_chunk_end, col_chunk_start:col_chunk_end, :]
                
                # Convert chunk data to CuPy array for GPU processing
                data_gpu = cp.array(data_chunk, dtype=cp.float32)

                # Compute FFT for the current chunk
                fft_result = cp.fft.fft(data_gpu, axis=-1)
                fft_result = fft_result[:, :, :N // 2 + 1]

                # Calculate power spectral density (PSD)
                psd = (1 / (2 * cp.pi * N)) * cp.abs(fft_result) ** 2
                psd[:, :, 1:-1] *= 2

                # Frequency vector
                frequencies = (fs * cp.arange(N // 2 + 1)) / N

                # Fetch the encountered frequency coefficient chunk
                encountered_frequency_coefficient_chunk = encountered_frequency_coefficient[current_idx * 256 + row_chunk_start:current_idx * 256 + row_chunk_end, col_chunk_start:col_chunk_end]

                # Check for NaN values in the encountered frequency coefficient
                nan_mask = cp.isnan(encountered_frequency_coefficient_chunk)
                
                # Adjust frequencies based on encountered frequency coefficient
                adjusted_frequencies = frequencies[None, :] * encountered_frequency_coefficient_chunk[:, :, None]

                # Initialize the result matrix for the current chunk
                chunk_result = cp.zeros((row_chunk_end - row_chunk_start, col_chunk_end - col_chunk_start))

                if not nan_mask.all():  # If not all are NaNs
                    # Process frequencies for the chunk
                    closest_indices = cp.empty(adjusted_frequencies.shape, dtype=cp.int32)
                    for freq_idx in range(adjusted_frequencies.shape[2]):
                        closest_indices[:, :, freq_idx] = cp.argmin(cp.abs(motionRAO_w[None, None, :] - adjusted_frequencies[:, :, freq_idx, None]), axis=-1)

                    # Get corresponding amplitudes from motionRAO_amp_at_0
                    resulting_amplitudes = motionRAO_amp_at_0[closest_indices]

                    # Compute the product of amplitude squared and PSD
                    amplitude_squared = resulting_amplitudes ** 2
                    product = amplitude_squared * psd

                    # Sum over the frequencies and store in chunk_result
                    chunk_result = cp.sum(product, axis=-1)

                    # Assign NaN to positions where encountered_frequency_coefficient is NaN
                    chunk_result[nan_mask] = cp.nan
                else:
                # 如果全是 NaN，则直接将 chunk_result 赋值为 NaN
                    chunk_result[:] = cp.nan
                # Accumulate the result into total_sum_matrix
                total_sum_matrix[current_idx * 256 + row_chunk_start:current_idx * 256 + row_chunk_end, col_chunk_start:col_chunk_end] = chunk_result

                # Free memory after processing each chunk
                cp.get_default_memory_pool().free_all_blocks()

                # Calculate and print the time taken for the chunk
                chunk_end_time = time.time()
                chunk_elapsed_time = chunk_end_time - chunk_start_time
                print(f'Chunk {row_chunk_start}-{row_chunk_end}, {col_chunk_start}-{col_chunk_end} processed in {chunk_elapsed_time:.2f} seconds')

end_time = time.time()  # End timing

# Compute and display total processing time
elapsed_time = end_time - start_time
print(f'Processing complete. Total sums calculated for each pixel in {elapsed_time:.2f} seconds.')

# Convert total_sum_matrix to a NumPy array and save as a .mat file
total_sum_matrix_np = cp.asnumpy(total_sum_matrix)
sio.savemat('total_sum_matrix.mat', {'total_sum_matrix': total_sum_matrix_np})

