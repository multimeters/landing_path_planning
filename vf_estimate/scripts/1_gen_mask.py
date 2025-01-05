import numpy as np
import os
from scipy.io import savemat
from PIL import Image

def generate_mask_files_and_image(
    input_txt='path/to/your/input_file.txt',  # Configure your input file path here
    output_dir='../data/processing/mask',
    threshold=0.75,
    n_rows=256,   # Configure number of rows (e.g., 256 or 512)
    n_cols=512    # Configure number of columns (e.g., 512)
):
    """
    Reads data from a specified text file with dimensions n_rows x n_cols,
    applies a threshold to generate a binary mask,
    and saves the mask as a .mat file and a PNG image.
    
    Parameters:
    - input_txt: str, path to the input text file.
    - output_dir: str, directory to save output files.
    - threshold: float, threshold value to create the mask.
    - n_rows: int, number of rows in the input data.
    - n_cols: int, number of columns in the input data.
    """

    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Read the text file
    try:
        # Specify the expected shape to prevent errors
        data = np.loadtxt(input_txt).reshape(n_rows, n_cols)
        print(f"Loaded data from {input_txt} with shape {data.shape}")
    except Exception as e:
        print(f"Error reading the input file: {e}")
        return

    # Generate the binary mask based on the threshold
    mask = (data > threshold).astype(np.uint8)
    print(f"Generated mask with threshold {threshold}")

    # Save the mask as a .mat file
    mat_save_path = os.path.join(output_dir, 'mask.mat')
    savemat(mat_save_path, {'mask': mask})
    print(f"Saved mask MAT file to: {mat_save_path}")

    # Save the mask as a PNG image
    img = Image.fromarray(mask * 255).convert('L')
    img_save_path = os.path.join(output_dir, 'mask.png')
    img.save(img_save_path)
    print(f"Saved mask image to: {img_save_path}")

    # (Optional) If you want to save the mask as a PNG with different naming
    # png_filename = "mask_full.png"
    # png_save_path = os.path.join(output_dir, png_filename)
    # img.save(png_save_path)
    # print(f"Saved mask PNG to: {png_save_path}")

if __name__ == "__main__":
    # Configuration Section
    INPUT_FILE_PATH = '../data/source/depth_a15_512x512.txt'  # Update this path
    OUTPUT_DIRECTORY = '../data/processing/mask'               # Desired output directory
    THRESHOLD_VALUE = 0.75                         # Threshold for mask generation
    NUM_ROWS = 256                                  # Number of rows in input data
    NUM_COLS = 512                                  # Number of columns in input data

    generate_mask_files_and_image(
        input_txt=INPUT_FILE_PATH,
        output_dir=OUTPUT_DIRECTORY,
        threshold=THRESHOLD_VALUE,
        n_rows=NUM_ROWS,
        n_cols=NUM_COLS
    )
