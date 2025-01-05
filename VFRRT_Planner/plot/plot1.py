import numpy as np
import matplotlib
# Use the 'Agg' backend for environments without a display (e.g., servers)
matplotlib.use('Agg')  
import matplotlib.pyplot as plt
import os
import sys

def load_vector_field(filename):
    """
    Load vector field data from a file.
    
    Each line in the file should contain four values: x, y, u, v
    where (x, y) is the position and (u, v) is the vector at that position.
    
    Parameters:
        filename (str): Path to the vector field file.
    
    Returns:
        tuple: Arrays of x, y, u, v components or (None, None, None, None) on failure.
    """
    try:
        data = np.loadtxt(filename)
        if data.ndim == 1:
            data = data.reshape(1, -1)
        if data.shape[1] != 4:
            raise ValueError(f"Vector field file format error: Each line must contain 4 values (x y u v). Found {data.shape[1]} columns.")
        x, y, u, v = data[:,0], data[:,1], data[:,2], data[:,3]
        if len(x) == 0:
            raise ValueError("Vector field data is empty.")
        print(f"First 5 vector field points:\nX: {x[:5]}\nY: {y[:5]}\nU: {u[:5]}\nV: {v[:5]}")
        return x, y, u, v
    except Exception as e:
        print(f"Error loading vector field data: {e}")
        return None, None, None, None

def load_path(filename):
    """
    Load path data from a file.
    
    Each line in the file should contain two values: x, y
    where (x, y) is a point on the path.
    
    Parameters:
        filename (str): Path to the path file.
    
    Returns:
        tuple: Arrays of x, y coordinates or (None, None) on failure.
    """
    try:
        data = np.loadtxt(filename)
        if data.ndim == 1:
            data = data.reshape(1, -1)
        if data.shape[1] != 2:
            raise ValueError(f"Path file format error: Each line must contain 2 values (x y). Found {data.shape[1]} columns.")
        x, y = data[:,0], data[:,1]
        if len(x) == 0:
            raise ValueError("Path data is empty.")
        print(f"First 5 path points:\nX: {x[:5]}\nY: {y[:5]}")
        return x, y
    except Exception as e:
        print(f"Error loading path data: {e}")
        return None, None

def plot_vector_field_and_path(vector_field_file, path_file):
    """
    Plot the vector field and the planned path.
    
    Parameters:
        vector_field_file (str): Path to the vector field file.
        path_file (str): Path to the path file.
    """
    # Load vector field data
    x_vf, y_vf, u_vf, v_vf = load_vector_field(vector_field_file)
    if x_vf is None:
        print("Failed to load vector field data. Cannot plot.")
        return
    
    # Load path data
    x_path, y_path = load_path(path_file)
    if x_path is None:
        print("Failed to load path data. Cannot plot.")
        return
    
    # Check if data is empty
    if len(x_vf) == 0 or len(y_vf) == 0 or len(u_vf) == 0 or len(v_vf) == 0:
        print("Vector field data is empty. Cannot plot.")
        return
    if len(x_path) == 0 or len(y_path) == 0:
        print("Path data is empty. Cannot plot.")
        return
    
    print(f"Vector field data loaded successfully with {len(x_vf)} points.")
    print(f"Path data loaded successfully with {len(x_path)} points.")
    
    # Create a figure and axis
    plt.figure(figsize=(12, 8))
    
    # Plot the vector field using quiver
    plt.quiver(x_vf, y_vf, u_vf, v_vf, color='lightgray', alpha=0.5, scale=50, width=0.0025, label='Vector Field')
    
    # Plot the path
    plt.plot(x_path, y_path, 'r-', linewidth=2, label='Planned Path')
    
    # Mark the start and end points
    plt.plot(x_path[0], y_path[0], 'go', markersize=8, label='Start')
    plt.plot(x_path[-1], y_path[-1], 'bo', markersize=8, label='End')
    
    # Set axis labels and title
    plt.xlabel('X (meters)', fontsize=14)
    plt.ylabel('Y (meters)', fontsize=14)
    plt.title('Vector Field and Planned Path', fontsize=16)
    
    # Add legend
    plt.legend(fontsize=12)
    
    # Add grid
    plt.grid(True)
    
    # Ensure equal scaling on both axes
    plt.axis('equal')
    
    # Save the figure to a file
    output_image = "vector_field_and_path.png"
    plt.savefig(output_image, dpi=300)
    print(f"Figure saved as: {output_image}")
    
    # If you need to display the plot interactively, remove or comment out the 'Agg' backend and uncomment the following line:
    # plt.show()

if __name__ == "__main__":
    # Define the paths to the vector field and path files
    vector_field_file = "/home/ganzhi/VFRRT_Planner/vector_Field.txt"  # Vector field file
    path_file = "/home/ganzhi/VFRRT_Planner/path_exploration_0.100000_lambda_1.100000_samples_1000000_cost_55.894269.txt"  # Path file
    
    # Check if the vector field file exists
    if not os.path.isfile(vector_field_file):
        print(f"Vector field file does not exist: {vector_field_file}")
        sys.exit(1)
    
    # Check if the path file exists
    if not os.path.isfile(path_file):
        print(f"Path file does not exist: {path_file}")
        sys.exit(1)
    
    # Plot the vector field and path
    plot_vector_field_and_path(vector_field_file, path_file)
