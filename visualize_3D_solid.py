import os
import numpy as np
from skimage import measure
import cv2
import plotly.graph_objects as go


def load_images(folder_path):
    """
    Load and process image stack from a folder and put them into a numpy array.

    Args:
        folder_path (str): Path to the folder containing the image files.

    Returns:
        np.ndarray: A numpy array containing the loaded images.
    """
    file_paths = sorted(
        [os.path.join(folder_path, f) for f in os.listdir(folder_path) if f.endswith('.tif')]
    )
    images = []
    for file_path in file_paths:
        img = cv2.imread(file_path, cv2.IMREAD_UNCHANGED)
        if img is not None:
            images.append(img)
    return np.array(images)


def extract_mesh(volume, threshold=1):
    """
    Extracts isosurface using the Marching Cubes algorithm.

    Args:
        volume (np.ndarray): The 3D volume data from which to extract the isosurface.
        threshold (float, optional): The threshold value to use for the isosurface extraction. Defaults to 1.

    Returns:
        tuple: A tuple containing the vertices, faces, and values of the extracted isosurface.
    """
    verts, faces, normals, values = measure.marching_cubes(volume, level=threshold)
    return verts, faces, values


def create_figure(verts, faces, values, plane_position, plane_normal):
    """Create the Plotly figure with a toggle for the slicing plane."""
    x, y, z = verts[:, 0], verts[:, 1], verts[:, 2]
    i, j, k = faces.T
    min_value, max_value = values.min(), values.max()
    normalized_values = (values - min_value) / (max_value - min_value)

    a, b, c = plane_normal
    d = -np.dot(plane_normal, plane_position)
    xx, yy = np.meshgrid(
        np.linspace(min(x), max(x), 50), np.linspace(min(y), max(y), 50)
    )
    zz = (-a * xx - b * yy - d) / c

    fig = go.Figure()

    fig.add_trace(go.Mesh3d(
        x=x, y=y, z=z,
        i=i, j=j, k=k,
        intensity=normalized_values,
        colorscale='Viridis',
        showscale=False,
        opacity=1.0,
        name="Isosurface",
        visible=True
    ))

    # Add the slicing plane
    fig.add_trace(go.Surface(
        x=xx, y=yy, z=zz,
        colorscale=[[0, 'red'], [1, 'red']],
        showscale=False,
        opacity=1.0,
        name="Slicing Plane",
        visible=False
    ))

    # Add toggle buttons
    fig.update_layout(
        scene=dict(
            xaxis=dict(visible=False),
            yaxis=dict(visible=False),
            zaxis=dict(visible=False),
        ),
        margin=dict(l=0, r=0, t=0, b=0),
    )

    return fig


# Main Script
folder_path = "/Users/alexandergadin/Desktop/Group 1/MRI 1"
volume = load_images(folder_path)

# Extract the mesh
threshold_value = 1
verts, faces, values = extract_mesh(volume, threshold=threshold_value)

# Define plane properties
plane_position = [0, 0, 0]
plane_normal = [1, 0, 0]

# Figure
fig = create_figure(verts, faces, values, plane_position, plane_normal)
fig.show()
