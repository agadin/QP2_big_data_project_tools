import os
import numpy as np
import cv2
import trimesh
from skimage import measure

def load_image_stack(folder_path):
    image_files = sorted([f for f in os.listdir(folder_path) if f.endswith(('.png', '.jpg', '.tif'))])
    stack = []
    for image_file in image_files:
        image = cv2.imread(os.path.join(folder_path, image_file), cv2.IMREAD_GRAYSCALE)
        if image is not None:
            stack.append(image)
    return np.array(stack)

def create_mesh_from_stack(image_stack, folder_path):
    binary_stack = (image_stack > 0).astype(np.uint8)
    spacing = [5, 10 / 17.53, 10 / 17.53] if "CT" in folder_path else [1, 10 / 17.53, 10 / 17.53]
    verts, faces, _, _ = measure.marching_cubes(binary_stack, level=0.5, spacing=spacing)
    mesh = trimesh.Trimesh(vertices=verts, faces=faces)
    return mesh

def save_mesh(mesh, output_path):
    if os.path.exists(output_path):
        os.remove(output_path)
    mesh.export(output_path)

def load_mesh_from_obj(file_path):
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"The file {file_path} does not exist.")
    return trimesh.load(file_path)

# Step 2: Analyze Surface Properties
def analyze_surface_properties_from_obj(obj_path):
    mesh = load_mesh_from_obj(obj_path)

    # Calculate volume and convert to mL
    volume_ml = abs(mesh.volume) / 1000

    bounding_box = mesh.bounds
    bounding_box_dims = {
        'x_dim_mm': bounding_box[1][0] - bounding_box[0][0],
        'y_dim_mm': bounding_box[1][1] - bounding_box[0][1],
        'z_dim_mm': bounding_box[1][2] - bounding_box[0][2]
    }

    avg_edge_length = np.mean(mesh.edges_unique_length)
    max_edge_length = np.max(mesh.edges_unique_length)
    min_edge_length = np.min(mesh.edges_unique_length)

    euler_number = mesh.euler_number
    non_manifold_edges = mesh.edges_unique.shape[0] - mesh.edges.shape[0]

    properties = {
        'surface_area_mm2': mesh.area,  # Surface area in square mm
        'volume_ml': volume_ml,  # Volume in mL
        'bounding_box_dims_mm': bounding_box_dims,  # Bounding box dimensions in mm
        'num_triangles': len(mesh.faces),  # Number of triangles in the mesh
        'avg_edge_length_mm': avg_edge_length,  # Average edge length in mm
        'max_edge_length_mm': max_edge_length,  # Maximum edge length in mm
        'min_edge_length_mm': min_edge_length,  # Minimum edge length in mm
        'euler_number': euler_number,  # Topological measure
        'non_manifold_edges': non_manifold_edges  # Number of non-manifold edges
    }
    return properties

def main():
    folder_path = "./sample_group/CT_1"
    output_path = "output_mesh.obj"

    image_stack = load_image_stack(folder_path)

    mesh = create_mesh_from_stack(image_stack, folder_path)
    save_mesh(mesh, output_path)
    print(f"Mesh saved to {output_path}")

    properties = analyze_surface_properties_from_obj(output_path)
    print("Surface Properties from OBJ:", properties)

if __name__ == "__main__":
    main()
