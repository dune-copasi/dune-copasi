import numpy as np

def read_gmsh(file_path):
    nodes = {}
    elements = []

    with open(file_path, 'r') as file:
        lines = file.readlines()

        reading_nodes = False
        reading_elements = False

        for line in lines:
            if line.startswith('$Nodes'):
                print("Start readingNodes ---")
                reading_nodes = True
                reading_elements = False
                continue
            elif line.startswith('$EndNodes'):
                reading_nodes = False
                continue
            elif line.startswith('$Elements'):
                print("Start reading Elements ---")
                reading_nodes = False
                reading_elements = True
                continue
            elif line.startswith('$EndElements'):
                reading_elements = False
                continue

            if reading_nodes:
                parts = line.strip().split()
                if len(parts) > 3:
                    node_id = int(parts[0])
                    x, y, z = map(float, parts[1:4])
                    nodes[node_id] = (x, y, z)

            if reading_elements:
                parts = line.strip().split()
                if len(parts) > 3:
                    element_type = int(parts[1])
                    # Assuming 2D elements with 3 nodes (change as needed)
                    if element_type == 2:  # Change this number based on the element type
                        element_nodes = list(map(int, parts[5:8]))
                        elements.append(element_nodes)

    return nodes, elements

# Function to generate lines with integers and doubles
def generate_lines(nodes, elements):
    with open('disk_griddata.txt', 'w') as file:
        line = "# --- griddata dune-copasi -----\n"
        file.write(line)
        line = str(len(elements)) + "\n"
        file.write(line)

        for i in range(1, len(elements)+1):
            # Generating the line with integer i and double i*1.5
            p = np.array([0.,0.,0.])
            for j in range(0, 3):
                print(f"{elements[i-1][j]} : {np.array(nodes[elements[i-1][j]])}")
                p += np.array(nodes[elements[i-1][j]])
            val = p / 3;
            print("==== ")
            print(val)
            line = f"{i} {val[0]:.4f}\n"  # Adjust decimal places as needed
            file.write(line)

# Example usage:
file_path = 'test_cell.msh'  # Replace this with your Gmsh file path
nodes, elements = read_gmsh(file_path)

generate_lines(nodes, elements)
print(f"output written to disk_griddata.txt")
