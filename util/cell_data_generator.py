#!/usr/bin/env python3

import numpy as np
import ast
import argparse
import sys


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


def generate_file(filename, expr, nodes, elements):

    ast_node = ast.parse(expr, mode='eval')
    expr_comp = compile(ast_node, '<string>', mode='eval')
    with open(filename, 'w') as file:
        line = "# --- griddata dune-copasi -----\n"
        line += f"# Automatically generated using '{" ".join(sys.argv)}'\n"
        file.write(line)
        line = str(len(elements)) + "\n"
        file.write(line)

        for i in range(1, len(elements)+1):
            # Generating the line with integer i and double i*1.5
            p = np.array([0., 0., 0.])
            for j in range(0, 3):
                print(
                    f"{elements[i-1][j]} : {np.array(nodes[elements[i-1][j]])}")
                p += np.array(nodes[elements[i-1][j]])
            val = p / 3
            [x, y, z] = val
            result = eval(expr_comp)
            print(f"pos := {val}; eval := {result}")
            line = f"{i-1} {result:.4f}\n"  # Adjust decimal places as needed
            file.write(line)
            print("==== ")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generates analytical cell data for dune-copasi from a GMSH file")
    parser.add_argument('-i', '--input', required=True,
                        help="GMSH file to read")
    parser.add_argument('-o', '--output', required=True,
                        help="file to write results")
    parser.add_argument('-e', '--expression', required=True,
                        help="expression with 'x', 'y', and 'z' to evaluate in each cell")
    args = parser.parse_args()
    nodes, elements = read_gmsh(args.input)

    generate_file(args.output, args.expression, nodes, elements)
    print(f"output written to {args.output}")
