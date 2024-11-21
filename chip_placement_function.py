
import pya

import random



#random_number_y = random.randint(500000, 1500000) #in nm
#random_number_x = random.randint(-10000000, -6000000) #in nm


random_number_x = 0 #in nm
random_number_y = 0 #in nm



#========================================= initial/reset stage===============
# Load the layout
layout = pya.Layout()
layout.read("C:\\Users\\limyu\\Downloads\\amf_grating_replaced.gds")

all_cells = []
for cell in layout.each_cell():
    cell_name = cell.name
    all_cells.append(cell_name)
#all_cells.remove('TOP')

OPA_index = all_cells.index('OPA2')
Ring_index = all_cells.index('Ring')

# Get the cells (assuming you know which ones have the shapes)
cell1 = layout.cell(Ring_index)  # First cell
cell2 = layout.cell(OPA_index)  # Second cell

# Get the top cell (assuming 'TOP' is the top-level cell)
top_cell = layout.cell('TOP')

# Function to get the origin of a specific cell instance within the top cell
def get_cell_origin(parent_cell, target_cell):
    for inst in parent_cell.each_inst():
        # Check if the instance is of the target cell
        if inst.cell.name == target_cell.name:
            # Get the transformation (trans) of the instance
            trans = inst.trans
            x_origin = trans.disp.x
            y_origin = trans.disp.y
            return x_origin, y_origin
    return None, None

# Find the origin of cell1 (Ring) inside the top cell
x_origin, y_origin = get_cell_origin(top_cell, cell1)

#========================================= action stage===============

#----------------------moves cell----------------------



def move_cell(cell, delta_x, delta_y):
    """
    Moves the cell by `delta_x` microns in the X direction and `delta_y` microns in the Y direction.
    This function applies the transformation to all shapes across all layers in the cell.
    """
    current_position_x = 0
    current_position_y = 0
    
    # Apply transformation relative to the current position (X and Y)
    transformation = pya.Trans(0, False, delta_x - current_position_x, delta_y - current_position_y)

    # Update the current position
    current_position_x = delta_x
    current_position_y = delta_y
    
    # Iterate over all possible layer indices (assuming layers are indexed starting from 0)
    for layer_index in range(cell.layout().layers()):
        # Check if the cell has shapes in the current layer
        if not cell.shapes(layer_index).is_empty():
            # Transform all shapes in this layer by the new transformation
            for shape in cell.shapes(layer_index).each():
                shape.transform(transformation)


# Function to rotate the cell around its center
def rotate_cell_around_center(cell, angle_degrees):
    # Get the bounding box of the cell
    bbox = cell.bbox()
    
    # Calculate the center of the bounding box
    center_x = (bbox.p1.x + bbox.p2.x) / 2
    center_y = (bbox.p1.y + bbox.p2.y) / 2
    
    # Create a transformation to move the cell to the origin
    to_origin = pya.Trans(-center_x, -center_y)
    
    # Create a rotation transformation (angle in degrees, counterclockwise)
    rotation = pya.Trans(pya.Trans.R90 if angle_degrees == 90 else angle_degrees, False)
    
    # Create a transformation to move the cell back to its original center
    back_to_center = pya.Trans(center_x, center_y)
    
    # Apply the transformations in sequence
    for layer_index in range(cell.layout().layers()):
        if not cell.shapes(layer_index).is_empty():
            for shape in cell.shapes(layer_index).each():
                shape.transform(to_origin)      # Move to origin
                shape.transform(rotation)       # Rotate
                shape.transform(back_to_center) # Move back to center


# Example: Move the cell by `random_number_x` and `random_number_y`
move_cell(cell1, random_number_x, random_number_y)

# Rotate cell1 (Ring) within the top cell in place
rotate_cell_around_center(cell1, 90)

# Save the layout to see the transformation
layout.write("C:\\Users\\limyu\\Downloads\\amf_grating_replaced_transformed.gds")
print(f"Moved cell1 by {random_number_x/1000} µm in X and {random_number_y/1000} µm in Y.")
#---------------------- end of moves cell----------------------



#========================================= reward evaluation stage===============


#----------------------loading moved klayout file----------------------



# Function to get shapes from a cell
def get_shapes_from_cell(cell, layer_index):
    shapes = cell.shapes(layer_index)
    return shapes


# Update evaluation logic to be consistent with the movement part
def evaluate_cells(reward, layout, cell1, cell2, layer_index, layer_index2, min_distance=40000):

    
    # Function to get shapes from a cell
    def get_shapes_from_cell(cell, layer_index):
        shapes = cell.shapes(layer_index)
        return shapes

    # Function to convert shape to polygon and retrieve its points
    def get_shape_coordinates(shape):
        coordinates = []
        if shape.is_path():
            path = shape.path
            for point in path.each_point():
                coordinates.append((point.x, point.y))
        elif shape.is_polygon():
            polygon = shape.polygon
            for point in polygon.each_point_hull():  # Get polygon points
                coordinates.append((point.x, point.y))
        elif shape.is_simple_polygon():
            simple_polygon = shape.simple_polygon
            for point in simple_polygon.each_point_hull():
                coordinates.append((point.x, point.y))
        elif shape.is_box():
            box = shape.box
            coordinates.append((box.p1.x, box.p1.y))
            coordinates.append((box.p2.x, box.p2.y))
        return coordinates

    # Function to compute bounding box
    def get_bounding_box(coords):
        if coords:
            x_values = [x for x, y in coords]
            y_values = [y for x, y in coords]
            min_x = min(x_values)
            max_x = max(x_values)
            min_y = min(y_values)
            max_y = max(y_values)
            return min_x, max_x, min_y, max_y
        return None, None, None, None

    def ranges_overlap(x1, x2, x3, x4):
        return not (x2 < x3 or x4 < x1)

    # Function to get the bounding box of a shape
    def get_bounding_box1(shape):
        if shape.is_path():
            return shape.path.bbox()  # Returns the bounding box of the path
        elif shape.is_polygon() or shape.is_simple_polygon():
            return shape.polygon.bbox()  # Returns the bounding box of the polygon
        elif shape.is_box():
            return shape.box  # Returns the bounding box of the box
        return None

    # Function to calculate the distance between two bounding boxes
    def calculate_distance(bbox1, bbox2):
        if bbox1.overlaps(bbox2):
            return 0  # If they overlap, the distance is zero
        dx = max(bbox2.p1.x - bbox1.p2.x, bbox1.p1.x - bbox2.p2.x, 0)
        dy = max(bbox2.p1.y - bbox1.p2.y, bbox1.p1.y - bbox2.p2.y, 0)
        return (dx**2 + dy**2) ** 0.5  # Return Euclidean distance

    # Get shapes from both cells
    shapes1 = get_shapes_from_cell(cell1, layer_index)
    shapes2 = get_shapes_from_cell(cell2, layer_index2)

    # Check for overlapping shapes
    for shape1 in shapes1.each():
        for shape2 in shapes2.each():
            coords1 = get_shape_coordinates(shape1)
            min_x1, max_x1, min_y1, max_y1 = get_bounding_box(coords1)

            coords2 = get_shape_coordinates(shape2)
            min_x2, max_x2, min_y2, max_y2 = get_bounding_box(coords2)

            if ranges_overlap(min_x1, max_x1, min_x2, max_x2) and ranges_overlap(min_y1, max_y1, min_y2, max_y2):
                # Deduct reward if shapes overlap
                reward -= 100

    # Check if shapes are too close
    distance_list = []
    for shape1 in shapes1.each():
        bbox1 = get_bounding_box1(shape1)
        if bbox1 is None:
            continue

        for shape2 in shapes2.each():
            bbox2 = get_bounding_box1(shape2)
            if bbox2 is None:
                continue

            distance = calculate_distance(bbox1, bbox2)
            distance_list.append(distance)
    

    for distance in distance_list:
        if distance < min_distance:  # Deduct reward if too close
            reward -= 100
    #ignore this, calculate reward when the cells are packed 
    # Check if cells are packed
    """
    def get_bbox_coords(cell):
        bbox = cell.bbox()
        dbu = layout.dbu
        xmin = bbox.p1.x * dbu
        ymin = bbox.p1.y * dbu
        xmax = bbox.p2.x * dbu
        ymax = bbox.p2.y * dbu
        return xmin, xmax, ymin, ymax

    xmin1, xmax1, ymin1, ymax1 = get_bbox_coords(cell1)
    xmin2, xmax2, ymin2, ymax2 = get_bbox_coords(cell2)

    x_width = abs(max(xmax1, xmax2) - min(xmin1, xmin2))
    y_height = abs(max(ymax1, ymax2) - min(ymin1, ymin2))

    # Deduct reward if cells are too far apart (not packed)
    reward -= x_width
    reward -= y_height
    """

    return reward



# Function to get all unique layer numbers (without datatypes) that contain shapes in a specific cell
def get_layer_numbers_from_cell(cell):
    layer_numbers = set()
    layout = cell.layout()
    # Iterate over all layer indices in the layout
    for layer_index in range(layout.layers()):
        shapes = cell.shapes(layer_index)
        if not shapes.is_empty():
            # Get the layer information (number only) for the layer index
            layer_info = layout.get_info(layer_index)
            layer_numbers.add(layer_info.layer)
    return list(layer_numbers)

def evaluate_all_layers(cell1, cell2, reward):
    # Assuming layout and cell2 are already defined as in your example
    layers1 = get_layer_numbers_from_cell(cell1)
    
    # Assuming layout and cell2 are already defined as in your example
    layers2 = get_layer_numbers_from_cell(cell2)
    
    
    for l1 in layers1:
        for l2 in layers2:
    
            # Define the layer index (use the appropriate layer number and datatype)
            layer_index = layout.layer(l1, 0)  # Layer 54, datatype 0
            layer_index2 = layout.layer(l2, 0) 
            

            
            # Call the updated evaluation function
            reward = evaluate_cells(reward, layout, cell1, cell2, layer_index, layer_index2)
    return reward

# Load the layout
layout = pya.Layout()
layout.read("C:\\Users\\limyu\\Downloads\\amf_grating_replaced_transformed.gds")

all_cells = []
for cell in layout.each_cell():
    cell_name = cell.name
    all_cells.append(cell_name)


OPA_index = all_cells.index('OPA2')
Ring_index = all_cells.index('Ring')
top_cell = layout.cell("TOP")

# Get the cells (assuming you know which ones have the shapes)
cell1 = layout.cell(Ring_index)  # First cell
cell2 = layout.cell(OPA_index)  # Second cell

reward  = evaluate_all_layers(cell1, cell2, reward = 0)
print('=======================rewards======================')
print('reward', reward)
