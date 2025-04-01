import Rhino
import math
import System
from Rhino.Geometry import Rectangle3d, Plane, Point3d, Curve, Transform
from Rhino.Geometry.Intersect import Intersection

# Update component info
ghenv.Component.Name = "Grid Pattern Generator and Trimmer"
ghenv.Component.NickName = "GridPattern"
ghenv.Component.Message = ""

# Initialize debug messages list globally
debug_messages = []

def debug_print(msg):
    # Only add to out parameter, not to component message
    print(msg)
    debug_messages.append(msg)

class GridCell:
    def __init__(self):
        self.curve = None
        self.center = None
        self.corners = []
    
    def get_test_points(self):
        return self.corners + [self.center] if self.center else self.corners
    
    def transform(self, xform):
        if self.curve:
            self.curve.Transform(xform)
        if self.center:
            self.center.Transform(xform)
        for i in range(len(self.corners)):
            self.corners[i].Transform(xform)

class GridGenerator:
    def __init__(self, boundary, cell_width, cell_height, grid_type="rectangular", plane=Plane.WorldXY):
        self.boundary = boundary
        self.cell_width = cell_width
        self.cell_height = cell_height
        self.grid_type = grid_type
        self.plane = plane
        
        # Create transforms between world and local coordinates
        self.to_local = Transform.PlaneToPlane(self.plane, Plane.WorldXY)
        self.to_world = Transform.PlaneToPlane(Plane.WorldXY, self.plane)
        
        # Transform boundary to local coordinates for calculations
        self.local_boundary = boundary.DuplicateCurve()
        self.local_boundary.Transform(self.to_local)
        
        # Get bounding box in local coordinates
        self.bbox = self.local_boundary.GetBoundingBox(True)
        debug_print(f"Local boundary box: Min({self.bbox.Min}), Max({self.bbox.Max})")
    
    def create_rectangular_cell(self, x, y, offset=0):
        cell = GridCell()
        # Apply offset to x-coordinate
        adjusted_x = x + offset
        
        cell.curve = Rectangle3d(
            Plane.WorldXY, 
            Point3d(adjusted_x, y, 0), 
            Point3d(adjusted_x + self.cell_width, y + self.cell_height, 0)
        ).ToNurbsCurve()
        cell.center = Point3d(adjusted_x + self.cell_width/2, y + self.cell_height/2, 0)
        cell.corners = [
            Point3d(adjusted_x, y, 0),
            Point3d(adjusted_x + self.cell_width, y, 0),
            Point3d(adjusted_x + self.cell_width, y + self.cell_height, 0),
            Point3d(adjusted_x, y + self.cell_height, 0)
        ]
        return cell
    
    def create_hexagonal_cell(self, center_x, center_y):
        cell = GridCell()
        radius = self.cell_width / 2
        corners = []
        for i in range(6):
            angle = i * math.pi / 3
            x = center_x + radius * math.cos(angle)
            y = center_y + radius * math.sin(angle)
            corners.append(Point3d(x, y, 0))
        
        cell.corners = corners
        cell.center = Point3d(center_x, center_y, 0)
        cell.curve = Curve.CreateInterpolatedCurve(corners + [corners[0]], 1)
        return cell
    
    def create_triangular_cell(self, base_x, base_y, inverted=False):
        cell = GridCell()
        height = self.cell_width * math.sqrt(3)/2
        
        if not inverted:
            corners = [
                Point3d(base_x, base_y, 0),
                Point3d(base_x + self.cell_width, base_y, 0),
                Point3d(base_x + self.cell_width/2, base_y + height, 0)
            ]
        else:
            corners = [
                Point3d(base_x + self.cell_width/2, base_y, 0),
                Point3d(base_x + self.cell_width, base_y + height, 0),
                Point3d(base_x, base_y + height, 0)
            ]
        
        cell.corners = corners
        cell.center = Point3d(sum(p.X for p in corners)/3, sum(p.Y for p in corners)/3, 0)
        cell.curve = Curve.CreateInterpolatedCurve(corners + [corners[0]], 1)
        return cell
    
    def generate_cells(self):
        debug_print(f"Generating {self.grid_type} grid with cell size: {self.cell_width} x {self.cell_height}")
        
        cells = []
        if self.grid_type in ["rectangular", "offset_rectangular"]:
            x = math.floor(self.bbox.Min.X / self.cell_width) * self.cell_width
            while x < self.bbox.Max.X:
                y = math.floor(self.bbox.Min.Y / self.cell_height) * self.cell_height
                row_index = int((y - self.bbox.Min.Y) / self.cell_height)
                
                while y < self.bbox.Max.Y:
                    # For offset_rectangular, apply 0.5 * cell_width offset on alternate rows
                    offset = (0.5 * self.cell_width) if (self.grid_type == "offset_rectangular" and row_index % 2 == 1) else 0
                    cell = self.create_rectangular_cell(x, y, offset)
                    cells.append(cell)
                    y += self.cell_height
                    row_index += 1
                x += self.cell_width
                
        elif self.grid_type == "hexagonal":
            h_spacing = self.cell_width * 3/4
            v_spacing = self.cell_width * math.sqrt(3)/2
            x = math.floor(self.bbox.Min.X / h_spacing) * h_spacing
            row_count = 0
            while x < self.bbox.Max.X:
                y = math.floor(self.bbox.Min.Y / v_spacing) * v_spacing
                offset = (v_spacing/2) if row_count % 2 else 0
                while y < self.bbox.Max.Y:
                    cell = self.create_hexagonal_cell(x, y + offset)
                    cells.append(cell)
                    y += v_spacing
                row_count += 1
                x += h_spacing
                
        elif self.grid_type == "triangular":
            height = self.cell_width * math.sqrt(3)/2
            x = math.floor(self.bbox.Min.X / self.cell_width) * self.cell_width
            while x < self.bbox.Max.X:
                y = math.floor(self.bbox.Min.Y / height) * height
                while y < self.bbox.Max.Y:
                    cells.append(self.create_triangular_cell(x, y, False))
                    cells.append(self.create_triangular_cell(x, y, True))
                    y += height
                x += self.cell_width
                
        debug_print(f"Generated {len(cells)} cells")
        return cells

def point_containment_test(point, curve, plane, tolerance=0.001):
    containment = curve.Contains(point, plane, tolerance)
    return containment == Rhino.Geometry.PointContainment.Inside

try:
    # Clear previous debug messages
    debug_messages = []
    debug_print("Starting script execution")
    
    # Input validation
    if not boundary:
        raise ValueError("No boundary curve provided")
    
    if isinstance(boundary, list):
        boundaries = [b.ToNurbsCurve() if isinstance(b, Rhino.Geometry.PolylineCurve) else b for b in boundary]
    else:
        boundaries = [boundary.ToNurbsCurve() if isinstance(boundary, Rhino.Geometry.PolylineCurve) else boundary]
    
    for b in boundaries:
        if not b.IsClosed:
            raise ValueError("All boundaries must be closed curves")
    
    grid_type = grid_type if grid_type else "rectangular"
    x_dim = x_dim if x_dim else 1.0
    y_dim = y_dim if y_dim else 1.0
    origin_plane = origin_plane if origin_plane else Plane.WorldXY
    
    debug_print(f"Processing with parameters: grid_type={grid_type}, x_dim={x_dim}, y_dim={y_dim}")
    
    # Initialize outputs with better names
    all_grid_cells = []  # instead of a
    cell_trim_status = []  # instead of b
    complete_cell_count = 0  # instead of c
    trimmed_cell_count = 0  # instead of d
    complete_cells_only = []  # instead of e
    
    for boundary_curve in boundaries:
        debug_print(f"Processing boundary: {boundary_curve}")
        
        # Generate grid
        grid = GridGenerator(boundary_curve, x_dim, y_dim, grid_type, origin_plane)
        cells = grid.generate_cells()
        
        # Process cells
        for cell in cells:
            # Transform cell to world coordinates for intersection testing
            cell.transform(grid.to_world)
            test_points = cell.get_test_points()
            points_inside = [point_containment_test(pt, boundary_curve, origin_plane) for pt in test_points]
            
            intersection_events = Intersection.CurveCurve(cell.curve, boundary_curve, 0.001, 0.001)
            has_intersection = intersection_events is not None and intersection_events.Count > 0
            
            if all(points_inside):
                all_grid_cells.append(cell.curve)
                complete_cells_only.append(cell.curve)
                cell_trim_status.append(0)
                complete_cell_count += 1
            elif any(points_inside) or has_intersection:
                intersect_curves = Rhino.Geometry.Curve.CreateBooleanIntersection(cell.curve, boundary_curve)
                if intersect_curves and len(intersect_curves) > 0:
                    for trimmed_curve in intersect_curves:
                        all_grid_cells.append(trimmed_curve)
                        cell_trim_status.append(1)
                        trimmed_cell_count += 1
    
    debug_print(f"Final counts - Complete cells: {complete_cell_count}, Trimmed cells: {trimmed_cell_count}")
    
    # Update component message to be cleaner
    ghenv.Component.Message = f"{grid_type.capitalize()} Grid\n{complete_cell_count} complete | {trimmed_cell_count} trimmed"
    
    # Set outputs with your preferred names
    allTiles = all_grid_cells
    panelStatus = cell_trim_status
    fullCount = complete_cell_count
    trimCount = trimmed_cell_count
    fullTiles = complete_cells_only
    
    # Assign to component outputs
    a = allTiles
    b = panelStatus
    c = fullCount
    d = trimCount
    e = fullTiles
    
    # Combine all debug messages into the out parameter
    out = "\n".join(debug_messages)

except Exception as e:
    error_msg = f"Error: {str(e)}\nType: {type(e)}"
    print(error_msg)
    ghenv.Component.Message = "Error occurred"
    a, b, c, d, e = None, None, None, None, None
    out = error_msg