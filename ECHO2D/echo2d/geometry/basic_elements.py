from __future__ import annotations

from pathlib import Path

from typing import Literal, TypeAlias


ELLIPSE_DIRECTION: dict[str, int] = {
    'clockwise' : 0,
    'counterclockwise': 1
}
EllipseDirection: TypeAlias = Literal['clockwise', 'counterclockwise']


LENGTH_UNIT: dict[str, float] = {
    'mm': 1e-3,
    'cm': 1e-2,
    'm' : 1
}
LengthUnit: TypeAlias = Literal['mm', 'cm', 'm']


class Coordinate:


    def __init__(self, z: float, y: float) -> None:
        self.z = z
        self.y = y


    def shift_y(self, delta_y: float):
        return Coordinate(self.z, self.y + delta_y)
    

    def shift_z(self, delta_z: float):
        return Coordinate(self.z + delta_z, self.y)
    

    def __add__(self, other: Coordinate) -> Coordinate:
        return Coordinate(self.z + other.z, self.y + other.y)
    

    def __mul__(self, scalar: float) -> Coordinate:
        return Coordinate(scalar * self.z, scalar * self.y)
    

    def __truediv__(self, scalar: float) -> Coordinate:
        return self * (1 / scalar)


class _Segment:


    def __init__(
        self,
        start: Coordinate,
        end: Coordinate,
        top_left: Coordinate | None = None,
        bottom_right: Coordinate | None = None,
        direction: EllipseDirection | None = None,
        wall_conductivity: float | None = None
    ) -> None:  
        
        self.start = start
        self.end = end
        self.top_left = top_left
        self.bottom_right = bottom_right
        self.direction = direction
        self.wall_conductivity = wall_conductivity


    def min_y(self) -> float:
        raise NotImplementedError


    def max_y(self) -> float:
        raise NotImplementedError


    def generate_text_line(self, output_unit: Literal['mm', 'cm', 'm']) -> str:

        parameters: list[float] = [
            self.start.z / LENGTH_UNIT[output_unit],
            self.start.y / LENGTH_UNIT[output_unit],
            self.end.z / LENGTH_UNIT[output_unit],
            self.end.y / LENGTH_UNIT[output_unit],
            self.top_left.z / LENGTH_UNIT[output_unit] if self.top_left is not None else 0,
            self.top_left.y / LENGTH_UNIT[output_unit] if self.top_left is not None else 0,
            self.bottom_right.z / LENGTH_UNIT[output_unit] if self.bottom_right is not None else 0,
            self.bottom_right.y / LENGTH_UNIT[output_unit] if self.bottom_right is not None else 0,
            ELLIPSE_DIRECTION[self.direction] if self.direction is not None else 0,
            self.wall_conductivity if self.wall_conductivity is not None else 0
        ]

        return ' '.join(str(p) for p in parameters)


class LineSegment(_Segment):


    def __init__(
        self, 
        start: Coordinate, end: Coordinate,
        wall_conductivity: float | None = None    
    ) -> None:   
             
        super().__init__(
            start, end,
            None, None,
            None,
            wall_conductivity
        )


    def min_y(self) -> float:
        return min(self.start.y, self.end.y)


    def max_y(self) -> float:
        return max(self.start.y, self.end.y)


class EllipseSegment(_Segment):


    def __init__(
        self, 
        start: Coordinate, end: Coordinate,
        top_left: Coordinate, bottom_right: Coordinate,
        direction: EllipseDirection,
        wall_conductivity: float | None = None    
    ) -> None:  
              
        super().__init__(
            start, end,
            top_left, bottom_right,
            direction,
            wall_conductivity
        )

        self.ellipse_center = (self.start + self.end) / 2


class _Material:


    def __init__(
        self,
        relative_permeability: float,
        relative_permittivity: float,
        conductivity: float,
        segments: list[_Segment] = [],
    ) -> None:
        
        self.segments = segments
        self.relative_permeability = relative_permeability
        self.relative_permittivity = relative_permittivity
        self.conductivity = conductivity


    def _generate_title_comment(self) -> str:
        raise NotImplementedError
    

    def _process_conductivity(self) -> None:
        raise NotImplementedError
    

    def _check_segment_connectivity(self, rel_tol: float = 1e-6) -> None:

        for this_seg, next_seg in zip(self.segments, self.segments[1:]):
            
            y0 = this_seg.end.y
            y1 = next_seg.start.y
            if abs(y1/y0 - 1) > rel_tol:
                raise ValueError(f'Disconnected segments (Y coordinate)! y end = {y0:f} vs. y start = {y1:f}')
            
            z0 = this_seg.end.z
            z1 = next_seg.start.z
            if abs(z1/z0 - 1) > rel_tol:
                raise ValueError(f'Disconnected segments (Z coordinate)! z end = {z0:f} vs. z start = {z1:f}')
        

    def start_z(self) -> float:
        try:
            return self.segments[0].start.z
        except IndexError: # no elements defined
            return 0
    

    def end_z(self) -> float:
        try:
            return self.segments[-1].end.z
        except IndexError:
            return 0    


    def generate_text_lines(self, output_unit: LengthUnit) -> list[str]:

        if len(self.segments) == 0:
            raise ValueError('Material needs at least one segment')
        
        self._process_conductivity()
        self._check_segment_connectivity()
        
        lines: list[str] = []
        
        lines.append(self._generate_title_comment())
        lines.append(f'{len(self.segments):d} {self.relative_permeability:f} {self.relative_permittivity:f} {self.conductivity:f}')
        lines.append('% Segments: lines and ellipses')
        lines.extend(seg.generate_text_line(output_unit) for seg in self.segments)
        
        return lines


class VolumeMaterial(_Material):


    def __init__(
        self,
        relative_permeability: float,
        relative_permittivity: float,
        conductivity: float,
        segments: list[_Segment] = [],
    ) -> None:
        
        super().__init__(
            relative_permeability,
            relative_permittivity,
            conductivity,
            segments
        )


    def _process_conductivity(self) -> None:
        if any(seg.wall_conductivity is not None for seg in self.segments):
            raise ValueError('For defining volume material, wall conductivity of individual segments must be None')
        

    def _generate_title_comment(self) -> str:
        return '% Definition of volume material (dielectric or lossy): number of segments, rel. permeability, rel. permittivity, conductivity'
        

class ConductiveWall(_Material):
    

    def __init__(
        self,
        wall_conductivity: float | None = None,
        segments: list[_Segment] = [],
    ) -> None:
        
        super().__init__(1, 1, 0, segments)
        self.wall_conductivity = wall_conductivity


    def _process_conductivity(self) -> None:
        if self.wall_conductivity is None:
            if any(seg.wall_conductivity is None for seg in self.segments):
                raise ValueError('For defining the conductive wall, wall conductivity (float > 0) must be provided either for the wall or for each segment individually')
        else:
            for seg in self.segments:
                seg.wall_conductivity = self.wall_conductivity
        

    def _generate_title_comment(self) -> str:
        return '% Definition of conductive wall: number of segments, rel. permeability, rel. permittivity, conductivity'
    

    def is_convex(self) -> bool:
        for this_seg, next_seg in zip(self.segments, self.segments[1:]):
            if this_seg.end.z > next_seg.start.z:
                return False
        return True
    

class Geometry:


    def __init__(
        self,
        materials: list[_Material] = []
    ) -> None:
        
        self.materials = materials


    def _check_material_definition(self) -> None:

        if len(self.materials) == 0:
            raise ValueError('Geometry needs at least one material')
        if not isinstance(self.materials[0], ConductiveWall):
            raise ValueError('First material must be a conductive wall')


    def generate_text_lines(self, output_length_unit: LengthUnit) -> list[str]:

        self._check_material_definition()
        
        lines: list[str] = []
        lines.append('% Number of materials')
        lines.append(f'{len(self.materials):d}')
        for material in self.materials:
            lines.extend(material.generate_text_lines(output_length_unit))

        return lines
    
    
    def write_geometry_file(self, file_path: str | Path, output_unit: LengthUnit) -> None:
        
        lines = self.generate_text_lines(output_unit)
        
        file_path = Path(file_path)
        file_path.write_text('\n'.join(lines))

        print(f'Geometry file written to "{file_path.resolve()}"')


    def write_geometry_file_echoz2(self, file_path: str | Path) -> None:
        self.write_geometry_file(file_path, 'cm') # EchoZ2 only accepts cm