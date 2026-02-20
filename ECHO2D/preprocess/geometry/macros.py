from .basic_elements import _Segment, Coordinate, LineSegment, EllipseSegment

def generate_beam_pipe(
        start_z: float,
        radius: float,
        length: float
    ) -> list[_Segment]:

    start = Coordinate(start_z, radius)
    return [LineSegment(start, start.shift_z(length))]


def generate_bellow_convolution(
    start_z: float,
    beam_pipe_radius: float,
    depth: float,
    gap: float,
    rounding_radius: float
) -> list[_Segment]:
    
    if depth < 2*rounding_radius:
        ValueError('Depth must not be smaller than 2 rounding radii')

    if gap < 2*rounding_radius:
        ValueError('Gap must not be smaller than 2 rounding radii')

    segments: list[_Segment] = []
    start = Coordinate(start_z, beam_pipe_radius)

    segments.append(
        EllipseSegment(
            start,
            start.shift_z(rounding_radius).shift_y(rounding_radius),
            start.shift_z(-rounding_radius).shift_y(2*rounding_radius),
            start.shift_z(rounding_radius),
            'counterclockwise'
        )
    )

    if depth > 2*rounding_radius:
        segments.append(
            LineSegment(
                segments[-1].end,
                segments[-1].end.shift_y(depth-2*rounding_radius)
            )
        )
    
    segments.append(
        EllipseSegment(
            segments[-1].end,
            segments[-1].end.shift_z(rounding_radius).shift_y(rounding_radius),
            segments[-1].end.shift_y(rounding_radius),
            segments[-1].end.shift_z(2*rounding_radius).shift_y(-rounding_radius),
            'clockwise'
        )
    )

    if gap > 2*rounding_radius:
        segments.append(
            LineSegment(
                segments[-1].end,
                segments[-1].end.shift_z(gap - 2*rounding_radius)
            )
        )

    segments.append(
        EllipseSegment(
            segments[-1].end,
            segments[-1].end.shift_z(rounding_radius).shift_y(-rounding_radius),
            segments[-1].end.shift_z(-rounding_radius),
            segments[-1].end.shift_z(rounding_radius).shift_y(-2*rounding_radius),
            'clockwise'
        )
    )

    if depth > 2*rounding_radius:
        segments.append(
            LineSegment(
                segments[-1].end,
                segments[-1].end.shift_y(-(depth-2*rounding_radius))
            )
        )

    segments.append(
        EllipseSegment(
            segments[-1].end,
            segments[-1].end.shift_z(rounding_radius).shift_y(-rounding_radius),
            segments[-1].end.shift_y(rounding_radius),
            segments[-1].end.shift_z(2*rounding_radius).shift_y(-rounding_radius),
            'counterclockwise'
        )
    )

    return segments


def generate_bellow(
    start_z: float,
    beam_pipe_radius: float,
    depth: float,
    gap: float,
    rounding_radius: float,
    period: float,
    number_of_convolutions: int
) -> list[_Segment]:
    
    if period < gap + 2*rounding_radius:
        raise ValueError('Bellow convolution period must not be smaller than (gap + 2 rounding radii)')

    if number_of_convolutions < 1:
        raise ValueError('Bellow must have at least one convolution')

    segments: list[_Segment] = []
    for n in range(number_of_convolutions):
        
        segments.extend(
            generate_bellow_convolution(
                start_z + n * period,
                beam_pipe_radius,
                depth,
                gap,
                rounding_radius
            )
        )

        if (period > gap + 2*rounding_radius) and (n < number_of_convolutions - 1):
            segments.append(
                LineSegment(
                    segments[-1].end,
                    segments[-1].end.shift_z(period - (gap + 2*rounding_radius))
                )
            )

    return segments