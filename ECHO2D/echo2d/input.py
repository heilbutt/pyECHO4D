from math import log, pi, sqrt
from pathlib import Path
import textwrap

from scipy.constants import speed_of_light as c

from typing import Literal

from ...preprocess.geometry.basic_elements import LengthUnit


def write_input_file(
    directory: str | Path,
    geometry_file_name: str,
    geometry_file_unit: LengthUnit,
    geometry_type: Literal['recta', 'round'],
    modes: list[int],
    wake_length: float,
    bunch_sigma: float,
    bunch_offset: float | None = None,
    bunch_offset_mesh_steps: int | None = None,
    z_mesh_steps_per_sigma: float = 5,
    y_mesh_steps_per_sigma: float | None = None,
    recta_geometry_width: float | None = None,
    recta_symmetry_condition: Literal['elec', 'magn'] | None = None,
    geometry_is_convex: bool = False,
    conductive_material_mesh_steps: int = 0
) -> None:
    
    # Consistency check for round/recta and geometry width

    if geometry_type == 'recta' and recta_geometry_width is None:
        raise ValueError('Must specify `recta_geometry_width` in meters for rectangular geometry')
    
    if geometry_type == 'round':
        if recta_geometry_width is None:
            recta_geometry_width = 0
        else:
            print('WARNING: Geometry is round but `recta_geometry_width` has been set, it will be ignored')

    # Consistency check for round/recta and symmetry condition

    if geometry_type == 'recta' and recta_symmetry_condition not in ['magn', 'elec']:
        raise ValueError('Must specify `recta_symmetry_condition` as `magn` or `elec` for rectangular geometry')
    
    if geometry_type == 'round':
        if recta_symmetry_condition is None:
            recta_symmetry_condition = 'magn'
        else:
            print('WARNING: Geometry is round but `recta_symmetry_condition` has been set, it will be ignored')

    # mesh steps from bunch length

    if y_mesh_steps_per_sigma is None:
        y_mesh_steps_per_sigma = z_mesh_steps_per_sigma

    z_mesh_step: float = bunch_sigma / z_mesh_steps_per_sigma
    y_mesh_step: float = bunch_sigma / y_mesh_steps_per_sigma

    # bunch offset

    if (bunch_offset is not None) and (bunch_offset_mesh_steps is None):
        bunch_offset_mesh_steps = int(bunch_offset / y_mesh_step)
    elif (bunch_offset is None) and (bunch_offset_mesh_steps is not None):
        pass
    else:
        raise ValueError('Specify one of: `bunch_offset` or `bunch_offset_mesh_steps`')
    
    if (bunch_offset_mesh_steps < 3) and any(mode > 0 for mode in modes):
        print('WARNING: Non-monopole modes have been requested, but input values would result in bunch offset of less than 3 mesh steps')

    # Calculate some figures of merit

    bunch_max_freq = sqrt(-(c**2)/(2*pi**2*bunch_sigma**2) * (-20/20) * log(10))
    freq_resolution = c/wake_length
    max_resolvable_Q = bunch_max_freq / freq_resolution
    
    print(f'Creating Echo2D input file. Relevant parameters:')
    print(f'  Bunch length     : {bunch_sigma*1000:.1f} mm (1 sigma)')
    print(f'                     {bunch_sigma/c*1e9:.3f} ns (1 sigma)')
    print(f'  Max freq (-20 dB): {bunch_max_freq/1e9:.1f} GHz')
    print(f'  Wake length      : {int(wake_length/z_mesh_step):d} mesh steps')
    print(f'                     {wake_length*1000:.1f} mm')
    print(f'                     {wake_length/c*1e9:.1f} ns')
    print(f'  Freq. resolution : {c/wake_length/1e6:.1f} MHz')
    print(f'  Max. resolvable Q: {max_resolvable_Q:.0f} (approx)')
    print(f'  Bunch offset     : {bunch_offset_mesh_steps:d} mesh steps')
    print(f'                     {bunch_offset_mesh_steps*y_mesh_step*1000:.1f} mm')
    print(f'  Mesh step Z      : {z_mesh_step*1000:.1f} mm')
    print(f'                     {z_mesh_steps_per_sigma:.1f} per sigma')
    print(f'  Mesh step Y      : {y_mesh_step*1000:.1f} mm')
    print(f'                     {y_mesh_steps_per_sigma:.1f} per sigma')

    # Write input file

    text = textwrap.dedent(f'''
        %%%%%%%%%%%%%% geometry %%%%%%%%%%%%%%%%%%%

        GeometryFile={geometry_file_name}
        Units={geometry_file_unit} % m/cm/mm, unit used by geometry file
        GeometryType={geometry_type} % recta / round
        Width={recta_geometry_width:f} % in meters, irrelevant for round geometry
        SymmetryCondition={recta_symmetry_condition} % magn/elec, irrelevant for round geometry
        Convex={int(geometry_is_convex):d} % 0/1, wether geometry is convex, faster computation in that case

        %%%%%%%%%%%%%% beam %%%%%%%%%%%%%%%%%%%%%%%

        InPartFile=- % - for Gaussian beam with BunchSigma
        BunchSigma={bunch_sigma:f} % rms bunch length in m
        Offset={bunch_offset_mesh_steps:d} % y0 in transverse mesh steps
        InjectionTimeStep=0

        %%%%%%%%%%%%%% field %%%%%%%%%%%%%%%%%%%%%%

        InFieldDir=-
        PortDir=-
        PortPosition=-1

        %%%%%%%%%%%%%% model %%%%%%%%%%%%%%%%%%%%%%

        WakeIntMethod=ind
        Modes={' '.join([str(mode) for mode in modes])} % modes to consider, list of ints separated by spaces
        ParticleMotion=0
        ParticleField=1
        CurrentFilter=0
        ParticleLoss=0

        %%%%%%%%%%%%%% mesh %%%%%%%%%%%%%%%%%%%%%%%

        MeshLength={int(wake_length/z_mesh_step):d} % length of moving mesh/calculation window in mesh steps
        StartPosition=0
        TimeSteps=-1
        StepY={y_mesh_step:f} % mesh step size in m in transverse direction
        StepZ={z_mesh_step:f} % mesh step size in m in longitudinal direction
        NStepsInConductive={conductive_material_mesh_steps:d}
        AdjustMesh=0
        MeshMotionFile=-

        %%%%%%%%%%%%%% monitors %%%%%%%%%%%%%%%%%%%%%%%

        DumpField=0
        DumpParticles=0
        DumpCurrent=0
        DumpMesh=0
        ''')
    
    file_path = Path(directory) / 'input_in.txt'
    file_path.write_text(text)
    print(f'Input written to "{file_path.resolve()}"')