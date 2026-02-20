from pathlib import Path
from math import pi

import numpy as np
from scipy.constants import speed_of_light as c
from scipy.fft import rfft, rfftfreq

from numpy.typing import NDArray


class PostProcessorRound:


    def __init__(self, working_directory: str | Path) -> None:

        self.work_dir = Path(working_directory) / 'round'
        
        self.wake_files = sorted(self.work_dir.glob('wakeL_*.txt'))
        if not self.wake_files:
            raise FileNotFoundError('No `wakeL_*.txt` files found. Paths correct and simulation complete?')
        
        print('Found wake files for round geometry for the following modes: ')
        for wake_file in self.wake_files:
            print(f'  {wake_file.name[6:8]}')

        # Load wake parameters
        wake_data = np.loadtxt(self.wake_files[0], comments='%')
        self.y_mesh_step = float(wake_data[0,0]) # m
        self.bunch_offset_mesh_steps = int(wake_data[0,1]) # dimensionless
        self.bunch_offset = (self.bunch_offset_mesh_steps + 0.5) * self.y_mesh_step # m
        self.bunch_sigma = float(wake_data[1,1]) # m
        self.s = wake_data[2:,0] # m
        self.ds = self.s[1] - self.s[0] # m

        # Load bunch parameters
        bunch_data = np.loadtxt(self.work_dir/'Iz0.txt', comments='%')
        s_bunch = bunch_data[:, 0] # m
        i_bunch = bunch_data[:, self.bunch_offset_mesh_steps+2] # C/m
        self.b = np.interp(self.s, s_bunch, i_bunch) # C/m, interpolate bunch on wake axis
        self.bunch_charge = np.trapezoid(self.b, self.s) # C

        # shift so bunch center is at s = 0
        self.s -= 5*self.bunch_sigma - self.ds/2

        print(f'Wake data loaded:')
        print(f'  Bunch length  : {self.bunch_sigma*1e3} mm (1 sigma)')
        print(f'                  {self.bunch_sigma/c*1e9:.3f} ns (1 sigma)')
        print(f'  Bunch charge  : {self.bunch_charge*1e9} nC')
        print(f'  Wake length   : {self.s[-1]*1e3:.1f} mm')
        print(f'                  {self.s[-1]/c*1e9:.1f} ns')
        print(f'  Mesh step Z   : {self.ds*1e3:.1f} mm')
        print(f'                  {self.bunch_sigma/self.ds:.1f} per sigma')


    def get_bunch(self) -> tuple[NDArray, NDArray]:
        return self.s, self.b # m, C
    

    def get_bunch_spectrum(
            self,
            oversampling: float = 1,
            cutoff_by_bunch_sigma: float | None = 3,
            normalize: bool = False,
        ) -> tuple[NDArray, NDArray]:

        n_samples = int(len(self.s) * oversampling + 0.5)

        f = rfftfreq(n_samples, d=self.ds/c) # Hz
        B = np.asarray(rfft(self.b, n=n_samples)) * self.ds # C
        B *= np.exp(-2j*pi*f * self.s[0]/c) # compensate phase offset from time shift

        if normalize:
            B /= np.max(np.abs(B))
        else:
            pass

        if cutoff_by_bunch_sigma is None:
            sample_cutoff = len(f)
        else:
            df = f[1] - f[0]
            sample_cutoff = int(cutoff_by_bunch_sigma  / (2*pi * self.bunch_sigma/c) / df)
        
        return f[:sample_cutoff], B[:sample_cutoff] # Hz, C
    

    def get_longitudinal_wake(self, mode: int) -> tuple[NDArray, NDArray]:
        
        wake_data = np.loadtxt(self.work_dir / f'wakeL_{mode:02d}.txt', comments='%')
        w = wake_data[2:,1]*1e9 # V/C/m^p, where p is mode index
        w /= self.bunch_offset ** (2*mode)

        return self.s, wake_data[2:,1]*1e9 # m, V/C/m^(2p)
    

    def get_transverse_wake(self, mode: int) -> tuple[NDArray, NDArray]:
        raise NotImplementedError
        

    def get_longitudinal_impedance(
            self,
            mode: int,
            oversampling: float = 1,
            deconvolution: bool = True,
            cutoff_by_bunch_sigma: float | None = 3,
        ) -> tuple[NDArray, NDArray]:

        s, w = self.get_longitudinal_wake(mode) # m, V/C/m^(2p)

        n_samples = int(len(self.s) * oversampling + 0.5)

        f = rfftfreq(n_samples, d=self.ds/c) # Hz
        Z = -np.asarray(rfft(w, n=n_samples)) * self.ds / c # Ohm/m^p
        Z *= np.exp(-2j*pi*f * self.s[0]/c) # compensate phase offset from time shift

        # divide by bunch spectrum to deconvolve
        if deconvolution:
            Z /= self.get_bunch_spectrum(
                normalize=True, cutoff_by_bunch_sigma=None, oversampling=oversampling
            )[1]

        # cutoff based on bunch spectrum
        if cutoff_by_bunch_sigma is None:
            sample_cutoff = len(f)
        else:
            df = f[1] - f[0]
            sample_cutoff = int(cutoff_by_bunch_sigma  / (2*pi * self.bunch_sigma/c) / df)

        return f[:sample_cutoff], Z[:sample_cutoff], # Hz, Ohm/m^p


    def get_transverse_impedance(
            self,
            mode: int,
            oversampling: float = 1,
            deconvolution: bool = True,
            cutoff_by_bunch_sigma: float | None = 3,
        ) -> tuple[NDArray, NDArray]:

        raise NotImplementedError
