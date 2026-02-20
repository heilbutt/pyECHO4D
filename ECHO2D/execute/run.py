from pathlib import Path
import subprocess

def run(executable_path: str | Path, working_directory: str | Path) -> None:
    print(f'Running Echo2D in {Path(working_directory).resolve()}')
    subprocess.run(executable_path, cwd=working_directory)
    print('Echo2D execution complete')