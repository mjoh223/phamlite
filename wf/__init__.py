from pathlib import Path
from typing import Annotated, List, Optional, Tuple
import sys
sys.path.insert(0, '/root/wf/')
import phamlite
from latch import small_task, workflow
from latch.types import LatchDir, LatchFile

@small_task
def launch_phamlite_app_tsk(genbank_files: List[LatchFile], run_name: str) -> (LatchFile):
    local_paths = [x.local_path for x in genbank_files]
    output_file = phamlite.phamlite_runner(local_paths, run_name)
    #import and run the phamlite scripts
    return LatchFile(output_file, "latch:///{}.html".format(run_name))

@workflow
def phamlite_wf(genbank_files: List[LatchFile], run_name: str) -> (LatchFile):
    """
    __metadata__:
        display_name: Phamlite
        author: Matthew Johnson
            name: Phamlite
            email: matthew.johnson@ucsf.edu
            github: https://github.com/mjoh223
        repository: https://github.com/mjoh223/phamlite
        license:
            id: MIT
    Args:
        genbank_files:
          input gb file.
          __metadata__:
            display_name: genbank file
        run_name:
          run name
          __metadata__:
            display_name: output filename
    """
    return launch_phamlite_app_tsk(genbank_files = genbank_files, run_name = run_name)

# if __name__ == "__main__":
#     phamlite_wf(genbank_file=LatchFile("/Users/matt/Desktop/PhiKZ.gb"))
