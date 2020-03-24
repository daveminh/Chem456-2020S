from MDAnalysis.core import flags
from MDAnalysis.coordinates import base
from MDAnalysis.lib import util

class INPCRDWriter(base.WriterBase):
    """INPCRD format writer

    Write a file in the AMBER INPCRD format.

    """
    format = 'INPCRD'
    multiframe = False
    units = {'time': 'ps', 'length': 'Angstrom'}

    def __init__(self, filename, n_atoms=None, convert_units=None):
        """Create a new INPCRDWriter

        Parameters
        ----------
        filename: str
            name of output file
        convert_units: bool (optional)
            units are converted to the MDAnalysis base format; ``None`` selects
            the value of :data:`MDAnalysis.core.flags` ['convert_lengths']
        """
        self.filename = filename
        if convert_units is None:
            convert_units = flags['convert_lengths']
        self.convert_units = convert_units  # convert length and time to base units

        self.frames_written = 0

        self.file = util.anyopen(self.filename, 'w')  # open file on init

    def close(self):
        self.file.close()

    def encode_block(self, obj):
        """
        Parameters
        ----------
        obj : AtomGroup or Universe
        """
        traj = obj.universe.trajectory
        ts = traj.ts

        block = f'works with {obj.filename[-80:]:<80}\n'
        block += f'{obj.atoms.n_atoms:5}{ts.time/1000.:15.7e}\n'

        # Positions
        for n in range(ts.positions.shape[0]):
          for dim in range(3):
            block += f'{ts.positions[n,dim]:12.7f}'
          if n%2==1:
            block += '\n'
        if n%2==0:
          block += '\n'

        # Box dimensions
        for c in range(6):
          block += f'{ts.dimensions[c]:12.7f}'
        block += '\n'

        return block

    def write(self, obj):
        """Write a new frame to the MOL2 file.

        Parameters
        ----------
        obj : AtomGroup or Universe
        """
        self.write_next_timestep(obj)

    def write_next_timestep(self, obj):
        """Write a new frame to the MOL2 file.

        Parameters
        ----------
        obj : AtomGroup or Universe
        """
        block = self.encode_block(obj)
        self.file.writelines(block)
