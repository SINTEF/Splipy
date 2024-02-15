from pathlib import Path
from itertools import groupby
from operator import itemgetter
from os.path import exists, isdir, join
from os import makedirs
from typing import Union, Optional, Type, Sequence
from types import TracebackType
from typing_extensions import Self

import numpy as np

from .master import MasterIO
from ..splineobject import SplineObject
from ..splinemodel import SplineModel


class OpenFOAM(MasterIO):
    target: str

    def __init__(self, target: Union[Path, str]) -> None:
        self.target = str(target)

    def __enter__(self) -> Self:
        # Create the target directory if it does not exist
        if not exists(self.target):
            makedirs(self.target)

        # If it does, ensure that it's a directory
        elif not isdir(self.target):
            raise FileExistsError('{} exists and is not a directory'.format(self.target))

        return self

    def __exit__(
        self,
        exc_type: Optional[Type[BaseException]],
        exc_val: Optional[BaseException],
        exc_tb: Optional[TracebackType]
    ) -> None:
        return

    def _header(self, cls: str, obj: str, note: Optional[str] = None) -> str:
        s = 'FoamFile\n{\n'
        s += '    version     2.0;\n'
        s += '    format      ascii;\n'
        s += f'    class       {cls};\n'
        if note:
            s += f'    note        "{note}";\n'
        s += f'    object      {obj};\n'
        s += '}\n'
        return s

    def write(self, model: Union[SplineObject, Sequence[SplineObject], SplineModel]) -> None:
        assert isinstance(model, SplineModel), "OpenFOAM.write only supports SplineModel objects"

        # Only linear volumes in 3D, please
        assert model.pardim == 3
        assert model.dimension == 3
        assert all(patch.obj.order() == (2,2,2) for patch in model.catalogue.top_nodes())

        # Generate all the information we need
        model.generate_cp_numbers()
        model.generate_cell_numbers()
        faces = model.faces()
        ninternal = sum(faces['name'] == None)  # noqa: E711
        note = 'nPoints: {} nCells: {} nFaces: {} nInternalFaces: {}'.format(
            model.ncps, model.ncells, len(faces), ninternal
        )

        # OpenFOAM is very particular about face ordering
        # In order of importance:
        # - Internal faces before boundary faces
        # - All faces in the same boundary must be contiguous
        # - Low number owners before high number owners
        # - Low number neighbors before high number neighbors
        faces_list = list(faces)
        faces_list = sorted(faces_list, key=itemgetter('neighbor'))
        faces_list = sorted(faces_list, key=itemgetter('owner'))
        faces_list = sorted(faces_list, key=lambda x: (x['name'] is not None, x['name']))
        faces = np.array(faces_list)

        # Write the points file (vertex coordinates)
        with open(join(self.target, 'points'), 'w') as f:
            f.write(self._header('vectorField', 'points'))
            f.write(str(model.ncps) + '\n(\n')
            for pt in model.cps():
                f.write('({})\n'.format(' '.join(str(p) for p in pt)))
            f.write(')\n')

        # Write the faces file (four vertex indices for each face)
        with open(join(self.target, 'faces'), 'w') as f:
            f.write(self._header('faceList', 'faces'))
            f.write(str(len(faces)) + '\n(\n')
            for face in faces['nodes']:
                f.write('({})\n'.format(' '.join(str(f) for f in face)))
            f.write(')\n')

        # Write the owner and neighbour files (cell indices for each face)
        with open(join(self.target, 'owner'), 'w') as f:
            f.write(self._header('labelList', 'owner', note=note))
            f.write(str(len(faces)) + '\n(\n')
            for owner in faces['owner']:
                f.write(str(owner) + '\n')
            f.write(')\n')
        with open(join(self.target, 'neighbour'), 'w') as f:
            f.write(self._header('labelList', 'neighbour', note=note))
            f.write(str(len(faces)) + '\n(\n')
            for neighbor in faces['neighbor']:
                f.write(str(neighbor) + '\n')
            f.write(')\n')

        # Write the boundary file
        with open(join(self.target, 'boundary'), 'w') as f:
            f.write(self._header('polyBoundaryMesh', 'boundary'))
            f.write(str(len(set(faces['name'])) - 1) + '\n(\n')
            start = 0
            for name, it in groupby(faces, key=itemgetter('name')):
                nfaces = len(list(it))
                if name is None:
                    start += nfaces
                    continue
                f.write(name + '\n{\n')
                f.write('    type patch;\n')
                f.write('    nFaces {};\n'.format(nfaces))
                f.write('    startFace {};\n'.format(start))
                f.write('}\n')
                start += nfaces
            f.write(')\n')
