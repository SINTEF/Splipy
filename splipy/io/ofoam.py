from __future__ import annotations

from itertools import groupby
from operator import itemgetter
from pathlib import Path
from typing import TYPE_CHECKING, Optional, Sequence, Union

import numpy as np
from typing_extensions import Self

from splipy.splinemodel import SplineModel

from .master import MasterIO

if TYPE_CHECKING:
    from types import TracebackType

    from splipy.splineobject import SplineObject


class OpenFOAM(MasterIO):
    target: Path

    def __init__(self, target: Union[Path, str]) -> None:
        self.target = Path(target)

    def __enter__(self) -> Self:
        # Create the target directory if it does not exist
        if not self.target.exists():
            self.target.mkdir(parents=True)

        # If it does, ensure that it's a directory
        elif not self.target.is_dir():
            raise FileExistsError(f"{self.target} exists and is not a directory")

        return self

    def __exit__(
        self,
        exc_type: Optional[type[BaseException]],
        exc_val: Optional[BaseException],
        exc_tb: Optional[TracebackType],
    ) -> None:
        return

    def _header(self, cls: str, obj: str, note: Optional[str] = None) -> str:
        s = "FoamFile\n{\n"
        s += "    version     2.0;\n"
        s += "    format      ascii;\n"
        s += f"    class       {cls};\n"
        if note:
            s += f'    note        "{note}";\n'
        s += f"    object      {obj};\n"
        s += "}\n"
        return s

    def write(self, model: Union[SplineObject, Sequence[SplineObject], SplineModel]) -> None:
        assert isinstance(model, SplineModel), "OpenFOAM.write only supports SplineModel objects"

        # Only linear volumes in 3D, please
        assert model.pardim == 3
        assert model.dimension == 3
        assert all(patch.obj.order() == (2, 2, 2) for patch in model.catalogue.top_nodes())

        # Generate all the information we need
        model.generate_cp_numbers()
        model.generate_cell_numbers()
        faces = model.faces()
        ninternal = sum(faces["name"] == None)  # noqa: E711
        note = "nPoints: {} nCells: {} nFaces: {} nInternalFaces: {}".format(
            model.ncps, model.ncells, len(faces), ninternal
        )

        # OpenFOAM is very particular about face ordering
        # In order of importance:
        # - Internal faces before boundary faces
        # - All faces in the same boundary must be contiguous
        # - Low number owners before high number owners
        # - Low number neighbors before high number neighbors
        faces_list = list(faces)
        faces_list = sorted(faces_list, key=itemgetter("neighbor"))
        faces_list = sorted(faces_list, key=itemgetter("owner"))
        faces_list = sorted(faces_list, key=lambda x: (x["name"] is not None, x["name"]))
        faces = np.array(faces_list)

        # Write the points file (vertex coordinates)
        with (self.target / "points").open("w") as f:
            f.write(self._header("vectorField", "points"))
            f.write(str(model.ncps) + "\n(\n")
            for pt in model.cps():
                f.write("({})\n".format(" ".join(str(p) for p in pt)))
            f.write(")\n")

        # Write the faces file (four vertex indices for each face)
        with (self.target / "faces").open("w") as f:
            f.write(self._header("faceList", "faces"))
            f.write(str(len(faces)) + "\n(\n")
            for face in faces["nodes"]:
                f.write("({})\n".format(" ".join(str(f) for f in face)))
            f.write(")\n")

        # Write the owner and neighbour files (cell indices for each face)
        with (self.target / "owner").open("w") as f:
            f.write(self._header("labelList", "owner", note=note))
            f.write(str(len(faces)) + "\n(\n")
            for owner in faces["owner"]:
                f.write(str(owner) + "\n")
            f.write(")\n")

        with (self.target / "neighbour").open("w") as f:
            f.write(self._header("labelList", "neighbour", note=note))
            f.write(str(len(faces)) + "\n(\n")
            for neighbor in faces["neighbor"]:
                f.write(str(neighbor) + "\n")
            f.write(")\n")

        # Write the boundary file
        with (self.target / "boundary").open("w") as f:
            f.write(self._header("polyBoundaryMesh", "boundary"))
            f.write(str(len(set(faces["name"])) - 1) + "\n(\n")
            start = 0
            for name, it in groupby(faces, key=itemgetter("name")):
                nfaces = len(list(it))
                if name is None:
                    start += nfaces
                    continue
                f.write(name + "\n{\n")
                f.write("    type patch;\n")
                f.write(f"    nFaces {nfaces};\n")
                f.write(f"    startFace {start};\n")
                f.write("}\n")
                start += nfaces
            f.write(")\n")
