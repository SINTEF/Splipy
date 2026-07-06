from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, ClassVar, Literal, Self, TextIO, cast, overload

import numpy as np
from numpy import pi, savetxt

from splipy import curve_factory, state, surface_factory
from splipy.basis import BSplineBasis
from splipy.splinemodel import SplineModel
from splipy.splineobject import SplineObject
from splipy.surface import Surface
from splipy.trimmedsurface import TrimmedSurface
from splipy.utils import flip_and_move_plane_geometry, rotate_local_x_axis

from .master import MasterIO

if TYPE_CHECKING:
    from collections.abc import Sequence
    from types import TracebackType

    from splipy.curve import Curve
    from splipy.typing import FloatArray
    from splipy.volume import Volume


class G2(MasterIO):
    fstream: TextIO
    filename: str
    trimming_curves: list[Curve]
    onlywrite: bool

    g2_type: ClassVar[list[int]] = [100, 200, 700]  # curve, surface, volume identifiers

    def read_next_non_whitespace(self) -> str:
        line = next(self.fstream).strip()
        while not line:
            line = next(self.fstream).strip()
        return line

    def read_next_param_range(self) -> tuple[float, float]:
        start, end = map(float, next(self.fstream).split())
        return start, end

    def read_next_bool(self) -> bool:
        return next(self.fstream).strip() != "0"

    def read_next_float(self) -> float:
        return float(next(self.fstream).strip())

    def read_next_array(self) -> FloatArray:
        return np.array(next(self.fstream).split(), dtype=float)

    def circle(self) -> Curve:
        int(self.read_next_non_whitespace().strip())
        r = self.read_next_float()
        center = self.read_next_array()
        normal = self.read_next_array()
        xaxis = self.read_next_array()
        param = self.read_next_param_range()
        reverse = self.read_next_bool()

        result = curve_factory.circle(r=r, center=center, normal=normal, xaxis=xaxis)
        result.reparam(param)
        if reverse:
            result.reverse()
        return result

    def ellipse(self) -> Curve:
        int(self.read_next_non_whitespace().strip())
        r1 = self.read_next_float()
        r2 = self.read_next_float()
        center = self.read_next_array()
        normal = self.read_next_array()
        xaxis = self.read_next_array()
        param = self.read_next_param_range()
        reverse = self.read_next_bool()

        result = curve_factory.ellipse(r1=r1, r2=r2, center=center, normal=normal, xaxis=xaxis)
        result.reparam(param)
        if reverse:
            result.reverse()
        return result

    def line(self) -> Curve:
        int(self.read_next_non_whitespace().strip())
        start = self.read_next_array()
        direction = self.read_next_array()
        finite = self.read_next_bool()
        param = self.read_next_param_range()
        reverse = self.read_next_bool()
        if not finite:
            param = (-state.unlimited, state.unlimited)

        result = curve_factory.line(start + direction * param[0], start + direction * param[1])
        if reverse:
            result.reverse()
        return result

    #   def cone(self):
    #       dim      = int(     self.read_next_non_whitespace().strip())
    #       r        = float(   next(self.fstream).strip())
    #       center   = np.array(next(self.fstream).split(' '), dtype=float)
    #       z_axis   = np.array(next(self.fstream).split(' '), dtype=float)
    #       x_axis   = np.array(next(self.fstream).split(' '), dtype=float)
    #       angle    = float(   next(self.fstream).strip())
    #       finite   =          next(self.fstream).strip() != '0'
    #       param_u  = np.array(next(self.fstream).split(' '), dtype=float)
    #       if finite:
    #           param_v=np.array(next(self.fstream).split(' '), dtype=float)

    def cylinder(self) -> Surface:
        int(self.read_next_non_whitespace().strip())
        r = self.read_next_float()
        center = self.read_next_array()
        z_axis = self.read_next_array()
        x_axis = self.read_next_array()
        finite = self.read_next_bool()
        param_u = self.read_next_param_range()
        param_v = self.read_next_param_range() if finite else (-state.unlimited, state.unlimited)
        swap = self.read_next_bool()

        center = center + z_axis * param_v[0]
        h = param_v[1] - param_v[0]
        result = surface_factory.cylinder(r=r, center=center, xaxis=x_axis, axis=z_axis, h=h)
        result.reparam(param_u, param_v)
        if swap:
            result.swap()
        return result

    def disc(self) -> Surface:
        int(self.read_next_non_whitespace().strip())
        center = self.read_next_array()
        r = self.read_next_float()
        z_axis = self.read_next_array()
        x_axis = self.read_next_array()
        degen = self.read_next_bool()
        angles = [self.read_next_float() for _ in range(4)]
        param_u = self.read_next_param_range()
        param_v = self.read_next_param_range()
        swap = self.read_next_bool()

        if degen:
            result = surface_factory.disc(r=r, center=center, xaxis=x_axis, normal=z_axis, type="radial")
        else:
            if not (np.allclose(np.diff(angles), pi / 2, atol=1e-10)):
                raise RuntimeError("Unknown square parametrization of disc elementary surface")
            result = surface_factory.disc(r=r, center=center, xaxis=x_axis, normal=z_axis, type="square")
        result.reparam(param_u, param_v)
        if swap:
            result.swap()
        return result

    def plane(self) -> Surface:
        int(self.read_next_non_whitespace().strip())
        center = self.read_next_array()
        normal = self.read_next_array()
        x_axis = self.read_next_array()
        finite = self.read_next_bool()
        if finite:
            param_u = self.read_next_param_range()
            param_v = self.read_next_param_range()
        else:
            param_u = (-state.unlimited, +state.unlimited)
            param_v = (-state.unlimited, +state.unlimited)
        swap = self.read_next_bool()

        result = Surface() * [param_u[1] - param_u[0], param_v[1] - param_v[0]] + [param_u[0], param_v[0]]
        result.rotate(rotate_local_x_axis(x_axis, normal))
        result = flip_and_move_plane_geometry(result, center, normal)
        result.reparam(param_u, param_v)
        if swap:
            result.swap()
        return result

    def torus(self) -> Surface:
        int(self.read_next_non_whitespace().strip())
        r2 = self.read_next_float()
        r1 = self.read_next_float()
        center = self.read_next_array()
        z_axis = self.read_next_array()
        x_axis = self.read_next_array()
        self.read_next_bool()  # I have no idea what this does :(
        param_u = self.read_next_param_range()
        param_v = self.read_next_param_range()
        swap = self.read_next_bool()

        result = surface_factory.torus(minor_r=r1, major_r=r2, center=center, normal=z_axis, xaxis=x_axis)
        result.reparam(param_u, param_v)
        if swap:
            result.swap()
        return result

    def sphere(self) -> Surface:
        int(self.read_next_non_whitespace().strip())
        r = self.read_next_float()
        center = self.read_next_array()
        z_axis = self.read_next_array()
        x_axis = self.read_next_array()
        param_u = self.read_next_param_range()
        param_v = self.read_next_param_range()
        swap = self.read_next_bool()

        result = surface_factory.sphere(r=r, center=center, xaxis=x_axis, zaxis=z_axis).swap()
        if swap:
            result.swap()
        result.reparam(param_u, param_v)
        return result

    @overload
    def splines(self, pardim: Literal[1]) -> Curve: ...

    @overload
    def splines(self, pardim: Literal[2]) -> Surface: ...

    @overload
    def splines(self, pardim: Literal[3]) -> Volume: ...

    @overload
    def splines(self, pardim: int) -> SplineObject: ...

    def splines(self, pardim: int) -> SplineObject:
        _, rational = self.read_next_non_whitespace().strip().split()

        bases = [self.read_basis() for _ in range(pardim)]
        ncps = 1
        for b in bases:
            ncps *= b.num_functions()

        cps = [tuple(map(float, next(self.fstream).split())) for _ in range(ncps)]
        return SplineObject.construct_subclass(bases, cps, bool(int(rational)), raw=False)

    def surface_of_linear_extrusion(self) -> Surface:
        int(self.read_next_non_whitespace().strip())
        crv = self.splines(1)
        normal = np.array(self.read_next_non_whitespace().split(), dtype=float)
        finite = self.read_next_bool()
        param_u = self.read_next_param_range()
        param_v = self.read_next_param_range() if finite else (-state.unlimited, +state.unlimited)
        swap = self.read_next_bool()

        result = surface_factory.extrude(crv + normal * param_v[0], normal * (param_v[1] - param_v[0]))
        result.reparam(param_u, param_v)

        if swap:
            result.swap()
        return result

    def bounded_surface(self) -> TrimmedSurface:
        objtype = int(next(self.fstream).strip())

        # create the underlying surface which all trimming curves are to be applied
        if objtype in G2.g2_generators:
            constructor = getattr(self, G2.g2_generators[objtype].__name__)
            surface = constructor()
        elif objtype == 200:
            surface = self.splines(2)
        else:
            raise OSError("Unsopported trimmed surface or malformed input file")

        # for all trimming loops
        numb_loops = int(self.read_next_non_whitespace())
        all_loops = []
        for i in range(numb_loops):
            # for all cuve pieces of that loop
            numb_crvs, space_epsilon = next(self.fstream).split()
            state.parametric_absolute_tolerance = float(space_epsilon)
            one_loop = []
            for j in range(int(numb_crvs)):
                # read a physical and parametric representation of the same curve
                _, parameter_curve_type, space_curve_type = map(int, self.read_next_non_whitespace().split())
                two_curves: list[Curve] = []
                for crv_type in [parameter_curve_type, space_curve_type]:
                    if crv_type in G2.g2_generators:
                        constructor = getattr(self, G2.g2_generators[crv_type].__name__)
                        crv = constructor()
                    elif crv_type == 100:
                        crv = self.splines(1)
                    else:
                        raise OSError("Unsopported trimming curve or malformed input file")
                    two_curves.append(crv)

                # only keep the parametric version (re-generate physical one if we need it)
                one_loop.append(two_curves[0])
                self.trimming_curves.append(two_curves[1])
            all_loops.append(one_loop)

        return TrimmedSurface(
            surface.bases[0], surface.bases[1], surface.controlpoints, surface.rational, all_loops, raw=True
        )

    g2_generators = {
        120: line,
        130: circle,
        140: ellipse,
        260: cylinder,
        292: disc,
        270: sphere,
        290: torus,
        250: plane,
        210: bounded_surface,
        261: surface_of_linear_extrusion,
    }  # , 280:cone

    def __init__(self, filename: str) -> None:
        if filename[-3:] != ".g2":
            filename += ".g2"
        self.filename = filename
        self.trimming_curves = []

    def __enter__(self) -> Self:
        return self

    def write(self, obj: Sequence[SplineObject] | SplineObject | SplineModel) -> None:
        """Write the object in GoTools format."""
        if not hasattr(self, "fstream"):
            self.onlywrite = True
            self.fstream = Path(self.filename).open("w")
        if not self.onlywrite:
            raise OSError(f"Could not write to file {self.filename}")

        if isinstance(obj, SplineModel):
            for o in obj.objects():
                self.write(o)
            return

        if not isinstance(obj, SplineObject):
            for o in obj:
                self.write(o)
            return

        assert isinstance(obj, SplineObject)

        for i in range(obj.pardim):
            if obj.periodic(i):
                # TODO(Eivind): Use a type-safe version of split that returns a SplineObject.
                obj = cast("SplineObject", obj.split(obj.start(i), i))

        self.fstream.write(f"{G2.g2_type[obj.pardim - 1]} 1 0 0\n")
        self.fstream.write(f"{obj.dimension} {int(obj.rational)}\n")
        for b in obj.bases:
            self.fstream.write(f"{len(b.knots) - b.order} {b.order}\n")
            self.fstream.write(" ".join(f"{k:.16g}" for k in b.knots))
            self.fstream.write("\n")

        savetxt(
            self.fstream,
            obj.controlpoints.reshape(-1, obj.dimension + obj.rational, order="F"),
            fmt="%.16g",
            delimiter=" ",
            newline="\n",
        )

    def read(self) -> list[SplineObject]:
        if not hasattr(self, "fstream"):
            self.onlywrite = False
            self.fstream = Path(self.filename).open()

        if self.onlywrite:
            raise OSError(f"Could not read from file {self.filename}")

        result: list[SplineObject] = []

        for line in self.fstream:
            line = line.strip()
            if not line:
                continue

            # read object type
            objtype, major, minor, patch = map(int, line.split())
            if (major, minor, patch) != (1, 0, 0):
                raise OSError("Unknown G2 format")

            # if obj type is in factory methods (cicle, torus etc), create it now
            if objtype in G2.g2_generators:
                constructor = getattr(self, G2.g2_generators[objtype].__name__)
                result.append(constructor())
                continue

            # for "normal" splines (Curves, Surfaces, Volumes) create it now
            pardim = [i for i in range(len(G2.g2_type)) if G2.g2_type[i] == objtype]
            if not pardim:
                raise OSError(f"Unknown G2 object type {objtype}")
            result.append(self.splines(pardim[0] + 1))

        return result

    def read_basis(self) -> BSplineBasis:
        ncps, order = map(int, next(self.fstream).split())
        kts = list(map(float, next(self.fstream).split()))
        return BSplineBasis(order, kts, -1)

    def __exit__(
        self,
        exc_type: type[BaseException],
        exc_value: BaseException,
        traceback: TracebackType,
    ) -> None:
        if hasattr(self, "fstream"):
            self.fstream.close()
