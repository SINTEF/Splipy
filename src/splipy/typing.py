from __future__ import annotations

from typing import SupportsFloat, TypeAlias, TypeVar

import numpy as np
import numpy.typing as npt

B = TypeVar("B", bound=npt.NBitBase)

ArrayLike: TypeAlias = npt.ArrayLike
ScalarLike: TypeAlias = SupportsFloat
FloatArray: TypeAlias = npt.NDArray[np.floating]
