from __future__ import annotations

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING, Optional, Sequence, Union

from typing_extensions import Self

if TYPE_CHECKING:
    from pathlib import Path
    from types import TracebackType

    from splipy.splinemodel import SplineModel
    from splipy.splineobject import SplineObject


class MasterIO(ABC):
    def __init__(self, filename: Union[str, Path]) -> None:
        """Create an IO object attached to a file.

        :param str filename: The file to read from or write to
        """
        raise NotImplementedError()

    @abstractmethod
    def __enter__(self) -> Self:
        ...

    @abstractmethod
    def __exit__(
        self,
        exc_type: Optional[type[BaseException]],
        exc_val: Optional[BaseException],
        exc_tb: Optional[TracebackType],
    ) -> None:
        ...

    @abstractmethod
    def write(self, obj: Union[SplineObject, Sequence[SplineObject], SplineModel]) -> None:
        """Write one or more objects to the file.

        :param obj: The object(s) to write
        :type obj: [:class:`splipy.SplineObject`] or :class:`splipy.SplineObject`
        """

    @abstractmethod
    def read(self) -> list[SplineObject]:
        """Reads all the objects from the file.

        :return: Objects
        :rtype: [:class:`splipy.SplineObject`]
        """
