from abc import ABC, abstractmethod
from typing import Union, Sequence, Type, Optional
from types import TracebackType
from pathlib import Path

from typing_extensions import Self

from ..splineobject import SplineObject


class MasterIO:

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
        exc_type: Optional[Type[BaseException]],
        exc_val: Optional[BaseException],
        exc_tb: Optional[TracebackType]
    ) -> None:
        ...

    @abstractmethod
    def write(self, obj: Union[SplineObject, Sequence[SplineObject]]) -> None:
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
