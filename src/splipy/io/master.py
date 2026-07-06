from __future__ import annotations

from typing import TYPE_CHECKING, Self

if TYPE_CHECKING:
    from collections.abc import Sequence
    from types import TracebackType

    from splipy.splinemodel import SplineModel
    from splipy.splineobject import SplineObject


class MasterIO:
    def __init__(self, filename: str) -> None:
        """Create an IO object attached to a file.

        :param str filename: The file to read from or write to
        """
        raise NotImplementedError()

    def __enter__(self) -> Self:
        raise NotImplementedError()

    def __exit__(
        self,
        exc_type: type[BaseException],
        exc_value: BaseException,
        traceback: TracebackType,
    ) -> None:
        pass

    def write(self, obj: SplineObject | Sequence[SplineObject] | SplineModel) -> None:
        """Write one or more objects to the file.

        :param obj: The object(s) to write
        :type obj: [:class:`splipy.SplineObject`] or :class:`splipy.SplineObject`
        """
        raise NotImplementedError()

    def read(self) -> list[SplineObject]:
        """Reads all the objects from the file.

        :return: Objects
        :rtype: [:class:`splipy.SplineObject`]
        """
        raise NotImplementedError()
