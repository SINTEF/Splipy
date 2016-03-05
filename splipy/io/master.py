class MasterIO(object):

    def __init__(self, filename):
        """Create an IO object attached to a file.

        :param str filename: The file to read from or write to
        """
        raise NotImplementedError()

    def __enter__(self):
        raise NotImplementedError()

    def write(self, obj):
        """Write one or more objects to the file.

        :param obj: The object(s) to write
        :type obj: [:class:`splipy.SplineObject`] or :class:`splipy.SplineObject`
        """
        raise NotImplementedError()

    def read(self):
        """Reads all the objects from the file.

        :return: Objects
        :rtype: [:class:`splipy.SplineObject`]
        """
        raise NotImplementedError()
