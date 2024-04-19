
class RGToolsInternalException(Exception):
    def __init__(self, *args: object) -> None:
        super().__init__(*args)

class GTFHandleFilterException(RGToolsInternalException):
    def __init__(self, *args: object) -> None:
        super().__init__(*args)

class GTFRecordNoFeatureException(RGToolsInternalException):
    def __init__(self, *args: object) -> None:
        super().__init__(*args)

class BedTableException(RGToolsInternalException):
    def __init__(self, *args: object) -> None:
        super().__init__(*args)

class BedTableLoadException(BedTableException):
    def __init__(self, *args: object) -> None:
        super().__init__(*args)