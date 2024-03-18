
class RGToolsInternalException(Exception):
    def __init__(self, *args: object) -> None:
        super().__init__(*args)

class GTFHandleFilterException(RGToolsInternalException):
    def __init__(self, *args: object) -> None:
        super().__init__(*args)

class GTFRecordNoFeatureException(RGToolsInternalException):
    def __init__(self, *args: object) -> None:
        super().__init__(*args)