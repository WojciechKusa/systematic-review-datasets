import abc


class ScreeningDataset(abc.ABC):
    output_directory: str

    @abc.abstractmethod
    def __init__(self, output_directory: str):
        self.output_directory = output_directory

    @abc.abstractmethod
    def is_prepared(self) -> bool:
        pass

    @abc.abstractmethod
    def prepare_dataset(self) -> None:
        pass
