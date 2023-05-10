import abc


class ScreeningDataset(abc.ABC):
    @abc.abstractmethod
    def is_prepared(self) -> bool:
        pass

    @abc.abstractmethod
    def prepare_dataset(self) -> None:
        pass
