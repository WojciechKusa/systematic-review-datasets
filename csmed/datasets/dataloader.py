import pathlib
from dataclasses import dataclass, field
from importlib import import_module
from types import ModuleType
from typing import Callable, Iterable, List, Optional

import datasets
from bigbio.dataloader import SCHEMA_TO_METADATA_CLS
from bigbio.utils.configs import BigBioConfig
from bigbio.utils.constants import SCHEMA_TO_TASKS
from datasets import load_dataset

# TODO: update this as fixes come in
_CURRENTLY_BROKEN_NAMES = {
    "nagel_source",
    "nagel_bigbio_kb",
    "pcr_source",
    "pcr_fixed_source",
    "pcr_bigbio_kb",
}

# large datasets take greater than ~ 10 minutes to load
_LARGE_CONFIG_NAMES = {
    "biomrc_large_A_source",
    "biomrc_large_B_source",
    "biomrc_large_A_bigbio_qa",
    "biomrc_large_B_bigbio_qa",
    "medal_source",
    "medal_bigbio_kb",
    "meddialog_zh_source",
    "meddialog_zh_bigbio_text",
    "pubtator_central_source",
    "pubtator_central_bigbio_kb",
}

# resource datasets are widely used but not annotated by human experts
# e.g. PubTator and MIMIC III
_RESOURCE_CONFIG_NAMES = {
    "pubtator_central_sample_source",
    "pubtator_central_sample_bigbio_kb",
    "pubtator_central_source",
    "pubtator_central_bigbio_kb",
}


@dataclass
class BigBioConfigHelper:
    """Metadata for one config of a dataset."""

    script: pathlib.Path
    dataset_name: str
    tasks: List[str]
    languages: List[str]
    config: BigBioConfig
    is_local: bool
    is_pubmed: bool
    is_bigbio_schema: bool
    bigbio_schema_caps: Optional[str]
    is_large: bool
    is_resource: bool
    is_default: bool
    is_broken: bool
    bigbio_version: str
    source_version: str
    citation: str
    description: str
    homepage: str
    display_name: str
    license: str

    _ds_module: datasets.load.DatasetModule = field(repr=False)
    _py_module: ModuleType = field(repr=False)
    _ds_cls: type = field(repr=False)

    def get_load_dataset_kwargs(
        self,
        from_hub=True,
        **extra_load_dataset_kwargs,
    ):
        if from_hub:
            path = f"bigbio/{self.dataset_name}"
        else:
            path = self.script
        return {
            "path": path,
            "name": self.config.name,
            **extra_load_dataset_kwargs,
        }

    def load_dataset(
        self,
        from_hub=False,  #  fixme: in bigbio is True
        **extra_load_dataset_kwargs,
    ):
        load_dataset_kwargs = self.get_load_dataset_kwargs(from_hub=from_hub)
        return load_dataset(
            **load_dataset_kwargs,
            **extra_load_dataset_kwargs,
        )

    def get_metadata(self, **extra_load_dataset_kwargs):
        if not self.is_bigbio_schema:
            raise ValueError("only supported for bigbio schemas")
        dsd = self.load_dataset(**extra_load_dataset_kwargs)
        split_metas = {}
        for split, ds in dsd.items():
            meta = SCHEMA_TO_METADATA_CLS[self.config.schema].from_dataset(ds)
            split_metas[split] = meta
        return split_metas


class BigBioConfigHelpers:
    """
    Handles creating and filtering BigBioDatasetConfigHelper instances.
    """

    def __init__(
        self,
        helpers: Optional[Iterable[BigBioConfigHelper]] = None,
        keep_broken: bool = False,
    ):

        path_to_here = pathlib.Path(__file__).parent.absolute()
        self.path_to_biodatasets = (
            path_to_here / "datasets"
        ).resolve()  # fixme hardcoded path
        self.dataloader_directories = sorted(
            [
                path
                for path in self.path_to_biodatasets.glob("*")
                if path.name != "__init__.py" and path.name != "__pycache__"
            ]
        )
        self.dataloader_scripts = [
            dpath / f"{dpath.name}.py" for dpath in self.dataloader_directories
        ]

        # if helpers are passed in, just attach and go
        if helpers is not None:
            if keep_broken:
                self._helpers = helpers
            else:
                self._helpers = [helper for helper in helpers if not helper.is_broken]
            return

        # print(self.dataloader_scripts)
        # otherwise, create all helpers available in package
        helpers = []
        for dataloader_script in self.dataloader_scripts:
            dataset_name = dataloader_script.stem
            py_module = import_module(
                f"csmed.datasets.datasets.{dataset_name}.{dataset_name}"
            )  # fixme hardcoded path
            ds_module = datasets.load.dataset_module_factory(
                dataloader_script.as_posix()
            )
            ds_cls = datasets.load.import_main_class(ds_module.module_path)

            for config in ds_cls.BUILDER_CONFIGS:

                is_bigbio_schema = config.schema.startswith("bigbio")
                if is_bigbio_schema:
                    # some bigbio datasets support tasks from multiple schemas
                    # here we just choose the tasks for the schema of this config
                    bigbio_schema_caps = config.schema.split("_")[1].upper()
                    schema_tasks = set(
                        [task.name for task in SCHEMA_TO_TASKS[bigbio_schema_caps]]
                    )
                    config_tasks = set(
                        [task.name for task in py_module._SUPPORTED_TASKS]
                    )
                    tasks = schema_tasks & config_tasks

                else:
                    tasks = set([task.name for task in py_module._SUPPORTED_TASKS])
                    bigbio_schema_caps = None

                helpers.append(
                    BigBioConfigHelper(
                        script=dataloader_script.as_posix(),
                        dataset_name=dataset_name,
                        tasks=tasks,
                        languages=py_module._LANGUAGES,
                        config=config,
                        is_local=py_module._LOCAL,
                        is_pubmed=py_module._PUBMED,
                        is_bigbio_schema=is_bigbio_schema,
                        bigbio_schema_caps=bigbio_schema_caps,
                        is_large=config.name in _LARGE_CONFIG_NAMES,
                        is_resource=config.name in _RESOURCE_CONFIG_NAMES,
                        is_default=config.name == ds_cls.DEFAULT_CONFIG_NAME,
                        is_broken=config.name in _CURRENTLY_BROKEN_NAMES,
                        bigbio_version=py_module._BIGBIO_VERSION,
                        source_version=py_module._SOURCE_VERSION,
                        citation=py_module._CITATION,
                        description=py_module._DESCRIPTION,
                        homepage=py_module._HOMEPAGE,
                        display_name=py_module._DISPLAYNAME,
                        license=py_module._LICENSE,
                        _ds_module=ds_module,
                        _py_module=py_module,
                        _ds_cls=ds_cls,
                    )
                )

        if keep_broken:
            self._helpers = helpers
        else:
            self._helpers = [helper for helper in helpers if not helper.is_broken]

    @property
    def available_dataset_names(self) -> List[str]:
        return sorted(list(set([helper.dataset_name for helper in self])))

    def for_dataset(self, dataset_name: str) -> "BigBioConfigHelpers":
        helpers = [helper for helper in self if helper.dataset_name == dataset_name]
        if len(helpers) == 0:
            raise ValueError(f"no helper with helper.dataset_name = {dataset_name}")
        return BigBioConfigHelpers(helpers=helpers)

    def for_config_name(self, config_name: str) -> "BigBioConfigHelper":
        helpers = [helper for helper in self if helper.config.name == config_name]
        if len(helpers) == 0:
            raise ValueError(f"no helper with helper.config.name = {config_name}")
        if len(helpers) > 1:
            raise ValueError(
                f"multiple helpers with helper.config.name = {config_name}"
            )
        return helpers[0]

    def default_for_dataset(self, dataset_name: str) -> BigBioConfigHelper:
        helpers = [
            helper
            for helper in self
            if helper.is_default and helper.dataset_name == dataset_name
        ]
        assert len(helpers) == 1
        return helpers[0]

    def filtered(
        self, is_keeper: Callable[[BigBioConfigHelper], bool]
    ) -> "BigBioConfigHelpers":
        """Return dataset config helpers that match is_keeper."""
        return BigBioConfigHelpers(
            helpers=[helper for helper in self if is_keeper(helper)]
        )

    def __repr__(self):
        return "\n\n".join([helper.__repr__() for helper in self])

    def __str__(self):
        return self.__repr__()

    def __iter__(self):
        for helper in self._helpers:
            yield helper

    def __len__(self):
        return len(self._helpers)

    def __getitem__(self, key):
        if isinstance(key, slice):
            start, stop, step = key.indices(len(self))
            return BigBioConfigHelpers(
                helpers=[self._helpers[ii] for ii in range(start, stop, step)]
            )
        elif isinstance(key, int):
            if key < 0:  # Handle negative indices
                key += len(self)
            if key < 0 or key >= len(self):
                raise IndexError(f"The index ({key}) is out of range.")
            return self._helpers[key]
        else:
            raise TypeError("Invalid argument type.")
