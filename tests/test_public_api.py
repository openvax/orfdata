"""Every public name is importable from ``orfdata``, although its submodule is imported on first use."""

import ast
from importlib import import_module
from pathlib import Path

import pytest

import orfdata
from orfdata import EXPORTED_FROM


def type_checking_imports():
    """Return {name: submodule} for the imports under ``if TYPE_CHECKING:`` in ``orfdata/__init__.py``."""
    tree = ast.parse(Path(orfdata.__file__).read_text())
    block = next(node for node in tree.body if isinstance(node, ast.If) and ast.unparse(node.test) == "TYPE_CHECKING")
    return {alias.name: node.module for node in block.body if isinstance(node, ast.ImportFrom) for alias in node.names}


@pytest.mark.parametrize("name", orfdata.__all__)
def test_every_public_name_resolves_to_its_submodule_object(name):
    value = getattr(orfdata, name)
    assert value is not None
    if name in EXPORTED_FROM:
        assert value is getattr(import_module(f"orfdata.{EXPORTED_FROM[name]}"), name)
    namespace = {}
    exec(f"from orfdata import {name}", namespace)
    assert namespace[name] is value


def test_all_lists_every_public_name_once_and_dir_shows_them():
    assert len(orfdata.__all__) == len(set(orfdata.__all__))
    assert set(EXPORTED_FROM) | {"__version__", "viral_subsets", "normalize_viral_subsets",
                                 "viral_subset_abbreviations"} == set(orfdata.__all__)
    assert set(orfdata.__all__) <= set(dir(orfdata))
    namespace = {}
    exec("from orfdata import *", namespace)
    assert set(orfdata.__all__) <= set(namespace)


def test_type_checkers_see_the_same_names_as_the_lazy_imports():
    assert type_checking_imports() == EXPORTED_FROM


def test_viral_subsets_stays_the_function_after_its_module_is_imported():
    import orfdata.viral_panel  # noqa: F401  Imports the orfdata.viral_subsets module.
    from orfdata.viral_subsets import viral_subsets

    assert orfdata.viral_subsets is viral_subsets and callable(orfdata.viral_subsets)


def test_submodules_are_attributes_and_unknown_names_are_not():
    from orfdata import liftover

    assert liftover is import_module("orfdata.liftover") is orfdata.liftover
    assert orfdata.build is import_module("orfdata.build")
    with pytest.raises(AttributeError, match="no attribute 'not_a_name'"):
        orfdata.not_a_name  # noqa: B018
    with pytest.raises(ImportError):
        exec("from orfdata import not_a_name", {})
    assert not hasattr(orfdata, "__not_a_dunder__")
