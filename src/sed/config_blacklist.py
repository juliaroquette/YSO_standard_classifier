# src/sed/config_blacklist.py
#
# This is a module that handles the list of tables from VizieR which have been blacklisted
#
import os
import yaml

def load_black_list(yaml_path):
    """
    Load a blacklist YAML file.

    Returns a dict where each top-level key maps to a list of values.
    If the file does not exist or is not valid YAML, an empty dict is returned.

    Example YAML:
        outdated_surveys:
          - II/349/apass9
        weird_surveys:
          - SomeSurvey
        custom_skip:
          - MySpecialCase

    This function does NOT enforce any specific key names.
    Your script decides which keys it wants to apply.
    """
    # if the path doesn’t exist → return empty
    if not os.path.exists(yaml_path):
        return {}

    try:
        with open(yaml_path, "r") as f:
            cfg = yaml.safe_load(f) or {}
    except Exception:
        # YAML parse error or permission issue → return empty
        return {}

    # ensure the result is a dict
    if not isinstance(cfg, dict):
        return {}

    # normalize single strings into lists
    # e.g. key: value → key: [value]
    normalized = {}
    for key, value in cfg.items():
        if value is None:
            normalized[key] = []
        elif isinstance(value, str):
            normalized[key] = [value]
        elif isinstance(value, (list, tuple)):
            normalized[key] = list(value)
        else:
            # unexpected value type → ignore
            normalized[key] = []

    return normalized


def tabnames_to_drop(cfg, keys=None):
    """
    From a blacklist config dict, build a set of tabnames to drop.

    Parameters
    ----------
    cfg : dict
        Dict returned by load_black_list()
    keys : list or tuple (optional)
        If provided, only collect values under these keys.
        If None, collect values from all keys in cfg.

    Returns
    -------
    set of strings
    """
    if not isinstance(cfg, dict):
        return set()

    tabs = set()

    if keys is None:
        key_iter = list(cfg.keys())
    else:
        # only use specified keys
        key_iter = list(keys)

    for key in key_iter:
        for tab in cfg.get(key, []):
            if tab is not None:
                tabs.add(str(tab))

    return tabs


def apply_blacklist(source, cfg, keys=None):
    """
    Apply blacklist to a Source object.

    Parameters
    ----------
    source :
        An object implementing dropCatalog(tabname)
    cfg : dict
        Loaded blacklist config
    keys : list/tuple of keys (optional)
        If provided, only apply those categories
        If None, all keys in cfg are applied.

    Example usage:
        keys_to_apply = ["outdated_surveys"]
        apply_blacklist(star, cfg, keys_to_apply)
    """
    tabs = tabnames_to_drop(cfg, keys)
    if not tabs:
        return

    # call dropCatalog once per tabname
    for tab in tabs:
        source.dropCatalog(tab)
