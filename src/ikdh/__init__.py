from ._ikdh import (
    DHTable, JointLimits, Solver, Robot,
    forward_kin, pose_from_xyzrpw, fk_error,
    load_robot as _load_robot,
)
import os as _os

def robots_dir():
    """Path to the bundled robots directory."""
    return _os.path.join(_os.path.dirname(__file__), "robots")

def load_robot(path):
    """Load a robot YAML file.

    Accepts an absolute path, a relative path (resolved from cwd), or a bare
    filename (e.g. "gofa5.yaml") which is looked up in the bundled robots dir.
    """
    if not _os.path.isabs(path) and not _os.path.exists(path):
        candidate = _os.path.join(robots_dir(), _os.path.basename(path))
        if _os.path.exists(candidate):
            path = candidate
    return _load_robot(path)
