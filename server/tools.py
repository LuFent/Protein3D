import os
import importlib.util
import sys
import inspect

def list_python_files(directory):
    """Returns a list of all .py files in the given directory and its subdirectories."""
    python_files = []
    for root, _, files in os.walk(directory):
        for file in files:
            if file.endswith(".py"):
                python_files.append(os.path.join(root, file))
    return python_files

def import_algorithm_subclasses(file_path):
    """Dynamically imports classes that subclass Algorithm from a given file if they exist."""
    module_name = os.path.splitext(os.path.basename(file_path))[0]
    spec = importlib.util.spec_from_file_location(module_name, file_path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    subclasses = []
    try:
        spec.loader.exec_module(module)
        for name, obj in inspect.getmembers(module, inspect.isclass):
            if issubclass(obj, module.Algorithm) and obj is not module.Algorithm:
                subclasses.append(obj)
    except Exception as e:
        print(f"Failed to import {module_name} from {file_path}: {e}")
    return subclasses

def import_all_algorithms(directory):
    """Imports classes that subclass Algorithm from all .py files in the specified directory."""
    python_files = list_python_files(directory)
    all_subclasses = []
    for file in python_files:
        subclasses = import_algorithm_subclasses(file)
        all_subclasses.extend(subclasses)
    return all_subclasses
